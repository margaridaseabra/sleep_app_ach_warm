function OUT = run_MA_fragmentation_from_overall(rows_overall, out_dir, varargin)
% RUN_MA_FRAGMENTATION_FROM_OVERALL
% -------------------------------------------------------------------------
% Compute a simple MA fragmentation index from the group OVERALL table:
%
%   FI_sleep = (# MA bouts) / (hours of SLEEP)
%
% where SLEEP = NREM + REM (baseline only).
%
% For each mouse (baseline):
%   - Count total MA bouts (n_bouts where state=="MA")
%   - Compute total SLEEP duration (NREM + REM, in hours)
%   - Fragmentation index = MA_bouts / sleep_hours
%
% Then:
%   - Compare WT vs APP (unpaired t-test, ranksum, Cohen's d)
%   - Plot a bar graph (mean ± SEM) with individual mice as dots
%
% INPUTS
%   rows_overall : OUT.overall from run_group_sleep_architecture
%   out_dir      : folder to save figure + CSV
%
% NAME–VALUE OPTIONS
%   'showIDs' : true/false, show mouse IDs next to dots (default: false)
%
% OUTPUT
%   OUT.data      : per-mouse table with FI_sleep
%   OUT.stats     : struct with WT vs APP stats
%   OUT.fig_file  : path to saved PNG
%   OUT.csv_file  : path to saved CSV
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

p = inputParser;
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
showIDs = p.Results.showIDs;

T = rows_overall;

% ---- 1) Baseline only ----
cond_str = lower(strtrim(T.condition));
T = T(cond_str == "baseline", :);

if isempty(T)
    warning('No baseline rows found in rows_overall. Nothing to do.');
    OUT = struct('success',false,'msg','no baseline data');
    return;
end

% Normalize to strings
T.state     = string(T.state);
T.mouse     = string(T.mouse);
T.genotype  = string(T.genotype);
T.condition = string(T.condition);

% We only care about states: MA, NREM, REM
hasMA   = any(T.state == "MA");
hasNREM = any(T.state == "NREM");
hasREM  = any(T.state == "REM");

if ~hasMA
    warning('No MA state found in baseline rows. Fragmentation index not defined.');
    OUT = struct('success',false,'msg','no MA state');
    return;
end
if ~(hasNREM || hasREM)
    warning('No NREM/REM states found in baseline rows. Fragmentation index not defined.');
    OUT = struct('success',false,'msg','no NREM/REM states');
    return;
end

% ---- 2) MA bouts per mouse (baseline) ----
T_MA = T(T.state == "MA", :);

if isempty(T_MA)
    warning('Baseline table has MA state label but no baseline MA rows?');
    OUT = struct('success',false,'msg','no MA rows');
    return;
end

G_MA = groupsummary(T_MA, {'mouse','genotype','condition'}, 'sum', 'n_bouts');
G_MA.Properties.VariableNames{end} = 'n_MA_bouts';

% ---- 3) Sleep duration (NREM + REM) per mouse ----
mask_sleep = (T.state == "NREM");
T_sleep = T(mask_sleep, :);

G_sleep = groupsummary(T_sleep, {'mouse','genotype','condition'}, 'sum', 'total_dur_s');
G_sleep.Properties.VariableNames{end} = 'sleep_dur_s';

% ---- 4) Join MA + sleep, compute FI_sleep ----
G = innerjoin(G_MA, G_sleep, ...
              'Keys', {'mouse','genotype','condition'});

G.sleep_h = G.sleep_dur_s / 3600;
G.FI_sleep = G.n_MA_bouts ./ G.sleep_h;   % MA bouts per hour of (NREM+REM) sleep

% ---- 5) Basic WT vs APP stats ----
geno = string(G.genotype);

FI_WT  = G.FI_sleep(geno=="WT");
FI_APP = G.FI_sleep(geno=="APP");

FI_WT  = FI_WT(~isnan(FI_WT));
FI_APP = FI_APP(~isnan(FI_APP));

nWT  = numel(FI_WT);
nAPP = numel(FI_APP);

p_t = NaN; p_rs = NaN; d = NaN;

if nWT >= 2 && nAPP >= 2
    [~, p_t] = ttest2(FI_WT, FI_APP, 'Vartype','unequal');
    p_rs     = ranksum(FI_WT, FI_APP);

    m1 = mean(FI_WT,'omitnan');  m2 = mean(FI_APP,'omitnan');
    s1 = std(FI_WT,'omitnan');   s2 = std(FI_APP,'omitnan');
    n1 = nWT;                    n2 = nAPP;
    sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d  = (m2 - m1) / sp;
end

fprintf('\n=== MA fragmentation index (baseline, MA bouts per hour of NREM+REM sleep) ===\n');
fprintf('WT  (n=%d): mean %.3f\n', nWT,  mean(FI_WT,'omitnan'));
fprintf('APP (n=%d): mean %.3f\n', nAPP, mean(FI_APP,'omitnan'));
fprintf('t-test p = %.4g, ranksum p = %.4g, Cohen d = %.2f (APP - WT)\n', ...
        p_t, p_rs, d);

% ---- 6) Plot bar + dots ----
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

figure('Color','w'); hold on;

Xpos = [1 2];
barWidth = 0.5;

mean_WT  = mean(FI_WT,'omitnan');
mean_APP = mean(FI_APP,'omitnan');
sem_WT   = std(FI_WT,'omitnan') / max(1,sqrt(nWT));
sem_APP  = std(FI_APP,'omitnan') / max(1,sqrt(nAPP));

if nWT>0
    bar(Xpos(1), mean_WT, barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
    errorbar(Xpos(1), mean_WT, sem_WT, 'k', 'LineStyle','none','LineWidth',1);
end
if nAPP>0
    bar(Xpos(2), mean_APP, barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
    errorbar(Xpos(2), mean_APP, sem_APP, 'k', 'LineStyle','none','LineWidth',1);
end

jit = 0.12;

% jittered dots WT
maskWT_rows  = (geno=="WT");
FI_WT_all    = G.FI_sleep(maskWT_rows);
mouse_WT_ids = G.mouse(maskWT_rows);

if ~isempty(FI_WT_all)
    xw = Xpos(1) + (rand(numel(FI_WT_all),1)-0.5)*2*jit;
    plot(xw, FI_WT_all, '.', 'Color', COL_WT_DOT, 'MarkerSize', 12);
    if showIDs
        for j = 1:numel(FI_WT_all)
            text(xw(j), FI_WT_all(j), char(mouse_WT_ids(j)), ...
                 'Rotation',45, ...
                 'HorizontalAlignment','left', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',8, ...
                 'Color',COL_WT_DOT);
        end
    end
end

% jittered dots APP
maskAPP_rows  = (geno=="APP");
FI_APP_all    = G.FI_sleep(maskAPP_rows);
mouse_APP_ids = G.mouse(maskAPP_rows);

if ~isempty(FI_APP_all)
    xa = Xpos(2) + (rand(numel(FI_APP_all),1)-0.5)*2*jit;
    plot(xa, FI_APP_all, '.', 'Color', COL_APP_DOT, 'MarkerSize', 12);
    if showIDs
        for j = 1:numel(FI_APP_all)
            text(xa(j), FI_APP_all(j), char(mouse_APP_ids(j)), ...
                 'Rotation',45, ...
                 'HorizontalAlignment','left', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',8, ...
                 'Color',COL_APP_DOT);
        end
    end
end

xlim([0.5 2.5]);
set(gca,'XTick',Xpos,'XTickLabel',{'WT','APP'},'FontSize',12);
ylabel('MA bouts per hour of (NREM+REM) sleep');
title('Baseline MA fragmentation index (WT vs APP)');
set(gca,'Box','off');

% annotate p-value, if available
% annotate p-value and significance stars, if available
y_max = max(G.FI_sleep,[],'omitnan');
if isempty(y_max) || isnan(y_max)
    y_max = 1;
end

% numeric text (p-values + Cohen's d)
if ~isnan(p_t)
    txt = sprintf('t-test p=%.3f; ranksum p=%.3f; d=%.2f', p_t, p_rs, d);
    text(1.5, y_max*1.05, txt, ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',10);
end

% significance stars based on t-test
starStr = '';
if ~isnan(p_t)
    if p_t < 0.001
        starStr = '***';
    elseif p_t < 0.01
        starStr = '**';
    elseif p_t < 0.05
        starStr = '*';
    end
end

if ~isempty(starStr)
    % horizontal line above the two bars
    y_star_line = y_max * 1.15;
    line([1 2], [y_star_line y_star_line], 'Color','k', 'LineWidth', 1.2);

    % star label slightly above the line
    text(1.5, y_max * 1.2, starStr, ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize', 14, ...
         'FontWeight','bold');
end

ylim([0, y_max*1.3]);

% ---- 7) Save figure + CSV ----
fig_file = fullfile(out_dir, 'baseline_MA_fragmentation_index_WTvsAPP.png');
saveas(gcf, fig_file);

csv_file = fullfile(out_dir, 'baseline_MA_fragmentation_index_per_mouse.csv');
writetable(G, csv_file);

% ---- 8) Pack output ----
STATS = struct();
STATS.nWT        = nWT;
STATS.nAPP       = nAPP;
STATS.mean_WT    = mean_WT;
STATS.mean_APP   = mean_APP;
STATS.sem_WT     = sem_WT;
STATS.sem_APP    = sem_APP;
STATS.p_ttest    = p_t;
STATS.p_ranksum  = p_rs;
STATS.Cohen_d    = d;

OUT = struct();
OUT.success   = true;
OUT.data      = G;
OUT.stats     = STATS;
OUT.fig_file  = fig_file;
OUT.csv_file  = csv_file;

fprintf('✅ MA fragmentation index figure saved to: %s\n', fig_file);
fprintf('✅ Per-mouse MA fragmentation index saved to: %s\n', csv_file);
end
