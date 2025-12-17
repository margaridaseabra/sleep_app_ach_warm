function OUT = make_baseline_boutduration_APPvsWT(PERHOUR, out_dir, states_to_use, varargin)
% make_baseline_boutduration_APPvsWT
% -------------------------------------------------------------------------
% Baseline-only plot of:
%   Y: mean bout duration (seconds) per state
%   X: states (e.g., WK, MA, NREM, REM)
%
% For each mouse and state:
%   mean_bout_dur = total_dur_s_across_hours / total_number_of_bouts
%   Here we approximate total number of bouts as sum(bouts_per_h) over hours.
%
% Genotypes:
%   - WT
%   - APP (all non-WT collapsed to APP)
%
% Plot (single figure):
%   - X: states
%   - Bars: WT vs APP (2 bars per state) with mean ± SEM
%   - Dots: individual mice (optional IDs)
%   - Stars: per-state WT vs APP difference (t-test, BH–FDR corrected)
%
% Stats per state:
%   - nWT, nAPP
%   - MeanWT, MeanAPP
%   - SD_WT, SD_APP
%   - CV_WT, CV_APP
%   - SD_ratio (APP / WT)
%   - Welch t-test (p_ttest)
%   - Mann–Whitney (p_ranksum)
%   - F-test on variances (p_var_Ftest, exploratory)
%   - Cohen's d (APP - WT)
%   - BH–FDR corrected t-test p (p_ttest_FDR) across states
%
% INPUT
%   PERHOUR       : OUT.per_hour table from run_group_sleep_architecture
%   out_dir       : folder to save the figure (default: pwd)
%   states_to_use : string/cell array, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'showIDs'             : true/false, show mouse IDs next to dots (default: true)
%   'useFDRforStars'      : true/false, use FDR p-values for stars (default: true)
%   'minNperGroupForStats': minimum n per genotype to draw a star (default: 3)
%
% OUTPUT
%   OUT.success      : logical
%   OUT.out_file     : PNG filename
%   OUT.state_order  : states actually plotted
%   OUT.stats        : stats table (one row per state)
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
if nargin < 3 || isempty(states_to_use)
    states_to_use = ["WK","MA","NREM","REM"];
end
states_to_use = string(states_to_use(:)).';

% ---------- optional args ----------
p = inputParser;
addParameter(p,'showIDs',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'useFDRforStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

showIDs              = p.Results.showIDs;
useFDRforStars       = p.Results.useFDRforStars;
minNperGroupForStats = p.Results.minNperGroupForStats;

PH = PERHOUR;

% ---------- 1) Keep only BASELINE ----------
if ~ismember('condition', PH.Properties.VariableNames)
    error('PERHOUR must contain a "condition" column.');
end
cond_str    = lower(strtrim(PH.condition));
is_baseline = cond_str == "baseline";
PH = PH(is_baseline, :);

if isempty(PH)
    warning('No baseline rows in PERHOUR. Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline data');
    return;
end

% ---------- 2) Check required columns ----------
needed = {'dur_s','bouts_per_h','state','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('PERHOUR table missing required columns for make_baseline_boutduration_APPvsWT.');
end

% Normalize text columns
PH.state    = string(PH.state);
PH.mouse    = string(PH.mouse);
PH.genotype = string(PH.genotype);

% Collapse any non-WT to APP
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% ---------- 3) Restrict to desired states ----------
avail_states = unique(PH.state);
state_order  = states_to_use(ismember(states_to_use, avail_states));

if isempty(state_order)
    warning('None of the requested states found in PERHOUR baseline data.');
    OUT = struct('success', false, 'msg', 'no states found');
    return;
end

% ---------- 4) Per-mouse mean bout duration per state ----------
[gid, m, gen, st] = findgroups(PH.mouse, PH.geno_group, PH.state);

tot_dur   = splitapply(@(x) sum(double(x),'omitnan'), PH.dur_s,       gid);
tot_bouts = splitapply(@(x) sum(double(x),'omitnan'), PH.bouts_per_h, gid);

mean_bout = tot_dur ./ max(tot_bouts, 1);
mean_bout(tot_bouts == 0) = NaN;

Tmouse = table(m, gen, st, mean_bout, ...
    'VariableNames', {'mouse','geno','state','mean_bout_dur_s'});

% Keep only states we care about
Tmouse = Tmouse(ismember(Tmouse.state, state_order), :);

if isempty(Tmouse)
    warning('No per-mouse bout durations for requested states.');
    OUT = struct('success', false, 'msg', 'no per-mouse data');
    return;
end

% ---------- Colors ----------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

% ---------- 5) Compute stats per state ----------
nStates = numel(state_order);

meanWT   = nan(1, nStates);
meanAPP  = nan(1, nStates);
sdWT     = nan(1, nStates);
sdAPP    = nan(1, nStates);
cvWT     = nan(1, nStates);
cvAPP    = nan(1, nStates);
sd_ratio = nan(1, nStates);
p_t      = nan(1, nStates);
p_rs     = nan(1, nStates);
p_var    = nan(1, nStates);
d_eff    = nan(1, nStates);
nWT_vec  = nan(1, nStates);
nAPP_vec = nan(1, nStates);
max_y    = nan(1, nStates);   % for star placement

for i = 1:nStates
    st_i = state_order(i);
    Ti   = Tmouse(Tmouse.state == st_i, :);
    if isempty(Ti), continue; end

    valsWT  = Ti.mean_bout_dur_s(Ti.geno=="WT");
    valsAPP = Ti.mean_bout_dur_s(Ti.geno=="APP");

    valsWT  = valsWT(~isnan(valsWT));
    valsAPP = valsAPP(~isnan(valsAPP));

    nWT_vec(i)  = numel(valsWT);
    nAPP_vec(i) = numel(valsAPP);

    if ~isempty(valsWT)
        meanWT(i) = mean(valsWT);
        sdWT(i)   = std(valsWT);
        cvWT(i)   = sdWT(i) / max(eps, meanWT(i));
    end
    if ~isempty(valsAPP)
        meanAPP(i) = mean(valsAPP);
        sdAPP(i)   = std(valsAPP);
        cvAPP(i)   = sdAPP(i) / max(eps, meanAPP(i));
    end

    sd_ratio(i) = sdAPP(i) / sdWT(i);

    if (nWT_vec(i) >= minNperGroupForStats) && (nAPP_vec(i) >= minNperGroupForStats)
        % Welch t-test
        [~, p_t(i)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
        % Mann–Whitney
        p_rs(i) = ranksum(valsWT, valsAPP);
        % F-test for variances (exploratory)
        try
            [~, pF] = vartest2(valsWT, valsAPP, 'Tail','both');
            p_var(i) = pF;
        catch
            p_var(i) = NaN;
        end
        % Cohen's d (APP - WT)
        m1 = mean(valsWT);  m2 = mean(valsAPP);
        s1 = std(valsWT);   s2 = std(valsAPP);
        n1 = numel(valsWT); n2 = numel(valsAPP);
        sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
        d_eff(i) = (m2 - m1) / sp;
    end

    if ~isempty([valsWT; valsAPP])
        max_y(i) = max([valsWT; valsAPP]);
    end
end

% ---------- 6) BH–FDR correction across STATES ----------
p_t_fdr = nan(size(p_t));
valid   = ~isnan(p_t);
pvals   = p_t(valid);
if ~isempty(pvals)
    [sorted_p, sort_idx] = sort(pvals(:));
    m = numel(sorted_p);
    adj = sorted_p .* (m ./ (1:m)');        % BH
    for k = m-1:-1:1
        adj(k) = min(adj(k), adj(k+1));
    end
    adj(adj>1) = 1;
    p_fdr_vals = nan(size(pvals));
    p_fdr_vals(sort_idx) = adj;
    p_t_fdr(valid) = p_fdr_vals;
end

% ---------- 7) Plot ----------
figure('Color','w'); hold on;
x = 1:nStates;
jitterFrac = 0.25;
y_offset   = 0.05;   % vertical offset as fraction of y-range (set later)

% Bars
barWidth = 0.4;
hasWT  = any(~isnan(meanWT));
hasAPP = any(~isnan(meanAPP));

if hasWT
    bar(x - barWidth/2, meanWT, barWidth, ...
        'FaceColor', COL_WT, 'EdgeColor', 'none');
end
if hasAPP
    bar(x + barWidth/2, meanAPP, barWidth, ...
        'FaceColor', COL_APP, 'EdgeColor', 'none');
end

% Error bars (SEM)
semWT  = sdWT  ./ sqrt(max(1,nWT_vec));
semAPP = sdAPP ./ sqrt(max(1,nAPP_vec));

if hasWT
    errorbar(x - barWidth/2, meanWT, semWT, 'k', ...
        'LineStyle','none','LineWidth',1);
end
if hasAPP
    errorbar(x + barWidth/2, meanAPP, semAPP, 'k', ...
        'LineStyle','none','LineWidth',1);
end

% Dots + optional IDs
for i = 1:nStates
    st_i = state_order(i);
    Ti   = Tmouse(Tmouse.state == st_i, :);
    if isempty(Ti), continue; end

    % WT
    maskWT = (Ti.geno == "WT");
    valsWT = Ti.mean_bout_dur_s(maskWT);
    mWT    = Ti.mouse(maskWT);
    if ~isempty(valsWT)
        xw = (x(i) - barWidth/2) + (rand(size(valsWT)) - 0.5) * barWidth * jitterFrac;
        plot(xw, valsWT, '.', 'Color', COL_WT_DOT, 'MarkerSize', 10);
        if showIDs
            for j = 1:numel(valsWT)
                thisID = char(mWT(j));
                text(xw(j), valsWT(j), thisID, ...
                    'Rotation', 45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, ...
                    'Color', COL_WT_DOT);
            end
        end
    end

    % APP
    maskAPP = (Ti.geno == "APP");
    valsAPP = Ti.mean_bout_dur_s(maskAPP);
    mAPP    = Ti.mouse(maskAPP);
    if ~isempty(valsAPP)
        xa = (x(i) + barWidth/2) + (rand(size(valsAPP)) - 0.5) * barWidth * jitterFrac;
        plot(xa, valsAPP, '.', 'Color', COL_APP_DOT, 'MarkerSize', 10);
        if showIDs
            for j = 1:numel(valsAPP)
                thisID = char(mAPP(j));
                text(xa(j), valsAPP(j), thisID, ...
                    'Rotation', 45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, ...
                    'Color', COL_APP_DOT);
            end
        end
    end
end

% Y-limits & star positioning
if any(~isnan(max_y))
    global_max_y = max(max_y(~isnan(max_y)));
else
    global_max_y = max(Tmouse.mean_bout_dur_s,[],'omitnan');
end
if isempty(global_max_y) || isnan(global_max_y)
    global_max_y = 1;
end
y_top = global_max_y * 1.3;
ylim([0, y_top]);

y_offset_abs = 0.03 * y_top;  % small offset for text

% Significance bars + stars
for i = 1:nStates
    if useFDRforStars
        p_here = p_t_fdr(i);
    else
        p_here = p_t(i);
    end
    if isnan(p_here) || p_here >= 0.05
        continue;
    end
    if nWT_vec(i) < minNperGroupForStats || nAPP_vec(i) < minNperGroupForStats
        continue;
    end

    % decide number of stars
    if p_here < 0.001
        stars = '***';
    elseif p_here < 0.01
        stars = '**';
    else
        stars = '*';
    end

    y_star = max_y(i);
    if isnan(y_star)
        y_star = global_max_y * 0.8;
    end
    y_star = y_star + y_offset_abs;

    line([x(i)-barWidth/2, x(i)+barWidth/2], [y_star, y_star], ...
         'Color','k','LineWidth',1.0);
    text(x(i), y_star + y_offset_abs, stars, ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',14, ...
         'FontWeight','bold');
end

% Axes & labels
xticks(x);
xticklabels(state_order);
xlabel('State');
ylabel('Mean bout duration (s)');
title('Baseline: mean bout duration per state (WT vs APP)');
set(gca,'Box','off','FontSize',12);

if hasWT && hasAPP
    legend({'WT mean \pm SEM','APP mean \pm SEM'}, ...
        'Location','northoutside','Orientation','horizontal');
end

% ---------- 8) Stats table ----------
stats_idx = ~isnan(nWT_vec + nAPP_vec);  % any state with at least some data
if any(stats_idx)
    stats_tbl = table( ...
        state_order(stats_idx)', ...
        nWT_vec(stats_idx)', ...
        nAPP_vec(stats_idx)', ...
        meanWT(stats_idx)', ...
        meanAPP(stats_idx)', ...
        sdWT(stats_idx)', ...
        sdAPP(stats_idx)', ...
        cvWT(stats_idx)', ...
        cvAPP(stats_idx)', ...
        sd_ratio(stats_idx)', ...
        p_t(stats_idx)', ...
        p_t_fdr(stats_idx)', ...
        p_rs(stats_idx)', ...
        p_var(stats_idx)', ...
        d_eff(stats_idx)', ...
        'VariableNames', {'State','nWT','nAPP', ...
                          'MeanWT','MeanAPP', ...
                          'SD_WT','SD_APP', ...
                          'CV_WT','CV_APP', ...
                          'SD_ratio_APP_vs_WT', ...
                          'p_ttest','p_ttest_FDR', ...
                          'p_ranksum','p_var_Ftest','Cohen_d'});
    fprintf('\nBaseline mean bout duration: WT vs APP stats (per state)\n');
    disp(stats_tbl);
else
    stats_tbl = table();
    fprintf('\n[Stats] No states with enough WT & APP data for tests.\n');
end

% ---------- 9) Save figure ----------
if showIDs
    id_suffix = '_withIDs';
else
    id_suffix = '_noIDs';
end
fname    = sprintf('baseline_boutduration_APPvsWT%s.png', id_suffix);
out_file = fullfile(out_dir, fname);
saveas(gcf, out_file);

% ---------- OUT ----------
OUT = struct();
OUT.success      = true;
OUT.out_file     = out_file;
OUT.state_order  = state_order;
OUT.stats        = stats_tbl;

fprintf('✅ Baseline bout duration state plot saved to: %s\n', out_file);
end
