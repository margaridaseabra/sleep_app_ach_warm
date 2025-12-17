function OUT = plot_REM_survival_3h_baseline_vs_ambtemp(ALL_REM_3h, out_dir)
% plot_REM_survival_3h_baseline_vs_ambtemp
% -------------------------------------------------------------------------
% Uses ALL_REM_3h (3 h baseline + 3 h ambtemp) to plot REM duration
% "survival" curves:
%
%   For each genotype (WT, APP):
%       - baseline vs ambtemp survival curves
%       - stats: t-test + ranksum on bout durations
%
% This is not a full censored survival analysis; bouts all end, so we use
% S(t) = 1 - ECDF(t), and compare durations with tests.
%
% Colour/linestyle coding:
%   - WT  = black
%   - APP = blue
%   - baseline = solid line
%   - ambtemp (warming) = dashed line
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

R = ALL_REM_3h;
R.cond  = lower(strtrim(string(R.cond)));
R.geno  = string(R.geno);
R.mouse = string(R.mouse);

R = R(R.cond=="baseline" | R.cond=="ambtemp", :);
if isempty(R)
    warning('No baseline/ambtemp REM bouts found in ALL_REM_3h.');
    OUT = struct('success',false); return;
end

genotypes = ["WT","APP"];

% New colour coding: WT = black, APP = blue
COL_WT  = [0 0 0];
COL_APP = [0.2 0.5 0.9];

OUT = struct();
OUT.surv_stats = struct();

% limit x-axis at 95th percentile across all durations
all_dur = R.dur_s(~isnan(R.dur_s));
if isempty(all_dur)
    warning('No REM durations found.'); OUT.success=false; return;
end
xmax = prctile(all_dur, 95);

figure('Color','w','Position',[200 200 700 450]); hold on;

legend_entries = {};
for gi = 1:numel(genotypes)
    gtype = genotypes(gi);
    G = R(R.geno==gtype,:);
    if isempty(G), continue; end

    d_base = G.dur_s(G.cond=="baseline");
    d_amb  = G.dur_s(G.cond=="ambtemp");

    d_base = d_base(~isnan(d_base));
    d_amb  = d_amb(~isnan(d_amb));

    if isempty(d_base) || isempty(d_amb)
        continue;
    end

    % ECDF -> survival
    [F_base,x_base] = ecdf(d_base);
    [F_amb, x_amb]  = ecdf(d_amb);

    S_base = 1 - F_base;
    S_amb  = 1 - F_amb;

    % keep only up to xmax
    keepB = x_base <= xmax;
    keepA = x_amb <= xmax;

    % Choose colour by genotype
    if gtype=="WT"
        col = COL_WT;   % black
    else
        col = COL_APP;  % blue
    end

    % Line style by condition:
    ls_base = '-';   % baseline solid
    ls_amb  = '--';  % warming (ambtemp) dashed

    % Plot baseline
    stairs(x_base(keepB), S_base(keepB), ...
        'LineStyle', ls_base, ...
        'Color',     col, ...
        'LineWidth', 1.5);
    legend_entries{end+1} = sprintf('%s baseline', gtype); %#ok<AGROW>

    % Plot ambtemp (warming)
    stairs(x_amb(keepA), S_amb(keepA), ...
        'LineStyle', ls_amb, ...
        'Color',     col, ...
        'LineWidth', 1.5);
    legend_entries{end+1} = sprintf('%s ambtemp', gtype); %#ok<AGROW>

    % simple stats on durations
    [~, p_t]  = ttest2(d_base, d_amb, 'Vartype','unequal');
    p_rs      = ranksum(d_base, d_amb);

    m1 = mean(d_base); m2 = mean(d_amb);
    s1 = std(d_base);  s2 = std(d_amb);
    n1 = numel(d_base); n2 = numel(d_amb);
    sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d  = (m2 - m1)/sp;

    OUT.surv_stats.(gtype).p_ttest   = p_t;
    OUT.surv_stats.(gtype).p_ranksum = p_rs;
    OUT.surv_stats.(gtype).Cohen_d   = d;
    OUT.surv_stats.(gtype).mean_base = m1;
    OUT.surv_stats.(gtype).mean_amb  = m2;

    fprintf('\n[REM duration 3 h] %s, baseline vs ambtemp:\n', gtype);
    fprintf('  mean base = %.2f s, mean amb = %.2f s\n', m1, m2);
    fprintf('  t-test p=%.3g, ranksum p=%.3g, d=%.2f\n', p_t, p_rs, d);
end

xlim([0 xmax]);
ylim([0 1]);
xlabel('REM bout duration (s)');
ylabel('Survival S(t) = 1 - F(t)');
title('REM duration 3 h windows: baseline vs ambtemp');

if ~isempty(legend_entries)
    legend(legend_entries, 'Location','southwest');
end
set(gca,'Box','off','FontSize',11);

fig_file = fullfile(out_dir, 'REM3h_survival_baseline_vs_ambtemp.png');
saveas(gcf, fig_file);

OUT.fig_file = fig_file;
OUT.success  = true;

% also write summary CSV
genos = fieldnames(OUT.surv_stats);
if ~isempty(genos)
    rows = cell(numel(genos), 7);
    for i = 1:numel(genos)
        g = genos{i};
        Sg = OUT.surv_stats.(g);
        rows(i,:) = {g, Sg.mean_base, Sg.mean_amb, ...
                     Sg.p_ttest, Sg.p_ranksum, Sg.Cohen_d, xmax};
    end
    tbl = cell2table(rows, 'VariableNames', ...
        {'Genotype','Mean_dur_base_s','Mean_dur_amb_s', ...
         'p_ttest','p_ranksum','Cohen_d','xmax_s'});
    writetable(tbl, fullfile(out_dir, 'REM3h_survival_stats.csv'));
end

fprintf('✅ REM 3 h survival-style plots + stats saved in %s\n', out_dir);
end
