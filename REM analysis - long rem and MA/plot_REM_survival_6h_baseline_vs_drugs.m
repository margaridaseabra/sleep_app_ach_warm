function OUT = plot_REM_survival_6h_baseline_vs_drugs(ALL_REM_6h, out_dir)
% plot_REM_survival_6h_baseline_vs_drugs
% -------------------------------------------------------------------------
% Uses ALL_REM_6h (6 h baseline + 6 h drugs) to plot REM duration
% "survival" curves:
%
%   For each genotype (WT, APP):
%       - baseline vs drugs survival curves
%       - stats: t-test + ranksum on bout durations
%
% This is not a full censored survival analysis; bouts all end, so we use
% S(t) = 1 - ECDF(t), and compare durations with tests.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

R = ALL_REM_6h;
R.cond  = lower(strtrim(string(R.cond)));
R.geno  = string(R.geno);
R.mouse = string(R.mouse);

% keep only baseline + drugs
R = R(R.cond=="baseline" | R.cond=="drugs", :);
if isempty(R)
    warning('No baseline/drugs REM bouts found in ALL_REM_6h.');
    OUT = struct('success',false); 
    return;
end

genotypes = ["WT","APP"];

% Colors & styles:
COL_WT        = [0 0 0];                     % black
COL_APP_BASE  = [100 149 237] / 255;        % cornflower blue
COL_APP_DRUG  = [ 38  72 128] / 255;        % dark cornflower blue (approx)

OUT = struct();
OUT.surv_stats = struct();

% limit x-axis at 95th percentile across all durations
all_dur = R.dur_s(~isnan(R.dur_s));
if isempty(all_dur)
    warning('No REM durations found.'); 
    OUT.success=false; 
    return;
end
xmax = prctile(all_dur, 95);

figure('Color','w','Position',[200 200 700 450]); hold on;

legend_entries = {};
for gi = 1:numel(genotypes)
    gtype = genotypes(gi);
    G = R(R.geno==gtype,:);
    if isempty(G), continue; end

    d_base = G.dur_s(G.cond=="baseline");
    d_drugs = G.dur_s(G.cond=="drugs");

    d_base = d_base(~isnan(d_base));
    d_drugs = d_drugs(~isnan(d_drugs));

    if isempty(d_base) || isempty(d_drugs)
        continue;
    end

    % ECDF -> survival
    [F_base,x_base] = ecdf(d_base);
    [F_drugs,x_drugs] = ecdf(d_drugs);

    S_base = 1 - F_base;
    S_drugs = 1 - F_drugs;

    % keep only up to xmax
    keepB = x_base <= xmax;
    keepD = x_drugs <= xmax;

    % Line colors and styles per genotype/condition
    if gtype=="WT"
        % WT baseline: black solid
        % WT drugs:   black dotted
        colB = COL_WT;
        colD = COL_WT;
        lsB  = '-';    % solid
        lsD  = ':';    % dotted

    elseif gtype=="APP"
        % APP baseline: cornflower blue solid
        % APP drugs:   dark cornflower blue dotted
        colB = COL_APP_BASE;
        colD = COL_APP_DRUG;
        lsB  = '-';    % solid
        lsD  = ':';    % dotted

    else
        % Fallback (in case of unexpected genotype)
        colB = [0.5 0.5 0.5];
        colD = [0.2 0.2 0.2];
        lsB  = '-';
        lsD  = ':';
    end

    stairs(x_base(keepB), S_base(keepB), 'LineStyle',lsB, ...
        'Color', colB, 'LineWidth',1.5);
    legend_entries{end+1} = sprintf('%s baseline', gtype); %#ok<AGROW>

    stairs(x_drugs(keepD), S_drugs(keepD), 'LineStyle',lsD, ...
        'Color', colD, 'LineWidth',1.5);
    legend_entries{end+1} = sprintf('%s drugs', gtype); %#ok<AGROW>

    % simple stats on durations
    [~, p_t]  = ttest2(d_base, d_drugs, 'Vartype','unequal');
    p_rs      = ranksum(d_base, d_drugs);

    m1 = mean(d_base); m2 = mean(d_drugs);
    s1 = std(d_base);  s2 = std(d_drugs);
    n1 = numel(d_base); n2 = numel(d_drugs);
    sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d  = (m2 - m1)/sp;

    OUT.surv_stats.(gtype).p_ttest    = p_t;
    OUT.surv_stats.(gtype).p_ranksum  = p_rs;
    OUT.surv_stats.(gtype).Cohen_d    = d;
    OUT.surv_stats.(gtype).mean_base  = m1;
    OUT.surv_stats.(gtype).mean_drugs = m2;

    fprintf('\n[REM duration 6 h] %s, baseline vs drugs:\n', gtype);
    fprintf('  mean base = %.2f s, mean drugs = %.2f s\n', m1, m2);
    fprintf('  t-test p=%.3g, ranksum p=%.3g, d=%.2f\n', p_t, p_rs, d);
end

xlim([0 xmax]);
ylim([0 1]);
xlabel('REM bout duration (s)');
ylabel('Survival S(t) = 1 - F(t)');
title('REM duration 6 h windows: baseline vs drugs');

if ~isempty(legend_entries)
    legend(legend_entries, 'Location','southwest');
end
set(gca,'Box','off','FontSize',11);

fig_file = fullfile(out_dir, 'REM6h_survival_baseline_vs_drugs.png');
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
        rows(i,:) = {g, Sg.mean_base, Sg.mean_drugs, ...
                     Sg.p_ttest, Sg.p_ranksum, Sg.Cohen_d, xmax};
    end
    tbl = cell2table(rows, 'VariableNames', ...
        {'Genotype','Mean_dur_base_s','Mean_dur_drugs_s', ...
         'p_ttest','p_ranksum','Cohen_d','xmax_s'});
    writetable(tbl, fullfile(out_dir, 'REM6h_survival_stats.csv'));
end

fprintf('✅ REM 6 h survival-style plots + stats saved in %s\n', out_dir);
end
