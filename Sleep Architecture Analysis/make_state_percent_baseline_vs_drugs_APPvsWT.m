function OUT = make_state_percent_baseline_vs_drugs_APPvsWT(rows_overall, out_dir, varargin)
% make_state_percent_baseline_vs_drugs_APPvsWT
%- -------------------------------------------------------------------------
% Compare 6 h BASELINE vs 6 h DRUGS, keeping WT vs APP:
%
%   Y: % of 6 h window spent in each state
%   Factors:
%       - Condition: baseline vs drugs  (within mouse, paired)
%       - Genotype : WT vs APP           (between mice)
%
% For each state (WK, MA, NREM, REM, SLEEP):
%   - Compute per-mouse %time in state for baseline and drugs
%   - Paired t-test (baseline vs drugs) separately for WT and APP
%   - Compare change Δ = drugs - baseline between WT and APP (unpaired t-test)
%
% Plot:
%   - 2 panels: top WT, bottom APP
%   - x-axis: states
%   - per mouse: line from baseline to drugs
%   - group means ± SEM as big markers
%   - stars for paired differences in each genotype
%
% INPUT
%   rows_overall : OVERALL table from run_group_sleep_architecture
%                  (but run on cropped 6h baseline + 6h drugs files)
%   out_dir      : folder to save figure + stats CSV
%
% OPTIONS
%   'minNperGroupForStats' : minimum n WT or APP to run stats (default: 3)
%   'useFDR'               : apply BH-FDR across states per genotype (default: true)
%
% OUTPUT
%   OUT.fig_file          : PNG path
%   OUT.within_stats      : table of paired test results (per state, per genotype)
%   OUT.delta_stats       : table of Δ comparisons WT vs APP
%   OUT.state_order       : order of states on x-axis
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

p = inputParser;
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});

minN = p.Results.minNperGroupForStats;
useFDR = p.Results.useFDR;

T = rows_overall;

% ---------- 1) Keep only baseline + drugs ----------
cond_str = lower(strtrim(string(T.condition)));
keep_cond = cond_str == "baseline" | cond_str == "drugs";
T = T(keep_cond, :);

if isempty(T)
    warning('No baseline/drugs rows found in rows_overall. Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline/drugs data');
    return;
end

% Normalise to string types
T.mouse     = string(T.mouse);
T.genotype  = string(T.genotype);
T.condition = lower(strtrim(string(T.condition)));
T.state     = string(T.state);

% ---------- 2) Base states and composite SLEEP ----------
base_desired = ["WK","MA","NREM","REM"];
avail_states = unique(T.state);
base_states  = base_desired(ismember(base_desired, avail_states));

if isempty(base_states)
    warning('No WK/MA/NREM/REM states found.');
    OUT = struct('success', false, 'msg', 'no base states');
    return;
end

% Per-mouse, per-cond, per-state durations
Gbase = groupsummary(T, {'mouse','genotype','condition','state'}, 'sum', 'total_dur_s');
Gbase.Properties.VariableNames{end} = 'dur_s';
Gbase = Gbase(ismember(Gbase.state, base_states), :);

% Total duration per mouse & condition across all base states
Gtot = groupsummary(Gbase, {'mouse','genotype','condition'}, 'sum', 'dur_s');
Gtot.Properties.VariableNames{end} = 'total_dur_allstates';

% Build SLEEP = NREM + REM per mouse & condition
sleep_components = ["NREM","REM"];
Gsleep_src = Gbase(ismember(Gbase.state, sleep_components), :);
if ~isempty(Gsleep_src)
    Gsleep = groupsummary(Gsleep_src, {'mouse','genotype','condition'}, 'sum', 'dur_s');
    Gsleep.Properties.VariableNames{end} = 'dur_s';
    Gsleep.state = repmat("SLEEP", height(Gsleep), 1);
else
    Gsleep = Gbase([],:);
end

% Combine base states + SLEEP
G = [Gbase; Gsleep];

% Attach total duration and compute %
G = innerjoin(G, Gtot, 'Keys', {'mouse','genotype','condition'});
G.pct = 100 * G.dur_s ./ G.total_dur_allstates;

% Final state order
desired_order = ["WK","MA","NREM","REM","SLEEP"];
avail_states2 = unique(G.state);
state_order   = desired_order(ismember(desired_order, avail_states2));

if isempty(state_order)
    warning('No usable states after adding SLEEP.');
    OUT = struct('success', false, 'msg', 'no states after SLEEP build');
    return;
end

% ---------- 3) Within-genotype paired stats (baseline vs drugs) ----------
genotypes = ["WT","APP"];
nStates = numel(state_order);
nGen    = numel(genotypes);

within_results = [];

% Store p-values for FDR
p_within = nan(nGen, nStates);

for gi = 1:nGen
    gtype = genotypes(gi);
    for si = 1:nStates
        st = state_order(si);

        Gst = G(G.genotype == gtype & G.state == st, :);
        if isempty(Gst), continue; end

        % pivot: one row per mouse, columns: baseline, drugs
        Twide = unstack(Gst(:,{'mouse','condition','pct'}), 'pct', 'condition');

        if ~all(ismember({'baseline','drugs'}, Twide.Properties.VariableNames))
            % some mice missing one condition
            continue;
        end

        x = Twide.baseline;
        y = Twide.drugs;
        ok = ~isnan(x) & ~isnan(y);
        x = x(ok);
        y = y(ok);
        mIDs = Twide.mouse(ok);

        nPairs = numel(x);
        if nPairs < minN
            continue;
        end

        [~, p_pair] = ttest(x, y);    % baseline vs drugs
        delta = y - x;

        mean_base  = mean(x, 'omitnan');
        mean_amb   = mean(y, 'omitnan');
        mean_delta = mean(delta, 'omitnan');
        sd_base    = std(x, 'omitnan');
        sd_amb     = std(y, 'omitnan');
        sd_delta   = std(delta, 'omitnan');

        p_within(gi, si) = p_pair;

        within_results = [within_results; ...
            {char(st), char(gtype), nPairs, ...
             mean_base, sd_base, ...
             mean_amb,  sd_amb, ...
             mean_delta, sd_delta, ...
             p_pair}];
    end
end

if isempty(within_results)
    warning('No within-genotype paired comparisons could be computed.');
    within_tbl = table();
else
    within_tbl = cell2table(within_results, ...
        'VariableNames', {'State','Genotype','n', ...
                          'Mean_baseline','SD_baseline', ...
                          'Mean_drugs','SD_drugs', ...
                          'Mean_delta','SD_delta', ...
                          'p_paired'});
end

% ---------- 4) FDR across states for each genotype ----------
p_within_fdr = nan(size(p_within));
if useFDR
    for gi = 1:nGen
        pv = p_within(gi,:);
        valid = ~isnan(pv);
        if ~any(valid), continue; end

        pvals = pv(valid);
        [sorted_p, sort_idx] = sort(pvals(:));
        m = numel(sorted_p);
        adj = sorted_p .* (m ./ (1:m)');
        for k = m-1:-1:1
            adj(k) = min(adj(k), adj(k+1));
        end
        adj(adj>1) = 1;
        p_adj = nan(size(pvals));
        p_adj(sort_idx) = adj;
        pv_fdr = nan(size(pv));
        pv_fdr(valid) = p_adj;

        p_within_fdr(gi,:) = pv_fdr;
    end
end

% attach FDR to within_tbl
if ~isempty(within_tbl)
    p_fdr_vec = nan(height(within_tbl),1);
    for r = 1:height(within_tbl)
        st = string(within_tbl.State(r));
        gt = string(within_tbl.Genotype(r));
        gi = find(genotypes == gt, 1);
        si = find(state_order == st, 1);
        if ~isempty(gi) && ~isempty(si)
            p_fdr_vec(r) = p_within_fdr(gi, si);
        end
    end
    within_tbl.p_paired_FDR = p_fdr_vec;
end

% ---------- 5) Δ comparison between genotypes (APP vs WT) ----------
delta_results = [];

for si = 1:nStates
    st = state_order(si);

    Gst = G(G.state == st, :);
    if isempty(Gst), continue; end

    Twide_all = unstack(Gst(:,{'mouse','genotype','condition','pct'}), ...
                        'pct','condition');
    if ~all(ismember({'baseline','drugs'}, Twide_all.Properties.VariableNames))
        continue;
    end

    Twide_all.delta = Twide_all.drugs - Twide_all.baseline;

    delta_WT  = Twide_all.delta(Twide_all.genotype=="WT");
    delta_APP = Twide_all.delta(Twide_all.genotype=="APP");

    delta_WT  = delta_WT(~isnan(delta_WT));
    delta_APP = delta_APP(~isnan(delta_APP));

    nWT = numel(delta_WT);
    nAPP = numel(delta_APP);

    if nWT >= minN && nAPP >= minN
        [~, p_delta] = ttest2(delta_WT, delta_APP, 'Vartype','unequal');

        mWT  = mean(delta_WT,'omitnan');
        mAPP = mean(delta_APP,'omitnan');
        sWT  = std(delta_WT,'omitnan');
        sAPP = std(delta_APP,'omitnan');

        % Cohen's d on delta (APP - WT)
        n1 = nWT; n2 = nAPP;
        sp = sqrt(((n1-1)*sWT^2 + (n2-1)*sAPP^2) / max(1,(n1+n2-2)));
        d  = (mAPP - mWT) / sp;
    else
        p_delta = NaN;
        mWT = NaN; mAPP = NaN; sWT = NaN; sAPP = NaN; d = NaN;
    end

    delta_results = [delta_results; ...
        {char(st), nWT, nAPP, mWT, sWT, mAPP, sAPP, p_delta, d}];
end

if isempty(delta_results)
    delta_tbl = table();
else
    delta_tbl = cell2table(delta_results, ...
        'VariableNames', {'State','nWT','nAPP', ...
                          'Mean_delta_WT','SD_delta_WT', ...
                          'Mean_delta_APP','SD_delta_APP', ...
                          'p_delta_WTvsAPP','Cohen_d_delta'});
end

% ---------- 6) Plot: paired lines per genotype ----------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.5 0.65]);

for gi = 1:nGen
    gtype = genotypes(gi);

    subplot(2,1,gi); hold on;
    title(sprintf('%s: 6 h baseline vs 6 h drugs', gtype));
    x = 1:nStates;

    y_all = [];

    for si = 1:nStates
        st = state_order(si);
        Gst = G(G.genotype == gtype & G.state == st, :);
        if isempty(Gst), continue; end

        Twide = unstack(Gst(:,{'mouse','condition','pct'}), 'pct','condition');
        if ~all(ismember({'baseline','drugs'}, Twide.Properties.VariableNames))
            continue;
        end
        base = Twide.baseline;
        amb  = Twide.drugs;
        ok   = ~isnan(base) & ~isnan(amb);
        base = base(ok);
        amb  = amb(ok);

        nPairs = numel(base);
        if nPairs == 0
            continue;
        end

        % store for y-limits
        y_all = [y_all; base; amb];

        % jitter each mouse slightly around x(si)
        jitter = 0.12;
        x_base = x(si) - 0.15 + (rand(size(base))-0.5)*2*jitter;
        x_amb  = x(si) + 0.15 + (rand(size(amb))-0.5)*2*jitter;

        col = (gtype=="WT") * COL_WT + (gtype=="APP") * COL_APP;

        for j = 1:numel(base)
            plot([x_base(j), x_amb(j)], [base(j), amb(j)], '-', ...
                 'Color', [col 0.4]);  % faint line
        end
        plot(x_base, base, 'o', 'MarkerFaceColor', col, ...
             'MarkerEdgeColor','k', 'MarkerSize',5);
        plot(x_amb,  amb,  'o', 'MarkerFaceColor', col, ...
             'MarkerEdgeColor','k', 'MarkerSize',5);

        % group means ± SEM
        m_base = mean(base,'omitnan');
        m_amb  = mean(amb,'omitnan');
        sem_base = std(base,'omitnan') / sqrt(numel(base));
        sem_amb  = std(amb,'omitnan')  / sqrt(numel(amb));

        errorbar(x(si)-0.18, m_base, sem_base, 'k','LineStyle','none','LineWidth',1.2);
        errorbar(x(si)+0.18, m_amb,  sem_amb,  'k','LineStyle','none','LineWidth',1.2);

        plot(x(si)-0.18, m_base, 's', 'MarkerSize',7, ...
             'MarkerFaceColor',col, 'MarkerEdgeColor','k');
        plot(x(si)+0.18, m_amb,  's', 'MarkerSize',7, ...
             'MarkerFaceColor',col, 'MarkerEdgeColor','k');

        % stars for this genotype/state
        if ~isempty(within_tbl)
            row = within_tbl.State == st & within_tbl.Genotype == gtype;
            if any(row)
                p_here  = within_tbl.p_paired(row);
                p_hereF = within_tbl.p_paired_FDR(row);

                if useFDR && ~isnan(p_hereF)
                    p_use = p_hereF;
                else
                    p_use = p_here;
                end

                if ~isnan(p_use) && p_use < 0.05 && nPairs >= minN
                    if p_use < 0.001
                        stars = '***';
                    elseif p_use < 0.01
                        stars = '**';
                    else
                        stars = '*';
                    end

                    y_star = max([base; amb]) + 5;
                    line([x(si)-0.2, x(si)+0.2], [y_star, y_star], 'Color','k','LineWidth',1.0);
                    text(x(si), y_star+1.5, stars, ...
                         'HorizontalAlignment','center', ...
                         'VerticalAlignment','bottom', ...
                         'FontSize',12, ...
                         'FontWeight','bold');
                end
            end
        end
    end

    if isempty(y_all)
        ylim([0 100]);
    else
        ymax = max(y_all);
        ylim([0, ymax+15]);
    end

    xticks(x);
    xticklabels(state_order);
    ylabel('% time in state');
    set(gca,'Box','off','FontSize',11);
end

xlabel('State');

fig_file = fullfile(out_dir, 'state_percent_6h_baseline_vs_drugs_APPvsWT.png');
saveas(gcf, fig_file);

OUT = struct();
OUT.success      = true;
OUT.fig_file     = fig_file;
OUT.within_stats = within_tbl;
OUT.delta_stats  = delta_tbl;
OUT.state_order  = state_order;

% Save stats to CSVs
if ~isempty(within_tbl)
    writetable(within_tbl, fullfile(out_dir, 'state_percent_within_cond_stats.csv'));
end
if ~isempty(delta_tbl)
    writetable(delta_tbl, fullfile(out_dir, 'state_percent_delta_genotype_stats.csv'));
end

fprintf('✅ 6h baseline vs drugs state-percent figure saved to: %s\n', fig_file);
end
