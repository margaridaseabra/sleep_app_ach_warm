function OUT = plot_bouts_per_hour_3h_baseline_vs_ambtemp_APPvsWT(rows_perhr_ambtemp, out_dir, varargin)
% plot_bouts_per_hour_3h_baseline_vs_ambtemp_APPvsWT
% -------------------------------------------------------------------------
% Bar plot comparing *mean bouts/hour over the 3 h windows* in:
%   - BASELINE vs AMBTEMP
%   - separately for WT vs APP
%
% Uses PER-HOUR style table from run_group_sleep_architecture_ambtemp:
%   hour_idx, ..., bouts_per_h, state, condition, mouse, genotype
%
% For each state (default: WK, MA, NREM, REM):
%   1) Compute per-mouse mean bouts/hour over 3h baseline and 3h ambtemp.
%   2) Make bars:
%        - WT baseline   (solid grey)
%        - WT ambtemp    (grey outline)
%        - APP baseline  (solid blue)
%        - APP ambtemp   (blue outline)
%
% Stats per state:
%   - Within-genotype paired t-test (baseline vs ambtemp)  -> stars.
%   - Δ = ambtemp - baseline per mouse, WT vs APP (t-test) -> "Δ p=..."
%   - Mixed model: Condition (within) x Genotype (between) -> "C p=..., G×C p=..."
%
% OUTPUT:
%   OUT.fig_file        : PNG path
%   OUT.within_stats    : table with paired t-tests
%   OUT.delta_stats     : table with Δ t-tests (WT vs APP)
%   OUT.anova           : table with Condition/Genotype/Interaction
%   OUT.states          : states used
%
% Optional name–value:
%   'states'               : which states to include (default ["WK","MA","NREM","REM"])
%   'minNperGroupForStats' : minimum #mice per group for stats (default 3)
%   'useFDR'               : BH–FDR for within-genotype tests (default true)
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'states',["WK","MA","NREM","REM"],@(x)isstring(x)||iscellstr(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
states_to_use = string(p.Results.states(:)).';
minN  = p.Results.minNperGroupForStats;
useFDR = p.Results.useFDR;

T = rows_perhr_ambtemp;

% ---------- 1) Keep only baseline + ambtemp ----------
T.condition = lower(strtrim(string(T.condition)));
keep_cond = T.condition=="baseline" | T.condition=="ambtemp";
T = T(keep_cond, :);

if isempty(T)
    warning('No baseline + ambtemp rows found in rows_perhr_ambtemp.');
    OUT = struct('success',false,'msg','no baseline/ambtemp data');
    return;
end

T.mouse    = string(T.mouse);
T.genotype = string(T.genotype);
T.state    = string(T.state);

% ---------- 2) Per-mouse mean bouts/hour over 3h ----------
% We aggregate across hours for each mouse × genotype × condition × state
G = groupsummary(T, {'mouse','genotype','condition','state'}, 'mean', 'bouts_per_h');
G.Properties.VariableNames{end} = 'mean_bouts_per_h';

% keep only desired states
avail_states = unique(G.state);
states = states_to_use(ismember(states_to_use, avail_states));
nStates = numel(states);
if nStates==0
    warning('None of the requested states found in data.');
    OUT = struct('success',false,'msg','no requested states');
    return;
end

genotypes = ["WT","APP"];
nGen = numel(genotypes);

% containers
within_results = {};
delta_results  = {};
anova_rows     = {};

p_within = nan(nGen, nStates);   % for FDR
meanMat  = nan(4, nStates);      % WT_b, WT_a, APP_b, APP_a
semMat   = nan(4, nStates);

F_genotype    = nan(1, nStates);
p_genotype    = nan(1, nStates);
F_condition   = nan(1, nStates);
p_condition   = nan(1, nStates);
F_interaction = nan(1, nStates);
p_interaction = nan(1, nStates);

% ---------- 3) Loop over states ----------
for si = 1:nStates
    st = states(si);
    Gst = G(G.state==st, :);

    % ----- (A) Within-genotype paired tests + bars -----
    for gi = 1:nGen
        gtype = genotypes(gi);
        Gg = Gst(Gst.genotype==gtype, :);
        if isempty(Gg), continue; end

        Twide_g = unstack(Gg(:,{'mouse','condition','mean_bouts_per_h'}), ...
                          'mean_bouts_per_h','condition');

        if ~all(ismember({'baseline','ambtemp'}, Twide_g.Properties.VariableNames))
            continue;
        end

        base = Twide_g.baseline;
        amb  = Twide_g.ambtemp;
        ok   = ~isnan(base) & ~isnan(amb);
        base = base(ok);
        amb  = amb(ok);
        nPairs = numel(base);
        if nPairs < minN, continue; end

        [~, p_pair] = ttest(base, amb);

        p_within(gi,si) = p_pair;
        mean_base = mean(base,'omitnan');
        mean_amb  = mean(amb,'omitnan');
        sem_base  = std(base,'omitnan')/sqrt(nPairs);
        sem_amb   = std(amb,'omitnan')/sqrt(nPairs);
        delta     = amb - base;

        within_results = [within_results; ...
            {char(st), char(gtype), nPairs, ...
             mean_base, mean_amb, mean(delta,'omitnan'), p_pair}];

        if gtype=="WT"
            meanMat(1,si) = mean_base;
            meanMat(2,si) = mean_amb;
            semMat(1,si)  = sem_base;
            semMat(2,si)  = sem_amb;
        else
            meanMat(3,si) = mean_base;
            meanMat(4,si)  = mean_amb;
            semMat(3,si)   = sem_base;
            semMat(4,si)   = sem_amb;
        end
    end

    % ----- (B) Δ WT vs APP (difference of differences) -----
    Twide_all = unstack(Gst(:,{'mouse','genotype','condition','mean_bouts_per_h'}), ...
                        'mean_bouts_per_h','condition');
    if ~all(ismember({'baseline','ambtemp'}, Twide_all.Properties.VariableNames))
        delta_results = [delta_results; {char(st), NaN, NaN, NaN, NaN}];
        anova_rows    = [anova_rows;    {char(st), 0, 0, NaN, NaN, NaN, NaN, NaN, NaN}];
        continue;
    end

    ok_all = ~isnan(Twide_all.baseline) & ~isnan(Twide_all.ambtemp);
    Twide_all = Twide_all(ok_all, :);

    nWT_all  = sum(Twide_all.genotype=="WT");
    nAPP_all = sum(Twide_all.genotype=="APP");

    Twide_all.delta = Twide_all.ambtemp - Twide_all.baseline;
    delta_WT  = Twide_all.delta(Twide_all.genotype=="WT");
    delta_APP = Twide_all.delta(Twide_all.genotype=="APP");
    delta_WT  = delta_WT(~isnan(delta_WT));
    delta_APP = delta_APP(~isnan(delta_APP));

    nWTd  = numel(delta_WT);
    nAPPd = numel(delta_APP);

    if nWTd>=minN && nAPPd>=minN
        [~, p_delta] = ttest2(delta_WT, delta_APP, 'Vartype','unequal');
        mWTd  = mean(delta_WT,'omitnan');
        mAPPd = mean(delta_APP,'omitnan');
        sWTd  = std(delta_WT,'omitnan');
        sAPPd = std(delta_APP,'omitnan');
        n1=nWTd; n2=nAPPd;
        sp = sqrt(((n1-1)*sWTd^2 + (n2-1)*sAPPd^2)/max(1,(n1+n2-2)));
        d  = (mAPPd - mWTd)/sp;
    else
        p_delta = NaN; mWTd = NaN; mAPPd = NaN; d = NaN;
    end

    delta_results = [delta_results; ...
        {char(st), mWTd, mAPPd, p_delta, d}];

    % ----- (C) Mixed-effects model: Condition x Genotype -----
    F_gen = NaN; p_gen = NaN;
    F_cond = NaN; p_cond = NaN;
    F_int  = NaN; p_int  = NaN;

    if height(Twide_all) >= 2 && nWT_all>0 && nAPP_all>0
        try
            L = stack(Twide_all, {'baseline','ambtemp'}, ...
                      'NewDataVariableName','bph', ...
                      'IndexVariableName','Condition');
            L.Condition = categorical(L.Condition, {'baseline','ambtemp'});
            L.genotype  = categorical(L.genotype);
            L.mouse     = categorical(L.mouse);

            lme = fitlme(L, 'bph ~ Condition * genotype + (1|mouse)', ...
                         'DummyVarCoding','effects');

            a = anova(lme);

            rowC = find(contains(lower(string(a.Term)),'condition') & ...
                        ~contains(lower(string(a.Term)),':'), 1);
            rowG = find(strcmpi(string(a.Term),'genotype'), 1);
            rowI = find(contains(lower(string(a.Term)),'condition:') | ...
                        (contains(lower(string(a.Term)),'condition') & ...
                         contains(lower(string(a.Term)),'genotype')), 1);

            if ~isempty(rowC)
                F_cond = a.FStat(rowC);
                p_cond = a.pValue(rowC);
            end
            if ~isempty(rowG)
                F_gen  = a.FStat(rowG);
                p_gen  = a.pValue(rowG);
            end
            if ~isempty(rowI)
                F_int  = a.FStat(rowI);
                p_int  = a.pValue(rowI);
            end
        catch ME
            warning('Mixed-effects ANOVA for state %s failed: %s', st, ME.message);
        end
    end

    F_genotype(si)    = F_gen;
    p_genotype(si)    = p_gen;
    F_condition(si)   = F_cond;
    p_condition(si)   = p_cond;
    F_interaction(si) = F_int;
    p_interaction(si) = p_int;

    anova_rows = [anova_rows; ...
        {char(st), nWT_all, nAPP_all, F_gen, p_gen, F_cond, p_cond, F_int, p_int}];
end

% ---------- 4) Stats tables ----------
if isempty(within_results)
    within_tbl = table();
else
    within_tbl = cell2table(within_results, ...
        'VariableNames', {'State','Genotype','n', ...
                          'Mean_baseline','Mean_ambtemp', ...
                          'Mean_delta','p_paired'});
end

if isempty(delta_results)
    delta_tbl = table();
else
    delta_tbl = cell2table(delta_results, ...
        'VariableNames', {'State','Mean_delta_WT','Mean_delta_APP', ...
                          'p_delta_WTvsAPP','Cohen_d_delta'});
end

if isempty(anova_rows)
    anova_tbl = table();
else
    anova_tbl = cell2table(anova_rows, ...
        'VariableNames', {'State','nWT','nAPP', ...
                          'F_genotype','p_genotype', ...
                          'F_condition','p_condition', ...
                          'F_interaction','p_interaction'});
end

% ---------- 5) FDR on within-genotype p's ----------
p_within_fdr = nan(size(p_within));
if useFDR
    for gi = 1:nGen
        pv = p_within(gi,:);
        valid = ~isnan(pv);
        if ~any(valid), continue; end
        pvals = pv(valid);
        [sorted_p, idx] = sort(pvals(:));
        m = numel(sorted_p);
        adj = sorted_p .* (m ./ (1:m)');
        for k = m-1:-1:1
            adj(k) = min(adj(k), adj(k+1));
        end
        adj(adj>1) = 1;
        p_adj = nan(size(pvals));
        p_adj(idx) = adj;
        pv_fdr = nan(size(pv));
        pv_fdr(valid) = p_adj;
        p_within_fdr(gi,:) = pv_fdr;
    end
end

if ~isempty(within_tbl)
    p_fdr_vec = nan(height(within_tbl),1);
    for r = 1:height(within_tbl)
        st = string(within_tbl.State(r));
        gt = string(within_tbl.Genotype(r));
        gi = find(genotypes==gt,1);
        si = find(states==st,1);
        if ~isempty(gi) && ~isempty(si)
            p_fdr_vec(r) = p_within_fdr(gi,si);
        end
    end
    within_tbl.p_paired_FDR = p_fdr_vec;
end

% ---------- 6) Plot bar figure ----------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w','Position',[200 200 900 450]); hold on;
x = 1:nStates;
barWidth = 0.16;
offsets  = [-0.24, -0.08, +0.08, +0.24];

for si = 1:nStates
    % WT baseline
    if ~isnan(meanMat(1,si))
        xb = x(si) + offsets(1);
        bar(xb, meanMat(1,si), barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
        errorbar(xb, meanMat(1,si), semMat(1,si), 'k','LineStyle','none');
    end
    % WT ambtemp
    if ~isnan(meanMat(2,si))
        xa = x(si) + offsets(2);
        bar(xa, meanMat(2,si), barWidth, 'FaceColor','none','EdgeColor',COL_WT,'LineWidth',1.2);
        errorbar(xa, meanMat(2,si), semMat(2,si), 'k','LineStyle','none');
    end
    % APP baseline
    if ~isnan(meanMat(3,si))
        xb = x(si) + offsets(3);
        bar(xb, meanMat(3,si), barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(xb, meanMat(3,si), semMat(3,si), 'k','LineStyle','none');
    end
    % APP ambtemp
    if ~isnan(meanMat(4,si))
        xa = x(si) + offsets(4);
        bar(xa, meanMat(4,si), barWidth, 'FaceColor','none','EdgeColor',COL_APP,'LineWidth',1.2);
        errorbar(xa, meanMat(4,si), semMat(4,si), 'k','LineStyle','none');
    end
end

all_means = meanMat(~isnan(meanMat));
if isempty(all_means)
    ymax = 10;
else
    ymax = max(all_means) + 3;
end
ylim([0 ymax]);

% ---- stars for paired baseline vs ambtemp ----
for si = 1:nStates
    st = states(si);
    for gi = 1:nGen
        gtype = genotypes(gi);
        if isempty(within_tbl), continue; end
        row = strcmp(string(within_tbl.State),string(st)) & ...
              strcmp(string(within_tbl.Genotype),string(gtype));
        if ~any(row), continue; end

        p_here  = within_tbl.p_paired(row);
        p_hereF = within_tbl.p_paired_FDR(row);
        p_use = p_here;
        if useFDR && ~isnan(p_hereF), p_use = p_hereF; end
        if isnan(p_use) || p_use >= 0.05, continue; end

        if p_use < 0.001, stars = '***';
        elseif p_use < 0.01, stars = '**';
        else, stars = '*'; end

        if gtype=="WT"
            x1 = x(si)+offsets(1);
            x2 = x(si)+offsets(2);
            y_star = ymax - 1.5;
        else
            x1 = x(si)+offsets(3);
            x2 = x(si)+offsets(4);
            y_star = ymax - 3.5;
        end

        line([x1 x2],[y_star y_star],'Color','k','LineWidth',1);
        text(mean([x1 x2]), y_star+0.3, stars, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
    end
end

% ---- text with Δ, Condition, Interaction p-values ----
for si = 1:nStates
    st = states(si);

    p_delta = NaN;
    if ~isempty(delta_tbl)
        rowD = strcmp(string(delta_tbl.State),string(st));
        if any(rowD)
            p_delta = delta_tbl.p_delta_WTvsAPP(find(rowD,1));
        end
    end
    pC = p_condition(si);
    pI = p_interaction(si);

    if isnan(p_delta) && isnan(pC) && isnan(pI), continue; end

    txt = sprintf('Δ p=%.3f | C p=%.3f | G×C p=%.3f', p_delta, pC, pI);
    text(x(si), ymax-0.5, txt, ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end

xticks(x);
xticklabels(states);
xlabel('State');
ylabel('Mean bouts per hour (3h window)');
title('3 h baseline vs 3 h ambtemp – bouts/hour (WT vs APP)');
set(gca,'Box','off','FontSize',11);

hb1 = bar(nan,nan,barWidth,'FaceColor',COL_WT,'EdgeColor','none');
hb2 = bar(nan,nan,barWidth,'FaceColor',COL_APP,'EdgeColor','none');
ha1 = bar(nan,nan,barWidth,'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1.2);
legend([hb1 hb2 ha1], {'WT','APP','ambtemp (outline vs baseline solid)'}, ...
    'Location','northoutside','Orientation','horizontal');

fig_file = fullfile(out_dir,'bouts_per_hour_3h_baseline_vs_ambtemp_APPvsWT_bar.png');
saveas(gcf, fig_file);

% ---------- 7) Save stats ----------
if ~isempty(within_tbl)
    writetable(within_tbl, fullfile(out_dir,'bouts_within_cond_stats_bar.csv'));
end
if ~isempty(delta_tbl)
    writetable(delta_tbl, fullfile(out_dir,'bouts_delta_genotype_stats_bar.csv'));
end
if ~isempty(anova_tbl)
    writetable(anova_tbl, fullfile(out_dir,'bouts_anova_stats_bar.csv'));
end

OUT = struct();
OUT.success      = true;
OUT.fig_file     = fig_file;
OUT.within_stats = within_tbl;
OUT.delta_stats  = delta_tbl;
OUT.anova        = anova_tbl;
OUT.states       = states;

fprintf('✅ Bouts/hour bar plot + stats (baseline vs ambtemp, WT vs APP) saved in %s\n', out_dir);
end
