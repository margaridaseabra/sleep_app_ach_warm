function OUT = plot_state_percent_bar_6h_baseline_vs_drugs_APPvsWT(rows_overall_drugs, out_dir, varargin)
% plot_state_percent_bar_6h_baseline_vs_drugs_APPvsWT
% -------------------------------------------------------------------------
% Bar plot comparing 6 h BASELINE vs 6 h drugs, separately in WT vs APP.
%
% OVERALL-style table columns required:
%   state, total_dur_s, condition, mouse, genotype
%
% States: WK, MA, NREM, REM, and composite SLEEP = NREM+REM.
%
% Bars per state (left -> right):
%   1) WT baseline   (solid grey)
%   2) WT drugs    (grey outline)
%   3) APP baseline  (solid blue)
%   4) APP drugs   (blue outline)
%
% Stats:
%   - Within-genotype paired t-test baseline vs drugs   => stars.
%   - Δ = drugs - baseline per mouse, WT vs APP (t-test)=> "Δ p=..."
%   - 2x2 mixed ANOVA (Condition x Genotype)              => "C p=..., G×C p=..."
%
% OUTPUT:
%   OUT.fig_file        : PNG path
%   OUT.within_stats    : table with paired t-tests
%   OUT.delta_stats     : table with Δ t-tests (WT vs APP)
%   OUT.anova           : table with 2x2 ANOVA per state
%   OUT.state_order     : order of states on x-axis
%
% Optional name–value:
%   'minNperGroupForStats' (default 3)
%   'useFDR'               (default true) -> BH-FDR on within-genotype tests
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
minN  = p.Results.minNperGroupForStats;
useFDR = p.Results.useFDR;

T = rows_overall_drugs;

% ---------- 1) Keep only baseline + drugs ----------
T.condition = lower(strtrim(string(T.condition)));
keep_cond = T.condition == "baseline" | T.condition == "drugs";
T = T(keep_cond, :);

if isempty(T)
    warning('No baseline + drugs rows found. Nothing to plot.');
    OUT = struct('success',false,'msg','no baseline/drugs data');
    return;
end

T.mouse    = string(T.mouse);
T.genotype = string(T.genotype);
T.state    = string(T.state);

% ---------- 2) Build % time per state (incl. SLEEP) ----------
base_desired = ["WK","MA","NREM","REM"];
avail_states = unique(T.state);
base_states  = base_desired(ismember(base_desired, avail_states));

if isempty(base_states)
    warning('No WK/MA/NREM/REM states found.');
    OUT = struct('success',false,'msg','no base states');
    return;
end

% per-mouse, per-genotype, per-condition, per-state
Gbase = groupsummary(T, {'mouse','genotype','condition','state'}, 'sum', 'total_dur_s');
Gbase.Properties.VariableNames{end} = 'dur_s';
Gbase = Gbase(ismember(Gbase.state, base_states), :);

% total duration per mouse & condition (across WK/MA/NREM/REM)
Gtot = groupsummary(Gbase, {'mouse','genotype','condition'}, 'sum', 'dur_s');
Gtot.Properties.VariableNames{end} = 'total_dur_allstates';

% SLEEP = NREM + REM
sleep_components = ["NREM","REM"];
Gsleep_src = Gbase(ismember(Gbase.state, sleep_components), :);
if ~isempty(Gsleep_src)
    Gsleep = groupsummary(Gsleep_src, {'mouse','genotype','condition'}, 'sum', 'dur_s');
    Gsleep.Properties.VariableNames{end} = 'dur_s';
    Gsleep.state = repmat("SLEEP", height(Gsleep), 1);
else
    Gsleep = Gbase([],:);
end

% combine and compute %
G = [Gbase; Gsleep];
G = innerjoin(G, Gtot, 'Keys', {'mouse','genotype','condition'});
G.pct = 100 * G.dur_s ./ G.total_dur_allstates;

desired_order = ["WK","MA","NREM","REM","SLEEP"];
avail_states2 = unique(G.state);
state_order   = desired_order(ismember(desired_order, avail_states2));
nStates       = numel(state_order);

if nStates == 0
    warning('No usable states after SLEEP build.');
    OUT = struct('success',false,'msg','no states');
    return;
end

genotypes = ["WT","APP"];
nGen      = numel(genotypes);

% ---------- 3) Containers for stats & means ----------
within_results = {};   % cell rows
delta_results  = {};
anova_rows     = {};

% within-genotype paired p-values (for FDR)
p_within = nan(nGen, nStates);

% bar means / SEMs: rows = [WT_base, WT_amb, APP_base, APP_amb]
meanMat = nan(4, nStates);
semMat  = nan(4, nStates);

% ANOVA vectors
F_genotype    = nan(1, nStates);
p_genotype    = nan(1, nStates);
F_condition   = nan(1, nStates);
p_condition   = nan(1, nStates);
F_interaction = nan(1, nStates);
p_interaction = nan(1, nStates);

% ---------- 4) Loop over states ----------
for si = 1:nStates
    st = state_order(si);
    Gst = G(G.state == st, :);

    % ----- (A) Within-genotype paired tests + means/SEMs -----
    for gi = 1:nGen
        gtype = genotypes(gi);
        Gg = Gst(Gst.genotype == gtype, :);
        if isempty(Gg), continue; end

        Twide_g = unstack(Gg(:,{'mouse','condition','pct'}), 'pct','condition');
        if ~all(ismember({'baseline','drugs'}, Twide_g.Properties.VariableNames))
            continue;
        end

        base = Twide_g.baseline;
        amb  = Twide_g.drugs;
        ok   = ~isnan(base) & ~isnan(amb);
        base = base(ok);
        amb  = amb(ok);
        nPairs = numel(base);

        if nPairs < minN
            continue;
        end

        [~, p_pair] = ttest(base, amb);
        p_within(gi,si) = p_pair;

        mean_base = mean(base,'omitnan');
        mean_amb  = mean(amb,'omitnan');
        sem_base  = std(base,'omitnan') / sqrt(nPairs);
        sem_amb   = std(amb,'omitnan') / sqrt(nPairs);
        delta     = amb - base;

        within_results = [within_results; ...
            {char(st), char(gtype), nPairs, ...
             mean_base, mean_amb, mean(delta,'omitnan'), p_pair}];

        % store for bars
        if gtype == "WT"
            meanMat(1,si) = mean_base;
            meanMat(2,si) = mean_amb;
            semMat(1,si)  = sem_base;
            semMat(2,si)  = sem_amb;
        else
            meanMat(3,si) = mean_base;
            meanMat(4,si) = mean_amb;
            semMat(3,si)  = sem_base;
            semMat(4,si)  = sem_amb;
        end
    end

    % ----- (B) Δ WT vs APP -----
    Twide_all = unstack(Gst(:,{'mouse','genotype','condition','pct'}), 'pct','condition');
    if ~all(ismember({'baseline','drugs'}, Twide_all.Properties.VariableNames))
        % still create empty delta & anova rows
        delta_results = [delta_results; {char(st), NaN, NaN, NaN, NaN}];
        anova_rows    = [anova_rows;    {char(st), 0, 0, NaN, NaN, NaN, NaN, NaN, NaN}];
        continue;
    end

    ok_all = ~isnan(Twide_all.baseline) & ~isnan(Twide_all.drugs);
    Twide_all = Twide_all(ok_all, :);

    nWT_all  = sum(Twide_all.genotype=="WT");
    nAPP_all = sum(Twide_all.genotype=="APP");

    Twide_all.delta = Twide_all.drugs - Twide_all.baseline;
    delta_WT  = Twide_all.delta(Twide_all.genotype=="WT");
    delta_APP = Twide_all.delta(Twide_all.genotype=="APP");
    delta_WT  = delta_WT(~isnan(delta_WT));
    delta_APP = delta_APP(~isnan(delta_APP));

    nWTd  = numel(delta_WT);
    nAPPd = numel(delta_APP);

    if nWTd >= minN && nAPPd >= minN
        [~, p_delta] = ttest2(delta_WT, delta_APP, 'Vartype','unequal');
        mWTd  = mean(delta_WT,'omitnan');
        mAPPd = mean(delta_APP,'omitnan');
        sWTd  = std(delta_WT,'omitnan');
        sAPPd = std(delta_APP,'omitnan');
        n1 = nWTd; n2 = nAPPd;
        sp = sqrt(((n1-1)*sWTd^2 + (n2-1)*sAPPd^2) / max(1,(n1+n2-2)));
        d  = (mAPPd - mWTd) / sp;
    else
        p_delta = NaN; mWTd = NaN; mAPPd = NaN; d = NaN;
    end

    delta_results = [delta_results; ...
        {char(st), mWTd, mAPPd, p_delta, d}];

    % ----- (C) 2x2 mixed ANOVA via mixed-effects (Condition x Genotype) -----
    F_gen = NaN; p_gen = NaN;
    F_cond = NaN; p_cond = NaN;
    F_int  = NaN; p_int  = NaN;

    if height(Twide_all) >= 2 && nWT_all > 0 && nAPP_all > 0
        try
            % Wide -> long for mixed model
            L = stack(Twide_all, {'baseline','drugs'}, ...
                    'NewDataVariableName','pct', 'IndexVariableName','Condition');
            L.Condition = categorical(L.Condition, {'baseline','drugs'}); % enforce order
            L.genotype  = categorical(L.genotype);                           % WT/APP
            L.mouse     = categorical(L.mouse);

            % Random intercept per mouse; fixed Condition, Genotype, and interaction
            lme = fitlme(L, 'pct ~ Condition * genotype + (1|mouse)', ...
                        'DummyVarCoding','effects');

            % Type III-style tests
            a = anova(lme);  % returns rows: (Intercept), Condition, genotype, Condition:genotype (labels may vary but are stable textually)

            % pull rows by name, case-insensitive
            rowC = find(contains(lower(string(a.Term)),'condition') & ~contains(lower(string(a.Term)),':'), 1);
            rowG = find(strcmpi(string(a.Term),'genotype'), 1);
            rowI = find(contains(lower(string(a.Term)),'condition:') | ...
                    (contains(lower(string(a.Term)),'condition') & contains(lower(string(a.Term)),'genotype')), 1);

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
% ---------- 5) Build stats tables ----------
if isempty(within_results)
    within_tbl = table();
else
    within_tbl = cell2table(within_results, ...
        'VariableNames', {'State','Genotype','n', ...
                          'Mean_baseline','Mean_drugs', ...
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

% ---------- 6) FDR for within-genotype paired tests ----------
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

% attach FDR column
if ~isempty(within_tbl)
    p_fdr_vec = nan(height(within_tbl),1);
    for r = 1:height(within_tbl)
        st = string(within_tbl.State(r));
        gt = string(within_tbl.Genotype(r));
        gi = find(genotypes == gt, 1);
        si = find(state_order == st, 1);
        if ~isempty(gi) && ~isempty(si)
            p_fdr_vec(r) = p_within_fdr(gi,si);
        end
    end
    within_tbl.p_paired_FDR = p_fdr_vec;
end

% ---------- 7) Plot bar figure ----------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w','Position',[200 200 900 450]);
hold on;

x = 1:nStates;
barWidth = 0.16;
offsets  = [-0.24, -0.08, +0.08, +0.24];   % WT_b, WT_a, APP_b, APP_a

for si = 1:nStates
    % WT baseline
    if ~isnan(meanMat(1,si))
        xb = x(si) + offsets(1);
        bar(xb, meanMat(1,si), barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
        errorbar(xb, meanMat(1,si), semMat(1,si), 'k', 'LineStyle','none');
    end
    % WT drugs
    if ~isnan(meanMat(2,si))
        xa = x(si) + offsets(2);
        bar(xa, meanMat(2,si), barWidth, 'FaceColor','none','EdgeColor',COL_WT,'LineWidth',1.2);
        errorbar(xa, meanMat(2,si), semMat(2,si), 'k', 'LineStyle','none');
    end
    % APP baseline
    if ~isnan(meanMat(3,si))
        xb = x(si) + offsets(3);
        bar(xb, meanMat(3,si), barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(xb, meanMat(3,si), semMat(3,si), 'k', 'LineStyle','none');
    end
    % APP drugs
    if ~isnan(meanMat(4,si))
        xa = x(si) + offsets(4);
        bar(xa, meanMat(4,si), barWidth, 'FaceColor','none','EdgeColor',COL_APP,'LineWidth',1.2);
        errorbar(xa, meanMat(4,si), semMat(4,si), 'k', 'LineStyle','none');
    end
end

all_means = meanMat(~isnan(meanMat));
if isempty(all_means)
    ymax = 100;
else
    ymax = max(all_means) + 15;
end
ylim([0 ymax]);

% --- stars for paired baseline vs drugs within genotype ---
for si = 1:nStates
    st = state_order(si);
    for gi = 1:nGen
        gtype = genotypes(gi);
        if isempty(within_tbl), continue; end
        row = within_tbl.State == st & within_tbl.Genotype == gtype;
        if ~any(row), continue; end
        p_here  = within_tbl.p_paired(row);
        p_hereF = within_tbl.p_paired_FDR(row);
        p_use = p_here;
        if useFDR && ~isnan(p_hereF), p_use = p_hereF; end
        if isnan(p_use) || p_use >= 0.05, continue; end

        if p_use < 0.001, stars = '***';
        elseif p_use < 0.01, stars = '**';
        else, stars = '*'; end

        if gtype == "WT"
            x1 = x(si) + offsets(1);
            x2 = x(si) + offsets(2);
            y_star = ymax - 5;
        else
            x1 = x(si) + offsets(3);
            x2 = x(si) + offsets(4);
            y_star = ymax - 10;
        end

        line([x1 x2], [y_star y_star], 'Color','k','LineWidth',1.0);
        text(mean([x1 x2]), y_star+1.5, stars, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
    end
end

% --- text with Δ, Condition, Interaction p-values ---
% --- text with Δ, Condition, Interaction p-values + ANOVA stars ---
for si = 1:nStates
    st = state_order(si);

    % ----- grab Δ p-value for this state -----
    p_delta = NaN;
    if ~isempty(delta_tbl)
        rowD = strcmp(string(delta_tbl.State), string(st));
        if any(rowD)
            p_delta = delta_tbl.p_delta_WTvsAPP(find(rowD,1));
        end
    end

    % ----- ANOVA p-values -----
    pC = p_condition(si);      % condition main effect
    pI = p_interaction(si);    % genotype × condition interaction

    % ----- text line with numeric p-values (as before) -----
    if ~(isnan(p_delta) && isnan(pC) && isnan(pI))
        txt = sprintf('Δ p=%.3f | C p=%.3f | G×C p=%.3f', p_delta, pC, pI);
        text(x(si), ymax-12, txt, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontSize',9);
    end

    % ----- decide ANOVA stars (priority: interaction, then condition) -----
    starLabel = '';
    % interaction first
    if ~isnan(pI) && pI < 0.05
        if     pI < 0.001, starsI = '***';
        elseif pI < 0.01,  starsI = '**';
        else               starsI = '*';
        end
        starLabel = sprintf('G×C %s', starsI);
    % if no significant interaction, show condition main effect
    elseif ~isnan(pC) && pC < 0.05
        if     pC < 0.001, starsC = '***';
        elseif pC < 0.01,  starsC = '**';
        else               starsC = '*';
        end
        starLabel = sprintf('C %s', starsC);
    end

    if ~isempty(starLabel)
        % line above whole state cluster
        y_line = ymax - 3;  % a bit below top of axis
        line([x(si)-0.35, x(si)+0.35], [y_line y_line], 'Color','k', 'LineWidth',1.0);
        text(x(si), y_line + 1.5, starLabel, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',10, ...
            'FontWeight','bold');
    end
end

xticks(x);
xticklabels(state_order);
xlabel('State');
ylabel('% time in state');
title('6 h baseline vs 6 h drugs – WT vs APP');

% legend
hb1 = bar(nan, nan, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
hb2 = bar(nan, nan, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
ha1 = bar(nan, nan, barWidth, 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1.2);
legend([hb1 hb2 ha1], {'WT','APP','drugs (outline vs baseline solid)'}, ...
    'Location','northoutside','Orientation','horizontal');

set(gca,'Box','off','FontSize',11);

fig_file = fullfile(out_dir, 'state_percent_6h_baseline_vs_drugs_APPvsWT_bar.png');
saveas(gcf, fig_file);

% ---------- 8) Save stats ----------
if ~isempty(within_tbl)
    writetable(within_tbl, fullfile(out_dir, 'state_percent_within_cond_stats_bar.csv'));
end
if ~isempty(delta_tbl)
    writetable(delta_tbl, fullfile(out_dir, 'state_percent_delta_genotype_stats_bar.csv'));
end
if ~isempty(anova_tbl)
    writetable(anova_tbl, fullfile(out_dir, 'state_percent_anova_stats_bar.csv'));
end

OUT = struct();
OUT.success      = true;
OUT.fig_file     = fig_file;
OUT.within_stats = within_tbl;
OUT.delta_stats  = delta_tbl;
OUT.anova        = anova_tbl;
OUT.state_order  = state_order;

fprintf('✅ Bar plot + stats (incl. 2-way ANOVA) saved in %s\n', out_dir);

end % function

% ---------- local helpers ----------
function labels = local_pickLabelColumn(tbl, candidates)
% pick the first existing label column from a list
labels = strings(height(tbl),1);
for k = 1:numel(candidates)
    if ismember(candidates(k), tbl.Properties.VariableNames)
        labels = string(tbl.(candidates(k)));
        return;
    end
end
end

function labels = local_fillEmptyWithRowNames(tbl, labels)
% fill empty/missing labels from row names if available
emptyW = (labels=="" | ismissing(labels));
if any(emptyW)
    try
        rn = string(tbl.Properties.RowNames);
        if numel(rn) == height(tbl)
            labels(emptyW) = rn(emptyW);
        end
    catch
        % noop if row names not present
    end
end
end
