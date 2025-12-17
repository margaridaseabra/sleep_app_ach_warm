function OUT = plot_bout_rate_3h_baseline_vs_ambtemp_APPvsWT(rows_overall_ambtemp, out_dir, varargin)
% plot_bout_rate_3h_baseline_vs_ambtemp_APPvsWT
% -------------------------------------------------------------------------
% Uses the CROPPED 3 h baseline + 3 h ambtemp OVERALL table to compute:
%
%   bouts_per_h = n_bouts / (total_dur_s/3600)
%
% for states WK, MA, NREM, REM (+ composite SLEEP = NREM+REM).
%
% Makes ONE bar figure:
%   For each state, 4 bars:
%       WT baseline   (solid grey)
%       WT ambtemp    (grey outline)
%       APP baseline  (solid blue)
%       APP ambtemp   (blue outline)
%
% Stats:
%   - Within-genotype paired t-test (baseline vs ambtemp), WT and APP
%       • BH-FDR across states
%       • stars drawn between the 2 bars (per genotype)
%   - Δ = ambtemp - baseline per mouse, WT vs APP (unpaired t-test)
%   - 2x2 mixed ANOVA (Condition × Genotype) via mixed-effects model
%       • prints C, G, G×C p-values above each state
%
% INPUT
%   rows_overall_ambtemp : table from run_group_sleep_architecture_ambtemp_overall
%   out_dir              : folder for PNG + CSVs
%
% NAME–VALUE OPTIONS
%   'states'              : which states to include (default ["WK","MA","NREM","REM","SLEEP"])
%   'minNperGroupForStats': minimum n for stats (default 3)
%   'useFDR'              : use FDR for within-genotype stars (default true)
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'states',["WK","MA","NREM","REM","SLEEP"]);
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});

states_desired = string(p.Results.states(:)).';
minN  = p.Results.minNperGroupForStats;
useFDR = p.Results.useFDR;

T = rows_overall_ambtemp;
T.state     = string(T.state);
T.condition = lower(strtrim(string(T.condition)));
T.mouse     = string(T.mouse);
T.genotype  = string(T.genotype);

% keep only baseline + ambtemp
mask_cond = T.condition=="baseline" | T.condition=="ambtemp";
T = T(mask_cond,:);
if isempty(T)
    warning('No baseline + ambtemp rows in rows_overall_ambtemp.');
    OUT = struct('success',false,'msg','no data');
    return;
end

% --- aggregate just in case & compute bouts/hour ---
G = groupsummary(T, {'mouse','genotype','condition','state'}, 'sum', {'n_bouts','total_dur_s'});
G.Properties.VariableNames(end-1:end) = {'n_bouts','total_dur_s'};
G.bouts_per_h = G.n_bouts ./ max(G.total_dur_s/3600, eps);

% build composite SLEEP (NREM+REM)
sleep_components = ["NREM","REM"];
Gsleep_src = G(ismember(G.state,sleep_components),:);
if ~isempty(Gsleep_src)
    Gsleep = groupsummary(Gsleep_src, {'mouse','genotype','condition'}, 'sum', {'n_bouts','total_dur_s'});
    Gsleep.Properties.VariableNames(end-1:end) = {'n_bouts','total_dur_s'};
    Gsleep.state = repmat("SLEEP", height(Gsleep),1);
    Gsleep.bouts_per_h = Gsleep.n_bouts ./ max(Gsleep.total_dur_s/3600, eps);
    G = [G; Gsleep];
end

% final state order (only those present)
avail = unique(G.state);
state_order = states_desired(ismember(states_desired, avail));
nStates = numel(state_order);
if nStates==0
    warning('None of the requested states are present.');
    OUT = struct('success',false,'msg','no states');
    return;
end

genotypes = ["WT","APP"];
nGen = numel(genotypes);

% containers
meanMat = nan(4,nStates);   % WT_b, WT_a, APP_b, APP_a
semMat  = nan(4,nStates);

p_within = nan(nGen,nStates);    % paired baseline vs ambtemp per genotype
p_within_FDR = nan(nGen,nStates);
delta_tbl_rows = {};
anova_tbl_rows = {};

% ---------- loop states ----------
for si = 1:nStates
    st  = state_order(si);
    Gst = G(G.state==st,:);

    % ----- within-genotype paired tests -----
    for gi = 1:nGen
        gtype = genotypes(gi);
        Gg = Gst(Gst.genotype==gtype,:);
        if isempty(Gg), continue; end

        Twide_g = unstack(Gg(:,{'mouse','condition','bouts_per_h'}), ...
                          'bouts_per_h','condition');
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

        m_base = mean(base,'omitnan');
        m_amb  = mean(amb,'omitnan');
        se_base = std(base,'omitnan')/sqrt(nPairs);
        se_amb  = std(amb,'omitnan')/sqrt(nPairs);
        dlt     = amb - base;

        if gtype=="WT"
            meanMat(1,si) = m_base;
            meanMat(2,si) = m_amb;
            semMat(1,si)  = se_base;
            semMat(2,si)  = se_amb;
        else
            meanMat(3,si) = m_base;
            meanMat(4,si) = m_amb;
            semMat(3,si)  = se_base;
            semMat(4,si)  = se_amb;
        end

        delta_tbl_rows = [delta_tbl_rows; ...
            {char(st), char(gtype), nPairs, ...
             mean(dlt,'omitnan'), p_pair}]; %#ok<AGROW>
    end

    % ----- Δ WT vs APP (ambtemp-baseline) -----
    Twide_all = unstack(Gst(:,{'mouse','genotype','condition','bouts_per_h'}), ...
                        'bouts_per_h','condition');
    if all(ismember({'baseline','ambtemp'}, Twide_all.Properties.VariableNames))
        ok_all = ~isnan(Twide_all.baseline) & ~isnan(Twide_all.ambtemp);
        Twide_all = Twide_all(ok_all,:);
        Twide_all.delta = Twide_all.ambtemp - Twide_all.baseline;

        dWT  = Twide_all.delta(Twide_all.genotype=="WT");
        dAPP = Twide_all.delta(Twide_all.genotype=="APP");
        dWT  = dWT(~isnan(dWT)); dAPP = dAPP(~isnan(dAPP));
        nWT  = numel(dWT); nAPP = numel(dAPP);
        if nWT>=minN && nAPP>=minN
            [~, p_delta] = ttest2(dWT,dAPP,'Vartype','unequal');
            mWT  = mean(dWT,'omitnan'); mAPP = mean(dAPP,'omitnan');
            sWT  = std(dWT,'omitnan');  sAPP = std(dAPP,'omitnan');
            sp   = sqrt(((nWT-1)*sWT^2 + (nAPP-1)*sAPP^2)/max(1,(nWT+nAPP-2)));
            d    = (mAPP-mWT)/sp;
        else
            p_delta = NaN; d = NaN;
        end
    else
        nWT = 0; nAPP = 0; p_delta = NaN; d = NaN;
    end

    % ----- 2x2 mixed ANOVA (Condition × Genotype) using mixed effects -----
    Fg=NaN; pg=NaN; Fc=NaN; pc=NaN; Fi=NaN; pi=NaN;

    if exist('Twide_all','var') && ~isempty(Twide_all) && nWT>0 && nAPP>0
        L = stack(Twide_all, {'baseline','ambtemp'}, ...
                  'NewDataVariableName','bph','IndexVariableName','Condition');
        L.Condition = categorical(L.Condition, {'baseline','ambtemp'});
        L.genotype  = categorical(L.genotype);
        L.mouse     = categorical(L.mouse);

        try
            lme = fitlme(L,'bph ~ Condition * genotype + (1|mouse)', ...
                           'DummyVarCoding','effects');
            a     = anova(lme);
            terms = string(a.Term);

            % main effect of condition
            rowC = find(strcmpi(terms,'Condition'), 1);

            % main effect of genotype
            rowG = find(strcmpi(terms,'genotype'), 1);

            % interaction: any term that contains both words and a colon
            rowI = find( contains(lower(terms),'condition') & ...
                        contains(lower(terms),'genotype') & ...
                        contains(terms,':'), 1 );

            Fg=NaN; pg=NaN; Fc=NaN; pc=NaN; Fi=NaN; pi=NaN;
            if ~isempty(rowC), Fc = a.FStat(rowC); pc = a.pValue(rowC); end
            if ~isempty(rowG), Fg = a.FStat(rowG); pg = a.pValue(rowG); end
            if ~isempty(rowI), Fi = a.FStat(rowI); pi = a.pValue(rowI); end

        catch ME
            warning('ANOVA failed for state %s: %s', st, ME.message);
        end
    end

    anova_tbl_rows = [anova_tbl_rows; ...
        {char(st), nWT, nAPP, Fg, pg, Fc, pc, Fi, pi, p_delta, d}]; %#ok<AGROW>
end

% ---------- FDR for within-genotype tests ----------
for gi = 1:nGen
    pv = p_within(gi,:);
    valid = ~isnan(pv);
    if ~any(valid), continue; end
    pvals = pv(valid);
    [sp,idx] = sort(pvals(:));
    m = numel(sp);
    adj = sp .* (m./(1:m))';
    for k = m-1:-1:1
        adj(k) = min(adj(k),adj(k+1));
    end
    adj(adj>1) = 1;
    p_adj = nan(size(pvals)); p_adj(idx)=adj;
    pv_fdr = nan(size(pv));   pv_fdr(valid)=p_adj;
    p_within_FDR(gi,:) = pv_fdr;
end

% ---------- build stats tables ----------
if isempty(delta_tbl_rows)
    within_tbl = table();
else
    within_tbl = cell2table(delta_tbl_rows, ...
        'VariableNames', {'State','Genotype','n_pairs', ...
                          'Mean_delta_bph','p_paired'});
    % add FDR
    pF = nan(height(within_tbl),1);
    for r = 1:height(within_tbl)
        st = string(within_tbl.State(r));
        gt = string(within_tbl.Genotype(r));
        si = find(state_order==st,1);
        gi = find(genotypes==gt,1);
        if ~isempty(si)&&~isempty(gi)
            pF(r) = p_within_FDR(gi,si);
        end
    end
    within_tbl.p_paired_FDR = pF;
end

if isempty(anova_tbl_rows)
    anova_tbl = table();
else
    anova_tbl = cell2table(anova_tbl_rows, ...
        'VariableNames', {'State','nWT','nAPP', ...
                          'F_genotype','p_genotype', ...
                          'F_condition','p_condition', ...
                          'F_interaction','p_interaction', ...
                          'p_delta_WTvsAPP','Cohen_d_delta'});
end

% ---------- plotting ----------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w','Position',[200 200 900 450]); hold on;
x = 1:nStates;
barWidth = 0.16;
offsets  = [-0.24,-0.08,+0.08,+0.24]; % WT_b, WT_a, APP_b, APP_a

for si = 1:nStates
    % WT baseline
    if ~isnan(meanMat(1,si))
        xb = x(si)+offsets(1);
        bar(xb, meanMat(1,si), barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
        errorbar(xb, meanMat(1,si), semMat(1,si),'k','LineStyle','none');
    end
    % WT ambtemp
    if ~isnan(meanMat(2,si))
        xa = x(si)+offsets(2);
        bar(xa, meanMat(2,si), barWidth, 'FaceColor','none','EdgeColor',COL_WT,'LineWidth',1.2);
        errorbar(xa, meanMat(2,si), semMat(2,si),'k','LineStyle','none');
    end
    % APP baseline
    if ~isnan(meanMat(3,si))
        xb = x(si)+offsets(3);
        bar(xb, meanMat(3,si), barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(xb, meanMat(3,si), semMat(3,si),'k','LineStyle','none');
    end
    % APP ambtemp
    if ~isnan(meanMat(4,si))
        xa = x(si)+offsets(4);
        bar(xa, meanMat(4,si), barWidth, 'FaceColor','none','EdgeColor',COL_APP,'LineWidth',1.2);
        errorbar(xa, meanMat(4,si), semMat(4,si),'k','LineStyle','none');
    end
end

all_means = meanMat(~isnan(meanMat));
if isempty(all_means), ymax = 1; else, ymax = max(all_means)+max(all_means)*0.3; end
ylim([0 ymax]);

% stars for paired baseline vs ambtemp
for si = 1:nStates
    st = state_order(si);
    for gi = 1:nGen
        gt = genotypes(gi);
        if isempty(within_tbl), continue; end
        row = strcmp(within_tbl.State, st) & strcmp(within_tbl.Genotype, gt);
        if ~any(row), continue; end
        p_raw  = within_tbl.p_paired(row);
        p_fdr  = within_tbl.p_paired_FDR(row);
        p_use  = p_raw;
        if useFDR && ~isnan(p_fdr), p_use = p_fdr; end
        if isnan(p_use) || p_use>=0.05, continue; end

        if p_use<0.001, stars='***';
        elseif p_use<0.01, stars='**';
        else, stars='*';
        end

        if gt=="WT"
            x1 = x(si)+offsets(1);
            x2 = x(si)+offsets(2);
            y_star = ymax*0.93;
        else
            x1 = x(si)+offsets(3);
            x2 = x(si)+offsets(4);
            y_star = ymax*0.80;
        end
        line([x1 x2],[y_star y_star],'Color','k','LineWidth',1);
        text(mean([x1 x2]), y_star+0.02*ymax, stars, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
    end
end

% text with Δ and ANOVA p-values
% --- text with Δ and ANOVA p-values + ANOVA stars ---
for si = 1:nStates
    st = state_order(si);
    rowA = strcmp(anova_tbl.State, st);
    if ~any(rowA), continue; end

    pD = anova_tbl.p_delta_WTvsAPP(rowA);
    pC = anova_tbl.p_condition(rowA);
    pI = anova_tbl.p_interaction(rowA);

    % If absolutely everything is NaN, skip this state
    if isnan(pD) && isnan(pC) && isnan(pI)
        continue;
    end

    % Pretty strings (avoid showing "NaN")
    if isnan(pD), sD = 'n/a'; else, sD = sprintf('%.3f', pD); end
    if isnan(pC), sC = 'n/a'; else, sC = sprintf('%.3f', pC); end
    if isnan(pI), sI = 'n/a'; else, sI = sprintf('%.3f', pI); end

    % Numeric p text (same layout as before)
    txt = sprintf('\\Delta p=%s | C p=%s | G\\timesC p=%s', sD, sC, sI);
    text(x(si), ymax*0.70, txt, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',9);

    % ---------- ANOVA stars ----------
    % Priority: interaction (G×C), then condition main effect (C)
    starLabel = '';

    if ~isnan(pI) && pI < 0.05
        if     pI < 0.001, star = '***';
        elseif pI < 0.01,  star = '**';
        else,              star = '*';  end
        starLabel = ['G×C ' star];
    elseif ~isnan(pC) && pC < 0.05
        if     pC < 0.001, star = '***';
        elseif pC < 0.01,  star = '**';
        else,              star = '*';  end
        starLabel = ['C ' star];
    end

    if ~isempty(starLabel)
        % bracket over all 4 bars for that state
        y_line = ymax*0.66;
        line([x(si)-0.35, x(si)+0.35], [y_line, y_line], ...
             'Color','k', 'LineWidth',1.0);
        text(x(si), y_line + 0.02*ymax, starLabel, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',10, 'FontWeight','bold');
    end
end


xticks(x);
xticklabels(state_order);
xlabel('State');
ylabel('Bouts per hour (3 h window)');
title('3 h baseline vs 3 h ambtemp – bouts/hour (WT vs APP)');

hb1 = bar(nan,nan,barWidth,'FaceColor',COL_WT,'EdgeColor','none');
hb2 = bar(nan,nan,barWidth,'FaceColor',COL_APP,'EdgeColor','none');
ha1 = bar(nan,nan,barWidth,'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',1.2);
legend([hb1 hb2 ha1], {'WT','APP','ambtemp (outline vs baseline solid)'}, ...
       'Location','northoutside','Orientation','horizontal');

set(gca,'Box','off','FontSize',11);

fig_file = fullfile(out_dir,'bout_rate_3h_baseline_vs_ambtemp_APPvsWT.png');
saveas(gcf, fig_file);

% save stats
if ~isempty(within_tbl)
    writetable(within_tbl, fullfile(out_dir,'bout_rate_3h_within_cond_stats.csv'));
end
if ~isempty(anova_tbl)
    writetable(anova_tbl, fullfile(out_dir,'bout_rate_3h_anova_stats.csv'));
end

OUT = struct();
OUT.success      = true;
OUT.fig_file     = fig_file;
OUT.within_stats = within_tbl;
OUT.anova        = anova_tbl;
OUT.state_order  = state_order;

fprintf('✅ 3 h bout-rate bar plot + stats saved in %s\n', out_dir);
end
