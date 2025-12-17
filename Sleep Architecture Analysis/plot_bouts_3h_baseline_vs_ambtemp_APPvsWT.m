function OUT = plot_bout_rate_3h_baseline_vs_ambtemp_APPvsWT(rows_overall_ambtemp, out_dir, varargin)
% plot_bout_rate_3h_baseline_vs_ambtemp_APPvsWT
% -------------------------------------------------------------------------
% Uses the cropped 3 h segments (baseline + ambtemp) OVERALL table:
%   state, n_bouts, total_dur_s, ..., condition, mouse, genotype
%
% For each state (WK, MA, NREM, REM, SLEEP):
%   - Compute bouts/hour in the 3 h window
%   - Bars:
%       WT baseline  (solid grey)
%       WT ambtemp   (grey outline)
%       APP baseline (solid blue)
%       APP ambtemp  (blue outline)
%   - Stats:
%       * Paired t-test baseline vs ambtemp within WT and within APP
%         (stars over the paired bars)
%       * Δ = ambtemp-baseline per mouse, WT vs APP (unpaired t-test)
%       * 2x2 mixed ANOVA (Condition x Genotype) via mixed-effects model
%
% INPUT
%   rows_overall_ambtemp : table from run_group_sleep_architecture_ambtemp_overall
%   out_dir              : folder for figure + CSVs
%
% OPTIONS
%   'states'               : which states to plot (default WK, MA, NREM, REM, SLEEP)
%   'minNperGroupForStats' : minimum n per group (default 3)
%
% OUTPUT struct OUT:
%   .fig_file
%   .within_stats   (paired tests)
%   .delta_stats    (Δ WT vs APP)
%   .anova          (Condition x Genotype)
%   .state_order
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'states',["WK","MA","NREM","REM","SLEEP"]);
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p,varargin{:});
state_order = string(p.Results.states(:))';
minN        = p.Results.minNperGroupForStats;

T = rows_overall_ambtemp;
T.condition = lower(strtrim(string(T.condition)));
T.state     = string(T.state);
T.mouse     = string(T.mouse);
T.genotype  = string(T.genotype);

% keep baseline + ambtemp only
keep = T.condition=="baseline" | T.condition=="ambtemp";
T    = T(keep,:);

if isempty(T)
    warning('No baseline/ambtemp rows in rows_overall_ambtemp.');
    OUT = struct('success',false);
    return;
end

% ------ 1) build bouts/hour (3h window) ------
T.window_h    = T.total_dur_s / 3600;
T.bouts_per_h = T.n_bouts ./ T.window_h;

% optionally add composite SLEEP (NREM+REM)
if any(state_order=="SLEEP")
    sleep_src = T(T.state=="NREM" | T.state=="REM",:);
    if ~isempty(sleep_src)
        Gsleep = groupsummary(sleep_src, ...
            {'mouse','genotype','condition'}, ...
            'sum',{'n_bouts','total_dur_s'});
        Gsleep.state        = repmat("SLEEP",height(Gsleep),1);
        Gsleep.bouts_per_h  = Gsleep.sum_n_bouts ./ (Gsleep.sum_total_dur_s/3600);
        Gsleep.total_dur_s  = Gsleep.sum_total_dur_s;
        Gsleep.n_bouts      = Gsleep.sum_n_bouts;
        Gsleep = removevars(Gsleep,{'sum_n_bouts','sum_total_dur_s'});
        T = [T; Gsleep];
    end
end

% keep only requested states
T = T(ismember(T.state,state_order),:);

genotypes = ["WT","APP"];
nStates   = numel(state_order);
nGen      = numel(genotypes);

% containers
within_rows = {};
delta_rows  = {};
anova_rows  = {};

% mean/SEM matrix for bars:
% rows: [WT_bas, WT_amb, APP_bas, APP_amb]
meanMat = nan(4,nStates);
semMat  = nan(4,nStates);

F_genotype    = nan(1,nStates);
p_genotype    = nan(1,nStates);
F_condition   = nan(1,nStates);
p_condition   = nan(1,nStates);
F_interaction = nan(1,nStates);
p_interaction = nan(1,nStates);

for si = 1:nStates
    st = state_order(si);
    Ts = T(T.state==st,:);
    if isempty(Ts), continue; end

    % ---------- A) within-genotype paired (baseline vs ambtemp) ----------
    for gi = 1:nGen
        g = genotypes(gi);
        Tg = Ts(Ts.genotype==g,:);
        if isempty(Tg), continue; end

        W = unstack(Tg(:,{'mouse','condition','bouts_per_h'}), ...
                    'bouts_per_h','condition');
        if ~all(ismember({'baseline','ambtemp'},W.Properties.VariableNames))
            continue;
        end
        base = W.baseline;
        amb  = W.ambtemp;
        ok   = ~isnan(base) & ~isnan(amb);
        base = base(ok);
        amb  = amb(ok);
        n    = numel(base);
        if n < minN, continue; end

        [~,p_pair] = ttest(base,amb);
        delta      = amb-base;

        within_rows = [within_rows; ...
            {char(st),char(g),n, ...
             mean(base,'omitnan'),mean(amb,'omitnan'), ...
             mean(delta,'omitnan'),p_pair}];

        mBase = mean(base,'omitnan');
        mAmb  = mean(amb,'omitnan');
        sBase = std(base,'omitnan')/sqrt(n);
        sAmb  = std(amb,'omitnan') /sqrt(n);

        if g=="WT"
            meanMat(1,si)=mBase; semMat(1,si)=sBase;
            meanMat(2,si)=mAmb;  semMat(2,si)=sAmb;
        else
            meanMat(3,si)=mBase; semMat(3,si)=sBase;
            meanMat(4,si)=mAmb;  semMat(4,si)=sAmb;
        end
    end

    % ---------- B) Δ WT vs APP ----------
    Wall = unstack(Ts(:,{'mouse','genotype','condition','bouts_per_h'}), ...
                   'bouts_per_h','condition');
    if ismember('baseline',Wall.Properties.VariableNames) && ...
       ismember('ambtemp',Wall.Properties.VariableNames)

        ok = ~isnan(Wall.baseline) & ~isnan(Wall.ambtemp);
        Wall = Wall(ok,:);
        Wall.delta = Wall.ambtemp - Wall.baseline;

        dWT  = Wall.delta(Wall.genotype=="WT");
        dAPP = Wall.delta(Wall.genotype=="APP");
        dWT  = dWT(~isnan(dWT));
        dAPP = dAPP(~isnan(dAPP));

        nWT  = numel(dWT);
        nAPP = numel(dAPP);
        if nWT>=minN && nAPP>=minN
            [~,p_delta] = ttest2(dWT,dAPP,'Vartype','unequal');
            mWT  = mean(dWT,'omitnan');
            mAPP = mean(dAPP,'omitnan');
            sWT  = std(dWT,'omitnan');
            sAPP = std(dAPP,'omitnan');
            sp   = sqrt(((nWT-1)*sWT^2 + (nAPP-1)*sAPP^2)/max(1,(nWT+nAPP-2)));
            d    = (mAPP-mWT)/sp;
        else
            p_delta = NaN; mWT=NaN; mAPP=NaN; d=NaN;
        end
    else
        p_delta=NaN; mWT=NaN; mAPP=NaN; d=NaN;
        nWT = sum(Wall.genotype=="WT");
        nAPP= sum(Wall.genotype=="APP");
    end

    delta_rows = [delta_rows; ...
        {char(st),mWT,mAPP,p_delta,d}];

    % ---------- C) 2x2 mixed ANOVA (Condition x Genotype) ----------
    Fg=NaN; pg=NaN; Fc=NaN; pc=NaN; Fi=NaN; pi=NaN;
    if exist('Wall','var') && height(Wall)>=2 && nWT>0 && nAPP>0
        try
            L = stack(Wall,{'baseline','ambtemp'}, ...
                      'NewDataVariableName','bouts', ...
                      'IndexVariableName','Condition');
            L.Condition = categorical(L.Condition,{'baseline','ambtemp'});
            L.genotype  = categorical(L.genotype);
            L.mouse     = categorical(L.mouse);

            lme = fitlme(L,'bouts ~ Condition*genotype + (1|mouse)', ...
                         'DummyVarCoding','effects');
            A = anova(lme);

            % rows
            rowC = find(strcmp(A.Term,'Condition'),1);
            rowG = find(strcmp(A.Term,'genotype'),1);
            rowI = find(contains(string(A.Term),'Condition:genotype'),1);

            if ~isempty(rowC), Fc=A.FStat(rowC); pc=A.pValue(rowC); end
            if ~isempty(rowG), Fg=A.FStat(rowG); pg=A.pValue(rowG); end
            if ~isempty(rowI), Fi=A.FStat(rowI); pi=A.pValue(rowI); end
        catch
            % leave NaNs
        end
    end

    F_genotype(si)    = Fg;  p_genotype(si)    = pg;
    F_condition(si)   = Fc;  p_condition(si)   = pc;
    F_interaction(si) = Fi;  p_interaction(si) = pi;

    anova_rows = [anova_rows; ...
        {char(st),nWT,nAPP,Fg,pg,Fc,pc,Fi,pi}];
end

% ------ build stats tables ------
if isempty(within_rows)
    within_tbl = table();
else
    within_tbl = cell2table(within_rows, ...
        'VariableNames',{'State','Genotype','n', ...
                         'Mean_baseline','Mean_ambtemp', ...
                         'Mean_delta','p_paired'});
end

if isempty(delta_rows)
    delta_tbl = table();
else
    delta_tbl = cell2table(delta_rows, ...
        'VariableNames',{'State','Mean_delta_WT','Mean_delta_APP', ...
                         'p_delta_WTvsAPP','Cohen_d_delta'});
end

if isempty(anova_rows)
    anova_tbl = table();
else
    anova_tbl = cell2table(anova_rows, ...
        'VariableNames',{'State','nWT','nAPP', ...
                         'F_genotype','p_genotype', ...
                         'F_condition','p_condition', ...
                         'F_interaction','p_interaction'});
end

% ------ 2) plotting ------
COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w','Position',[200 200 900 450]); hold on;
x = 1:nStates;
barWidth = 0.16;
offsets  = [-0.24 -0.08 0.08 0.24];

for si = 1:nStates
    if ~isnan(meanMat(1,si))
        xb = x(si)+offsets(1);
        bar(xb,meanMat(1,si),barWidth,'FaceColor',COL_WT,'EdgeColor','none');
        errorbar(xb,meanMat(1,si),semMat(1,si),'k','LineStyle','none');
    end
    if ~isnan(meanMat(2,si))
        xa = x(si)+offsets(2);
        bar(xa,meanMat(2,si),barWidth,'FaceColor','none','EdgeColor',COL_WT,'LineWidth',1.2);
        errorbar(xa,meanMat(2,si),semMat(2,si),'k','LineStyle','none');
    end
    if ~isnan(meanMat(3,si))
        xb = x(si)+offsets(3);
        bar(xb,meanMat(3,si),barWidth,'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(xb,meanMat(3,si),semMat(3,si),'k','LineStyle','none');
    end
    if ~isnan(meanMat(4,si))
        xa = x(si)+offsets(4);
        bar(xa,meanMat(4,si),barWidth,'FaceColor','none','EdgeColor',COL_APP,'LineWidth',1.2);
        errorbar(xa,meanMat(4,si),semMat(4,si),'k','LineStyle','none');
    end
end

all_means = meanMat(~isnan(meanMat));
if isempty(all_means), ymax = 10; else, ymax = max(all_means)+5; end
ylim([0 ymax]);

% stars for paired tests
for si = 1:nStates
    st = state_order(si);
    for gi = 1:nGen
        g = genotypes(gi);
        if isempty(within_tbl), continue; end
        row = strcmp(within_tbl.State,st) & strcmp(within_tbl.Genotype,g);
        if ~any(row), continue; end
        p_here = within_tbl.p_paired(row);
        if isnan(p_here) || p_here>=0.05, continue; end
        if p_here<0.001, stars='***';
        elseif p_here<0.01, stars='**';
        else, stars='*'; end

        if g=="WT"
            x1=x(si)+offsets(1); x2=x(si)+offsets(2); y_star=ymax-2;
        else
            x1=x(si)+offsets(3); x2=x(si)+offsets(4); y_star=ymax-5;
        end
        line([x1 x2],[y_star y_star],'Color','k','LineWidth',1);
        text(mean([x1 x2]),y_star+0.5,stars,'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom','FontSize',11,'FontWeight','bold');
    end
end

% text with Δ/Condition/Interaction p-values
for si = 1:nStates
    st = state_order(si);
    pD = NaN;
    if ~isempty(delta_tbl)
        rowD = strcmp(delta_tbl.State,st);
        if any(rowD), pD = delta_tbl.p_delta_WTvsAPP(find(rowD,1)); end
    end
    pC = p_condition(si);
    pI = p_interaction(si);
    if isnan(pD) && isnan(pC) && isnan(pI), continue; end
    txt = sprintf('Δ p=%.3f | C p=%.3f | G×C p=%.3f',pD,pC,pI);
    text(x(si),ymax-8,txt,'HorizontalAlignment','center','FontSize',8);
end

xticks(x); xticklabels(state_order);
xlabel('State');
ylabel('Bouts per hour (3 h window)');
title('3 h baseline vs 3 h ambtemp – bouts/hour (WT vs APP)');

hb1 = bar(nan,nan,barWidth,'FaceColor',COL_WT,'EdgeColor','none');
hb2 = bar(nan,nan,barWidth,'FaceColor',COL_APP,'EdgeColor','none');
ha1 = bar(nan,nan,barWidth,'FaceColor','none','EdgeColor','k','LineWidth',1.2);
legend([hb1 hb2 ha1],{'WT','APP','ambtemp (outline vs baseline solid)'}, ...
       'Location','northoutside','Orientation','horizontal');
set(gca,'Box','off','FontSize',11);

fig_file = fullfile(out_dir,'bout_rate_3h_baseline_vs_ambtemp_APPvsWT.png');
saveas(gcf,fig_file);

% save stats
if ~isempty(within_tbl)
    writetable(within_tbl,fullfile(out_dir,'bout_rate_within_cond_stats.csv'));
end
if ~isempty(delta_tbl)
    writetable(delta_tbl,fullfile(out_dir,'bout_rate_delta_genotype_stats.csv'));
end
if ~isempty(anova_tbl)
    writetable(anova_tbl,fullfile(out_dir,'bout_rate_anova_stats.csv'));
end

OUT = struct('success',true, ...
             'fig_file',fig_file, ...
             'within_stats',within_tbl, ...
             'delta_stats',delta_tbl, ...
             'anova',anova_tbl, ...
             'state_order',state_order);
fprintf('✅ 3 h bout-rate plot + stats saved in %s\n',out_dir);
end
