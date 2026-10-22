%% RUN_15_SPINDLE_PROBABILITY_STATS
% Mouse-level baseline/genotype and within-mouse intervention analysis.
% Uses run_14_spindle_probability.m exports; does not redetect spindles.
%
% Primary hypothesis tests:
%  (1) Cohort-adjusted genotype contrast at baseline, by NREM outcome.
%  (2) Cohort-adjusted genotype contrast of baseline-controlled responses
%      to warming and donepezil: (post-pre)_treatment - (late-early)_BL.
%
% Stratified permutations shuffle GENOTYPE within cohort, not bouts.
% q_BH controls exploratory FDR within each family of reported contrasts.
% This is a discovery screen, not confirmatory evidence of rescue.

clear; clc;
scriptFile=mfilename('fullpath');
if isempty(scriptFile), error('Save as run_15_spindle_probability_stats.m and run the saved script.'); end
repoDir=fileparts(fileparts(scriptFile));
addpath(fullfile(repoDir,'config'));
P=paths_config();
rootDir=fullfile(P.eeg_output,'spindles','probability');
figDir=fullfile(rootDir,'figures');
if ~isfolder(figDir), mkdir(figDir); end
inFile=fullfile(rootDir,'mouse_probability_summary_25s.csv');
assert(isfile(inFile),'Missing %s. Run run_14_spindle_probability.m first.',inFile);

M=readtable(inFile,'TextType','string');
required={'mouse_id','genotype','cohort','condition','phase', ...
    'transition_type','n_bouts','occupancy_probability_25_15', ...
    'any_spindle_25','occupancy_fraction_25'};
missing=setdiff(required,M.Properties.VariableNames);
assert(isempty(missing),'Missing columns in mouse-level summary: %s',strjoin(missing,', '));
for k={'mouse_id','genotype','cohort','condition','phase','transition_type'}
    field=k{1}; M.(field)=string(M.(field));
end

MIN_BOUTS_MOUSE=3;
N_PERM=4000; % reproducible Monte Carlo randomization p-values
rng(20261022,'twister');
metrics=["occupancy_probability_25_15", "any_spindle_25", ...
    "occupancy_fraction_25"];
transitions=["NREM_to_REM","NREM_to_WAKE"];
basePhases=["BASELINE_ZT3_9","BASELINE_ZT3_6","BASELINE_ZT6_9"];
interventions=["warming","donepezil"];
keys={'mouse_id','genotype','cohort','transition_type'};

BaseRows=struct('phase',{},'transition_type',{},'metric',{}, ...
    'n_WT',{},'n_APP',{},'mean_WT',{},'mean_APP',{}, ...
    'APP_minus_WT_adj_cohort',{},'p_stratified_permutation',{});
IntRows=struct('intervention',{},'transition_type',{},'metric',{}, ...
    'n_WT',{},'n_APP',{},'mean_DiD_WT',{},'mean_DiD_APP',{}, ...
    'APP_minus_WT_DiD_adj_cohort',{},'p_stratified_permutation',{});
CohortRows=struct('analysis',{},'phase_or_intervention',{}, ...
    'transition_type',{},'metric',{},'cohort',{},'n_WT',{},'n_APP',{}, ...
    'mean_WT',{},'mean_APP',{},'APP_minus_WT',{});
Responses=table();

fprintf('\n=== SPINDLE PROBABILITY / MOUSE-LEVEL EXPLORATORY STATS ===\n');
fprintf('Min eligible bouts/mouse/phase/outcome: %d\n',MIN_BOUTS_MOUSE);

%% 1) Genotype differences during baseline
for b=1:numel(basePhases)
  for t=1:numel(transitions)
    S=M(M.condition=="baseline" & M.phase==basePhases(b) & ...
      M.transition_type==transitions(t) & M.n_bouts>=MIN_BOUTS_MOUSE,:);
    for mm=1:numel(metrics)
      variable=char(metrics(mm));
      y=double(S.(variable));
      [beta,p,nWT,nAPP,mWT,mAPP]=strat_perm( ...
          y,S.genotype,S.cohort,N_PERM);
      BaseRows(end+1)=struct( ...
        'phase',basePhases(b), ...
        'transition_type',transitions(t),'metric',metrics(mm), ...
        'n_WT',nWT,'n_APP',nAPP,'mean_WT',mWT, ...
        'mean_APP',mAPP,'APP_minus_WT_adj_cohort',beta, ...
        'p_stratified_permutation',p); %#ok<SAGROW>
      CohortRows=append_cohort_rows(CohortRows,'baseline',basePhases(b), ...
          transitions(t),metrics(mm),S,y);
    end
  end
end

%% 2) Mouse-paired treatment deltas minus matched baseline ZT changes
for iv=1:numel(interventions)
  intervention=interventions(iv);
  if intervention=="warming"
    prePhase="PREWARM_ZT3_6"; postPhase="WARM_ACTIVE";
  else
    prePhase="PREDRUG_ZT3_6"; postPhase="DRUG_0_3H";
  end
  for t=1:numel(transitions)
    for mm=1:numel(metrics)
      variable=char(metrics(mm));
      Tpre=pick_window(M,intervention,prePhase,transitions(t), ...
          variable,'treatment_pre',keys,MIN_BOUTS_MOUSE);
      Tpost=pick_window(M,intervention,postPhase,transitions(t), ...
          variable,'treatment_post',keys,MIN_BOUTS_MOUSE);
      BLearly=pick_window(M,"baseline","BASELINE_ZT3_6",transitions(t), ...
          variable,'baseline_early',keys,MIN_BOUTS_MOUSE);
      BLlate=pick_window(M,"baseline","BASELINE_ZT6_9",transitions(t), ...
          variable,'baseline_late',keys,MIN_BOUTS_MOUSE);
      D=innerjoin(Tpre,Tpost,'Keys',keys);
      D=innerjoin(D,BLearly,'Keys',keys);
      D=innerjoin(D,BLlate,'Keys',keys);
      D.delta_treatment=D.treatment_post-D.treatment_pre;
      D.delta_baseline=D.baseline_late-D.baseline_early;
      D.difference_in_differences=D.delta_treatment-D.delta_baseline;
      D.metric=repmat(metrics(mm),height(D),1);
      D.intervention=repmat(intervention,height(D),1);
      if height(D)>0
          if isempty(Responses)
              Responses=D;
          else
              Responses=[Responses;D]; %#ok<AGROW>
          end
      end
      y=double(D.difference_in_differences);
      [beta,p,nWT,nAPP,mWT,mAPP]=strat_perm( ...
          y,D.genotype,D.cohort,N_PERM);
      IntRows(end+1)=struct( ...
        'intervention',intervention, ...
        'transition_type',transitions(t),'metric',metrics(mm), ...
        'n_WT',nWT,'n_APP',nAPP,'mean_DiD_WT',mWT, ...
        'mean_DiD_APP',mAPP,'APP_minus_WT_DiD_adj_cohort',beta, ...
        'p_stratified_permutation',p); %#ok<SAGROW>
      CohortRows=append_cohort_rows(CohortRows,'intervention_DiD',intervention, ...
          transitions(t),metrics(mm),D,y);
    end
  end
end

BaselineStats=struct2table(BaseRows);
InterventionStats=struct2table(IntRows);
CohortDescriptives=struct2table(CohortRows);
BaselineStats.q_BH=bh_qvalues(BaselineStats.p_stratified_permutation);
InterventionStats.q_BH=bh_qvalues(InterventionStats.p_stratified_permutation);
writetable(BaselineStats,fullfile(rootDir,'stats_baseline_genotype.csv'));
writetable(InterventionStats,fullfile(rootDir,'stats_intervention_DiD.csv'));
writetable(CohortDescriptives,fullfile(rootDir,'cohort_descriptives.csv'));
if ~isempty(Responses)
    writetable(Responses,fullfile(rootDir,'mouse_paired_intervention_responses.csv'));
else
    warning('No complete 4-window mouse pairs found for intervention DiD.');
end

%% 3) Plot baseline scalar effects (each dot = one mouse)
colors=[0.42 0.42 0.42; 0.16 0.46 0.78];
f=figure('Color','w','Position',[70 70 1360 690]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for t=1:2
 for mm=1:3
    nexttile; hold on;
    S=M(M.condition=="baseline" & M.phase=="BASELINE_ZT3_9" & ...
        M.transition_type==transitions(t) & M.n_bouts>=MIN_BOUTS_MOUSE,:);
    y=100*double(S.(char(metrics(mm))));
    mouse_points(S,y,colors);
    ylabel('Mouse mean probability (%)');
    title(sprintf('%s | %s',char(transitions(t)),char(metrics(mm))), ...
        'Interpreter','none');
    box off;
 end
end
export_three(f,figDir,'baseline_mouse_probabilities');

%% 4) Baseline-controlled responses (difference in percentage points)
for mm=1:3
 f=figure('Color','w','Position',[70 70 1020 680]);
 tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
 for t=1:2
   for iv=1:2
      nexttile; hold on;
      yline(0,':','Color',[.5 .5 .5],'HandleVisibility','off');
      if isempty(Responses), R=table();
      else
          R=Responses(Responses.transition_type==transitions(t) & ...
              Responses.intervention==interventions(iv) & ...
              Responses.metric==metrics(mm),:);
      end
      if ~isempty(R)
          mouse_points(R,100*R.difference_in_differences,colors);
      else
          xlim([.5 4.5]); xticks(1:4);
          xticklabels({'Oct WT','Oct APP','June WT','June APP'});
          text(2.5,0,'No complete pairs','HorizontalAlignment','center');
      end
      ylabel('Baseline-controlled change (percentage points)');
      title(sprintf('%s | %s',char(interventions(iv)),char(transitions(t))), ...
           'Interpreter','none');
      box off;
   end
 end
 export_three(f,figDir,'intervention_DiD_'+metrics(mm));
end

% Human-readable analysis notes, including important null model caveats.
notes=fullfile(rootDir,'interpretation_notes.txt');
fid=fopen(notes,'w');
fprintf(fid,'SPINDLE PROBABILITY: ANALYSIS AND INTERPRETATION\n');
fprintf(fid,'1. Each dot is ONE MOUSE; bouts are not independent group observations.\n');
fprintf(fid,'2. 1-s spindle occupancy: proportion of aligned seconds overlapped by a spindle.\n');
fprintf(fid,'3. occupancy_probability_25_15: ten bins [-25,-15) s before NREM end.\n');
fprintf(fid,'4. any_spindle_25: fraction of eligible transitions having a spindle in [-25,0) s.\n');
fprintf(fid,'5. Baseline: APP minus WT coefficient after adjusting for October/June cohort.\n');
fprintf(fid,'6. Treatment DiD: (post - pre)_treatment - (late - early)_baseline, per mouse.\n');
fprintf(fid,'7. DiD does NOT rule out carryover/injection/circadian confounding.\n');
fprintf(fid,'8. Stratified permutation exchanges genotype labels within cohorts only.\n');
fprintf(fid,'9. q_BH controls FDR per results family; many analyses remain exploratory.\n');
fprintf(fid,'10. Cohort or genotype-specific missing mice reduce power; inspect n_WT/n_APP.\n');
fprintf(fid,'11. MIN_BOUTS_MOUSE=%d, N_PERM=%d.\n',MIN_BOUTS_MOUSE,N_PERM);
fprintf(fid,'12. This screen cannot establish circuit mechanism or cognitive benefit.\n');
fclose(fid);

fprintf('\nDone. Baseline tests: %d | treatment tests: %d\n', ...
    height(BaselineStats),height(InterventionStats));
fprintf('Review %s\n',rootDir);
disp(BaselineStats(BaselineStats.phase=="BASELINE_ZT3_9", ...
    {'transition_type','metric','n_WT','n_APP', ...
     'APP_minus_WT_adj_cohort','p_stratified_permutation','q_BH'}));
disp(InterventionStats(:,{'intervention','transition_type','metric', ...
    'n_WT','n_APP','APP_minus_WT_DiD_adj_cohort', ...
    'p_stratified_permutation','q_BH'}));

%% Local helpers
function S=pick_window(M,condition,phase,transition,metric,newname,keys,minBouts)
    q=M.condition==condition & M.phase==phase & ...
       M.transition_type==transition & M.n_bouts>=minBouts;
    S=M(q,[keys,{metric,'n_bouts'}]);
    vars=S.Properties.VariableNames;
    vars{strcmp(vars,metric)}=newname;
    vars{strcmp(vars,'n_bouts')}=['n_' newname];
    S.Properties.VariableNames=vars;
end

function [effect,p,nWT,nAPP,wtMean,appMean]=strat_perm(y,genotype,cohort,nPerm)
    genotype=string(genotype); cohort=string(cohort);
    good=isfinite(y) & ismember(genotype,["WT","APP"]) & ~ismissing(cohort);
    y=y(good); genotype=genotype(good); cohort=cohort(good);
    nWT=nnz(genotype=="WT"); nAPP=nnz(genotype=="APP");
    wtMean=NaN; appMean=NaN; effect=NaN; p=NaN;
    if nWT>0, wtMean=mean(y(genotype=="WT")); end
    if nAPP>0, appMean=mean(y(genotype=="APP")); end
    if nWT<2 || nAPP<2, return; end
    g=double(genotype=="APP");
    strata=unique(cohort);
    gc=zeros(size(g)); yc=zeros(size(y));
    for s=1:numel(strata)
        ix=cohort==strata(s);
        gc(ix)=g(ix)-mean(g(ix));
        yc(ix)=y(ix)-mean(y(ix));
    end
    den=sum(gc.^2);
    if den<=eps, return; end
    effect=sum(gc.*yc)/den;
    nullStats=zeros(nPerm,1);
    for j=1:nPerm
        gp=g;
        for s=1:numel(strata)
            ii=find(cohort==strata(s));
            gp(ii)=g(ii(randperm(numel(ii))));
        end
        z=zeros(size(gp));
        for s=1:numel(strata)
            ii=cohort==strata(s);
            z(ii)=gp(ii)-mean(gp(ii));
        end
        nullStats(j)=sum(z.*yc)/den;
    end
    p=(1+nnz(abs(nullStats)>=abs(effect)-1e-12))/(nPerm+1);
end

function q=bh_qvalues(p)
    q=nan(size(p));
    good=find(isfinite(p));
    if isempty(good), return; end
    [a,ord]=sort(p(good));
    m=numel(a);
    adjusted=min(1,a.*m./(1:m)');
    adjusted=flipud(cummin(flipud(adjusted)));
    q(good(ord))=adjusted;
end

function rows=append_cohort_rows(rows,analysis,phase,transition,metric,S,y)
    if isempty(S), return; end
    for c=unique(S.cohort)'
        select=S.cohort==c;
        wt=select & S.genotype=="WT";
        app=select & S.genotype=="APP";
        nWT=nnz(wt); nAPP=nnz(app);
        mWT=NaN; mAPP=NaN;
        if nWT>0, mWT=mean(y(wt)); end
        if nAPP>0, mAPP=mean(y(app)); end
        rows(end+1)=struct('analysis',string(analysis), ...
            'phase_or_intervention',string(phase), ...
            'transition_type',string(transition),'metric',string(metric), ...
            'cohort',c,'n_WT',nWT,'n_APP',nAPP, ...
            'mean_WT',mWT,'mean_APP',mAPP, ...
            'APP_minus_WT',mAPP-mWT); %#ok<AGROW>
    end
end

function mouse_points(T,y,colors)
    cohortNames=["thesis","new"];
    genotypeNames=["WT","APP"];
    for c=1:2
        for g=1:2
            x=(c-1)*2+g;
            mask=T.cohort==cohortNames(c) & T.genotype==genotypeNames(g) & isfinite(y);
            vals=y(mask);
            if isempty(vals), continue; end
            jitter=linspace(-.10,.10,numel(vals))';
            plot(x+jitter,vals,'o','MarkerSize',6,'LineStyle','none', ...
                'MarkerEdgeColor',colors(g,:),'MarkerFaceColor',colors(g,:));
            mu=mean(vals);
            if numel(vals)>1, sem=std(vals)/sqrt(numel(vals));
            else, sem=0; end
            errorbar(x,mu,sem,'k','LineWidth',1.8,'CapSize',9);
        end
    end
    xlim([.5 4.5]); xticks(1:4);
    xticklabels({'Oct WT','Oct APP','June WT','June APP'});
    xtickangle(20);
end

function export_three(f,figDir,name)
    exportgraphics(f,fullfile(figDir,[char(name) '.png']),'Resolution',300);
    exportgraphics(f,fullfile(figDir,[char(name) '.pdf']),'ContentType','vector');
    savefig(f,fullfile(figDir,[char(name) '.fig']));
    close(f);
end
