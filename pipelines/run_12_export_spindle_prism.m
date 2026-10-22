%% RUN_12_EXPORT_SPINDLE_PRISM
% Export mouse-level spindle datasets for GraphPad Prism.
%
% Requires:
%   run('pipelines/run_10_spindle_analysis.m')
%
% Output:
%   results/prism/spindles/

clear;
clc;

pipelineDir=fileparts(mfilename('fullpath'));
repoDir=fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));
P=paths_config();

inDir=fullfile(P.eeg_output,'spindles');
outDir=fullfile(repoDir,'results','prism','spindles');

if ~isfolder(outDir)
    mkdir(outDir);
end

WIN=readtable( ...
    fullfile(inDir,'spindle_window_summary.csv'), ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

TM=readtable( ...
    fullfile(inDir,'spindle_transition_mouse_summary.csv'), ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

TC=readtable( ...
    fullfile(inDir,'spindle_transition_timecourse_mouse.csv'), ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

%% ------------------------------------------------------------------------
% Figure 1: baseline spindle phenotype
% -------------------------------------------------------------------------

B=WIN( ...
    WIN.condition=="baseline" & ...
    WIN.window=="BASELINE_ZT3_9",:);

writetable(B, ...
    fullfile(outDir,'baseline_spindle_density_ZT3_9.csv'));

QB=TM( ...
    TM.condition=="baseline" & ...
    TM.session_phase=="BASELINE_ZT3_9",:);

BASETRANS=make_transition_wide(QB);

writetable(BASETRANS, ...
    fullfile(outDir,'baseline_preREM_vs_preWake_spindle_rate.csv'));

TBASE=TC( ...
    TC.condition=="baseline" & ...
    TC.session_phase=="BASELINE_ZT3_9",:);

writetable(TBASE, ...
    fullfile(outDir,'baseline_transition_timecourse.csv'));

%% ------------------------------------------------------------------------
% Figure 3: warming spindle response
% -------------------------------------------------------------------------

WARM=make_window_response( ...
    WIN, ...
    "warming", ...
    "WARM_PRE", ...
    "WARM_ACTIVE", ...
    "WARM_MATCH_PRE", ...
    "WARM_MATCH_ACTIVE");

writetable(WARM, ...
    fullfile(outDir,'warming_spindle_density_response.csv'));

QW=TM(TM.condition=="warming",:);

writetable(QW, ...
    fullfile(outDir,'warming_transition_spindle_metrics.csv'));

TW=TC(TC.condition=="warming",:);

writetable(TW, ...
    fullfile(outDir,'warming_transition_timecourse.csv'));

%% ------------------------------------------------------------------------
% Donepezil spindle response
% -------------------------------------------------------------------------

DRUG=make_window_response( ...
    WIN, ...
    "donepezil", ...
    "DRUG_PRE", ...
    "DRUG_0_3H", ...
    "DRUG_MATCH_PRE", ...
    "DRUG_MATCH_0_3H");

writetable(DRUG, ...
    fullfile(outDir,'donepezil_spindle_density_response.csv'));

QD=TM(TM.condition=="donepezil",:);

writetable(QD, ...
    fullfile(outDir,'donepezil_transition_spindle_metrics.csv'));

TD=TC(TC.condition=="donepezil",:);

writetable(TD, ...
    fullfile(outDir,'donepezil_transition_timecourse.csv'));

%% Manifest

files=dir(fullfile(outDir,'*.csv'));

Manifest=table( ...
    string({files.name})', ...
    [files.bytes]', ...
    'VariableNames',{'file','bytes'});

writetable(Manifest, ...
    fullfile(outDir,'spindle_prism_manifest.csv'));

fprintf('\n============================================================\n');
fprintf('SPINDLE PRISM EXPORT COMPLETE\n');
fprintf('============================================================\n');
fprintf('%s\n',outDir);
fprintf('============================================================\n\n');

%% ========================================================================
% Local functions
% =========================================================================

function O = make_transition_wide(T)

O=table();

if isempty(T)
    return;
end

keys=[ ...
    "mouse_id","genotype","sex","cohort","cohort_display"];

[G,K]=findgroups(T(:,keys));

O=K;

O.preREM_rate=nan(height(K),1);
O.preWake_rate=nan(height(K),1);

O.preREM_n_bouts=zeros(height(K),1);
O.preWake_n_bouts=zeros(height(K),1);

O.preREM_last_spindle_distance_s=nan(height(K),1);
O.preWake_last_spindle_distance_s=nan(height(K),1);

for g=1:height(K)

    r=T(G==g & T.transition_type=="NREM_to_REM",:);
    w=T(G==g & T.transition_type=="NREM_to_WAKE",:);

    if ~isempty(r)
        O.preREM_rate(g)=mean( ...
            r.mean_spindle_rate_last25_per_min,'omitnan');

        O.preREM_n_bouts(g)=sum(r.n_transition_bouts,'omitnan');

        O.preREM_last_spindle_distance_s(g)=mean( ...
            r.mean_last_spindle_peak_distance_s,'omitnan');
    end

    if ~isempty(w)
        O.preWake_rate(g)=mean( ...
            w.mean_spindle_rate_last25_per_min,'omitnan');

        O.preWake_n_bouts(g)=sum(w.n_transition_bouts,'omitnan');

        O.preWake_last_spindle_distance_s(g)=mean( ...
            w.mean_last_spindle_peak_distance_s,'omitnan');
    end
end

O.preREM_minus_preWake_rate = ...
    O.preREM_rate-O.preWake_rate;
end


function O = make_window_response( ...
    W,condition,preWindow,activeWindow,basePreWindow,baseActiveWindow)

mice=unique(W.mouse_id(W.condition==condition));

O=table();

for i=1:numel(mice)

    m=mice(i);

    C=W(W.mouse_id==m & W.condition==condition,:);
    B=W(W.mouse_id==m & W.condition=="baseline",:);

    if isempty(C)
        continue;
    end

    pre=get_density(C,preWindow);
    active=get_density(C,activeWindow);

    basePre=get_density(B,basePreWindow);
    baseActive=get_density(B,baseActiveWindow);

    response=active-pre;
    baselineChange=baseActive-basePre;
    did=response-baselineChange;

    row=table( ...
        m, ...
        string(C.genotype(1)), ...
        string(C.sex(1)), ...
        string(C.cohort(1)), ...
        string(C.cohort_display(1)), ...
        basePre,baseActive,pre,active, ...
        baselineChange,response,did, ...
        'VariableNames',{ ...
        'mouse_id','genotype','sex','cohort','cohort_display', ...
        'baseline_pre','baseline_active_matched', ...
        'condition_pre','condition_active', ...
        'baseline_change','condition_change', ...
        'difference_in_differences'});

    O=[O;row]; %#ok<AGROW>
end
end


function x = get_density(T,windowName)

x=NaN;

Q=T(T.window==windowName,:);

if height(Q)==1
    x=double(Q.spindles_per_min_nrem);
end
end
