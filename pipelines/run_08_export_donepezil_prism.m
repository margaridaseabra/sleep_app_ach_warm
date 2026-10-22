%% RUN_08_EXPORT_DONEPEZIL_PRISM
%
% Export the injection-aligned donepezil analysis into lightweight Prism-ready CSV files.
%
% Run AFTER:
%   run('pipelines/run_donepezil_complete.m')
%
% Output:
%   results/prism/donepezil/
%
% Creates:
%   source/   = transparent mousewise tables
%   grouped/  = Prism-style grouped tables for raw and delta plots
%   survival/ = exact mouse-balanced KM coordinates
%   stats/    = MATLAB statistics
%
% -------------------------------------------------------------------------

clear;
clc;

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));

P = paths_config();

Mice = readtable(P.mice, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');
Mice.cohort(Mice.cohort == "thesis") = "original";

inDir = fullfile(P.sleep_output,'donepezil_complete');
survInDir = fullfile(inDir,'survival');

outDir = fullfile(repoDir,'results','prism','donepezil');
sourceDir = fullfile(outDir,'source');
groupDir = fullfile(outDir,'grouped');
survDir = fullfile(outDir,'survival');
statsDir = fullfile(outDir,'stats');

dirs = {outDir,sourceDir,groupDir,survDir,statsDir};

for i = 1:numel(dirs)
    if ~isfolder(dirs{i}), mkdir(dirs{i}); end
end

RAW = readtable( ...
    fullfile(inDir,'donepezil_raw_metrics_long.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

DELTA = readtable( ...
    fullfile(inDir,'donepezil_delta_from_pre.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

MATCHED = readtable( ...
    fullfile(inDir,'donepezil_matched_baseline_vs_donepezil.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

DID = readtable( ...
    fullfile(inDir,'donepezil_difference_in_differences.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

SURV = readtable( ...
    fullfile(survInDir,'donepezil_survival_dataset.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

%% Normalize historical cohort name for figure readability

RAW.cohort(RAW.cohort == "thesis") = "original";
DELTA.cohort(DELTA.cohort == "thesis") = "original";
MATCHED.cohort(MATCHED.cohort == "thesis") = "original";
DID.cohort(DID.cohort == "thesis") = "original";
SURV.cohort(SURV.cohort == "thesis") = "original";

%% Metric definitions

defs = {
    "state_percent","percent_time","WK","state_percent_WK";
    "state_percent","percent_time","NREM","state_percent_NREM";
    "state_percent","percent_time","REM","state_percent_REM";
    "state_percent","percent_time","MA","state_percent_MA";
    "bout_metric","bouts_per_h","WK","bout_rate_per_h_WK";
    "bout_metric","bouts_per_h","NREM","bout_rate_per_h_NREM";
    "bout_metric","bouts_per_h","REM","bout_rate_per_h_REM";
    "bout_metric","bouts_per_h","MA","bout_rate_per_h_MA";
    "bout_metric","mean_bout_duration_s","WK","bout_mean_duration_WK";
    "bout_metric","mean_bout_duration_s","NREM","bout_mean_duration_NREM";
    "bout_metric","mean_bout_duration_s","REM","bout_mean_duration_REM";
    "bout_metric","mean_bout_duration_s","MA","bout_mean_duration_MA";
    "bout_metric","median_bout_duration_s","WK","bout_median_duration_WK";
    "bout_metric","median_bout_duration_s","NREM","bout_median_duration_NREM";
    "bout_metric","median_bout_duration_s","REM","bout_median_duration_REM";
    "bout_metric","median_bout_duration_s","MA","bout_median_duration_MA";
    "fragmentation","MA_per_h_state","NREM","fragmentation_NREM";
    "fragmentation","MA_per_h_state","REM","fragmentation_REM";
    };

MANIFEST = table();

%% Source + grouped Prism tables

for i = 1:size(defs,1)

    domain = string(defs{i,1});
    metric = string(defs{i,2});
    state = string(defs{i,3});
    stem = string(defs{i,4});

    %% Source table: one mouse per row

    mouseIDs = unique( ...
        DELTA.mouse_id( ...
            DELTA.domain == domain & ...
            DELTA.metric == metric & ...
            DELTA.state == state), ...
        'stable');

    if isempty(mouseIDs)
        continue;
    end

    % Metadata comes from mice.csv, not from repeated analysis rows.
    Tsource = table();
    for mm = 1:numel(mouseIDs)
        idxM = find(Mice.mouse_id == mouseIDs(mm),1);
        if isempty(idxM)
            error('Mouse %s not found in mice.csv.',mouseIDs(mm));
        end
        Tsource = [Tsource; Mice(idxM,["mouse_id","genotype","cohort","sex"])]; %#ok<AGROW>
    end

    sourceCols = [ ...
        "donepezil_pre","donepezil_0_3h","donepezil_3_6h", ...
        "delta_pre","delta_0_3h","delta_3_6h", ...
        "baseline_pre","baseline_0_3h","baseline_3_6h", ...
        "donepezil_minus_baseline_pre", ...
        "donepezil_minus_baseline_0_3h", ...
        "donepezil_minus_baseline_3_6h", ...
        "DiD_0_3h","DiD_3_6h"];

    for c = sourceCols
        Tsource.(c) = nan(height(Tsource),1);
    end

    for r = 1:height(Tsource)

        m = Tsource.mouse_id(r);

        Dd = DELTA( ...
            DELTA.mouse_id == m & ...
            DELTA.domain == domain & ...
            DELTA.metric == metric & ...
            DELTA.state == state,:);

        Tsource.donepezil_pre(r) = pick_value(Dd,"PRE","value");
        Tsource.donepezil_0_3h(r) = pick_value(Dd,"DRUG_0_3H","value");
        Tsource.donepezil_3_6h(r) = pick_value(Dd,"DRUG_3_6H","value");

        Tsource.delta_pre(r) = pick_value(Dd,"PRE","delta_from_pre");
        Tsource.delta_0_3h(r) = pick_value(Dd,"DRUG_0_3H","delta_from_pre");
        Tsource.delta_3_6h(r) = pick_value(Dd,"DRUG_3_6H","delta_from_pre");

        Dm = MATCHED( ...
            MATCHED.mouse_id == m & ...
            MATCHED.domain == domain & ...
            MATCHED.metric == metric & ...
            MATCHED.state == state,:);

        Tsource.baseline_pre(r) = pick_value(Dm,"PRE","baseline_value");
        Tsource.baseline_0_3h(r) = pick_value(Dm,"DRUG_0_3H","baseline_value");
        Tsource.baseline_3_6h(r) = pick_value(Dm,"DRUG_3_6H","baseline_value");

        Tsource.donepezil_minus_baseline_pre(r) = ...
            pick_value(Dm,"PRE","donepezil_minus_baseline");

        Tsource.donepezil_minus_baseline_0_3h(r) = ...
            pick_value(Dm,"DRUG_0_3H","donepezil_minus_baseline");

        Tsource.donepezil_minus_baseline_3_6h(r) = ...
            pick_value(Dm,"DRUG_3_6H","donepezil_minus_baseline");

        Di = DID( ...
            DID.mouse_id == m & ...
            DID.domain == domain & ...
            DID.metric == metric & ...
            DID.state == state,:);

        Tsource.DiD_0_3h(r) = pick_value(Di,"DRUG_0_3H","difference_in_differences");
        Tsource.DiD_3_6h(r) = pick_value(Di,"DRUG_3_6H","difference_in_differences");
    end

    sourceFile = fullfile(sourceDir,char(stem + ".csv"));
    writetable(Tsource,sourceFile);

    %% Grouped raw donepezil-session table

    draw = DELTA( ...
        DELTA.domain == domain & ...
        DELTA.metric == metric & ...
        DELTA.state == state,:);

    rawFile = fullfile(groupDir,char(stem + "_donepezil_raw_grouped.csv"));
    write_grouped_three_phase(draw,"value",rawFile);

    %% Grouped delta-from-pre table

    deltaFile = fullfile(groupDir,char(stem + "_delta_from_pre_grouped.csv"));
    write_grouped_three_phase(draw,"delta_from_pre",deltaFile);

    %% Grouped matched baseline-vs-donepezil at donepezil window

    Dmatch = MATCHED( ...
        MATCHED.domain == domain & ...
        MATCHED.metric == metric & ...
        MATCHED.state == state & ...
        MATCHED.time_bin == "DRUG_0_3H",:);

    conditionFile = fullfile(groupDir,char(stem + "_DRUG_0_3H_baseline_vs_donepezil_grouped.csv"));
    write_grouped_condition(Dmatch,conditionFile);

    %% Manifest

    MANIFEST = [MANIFEST; table( ...
        stem,domain,state,metric, ...
        string(sourceFile),string(rawFile),string(deltaFile),string(conditionFile), ...
        'VariableNames',{ ...
            'figure_key','domain','state','metric', ...
            'source_csv','raw_grouped_csv','delta_grouped_csv', ...
            'DRUG_0_3H_condition_grouped_csv'})]; %#ok<AGROW>
end

%% Survival coordinates

views = ["original","new","all"];
bins = ["PRE","DRUG_0_3H","DRUG_3_6H"];
states = ["NREM","REM"];

for v = views

    if v == "all"
        V = SURV;
    else
        V = SURV(SURV.cohort == v,:);
    end

    for b = bins
        for st = states

            %% Donepezil-only WT vs APP

            D = V( ...
                V.condition == "donepezil" & ...
                V.time_bin == b & ...
                V.state == st,:);

            if ~isempty(D)
                T = two_group_km_table(D);
                writetable(T, ...
                    fullfile(survDir, ...
                        sprintf('donepezil_%s_%s_%s.csv',char(v),char(b),char(st))));
            end

            %% Four curves: WT/APP x baseline/donepezil

            D = V(V.time_bin == b & V.state == st,:);

            if ~isempty(D)
                T = four_group_km_table(D);
                writetable(T, ...
                    fullfile(survDir, ...
                        sprintf('condition_%s_%s_%s.csv',char(v),char(b),char(st))));
            end
        end
    end
end

%% Copy statistics

statFiles = {
    'donepezil_delta_genotype_stats.csv'
    'donepezil_DiD_genotype_stats.csv'
    'donepezil_LME_time_within_session.csv'
    'donepezil_LME_baseline_vs_donepezil.csv'
    };

for i = 1:numel(statFiles)
    src = fullfile(inDir,statFiles{i});
    if isfile(src)
        copyfile(src,fullfile(statsDir,statFiles{i}));
    end
end

survStatFiles = {
    'donepezil_only_survival_cox.csv'
    'donepezil_only_survival_km_summary.csv'
    'donepezil_within_session_pre_vs_drug_cox_interaction.csv'
    'donepezil_baseline_vs_donepezil_cox_interaction.csv'
    'donepezil_survival_interaction_leave_one_mouse_out.csv'
    };

for i = 1:numel(survStatFiles)
    src = fullfile(survInDir,survStatFiles{i});
    if isfile(src)
        copyfile(src,fullfile(statsDir,survStatFiles{i}));
    end
end

writetable(MANIFEST,fullfile(outDir,'donepezil_prism_manifest.csv'));

fprintf('\n============================================================\n');
fprintf('DONEPEZIL PRISM EXPORT COMPLETE\n');
fprintf('============================================================\n');
fprintf('Output:\n%s\n\n',outDir);

%% ========================================================================
% Local functions
% =========================================================================

function value = pick_value(T,timeBin,varName)

value = NaN;

if isempty(T)
    return;
end

x = double(T.(varName)(T.time_bin == timeBin));
x = x(isfinite(x));

if isempty(x)
    return;
end

xu = unique(x);

if numel(xu) == 1
    value = xu(1);
else
    error('Conflicting duplicate values for %s at %s.',varName,timeBin);
end
end


function write_grouped_three_phase(D,varName,outFile)

WTmice = unique(D.mouse_id(D.genotype == "WT"));
APPmice = unique(D.mouse_id(D.genotype == "APP"));

nWT = numel(WTmice);
nAPP = numel(APPmice);

C = cell(5,1+nWT+nAPP);

C{1,1} = '';
C{2,1} = 'Phase';

for j = 1:nWT
    C{1,1+j} = 'WT';
    C{2,1+j} = char(WTmice(j));
end

for j = 1:nAPP
    col = 1+nWT+j;
    C{1,col} = 'APP';
    C{2,col} = char(APPmice(j));
end

bins = ["PRE","DRUG_0_3H","DRUG_3_6H"];
labels = {'Pre','Drug 0-3 h','Drug 3-6 h'};

for b = 1:3

    C{2+b,1} = labels{b};

    for j = 1:nWT
        x = double(D.(varName)(D.mouse_id == WTmice(j) & D.time_bin == bins(b)));
        x = unique(x(isfinite(x)));
        if numel(x) == 1
            C{2+b,1+j} = x;
        elseif numel(x) > 1
            error('Conflicting duplicate WT values for %s.',WTmice(j));
        end
    end

    for j = 1:nAPP
        col = 1+nWT+j;
        x = double(D.(varName)(D.mouse_id == APPmice(j) & D.time_bin == bins(b)));
        x = unique(x(isfinite(x)));
        if numel(x) == 1
            C{2+b,col} = x;
        elseif numel(x) > 1
            error('Conflicting duplicate APP values for %s.',APPmice(j));
        end
    end
end

writecell(C,outFile);
end


function write_grouped_condition(D,outFile)

WTmice = unique(D.mouse_id(D.genotype == "WT"));
APPmice = unique(D.mouse_id(D.genotype == "APP"));

nWT = numel(WTmice);
nAPP = numel(APPmice);

C = cell(4,1+nWT+nAPP);

C{1,1} = '';
C{2,1} = 'Condition';

for j = 1:nWT
    C{1,1+j} = 'WT';
    C{2,1+j} = char(WTmice(j));
end

for j = 1:nAPP
    col = 1+nWT+j;
    C{1,col} = 'APP';
    C{2,col} = char(APPmice(j));
end

C{3,1} = 'Baseline';
C{4,1} = 'Donepezil';

for j = 1:nWT

    q = D(D.mouse_id == WTmice(j),:);

    if height(q) == 1
        if isfinite(q.baseline_value)
            C{3,1+j} = double(q.baseline_value);
        end
        if isfinite(q.donepezil_value)
            C{4,1+j} = double(q.donepezil_value);
        end
    end
end

for j = 1:nAPP

    col = 1+nWT+j;
    q = D(D.mouse_id == APPmice(j),:);

    if height(q) == 1
        if isfinite(q.baseline_value)
            C{3,col} = double(q.baseline_value);
        end
        if isfinite(q.donepezil_value)
            C{4,col} = double(q.donepezil_value);
        end
    end
end

writecell(C,outFile);
end


function T = two_group_km_table(D)

WT = D(D.genotype == "WT",:);
APP = D(D.genotype == "APP",:);

[tWT,sWT] = km_xy(WT);
[tAPP,sAPP] = km_xy(APP);

n = max(numel(tWT),numel(tAPP));

WT_time_s = nan(n,1);
WT_survival = nan(n,1);
APP_time_s = nan(n,1);
APP_survival = nan(n,1);

WT_time_s(1:numel(tWT)) = tWT;
WT_survival(1:numel(sWT)) = sWT;
APP_time_s(1:numel(tAPP)) = tAPP;
APP_survival(1:numel(sAPP)) = sAPP;

T = table(WT_time_s,WT_survival,APP_time_s,APP_survival);
end


function T = four_group_km_table(D)

labels = ["WT_baseline","WT_donepezil","APP_baseline","APP_donepezil"];

times = cell(4,1);
survs = cell(4,1);

for i = 1:4

    if labels(i) == "WT_baseline"
        G = D(D.genotype == "WT" & D.condition == "baseline",:);
    elseif labels(i) == "WT_donepezil"
        G = D(D.genotype == "WT" & D.condition == "donepezil",:);
    elseif labels(i) == "APP_baseline"
        G = D(D.genotype == "APP" & D.condition == "baseline",:);
    else
        G = D(D.genotype == "APP" & D.condition == "donepezil",:);
    end

    [times{i},survs{i}] = km_xy(G);
end

n = max(cellfun(@numel,times));

vars = struct();

for i = 1:4
    tt = nan(n,1);
    ss = nan(n,1);

    tt(1:numel(times{i})) = times{i};
    ss(1:numel(survs{i})) = survs{i};

    vars.(labels(i)+"_time_s") = tt;
    vars.(labels(i)+"_survival") = ss;
end

T = struct2table(vars);
end


function [txy,sxy] = km_xy(D)

if isempty(D)
    txy = [];
    sxy = [];
    return;
end

w = mouse_equal_weights(D.mouse_id);

[t,s] = weighted_km(double(D.time_s),logical(D.event),w);

txy = [0;t];
sxy = [1;s];
end


function w = mouse_equal_weights(mouseID)

mouseID = string(mouseID);
mice = unique(mouseID);

n = numel(mouseID);
nMice = numel(mice);

w = zeros(n,1);

for i = 1:nMice
    idx = mouseID == mice(i);
    nMouse = sum(idx);

    if nMouse > 0
        w(idx) = 1/nMouse;
    end
end

if nMice > 0
    w = w*(n/nMice);
end
end


function [tEvent,S] = weighted_km(time,event,w)

time = double(time(:));
event = logical(event(:));
w = double(w(:));

valid = isfinite(time) & time>0 & isfinite(w) & w>0;

time = time(valid);
event = event(valid);
w = w(valid);

eventTimes = sort(unique(time(event)));

S = nan(numel(eventTimes),1);
surv = 1;

for i = 1:numel(eventTimes)

    tt = eventTimes(i);

    atRisk = time >= tt;
    eventsNow = time == tt & event;

    riskWeight = sum(w(atRisk));
    eventWeight = sum(w(eventsNow));

    if riskWeight > 0
        surv = surv*(1-eventWeight/riskWeight);
    end

    S(i) = surv;
end

tEvent = eventTimes;
end
