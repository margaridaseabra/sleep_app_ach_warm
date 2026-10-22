%% RUN_06_EXPORT_BASELINE_PRISM
%
% Export the unified baseline analysis into Prism-friendly CSV files.
%
% Run AFTER:
%   run('pipelines/run_baseline_complete.m')
%
% INPUT:
%   outputs/sleep_architecture/baseline_complete/
%
% OUTPUT (tracked, lightweight figure-source data):
%   results/prism/baseline/
%       mousewise/
%       survival/
%       stats/
%
% Mousewise files preserve one row per mouse and one numeric column per
% analysis window. NaNs are written as blank cells for easy Prism import.
%
% Survival files contain the exact mouse-balanced Kaplan-Meier step
% coordinates plotted by the MATLAB analysis. Import each file into an XY
% table in Prism using the WT x/y and APP x/y column pairs.

clear;
clc;

%% Setup

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));
P = paths_config();

analysisDir = fullfile(P.sleep_output,'baseline_complete');
survAnalysisDir = fullfile(analysisDir,'survival');

outDir = fullfile(repoDir,'results','prism','baseline');
mouseDir = fullfile(outDir,'mousewise');
survDir = fullfile(outDir,'survival');
statsDir = fullfile(outDir,'stats');

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(mouseDir), mkdir(mouseDir); end
if ~isfolder(survDir), mkdir(survDir); end
if ~isfolder(statsDir), mkdir(statsDir); end

windows = ["WHOLE","ZT3_9","ZT3_6","ZT6_9","ZT9_12"];
states = ["WK","NREM","REM","MA"];
carrierStates = ["NREM","REM"];
views = ["thesis","new","all"];

%% Load baseline-complete outputs

STATE = readtable( ...
    fullfile(analysisDir,'baseline_state_percent_per_mouse.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

BOUT = readtable( ...
    fullfile(analysisDir,'baseline_bout_metrics_per_mouse.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

FRAG = readtable( ...
    fullfile(analysisDir,'baseline_fragmentation_per_mouse.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

SURV = readtable( ...
    fullfile(survAnalysisDir,'baseline_survival_dataset.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

%% Normalize types

STATE.mouse_id = string(STATE.mouse_id);
STATE.genotype = string(STATE.genotype);
STATE.cohort = string(STATE.cohort);
STATE.sex = string(STATE.sex);
STATE.time_bin = string(STATE.time_bin);
STATE.state = string(STATE.state);

BOUT.mouse_id = string(BOUT.mouse_id);
BOUT.genotype = string(BOUT.genotype);
BOUT.cohort = string(BOUT.cohort);
BOUT.sex = string(BOUT.sex);
BOUT.time_bin = string(BOUT.time_bin);
BOUT.state = string(BOUT.state);

FRAG.mouse_id = string(FRAG.mouse_id);
FRAG.genotype = string(FRAG.genotype);
FRAG.cohort = string(FRAG.cohort);
FRAG.sex = string(FRAG.sex);
FRAG.time_bin = string(FRAG.time_bin);
FRAG.carrier_state = string(FRAG.carrier_state);

SURV.mouse_id = string(SURV.mouse_id);
SURV.genotype = string(SURV.genotype);
SURV.cohort = string(SURV.cohort);
SURV.time_bin = string(SURV.time_bin);
SURV.state = string(SURV.state);

%% Mousewise state percentages

for s = 1:numel(states)

    stateName = states(s);

    D = STATE(STATE.state == stateName,:);

    T = make_mousewise_table( ...
        D,windows,"percent_time");

    write_prism_table( ...
        T, ...
        fullfile(mouseDir, ...
        sprintf('state_percent_%s.csv',char(stateName))));
end

%% Mousewise bout metrics

boutMetrics = [ ...
    "mean_bout_duration_s", ...
    "median_bout_duration_s", ...
    "bouts_per_h"];

boutStems = [ ...
    "bout_mean_duration", ...
    "bout_median_duration", ...
    "bout_rate_per_h"];

for s = 1:numel(states)

    stateName = states(s);
    Dstate = BOUT(BOUT.state == stateName,:);

    for m = 1:numel(boutMetrics)

        T = make_mousewise_table( ...
            Dstate,windows,boutMetrics(m));

        write_prism_table( ...
            T, ...
            fullfile(mouseDir, ...
            sprintf('%s_%s.csv', ...
            char(boutStems(m)),char(stateName))));
    end
end

%% Mousewise fragmentation

for s = 1:numel(carrierStates)

    stateName = carrierStates(s);

    D = FRAG(FRAG.carrier_state == stateName,:);

    T = make_mousewise_table( ...
        D,windows,"MA_per_h_state");

    write_prism_table( ...
        T, ...
        fullfile(mouseDir, ...
        sprintf('fragmentation_%s.csv',char(stateName))));
end

%% Survival curve coordinates
%
% One CSV per cohort view x window x state.
% Columns:
%   WT_time_s, WT_survival, APP_time_s, APP_survival

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = SURV;
    else
        V = SURV(SURV.cohort == viewName,:);
    end

    for w = 1:numel(windows)

        windowName = windows(w);

        for stateName = ["NREM","REM"]

            D = V( ...
                V.time_bin == windowName & ...
                V.state == stateName,:);

            if isempty(D)
                continue;
            end

            WT = D(D.genotype == "WT",:);
            APP = D(D.genotype == "APP",:);

            [xWT,yWT] = survival_step_coordinates(WT);
            [xAPP,yAPP] = survival_step_coordinates(APP);

            n = max(numel(xWT),numel(xAPP));

            wtX = nan(n,1);
            wtY = nan(n,1);
            appX = nan(n,1);
            appY = nan(n,1);

            wtX(1:numel(xWT)) = xWT;
            wtY(1:numel(yWT)) = yWT;
            appX(1:numel(xAPP)) = xAPP;
            appY(1:numel(yAPP)) = yAPP;

            Tcurve = table( ...
                wtX,wtY,appX,appY, ...
                'VariableNames',{ ...
                    'WT_time_s','WT_survival', ...
                    'APP_time_s','APP_survival'});

            write_prism_table( ...
                Tcurve, ...
                fullfile(survDir, ...
                sprintf('survival_%s_%s_%s.csv', ...
                char(viewName),char(windowName),char(stateName))));
        end
    end
end

%% Copy analysis/statistics summaries into the Prism result package

copy_csv_if_exists( ...
    fullfile(analysisDir,'baseline_welch_stats.csv'), ...
    fullfile(statsDir,'baseline_welch_stats.csv'));

copy_csv_if_exists( ...
    fullfile(analysisDir,'baseline_genotype_cohort_models.csv'), ...
    fullfile(statsDir,'baseline_genotype_cohort_models.csv'));

copy_csv_if_exists( ...
    fullfile(survAnalysisDir,'baseline_survival_cox.csv'), ...
    fullfile(statsDir,'baseline_survival_cox.csv'));

copy_csv_if_exists( ...
    fullfile(survAnalysisDir,'baseline_survival_cox_interaction.csv'), ...
    fullfile(statsDir,'baseline_survival_cox_interaction.csv'));

copy_csv_if_exists( ...
    fullfile(survAnalysisDir,'baseline_survival_km_summary.csv'), ...
    fullfile(statsDir,'baseline_survival_km_summary.csv'));

copy_csv_if_exists( ...
    fullfile(survAnalysisDir,'baseline_survival_leave_one_mouse_out.csv'), ...
    fullfile(statsDir,'baseline_survival_leave_one_mouse_out.csv'));

%% Export a compact manifest

manifest = table( ...
    ["WHOLE";"ZT3_9";"ZT3_6";"ZT6_9";"ZT9_12"], ...
    ["Entire available baseline recording"; ...
     "Common continuous first 6 h (primary overall window)"; ...
     "First 3 h"; ...
     "Second 3 h"; ...
     "Third 3 h (exploratory when coverage is incomplete)"], ...
    ["sensitivity";"primary";"primary_time_resolved"; ...
     "primary_time_resolved";"exploratory"], ...
    'VariableNames',{'time_bin','description','recommended_role'});

write_prism_table( ...
    manifest, ...
    fullfile(outDir,'baseline_prism_manifest.csv'));

fprintf('\n============================================================\n');
fprintf('BASELINE PRISM EXPORT COMPLETE\n');
fprintf('============================================================\n');
fprintf('Mousewise data: %s\n',mouseDir);
fprintf('Survival curves: %s\n',survDir);
fprintf('Statistics: %s\n',statsDir);
fprintf('\nThese files are ready to import into GraphPad Prism.\n\n');

%% ========================================================================
% Local functions
% =========================================================================

function T = make_mousewise_table(D,windows,metricName)

mice = unique(D.mouse_id,'stable');

n = numel(mice);

mouse_id = strings(n,1);
genotype = strings(n,1);
cohort = strings(n,1);
cohort_label = strings(n,1);
sex = strings(n,1);

X = nan(n,numel(windows));

for i = 1:n

    mouseID = mice(i);
    M = D(D.mouse_id == mouseID,:);

    mouse_id(i) = mouseID;
    genotype(i) = M.genotype(1);
    cohort(i) = M.cohort(1);
    sex(i) = M.sex(1);

    if cohort(i) == "thesis"
        cohort_label(i) = "original";
    else
        cohort_label(i) = cohort(i);
    end

    for w = 1:numel(windows)

        R = M(M.time_bin == windows(w),:);

        if isempty(R)
            continue;
        end

        % Do not export a value to Prism when that mouse/window fails the
        % primary >=95% coverage criterion.
        if ~logical(R.include_primary(1))
            continue;
        end

        value = double(R.(metricName)(1));

        if isfinite(value)
            X(i,w) = value;
        end
    end
end

T = table( ...
    mouse_id,genotype,cohort_label,sex, ...
    'VariableNames',{'mouse_id','genotype','cohort','sex'});

for w = 1:numel(windows)
    T.(char(windows(w))) = X(:,w);
end

% Stable, human-friendly order.
cohortRank = zeros(height(T),1);
genoRank = zeros(height(T),1);

cohortRank(T.cohort == "original") = 1;
cohortRank(T.cohort == "new") = 2;

genoRank(T.genotype == "WT") = 1;
genoRank(T.genotype == "APP") = 2;

T.cohort_order = cohortRank;
T.genotype_order = genoRank;

T = sortrows(T,{'cohort_order','genotype_order','mouse_id'});
T(:,{'cohort_order','genotype_order'}) = [];
end


function [xStep,yStep] = survival_step_coordinates(D)

if isempty(D)
    xStep = [];
    yStep = [];
    return;
end

w = mouse_equal_weights(D.mouse_id);

[tEvent,S] = weighted_km( ...
    double(D.time_s), ...
    logical(D.event), ...
    w);

if isempty(tEvent)
    xStep = 0;
    yStep = 1;
    return;
end

% Convert stairs() representation into plain XY coordinates so Prism can
% reproduce the same step curve with ordinary line connections.

n = numel(tEvent);

xStep = nan(2*n+1,1);
yStep = nan(2*n+1,1);

xStep(1) = 0;
yStep(1) = 1;

prevS = 1;

for i = 1:n

    j = 2*i;

    xStep(j) = tEvent(i);
    yStep(j) = prevS;

    xStep(j+1) = tEvent(i);
    yStep(j+1) = S(i);

    prevS = S(i);
end
end


function w = mouse_equal_weights(mouseID)

mouseID = string(mouseID);
mice = unique(mouseID);

n = numel(mouseID);
nMice = numel(mice);

w = zeros(n,1);

for i = 1:nMice
    idx = mouseID == mice(i);
    nBoutMouse = sum(idx);

    if nBoutMouse > 0
        w(idx) = 1/nBoutMouse;
    end
end

if nMice > 0
    w = w * (n/nMice);
end
end


function [tEvent,S] = weighted_km(time,event,w)

time = double(time(:));
event = logical(event(:));
w = double(w(:));

valid = ...
    isfinite(time) & ...
    time > 0 & ...
    isfinite(w) & ...
    w > 0;

time = time(valid);
event = event(valid);
w = w(valid);

eventTimes = unique(time(event));
eventTimes = sort(eventTimes);

S = nan(numel(eventTimes),1);
surv = 1;

for i = 1:numel(eventTimes)

    tt = eventTimes(i);

    atRisk = time >= tt;
    eventsNow = time == tt & event;

    riskWeight = sum(w(atRisk));
    eventWeight = sum(w(eventsNow));

    if riskWeight > 0
        surv = surv * (1-eventWeight/riskWeight);
    end

    S(i) = surv;
end

tEvent = eventTimes;
end


function write_prism_table(T,filePath)

% Write NaNs as blank cells instead of literal "NaN", which imports more
% cleanly into GraphPad Prism.

headers = T.Properties.VariableNames;
body = table2cell(T);

for r = 1:size(body,1)
    for c = 1:size(body,2)
        x = body{r,c};

        if isnumeric(x) && isscalar(x) && isnan(x)
            body{r,c} = '';
        elseif ismissing_scalar(x)
            body{r,c} = '';
        end
    end
end

C = [headers; body];
writecell(C,filePath);
end


function tf = ismissing_scalar(x)

tf = false;

try
    if isstring(x) && isscalar(x)
        tf = ismissing(x);
    elseif ischar(x)
        tf = false;
    end
catch
    tf = false;
end
end


function copy_csv_if_exists(src,dst)

if isfile(src)
    T = readtable(src,'TextType','string','VariableNamingRule','preserve');
    write_prism_table(T,dst);
else
    warning('Expected statistics file not found: %s',src);
end
end
