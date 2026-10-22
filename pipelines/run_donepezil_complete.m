%% RUN_DONEPEZIL_COMPLETE
%
% Unified donepezil sleep analysis for paired baseline/donepezil sessions.
%
% Requires standardized 1-s score CSVs from run_01_sleep_architecture.m.
%
% Experimental windows are aligned to each mouse's ACTUAL donepezil injection:
%   PRE        = up to 3 h immediately before injection
%   DRUG_0_3H  = injection to +3 h (primary acute drug response)
%   DRUG_3_6H  = +3 h to +6 h (secondary late response)
%
% Exact injection times are read from config/drug_times.csv. The same
% mouse-specific elapsed-time boundaries are applied to the paired baseline
% session, which controls the matched baseline comparison for time of day.
%
% Main analyses:
%   1) Raw donepezil-session metrics across the three phases.
%   2) Within-session delta from PRE.
%   3) Matched separate baseline vs donepezil at the same elapsed time.
%   4) Difference-in-differences (DiD):
%        (drug_target - drug_PRE) - (baseline_target - baseline_PRE)
%   5) NREM/REM bout survival:
%        - donepezil-only genotype comparisons
%        - within-donepezil-session PRE vs DRUG_0_3H genotype x phase
%        - matched baseline-vs-donepezil genotype x condition
%        - leave-one-mouse-out sensitivity for the matched interaction
%
% Shared conventions with run_baseline_complete.m:
%   - Wake bouts <=15 s -> micro-arousals (MA)
%   - windows require >=95% valid scored coverage
%   - bout number = reliable observed bout onsets
%   - bout duration = fully observed bouts contained within window
%   - fragmentation = MA onsets from NREM/REM per hour in carrier state
%   - survival is mouse-balanced for KM and uses mouse-clustered robust Cox SEs
%
% Saline:
%   The main analysis is aligned to donepezil injection. Saline timing is
%   intentionally NOT imputed where it was not recorded (October cohort).
%   A June-only saline sensitivity analysis can be added separately.
%
% Outputs:
%   outputs/sleep_architecture/donepezil_complete/
%
% NOTE:
%   "thesis" is the historical/original cohort internally.
%
% -------------------------------------------------------------------------

clear;
clc;

%% Setup

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));

P = paths_config();

Sessions = readtable(P.sessions, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

Mice = readtable(P.mice, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

standardDir = fullfile(P.sleep_output,'all','standardized_scores');

outDir = fullfile(P.sleep_output,'donepezil_complete');
figDir = fullfile(outDir,'figures');
survDir = fullfile(outDir,'survival');

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(survDir), mkdir(survDir); end

%% Definitions

conditions = ["baseline","donepezil"];

windowNames = ["PRE","DRUG_0_3H","DRUG_3_6H"];

% Default is only a safety fallback. All active drug sessions currently
% have a recorded injection time in config/drug_times.csv.
defaultDrugStartSec = 3*3600;
drugWindowSec = 3*3600;

if isfield(P,'drug_times')
    timingFile = P.drug_times;
else
    timingFile = fullfile(repoDir,'config','drug_times.csv');
end

if isfile(timingFile)
    DrugTimes = readtable( ...
        timingFile, ...
        'TextType','string', ...
        'VariableNamingRule','preserve');
else
    DrugTimes = table();
    warning(['No config/drug_times.csv found. ', ...
             'Using 03:00:00 injection time as a fallback for all mice.']);
end

stateNames = ["WK","NREM","REM","MA"];
stateCodes = [0,1,2,15];

carrierStates = ["NREM","REM"];
survivalStates = ["NREM","REM"];

views = ["thesis","new","all"];

minimumValidCoverage = 95;
maThresholdSec = 15;

fprintf('\n============================================================\n');
fprintf('DONEPEZIL ANALYSIS\n');
fprintf('============================================================\n\n');

%% Select paired baseline/donepezil mice

S = Sessions( ...
    ismember(Sessions.condition,conditions) & ...
    logical(Sessions.include_sleep),:);

allMice = unique(S.mouse_id);

pairedMice = strings(0,1);

for i = 1:numel(allMice)
    m = allMice(i);
    conds = unique(S.condition(S.mouse_id == m));

    if all(ismember(conditions,conds))
        pairedMice(end+1,1) = m; %#ok<SAGROW>
    end
end

S = S(ismember(S.mouse_id,pairedMice),:);

fprintf('Paired baseline + donepezil mice: %d\n',numel(pairedMice));

for g = ["WT","APP"]
    ids = Mice.mouse_id(Mice.genotype == g & ismember(Mice.mouse_id,pairedMice));
    fprintf('  %-3s: %d\n',g,numel(ids));
end
fprintf('\n');

%% Collectors

RAW = table();
ALLBOUTS = table();

%% Process baseline + donepezil sessions

for i = 1:height(S)

    sid = string(S.session_id(i));
    mouseID = string(S.mouse_id(i));
    condition = string(S.condition(i));
    cohort = string(S.cohort(i));

    mouseIdx = find(Mice.mouse_id == mouseID,1);

    if isempty(mouseIdx)
        warning('Mouse %s missing from mice.csv. Skipping.',mouseID);
        continue;
    end

    genotype = string(Mice.genotype(mouseIdx));
    sex = string(Mice.sex(mouseIdx));

    scoreFile = fullfile(standardDir,char(sid + ".csv"));

    if ~isfile(scoreFile)
        warning('Missing standardized score file: %s',scoreFile);
        continue;
    end

    T = readtable(scoreFile);

    t = double(T.time_s(:));
    score = double(T.score(:));

    if isempty(t)
        warning('Empty score file: %s',scoreFile);
        continue;
    end

    t = t - min(t);

    scoreMA = reclassify_short_wake(score,maThresholdSec);

    B = build_bout_table(scoreMA,t);

    if ~isempty(B)
        B.session_id = repmat(sid,height(B),1);
        B.mouse_id = repmat(mouseID,height(B),1);
        B.genotype = repmat(genotype,height(B),1);
        B.sex = repmat(sex,height(B),1);
        B.cohort = repmat(cohort,height(B),1);
        B.condition = repmat(condition,height(B),1);
        ALLBOUTS = [ALLBOUTS; B]; %#ok<AGROW>
    end

    [thisWindowStart,thisWindowEnd,timingStatus] = get_mouse_windows( ...
        mouseID,DrugTimes,defaultDrugStartSec,drugWindowSec);

    fprintf('[%02d/%02d] %-8s %-10s %-4s | drug %.2f h | windows %.2f-%.2f-%.2f-%.2f h (%s)\n', ...
        i,height(S),mouseID,condition,genotype, ...
        thisWindowStart(2)/3600, ...
        thisWindowStart(1)/3600,thisWindowStart(2)/3600, ...
        thisWindowEnd(2)/3600,thisWindowEnd(3)/3600,timingStatus);

    for w = 1:numel(windowNames)

        windowName = windowNames(w);
        wStart = thisWindowStart(w);
        wEnd = thisWindowEnd(w);

        inWindow = t >= wStart & t < wEnd;

        expectedSec = wEnd-wStart;
        nAvailable = sum(inWindow);

        validMask = inWindow & ismember(scoreMA,[0 1 2 15]);
        validSec = sum(validMask);

        recordingCoverage = 100*nAvailable/expectedSec;
        validCoverage = 100*validSec/expectedSec;
        includePrimary = validCoverage >= minimumValidCoverage;
        validHours = validSec/3600;

        %% State occupancy

        for s = 1:numel(stateNames)

            stateName = stateNames(s);
            stateCode = stateCodes(s);

            stateSec = sum(inWindow & scoreMA == stateCode);

            if validSec > 0
                pct = 100*stateSec/validSec;
            else
                pct = NaN;
            end

            RAW = [RAW; metric_row( ...
                sid,mouseID,genotype,sex,cohort,condition, ...
                windowName,wStart,wEnd,timingStatus, ...
                "state_percent",stateName,"percent_time",pct, ...
                validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>
        end

        %% Bout metrics

        for s = 1:numel(stateNames)

            stateName = stateNames(s);
            stateCode = stateCodes(s);

            if isempty(B)
                BR = B;
            else
                BR = B( ...
                    B.state_code == stateCode & ...
                    B.onset_s >= wStart & ...
                    B.onset_s < wEnd,:);
            end

            if isempty(BR)
                nBouts = 0;
                boutsPerHour = 0;
                meanDuration = NaN;
                medianDuration = NaN;
            else
                reliableStart = ~BR.left_censored;
                nBouts = sum(reliableStart);

                if validHours > 0
                    boutsPerHour = nBouts/validHours;
                else
                    boutsPerHour = NaN;
                end

                useDuration = ...
                    ~BR.left_censored & ...
                    ~BR.right_censored & ...
                    BR.end_exclusive_s <= wEnd;

                durations = double(BR.duration_s(useDuration));

                if isempty(durations)
                    meanDuration = NaN;
                    medianDuration = NaN;
                else
                    meanDuration = mean(durations);
                    medianDuration = median(durations);
                end
            end

            RAW = [RAW; metric_row( ...
                sid,mouseID,genotype,sex,cohort,condition, ...
                windowName,wStart,wEnd,timingStatus, ...
                "bout_metric",stateName,"bouts_per_h",boutsPerHour, ...
                validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>

            RAW = [RAW; metric_row( ...
                sid,mouseID,genotype,sex,cohort,condition, ...
                windowName,wStart,wEnd,timingStatus, ...
                "bout_metric",stateName,"mean_bout_duration_s",meanDuration, ...
                validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>

            RAW = [RAW; metric_row( ...
                sid,mouseID,genotype,sex,cohort,condition, ...
                windowName,wStart,wEnd,timingStatus, ...
                "bout_metric",stateName,"median_bout_duration_s",medianDuration, ...
                validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>
        end

        %% Fragmentation

        if isempty(B)
            maNREM = 0;
            maREM = 0;
        else
            MA = B( ...
                B.state_code == 15 & ...
                B.onset_s >= wStart & ...
                B.onset_s < wEnd & ...
                ~B.left_censored,:);

            maNREM = sum(MA.state_before == 1);
            maREM = sum(MA.state_before == 2);
        end

        nremSec = sum(inWindow & scoreMA == 1);
        remSec = sum(inWindow & scoreMA == 2);

        if nremSec > 0
            nremFrag = maNREM/(nremSec/3600);
        else
            nremFrag = NaN;
        end

        if remSec > 0
            remFrag = maREM/(remSec/3600);
        else
            remFrag = NaN;
        end

        RAW = [RAW; metric_row( ...
            sid,mouseID,genotype,sex,cohort,condition, ...
            windowName,wStart,wEnd,timingStatus, ...
            "fragmentation","NREM","MA_per_h_state",nremFrag, ...
            validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>

        RAW = [RAW; metric_row( ...
            sid,mouseID,genotype,sex,cohort,condition, ...
            windowName,wStart,wEnd,timingStatus, ...
            "fragmentation","REM","MA_per_h_state",remFrag, ...
            validCoverage,recordingCoverage,includePrimary)]; %#ok<AGROW>
    end
end

writetable(RAW,fullfile(outDir,'donepezil_raw_metrics_long.csv'));
writetable(ALLBOUTS,fullfile(outDir,'donepezil_all_bouts.csv'));

%% QC summary

QC = RAW( ...
    RAW.domain == "state_percent" & ...
    RAW.state == "REM" & ...
    RAW.metric == "percent_time",:);

fprintf('\nUsable sessions by condition/window:\n');

for c = conditions
    for w = windowNames
        D = QC(QC.condition == c & QC.time_bin == w,:);
        fprintf('  %-8s %-8s : %d/%d\n', ...
            c,w,sum(D.include_primary),height(D));
    end
end

%% Primary valid rows

R = RAW(logical(RAW.include_primary) & isfinite(RAW.value),:);

%% Build donepezil-session delta-from-PRE table

DELTA = table();

donepezilMetrics = R(R.condition == "donepezil",:);

% IMPORTANT: use only stable, non-missing identifiers as unique keys.
% Including sex here caused repeated keys for the historical cohort because
% sex is unknown/missing for those mice. Metadata are recovered from D below.
keys = unique( ...
    donepezilMetrics(:,["mouse_id","domain","state","metric"]), ...
    'rows','stable');

for i = 1:height(keys)

    k = keys(i,:);

    D = donepezilMetrics( ...
        donepezilMetrics.mouse_id == k.mouse_id & ...
        donepezilMetrics.domain == k.domain & ...
        donepezilMetrics.state == k.state & ...
        donepezilMetrics.metric == k.metric,:);

    pre = D.value(D.time_bin == "PRE");

    if numel(pre) ~= 1
        continue;
    end

    for w = 1:numel(windowNames)

        target = D.value(D.time_bin == windowNames(w));

        if numel(target) ~= 1
            continue;
        end

        deltaValue = target-pre;

        genotype = string(D.genotype(1));
        sex = string(D.sex(1));
        cohort = string(D.cohort(1));

        row = table( ...
            k.mouse_id,genotype,sex,cohort, ...
            k.domain,k.state,k.metric, ...
            windowNames(w), ...
            pre,target,deltaValue, ...
            'VariableNames',{ ...
                'mouse_id','genotype','sex','cohort', ...
                'domain','state','metric','time_bin', ...
                'pre_value','value','delta_from_pre'});

        DELTA = [DELTA; row]; %#ok<AGROW>
    end
end

assert_unique_mouse_metric_rows( ...
    DELTA,["mouse_id","domain","state","metric","time_bin"], ...
    'DELTA');

writetable(DELTA,fullfile(outDir,'donepezil_delta_from_pre.csv'));

%% Matched baseline-vs-donepezil at same elapsed time

MATCHED = table();

metricKeys = unique( ...
    R(:,["mouse_id","domain","state","metric","time_bin"]), ...
    'rows','stable');

for i = 1:height(metricKeys)

    k = metricKeys(i,:);

    D = R( ...
        R.mouse_id == k.mouse_id & ...
        R.domain == k.domain & ...
        R.state == k.state & ...
        R.metric == k.metric & ...
        R.time_bin == k.time_bin,:);

    b = D.value(D.condition == "baseline");
    a = D.value(D.condition == "donepezil");

    if numel(b) ~= 1 || numel(a) ~= 1
        continue;
    end

    genotype = string(D.genotype(1));
    sex = string(D.sex(1));
    cohort = string(D.cohort(1));

    row = table( ...
        k.mouse_id,genotype,sex,cohort, ...
        k.domain,k.state,k.metric,k.time_bin, ...
        b,a,a-b, ...
        'VariableNames',{ ...
            'mouse_id','genotype','sex','cohort', ...
            'domain','state','metric','time_bin', ...
            'baseline_value','donepezil_value','donepezil_minus_baseline'});

    MATCHED = [MATCHED; row]; %#ok<AGROW>
end

assert_unique_mouse_metric_rows( ...
    MATCHED,["mouse_id","domain","state","metric","time_bin"], ...
    'MATCHED');

writetable(MATCHED,fullfile(outDir,'donepezil_matched_baseline_vs_donepezil.csv'));

%% Difference-in-differences (matched baseline day control)

DID = table();

targets = ["DRUG_0_3H","DRUG_3_6H"];

metricKeys2 = unique( ...
    R(:,["mouse_id","domain","state","metric"]), ...
    'rows','stable');

for i = 1:height(metricKeys2)

    k = metricKeys2(i,:);

    D = R( ...
        R.mouse_id == k.mouse_id & ...
        R.domain == k.domain & ...
        R.state == k.state & ...
        R.metric == k.metric,:);

    bPre = D.value(D.condition == "baseline" & D.time_bin == "PRE");
    wPre = D.value(D.condition == "donepezil"  & D.time_bin == "PRE");

    if numel(bPre) ~= 1 || numel(wPre) ~= 1
        continue;
    end

    for tIdx = 1:numel(targets)

        targetBin = targets(tIdx);

        bTarget = D.value(D.condition == "baseline" & D.time_bin == targetBin);
        wTarget = D.value(D.condition == "donepezil" & D.time_bin == targetBin);

        if numel(bTarget) ~= 1 || numel(wTarget) ~= 1
            continue;
        end

        baselineChange = bTarget-bPre;
        donepezilChange = wTarget-wPre;
        didValue = donepezilChange-baselineChange;

        genotype = string(D.genotype(1));
        sex = string(D.sex(1));
        cohort = string(D.cohort(1));

        row = table( ...
            k.mouse_id,genotype,sex,cohort, ...
            k.domain,k.state,k.metric,targetBin, ...
            bPre,bTarget,wPre,wTarget, ...
            baselineChange,donepezilChange,didValue, ...
            'VariableNames',{ ...
                'mouse_id','genotype','sex','cohort', ...
                'domain','state','metric','time_bin', ...
                'baseline_pre','baseline_target', ...
                'donepezil_pre','donepezil_target', ...
                'baseline_change','donepezil_change','difference_in_differences'});

        DID = [DID; row]; %#ok<AGROW>
    end
end

assert_unique_mouse_metric_rows( ...
    DID,["mouse_id","domain","state","metric","time_bin"], ...
    'DID');

writetable(DID,fullfile(outDir,'donepezil_difference_in_differences.csv'));

%% Delta genotype comparisons

DELTA_STATS = table();

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DELTA;
    else
        V = DELTA(DELTA.cohort == viewName,:);
    end

    combinations = unique(V(:,["domain","state","metric"]), 'rows');

    for j = 1:height(combinations)

        q = combinations(j,:);

        for targetBin = targets

            D = V( ...
                V.domain == q.domain & ...
                V.state == q.state & ...
                V.metric == q.metric & ...
                V.time_bin == targetBin,:);

            WT = double(D.delta_from_pre(D.genotype == "WT"));
            APP = double(D.delta_from_pre(D.genotype == "APP"));

            [p,tstat,df,d] = welch_stats(WT,APP);

            row = table( ...
                viewName,targetBin,q.domain,q.state,q.metric, ...
                numel(WT),numel(APP), ...
                mean(WT,'omitnan'),mean(APP,'omitnan'), ...
                mean(APP,'omitnan')-mean(WT,'omitnan'), ...
                p,tstat,df,d, ...
                'VariableNames',{ ...
                    'cohort_view','time_bin','domain','state','metric', ...
                    'n_WT','n_APP','WT_mean_delta','APP_mean_delta', ...
                    'APP_minus_WT_delta','p_welch','t_stat','df','Cohen_d'});

            DELTA_STATS = [DELTA_STATS; row]; %#ok<AGROW>
        end
    end
end

DELTA_STATS = add_sidak_two_phase(DELTA_STATS);
writetable(DELTA_STATS,fullfile(outDir,'donepezil_delta_genotype_stats.csv'));

%% Difference-in-differences (matched baseline day control) genotype comparisons

DID_STATS = table();

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DID;
    else
        V = DID(DID.cohort == viewName,:);
    end

    combinations = unique(V(:,["domain","state","metric"]), 'rows');

    for j = 1:height(combinations)

        q = combinations(j,:);

        for targetBin = targets

            D = V( ...
                V.domain == q.domain & ...
                V.state == q.state & ...
                V.metric == q.metric & ...
                V.time_bin == targetBin,:);

            WT = double(D.difference_in_differences(D.genotype == "WT"));
            APP = double(D.difference_in_differences(D.genotype == "APP"));

            [p,tstat,df,d] = welch_stats(WT,APP);

            row = table( ...
                viewName,targetBin,q.domain,q.state,q.metric, ...
                numel(WT),numel(APP), ...
                mean(WT,'omitnan'),mean(APP,'omitnan'), ...
                mean(APP,'omitnan')-mean(WT,'omitnan'), ...
                p,tstat,df,d, ...
                'VariableNames',{ ...
                    'cohort_view','time_bin','domain','state','metric', ...
                    'n_WT','n_APP','WT_mean_DiD','APP_mean_DiD', ...
                    'APP_minus_WT_DiD','p_welch','t_stat','df','Cohen_d'});

            DID_STATS = [DID_STATS; row]; %#ok<AGROW>
        end
    end
end

DID_STATS = add_sidak_two_phase(DID_STATS);
writetable(DID_STATS,fullfile(outDir,'donepezil_DiD_genotype_stats.csv'));

%% Mixed-effects models: donepezil session across PRE / acute / late windows

LME_TIME = table();

metricCombos = unique(R(:,["domain","state","metric"]), 'rows');

for v = 1:numel(views)

    viewName = views(v);

    W = R(R.condition == "donepezil",:);

    if viewName ~= "all"
        W = W(W.cohort == viewName,:);
    end

    for j = 1:height(metricCombos)

        q = metricCombos(j,:);

        D = W( ...
            W.domain == q.domain & ...
            W.state == q.state & ...
            W.metric == q.metric,:);

        if height(D) < 6
            continue;
        end

        Tm = table( ...
            double(D.value), ...
            categorical(D.genotype), ...
            categorical(D.time_bin,windowNames,{'Pre','Warm','Cool'}), ...
            categorical(D.cohort), ...
            categorical(D.mouse_id), ...
            'VariableNames',{'Value','Genotype','TimeBin','Cohort','Mouse'});

        if viewName == "all"
            formula = 'Value ~ Genotype*TimeBin + Cohort + (1|Mouse)';
        else
            formula = 'Value ~ Genotype*TimeBin + (1|Mouse)';
        end

        A = safe_lme_anova(Tm,formula);

        if ~isempty(A)
            A.cohort_view = repmat(viewName,height(A),1);
            A.domain = repmat(q.domain,height(A),1);
            A.state = repmat(q.state,height(A),1);
            A.metric = repmat(q.metric,height(A),1);

            LME_TIME = [LME_TIME; ...
                movevars(A,["cohort_view","domain","state","metric"],'Before',1)]; %#ok<AGROW>
        end
    end
end

writetable(LME_TIME,fullfile(outDir,'donepezil_LME_time_within_session.csv'));

%% Mixed-effects models: baseline vs donepezil within each matched elapsed-time window

LME_CONDITION = table();

for v = 1:numel(views)

    viewName = views(v);

    V = R;

    if viewName ~= "all"
        V = V(V.cohort == viewName,:);
    end

    for w = 1:numel(windowNames)

        for j = 1:height(metricCombos)

            q = metricCombos(j,:);

            D = V( ...
                V.time_bin == windowNames(w) & ...
                V.domain == q.domain & ...
                V.state == q.state & ...
                V.metric == q.metric,:);

            % Keep only mice represented in both conditions for this metric/bin.
            paired = paired_ids_for_metric(D);
            D = D(ismember(D.mouse_id,paired),:);

            if height(D) < 6 || numel(unique(D.mouse_id)) < 3
                continue;
            end

            Tm = table( ...
                double(D.value), ...
                categorical(D.genotype), ...
                categorical(D.condition,["baseline","donepezil"],{'Baseline','Donepezil'}), ...
                categorical(D.cohort), ...
                categorical(D.mouse_id), ...
                'VariableNames',{'Value','Genotype','Condition','Cohort','Mouse'});

            if viewName == "all"
                formula = 'Value ~ Genotype*Condition + Cohort + (1|Mouse)';
            else
                formula = 'Value ~ Genotype*Condition + (1|Mouse)';
            end

            A = safe_lme_anova(Tm,formula);

            if ~isempty(A)
                A.cohort_view = repmat(viewName,height(A),1);
                A.time_bin = repmat(windowNames(w),height(A),1);
                A.domain = repmat(q.domain,height(A),1);
                A.state = repmat(q.state,height(A),1);
                A.metric = repmat(q.metric,height(A),1);

                LME_CONDITION = [LME_CONDITION; ...
                    movevars(A,["cohort_view","time_bin","domain","state","metric"],'Before',1)]; %#ok<AGROW>
            end
        end
    end
end

writetable(LME_CONDITION,fullfile(outDir,'donepezil_LME_baseline_vs_donepezil.csv'));

%% Build survival dataset for baseline + donepezil sessions

SURV = table();

for c = conditions

    sessionsC = unique(RAW.session_id(RAW.condition == c));

    for w = 1:numel(windowNames)

        binName = windowNames(w);

        QCbin = RAW( ...
            RAW.condition == c & ...
            RAW.time_bin == binName & ...
            RAW.domain == "state_percent" & ...
            RAW.state == "NREM" & ...
            RAW.metric == "percent_time" & ...
            logical(RAW.include_primary),:);

        primarySessions = unique(QCbin.session_id);

        for st = survivalStates

            for q = 1:numel(primarySessions)

                sid = primarySessions(q);

                QCs = QCbin(QCbin.session_id == sid,:);

                if isempty(QCs)
                    continue;
                end

                bStart = double(QCs.window_start_s(1));
                bEnd = double(QCs.window_end_s(1));

                D = ALLBOUTS( ...
                    ALLBOUTS.session_id == sid & ...
                    ALLBOUTS.state == st & ...
                    ALLBOUTS.onset_s >= bStart & ...
                    ALLBOUTS.onset_s < bEnd,:);

                if isempty(D)
                    continue;
                end

                D = D(~logical(D.left_censored),:);

                if isempty(D)
                    continue;
                end

                remaining = bEnd-double(D.onset_s);

                observedTime = min(double(D.duration_s),remaining);

                event = ...
                    ~logical(D.right_censored) & ...
                    double(D.end_exclusive_s) <= bEnd;

                censored = ~event;

                valid = isfinite(observedTime) & observedTime > 0;

                D = D(valid,:);
                observedTime = observedTime(valid);
                event = event(valid);
                censored = censored(valid);

                if isempty(D)
                    continue;
                end

                row = table( ...
                    D.session_id,D.mouse_id,D.genotype,D.cohort,D.condition, ...
                    repmat(binName,height(D),1), ...
                    repmat(st,height(D),1), ...
                    observedTime,event,censored, ...
                    'VariableNames',{ ...
                        'session_id','mouse_id','genotype','cohort','condition', ...
                        'time_bin','state','time_s','event','censored'});

                SURV = [SURV; row]; %#ok<AGROW>
            end
        end
    end
end

writetable(SURV,fullfile(survDir,'donepezil_survival_dataset.csv'));

%% Donepezil-only genotype Cox models + KM medians

COX_DONEPEZIL = table();
KM_DONEPEZIL = table();

DRUG_SURV = SURV(SURV.condition == "donepezil",:);

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DRUG_SURV;
    else
        V = DRUG_SURV(DRUG_SURV.cohort == viewName,:);
    end

    for w = 1:numel(windowNames)
        for st = survivalStates

            D = V(V.time_bin == windowNames(w) & V.state == st,:);

            if isempty(D) || numel(unique(D.genotype)) < 2
                continue;
            end

            for g = ["WT","APP"]

                G = D(D.genotype == g,:);

                if isempty(G), continue; end

                weights = mouse_equal_weights(G.mouse_id);

                [~,~,medianSurvival] = weighted_km( ...
                    double(G.time_s),logical(G.event),weights);

                KM_DONEPEZIL = [KM_DONEPEZIL; table( ...
                    viewName,windowNames(w),st,g, ...
                    numel(unique(G.mouse_id)),height(G),sum(G.event), ...
                    medianSurvival, ...
                    'VariableNames',{ ...
                        'cohort_view','time_bin','state','genotype', ...
                        'n_mice','n_bouts','n_events','KM_median_s'})]; %#ok<AGROW>
            end

            genotype = double(D.genotype == "APP");

            if viewName == "all"
                cohort = double(D.cohort == "new");
                FIT = clustered_cox( ...
                    [genotype cohort],double(D.time_s),logical(D.censored), ...
                    D.mouse_id,["genotype_APP","cohort_new"]);
            else
                FIT = clustered_cox( ...
                    genotype,double(D.time_s),logical(D.censored), ...
                    D.mouse_id,"genotype_APP");
            end

            if FIT.success

                Rfit = FIT.results(FIT.results.predictor == "genotype_APP",:);

                WT = D(D.genotype == "WT",:);
                APP = D(D.genotype == "APP",:);

                row = table( ...
                    viewName,windowNames(w),st, ...
                    numel(unique(WT.mouse_id)),numel(unique(APP.mouse_id)), ...
                    height(WT),height(APP),sum(WT.event),sum(APP.event), ...
                    Rfit.HR,Rfit.CI_low,Rfit.CI_high,Rfit.p_robust_t, ...
                    'VariableNames',{ ...
                        'cohort_view','time_bin','state', ...
                        'n_mice_WT','n_mice_APP','n_bouts_WT','n_bouts_APP', ...
                        'n_events_WT','n_events_APP', ...
                        'HR_APP_vs_WT','CI_low','CI_high','p_robust_t'});

                COX_DONEPEZIL = [COX_DONEPEZIL; row]; %#ok<AGROW>
            end
        end
    end
end

writetable(COX_DONEPEZIL,fullfile(survDir,'donepezil_only_survival_cox.csv'));
writetable(KM_DONEPEZIL,fullfile(survDir,'donepezil_only_survival_km_summary.csv'));

%% Within-donepezil-session PRE vs acute-drug Cox interaction
%
% This is the preferred survival test of an acute drug response because
% PRE and DRUG_0_3H come from the same recording day.
%
% HR_phase < 1 means lower bout-termination hazard during drug than PRE.
% HR_genotype_x_phase tests whether that change differs between APP and WT.

COX_WITHIN_SESSION = table();
DRUGSURV = SURV(SURV.condition == "donepezil" & ...
    ismember(SURV.time_bin,["PRE","DRUG_0_3H"]),:);

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DRUGSURV;
    else
        V = DRUGSURV(DRUGSURV.cohort == viewName,:);
    end

    for st = survivalStates

        D = V(V.state == st,:);

        paired = paired_ids_phase_survival(D,"PRE","DRUG_0_3H");
        D = D(ismember(D.mouse_id,paired),:);

        if isempty(D) || ...
                numel(unique(D.genotype)) < 2 || ...
                numel(unique(D.time_bin)) < 2
            continue;
        end

        genotype = double(D.genotype == "APP");
        phase = double(D.time_bin == "DRUG_0_3H");
        interaction = genotype.*phase;

        if viewName == "all"
            cohort = double(D.cohort == "new");
            X = [genotype phase interaction cohort];
            names = ["genotype_APP","drug_phase","genotype_x_drug_phase","cohort_new"];
        else
            X = [genotype phase interaction];
            names = ["genotype_APP","drug_phase","genotype_x_drug_phase"];
        end

        FIT = clustered_cox( ...
            X,double(D.time_s),logical(D.censored),D.mouse_id,names);

        if ~FIT.success
            continue;
        end

        Rfit = FIT.results;

        iG = find(Rfit.predictor == "genotype_APP");
        iP = find(Rfit.predictor == "drug_phase");
        iI = find(Rfit.predictor == "genotype_x_drug_phase");

        gPre = linear_combo_hr(FIT,iG,1);
        gDrug = linear_combo_hr(FIT,[iG iI],[1 1]);
        phaseWT = linear_combo_hr(FIT,iP,1);
        phaseAPP = linear_combo_hr(FIT,[iP iI],[1 1]);
        intHR = linear_combo_hr(FIT,iI,1);

        row = table( ...
            viewName,st,numel(unique(D.mouse_id)),height(D),sum(D.event), ...
            gPre.HR,gPre.CI_low,gPre.CI_high,gPre.p, ...
            gDrug.HR,gDrug.CI_low,gDrug.CI_high,gDrug.p, ...
            phaseWT.HR,phaseWT.CI_low,phaseWT.CI_high,phaseWT.p, ...
            phaseAPP.HR,phaseAPP.CI_low,phaseAPP.CI_high,phaseAPP.p, ...
            intHR.HR,intHR.CI_low,intHR.CI_high,intHR.p, ...
            'VariableNames',{ ...
                'cohort_view','state','n_mice','n_bouts','n_events', ...
                'HR_APP_vs_WT_PRE','CI_gpre_low','CI_gpre_high','p_gpre', ...
                'HR_APP_vs_WT_DRUG_0_3H','CI_gdrug_low','CI_gdrug_high','p_gdrug', ...
                'HR_DRUG_vs_PRE_WT','CI_wt_low','CI_wt_high','p_wt', ...
                'HR_DRUG_vs_PRE_APP','CI_app_low','CI_app_high','p_app', ...
                'HR_genotype_x_drug_phase','CI_int_low','CI_int_high','p_interaction'});

        COX_WITHIN_SESSION = [COX_WITHIN_SESSION; row]; %#ok<AGROW>
    end
end

writetable(COX_WITHIN_SESSION, ...
    fullfile(survDir,'donepezil_within_session_pre_vs_drug_cox_interaction.csv'));

%% Baseline-vs-donepezil Cox interaction

COX_INTERACTION = table();

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = SURV;
    else
        V = SURV(SURV.cohort == viewName,:);
    end

    for w = 1:numel(windowNames)

        for st = survivalStates

            D = V(V.time_bin == windowNames(w) & V.state == st,:);

            paired = paired_ids_survival(D);
            D = D(ismember(D.mouse_id,paired),:);

            if isempty(D) || ...
                    numel(unique(D.genotype)) < 2 || ...
                    numel(unique(D.condition)) < 2
                continue;
            end

            genotype = double(D.genotype == "APP");
            donepezil = double(D.condition == "donepezil");
            interaction = genotype.*donepezil;

            if viewName == "all"
                cohort = double(D.cohort == "new");
                X = [genotype donepezil interaction cohort];
                names = ["genotype_APP","donepezil","genotype_x_donepezil","cohort_new"];
            else
                X = [genotype donepezil interaction];
                names = ["genotype_APP","donepezil","genotype_x_donepezil"];
            end

            FIT = clustered_cox( ...
                X,double(D.time_s),logical(D.censored),D.mouse_id,names);

            if ~FIT.success
                continue;
            end

            Rfit = FIT.results;

            iG = find(Rfit.predictor == "genotype_APP");
            iW = find(Rfit.predictor == "donepezil");
            iI = find(Rfit.predictor == "genotype_x_donepezil");

            gBase = linear_combo_hr(FIT,[iG],[1]);
            gDrug = linear_combo_hr(FIT,[iG iI],[1 1]);
            drugWT = linear_combo_hr(FIT,[iW],[1]);
            drugAPP = linear_combo_hr(FIT,[iW iI],[1 1]);
            intHR = linear_combo_hr(FIT,[iI],[1]);

            row = table( ...
                viewName,windowNames(w),st, ...
                numel(unique(D.mouse_id)),height(D),sum(D.event), ...
                gBase.HR,gBase.CI_low,gBase.CI_high,gBase.p, ...
                gDrug.HR,gDrug.CI_low,gDrug.CI_high,gDrug.p, ...
                drugWT.HR,drugWT.CI_low,drugWT.CI_high,drugWT.p, ...
                drugAPP.HR,drugAPP.CI_low,drugAPP.CI_high,drugAPP.p, ...
                intHR.HR,intHR.CI_low,intHR.CI_high,intHR.p, ...
                'VariableNames',{ ...
                    'cohort_view','time_bin','state','n_mice','n_bouts','n_events', ...
                    'HR_APP_vs_WT_baseline','CI_gbase_low','CI_gbase_high','p_gbase', ...
                    'HR_APP_vs_WT_donepezil','CI_gdrug_low','CI_gdrug_high','p_gdrug', ...
                    'HR_donepezil_vs_baseline_WT','CI_wt_low','CI_wt_high','p_wt', ...
                    'HR_donepezil_vs_baseline_APP','CI_app_low','CI_app_high','p_app', ...
                    'HR_genotype_x_donepezil','CI_int_low','CI_int_high','p_interaction'});

            COX_INTERACTION = [COX_INTERACTION; row]; %#ok<AGROW>
        end
    end
end

writetable(COX_INTERACTION, ...
    fullfile(survDir,'donepezil_baseline_vs_donepezil_cox_interaction.csv'));

%% Leave-one-mouse-out interaction sensitivity (combined cohort)

LOO = table();

for w = 1:numel(windowNames)

    for st = survivalStates

        D0 = SURV(SURV.time_bin == windowNames(w) & SURV.state == st,:);

        paired = paired_ids_survival(D0);
        D0 = D0(ismember(D0.mouse_id,paired),:);

        fitFull = fit_condition_interaction_cox(D0,true);

        if ~fitFull.success
            continue;
        end

        miceHere = unique(D0.mouse_id);

        for m = 1:numel(miceHere)

            removedMouse = miceHere(m);
            D = D0(D0.mouse_id ~= removedMouse,:);

            fitOne = fit_condition_interaction_cox(D,true);

            g = unique(D0.genotype(D0.mouse_id == removedMouse));
            c = unique(D0.cohort(D0.mouse_id == removedMouse));

            if fitOne.success
                hr = fitOne.HR;
                p = fitOne.p;
            else
                hr = NaN;
                p = NaN;
            end

            row = table( ...
                windowNames(w),st,removedMouse,g(1),c(1), ...
                fitFull.HR,fitFull.p,hr,p, ...
                hr-fitFull.HR,abs(hr-fitFull.HR), ...
                'VariableNames',{ ...
                    'time_bin','state','removed_mouse','removed_genotype','removed_cohort', ...
                    'full_interaction_HR','full_interaction_p', ...
                    'interaction_HR_without_mouse','interaction_p_without_mouse', ...
                    'change_HR','abs_change_HR'});

            LOO = [LOO; row]; %#ok<AGROW>
        end
    end
end

writetable(LOO, ...
    fullfile(survDir,'donepezil_survival_interaction_leave_one_mouse_out.csv'));

%% Donepezil-only survival figures

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DRUG_SURV;
        displayName = "Combined cohorts";
    elseif viewName == "thesis"
        V = DRUG_SURV(DRUG_SURV.cohort == "thesis",:);
        displayName = "Original cohort";
    else
        V = DRUG_SURV(DRUG_SURV.cohort == "new",:);
        displayName = "New cohort";
    end

    fig = figure('Color','w','Position',[80 80 1250 720]);

    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    title(tl,displayName + " - donepezil session bout survival");

    for s = 1:numel(survivalStates)
        for w = 1:numel(windowNames)

            ax = nexttile;

            D = V( ...
                V.state == survivalStates(s) & ...
                V.time_bin == windowNames(w),:);

            plot_survival_panel(ax,D,viewName,windowNames(w),survivalStates(s));
        end
    end

    exportgraphics( ...
        fig, ...
        fullfile(figDir,sprintf('donepezil_survival_%s.png',char(viewName))), ...
        'Resolution',300);
end

%% Key response figures: delta from pre

plot_delta_metric(DELTA,views,"state_percent","percent_time", ...
    ["WK","NREM","REM"], ...
    "Change from pre in state time (%)", ...
    figDir,"donepezil_delta_state_percent");

plot_delta_metric(DELTA,views,"bout_metric","bouts_per_h", ...
    ["WK","NREM","REM"], ...
    "Change from pre in bout onsets / h", ...
    figDir,"donepezil_delta_bout_rate");

plot_delta_metric(DELTA,views,"bout_metric","mean_bout_duration_s", ...
    ["WK","NREM","REM"], ...
    "Change from pre in mean bout duration (s)", ...
    figDir,"donepezil_delta_mean_bout_duration");

plot_delta_metric(DELTA,views,"fragmentation","MA_per_h_state", ...
    ["NREM","REM"], ...
    "Change from pre in fragmentation (MA / h state)", ...
    figDir,"donepezil_delta_fragmentation");

%% Console summary

fprintf('\n============================================================\n');
fprintf('KEY DONEPEZIL DELTA RESULTS\n');
fprintf('============================================================\n\n');

Dshow = DELTA_STATS( ...
    DELTA_STATS.cohort_view == "all" & ...
    ismember(DELTA_STATS.state,["NREM","REM"]) & ...
    ( ...
      (DELTA_STATS.domain == "state_percent") | ...
      (DELTA_STATS.metric == "mean_bout_duration_s") | ...
      (DELTA_STATS.domain == "fragmentation") ...
    ),:);

disp(Dshow(:, ...
    ["time_bin","domain","state","metric", ...
     "n_WT","n_APP","WT_mean_delta","APP_mean_delta", ...
     "p_welch","p_sidak_2phase","Cohen_d"]));

fprintf('\n============================================================\n');
fprintf('KEY DIFFERENCE-IN-DIFFERENCES RESULTS\n');
fprintf('============================================================\n\n');

Dshow = DID_STATS( ...
    DID_STATS.cohort_view == "all" & ...
    ismember(DID_STATS.state,["NREM","REM"]) & ...
    ( ...
      (DID_STATS.domain == "state_percent") | ...
      (DID_STATS.metric == "mean_bout_duration_s") | ...
      (DID_STATS.domain == "fragmentation") ...
    ),:);

disp(Dshow(:, ...
    ["time_bin","domain","state","metric", ...
     "n_WT","n_APP","WT_mean_DiD","APP_mean_DiD", ...
     "p_welch","p_sidak_2phase","Cohen_d"]));

fprintf('\n============================================================\n');
fprintf('DONEPEZIL-ONLY CLUSTERED COX\n');
fprintf('============================================================\n\n');

disp(COX_DONEPEZIL);

fprintf('\n============================================================\n');
fprintf('WITHIN-SESSION PRE vs DRUG_0_3H COX INTERACTION\n');
fprintf('============================================================\n\n');

disp(COX_WITHIN_SESSION);

fprintf('\n============================================================\n');
fprintf('BASELINE vs DONEPEZIL COX INTERACTION\n');
fprintf('============================================================\n\n');

disp(COX_INTERACTION);

fprintf('\n============================================================\n');
fprintf('MOST INFLUENTIAL MICE FOR SURVIVAL INTERACTION\n');
fprintf('============================================================\n\n');

if ~isempty(LOO)
    for st = survivalStates
        D = LOO(LOO.state == st & LOO.time_bin == "DRUG_0_3H",:);

        if isempty(D), continue; end

        D = sortrows(D,'abs_change_HR','descend');

        fprintf('\n%s DRUG_0_3H top 5:\n',st);
        disp(D(1:min(5,height(D)), ...
            ["removed_mouse","removed_genotype","removed_cohort", ...
             "full_interaction_HR","interaction_HR_without_mouse", ...
             "interaction_p_without_mouse","abs_change_HR"]));
    end
end

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');
fprintf('Outputs:  %s\n',outDir);
fprintf('Figures:  %s\n',figDir);
fprintf('Survival: %s\n\n',survDir);

%% ========================================================================
% Local functions
% =========================================================================


function [windowStart,windowEnd,status] = get_mouse_windows( ...
    mouseID,DrugTimes,defaultDrugStartSec,drugWindowSec)

% Define manipulation-aligned windows from the ACTUAL donepezil injection.
%
% PRE        = up to 3 h immediately before injection
% DRUG_0_3H  = injection -> +3 h        (primary acute drug window)
% DRUG_3_6H  = +3 h -> +6 h             (secondary late window)
%
% The same elapsed-time boundaries are used in the paired baseline session.

drugStart = defaultDrugStartSec;
status = "default-drug-start-03:00:00";

if ~isempty(DrugTimes)

    vars = string(DrugTimes.Properties.VariableNames);

    if ~ismember("mouse_id",vars)
        warning('drug_times.csv has no mouse_id column. Using default timing.');
    else
        enabled = true(height(DrugTimes),1);

        if ismember("include_analysis",vars)
            x = DrugTimes.include_analysis;

            if isnumeric(x) || islogical(x)
                enabled = logical(x);
            else
                sx = lower(strtrim(string(x)));
                enabled = ismember(sx,["1","true","yes","y"]);
            end
        end

        idx = find( ...
            string(DrugTimes.mouse_id) == string(mouseID) & ...
            enabled);

        if numel(idx) > 1
            error('drug_times.csv contains multiple active rows for mouse %s.',mouseID);
        end

        if numel(idx) == 1

            source = "timing-registry";
            if ismember("timing_source",vars)
                q = string(DrugTimes.timing_source(idx));
                if ~ismissing(q) && strlength(strtrim(q)) > 0
                    source = strtrim(q);
                end
            end

            startSec = NaN;

            % Prefer numeric seconds: October timings have second precision.
            if ismember("drug_elapsed_s",vars)
                startSec = scalar_to_double(DrugTimes.drug_elapsed_s(idx));
            end

            if ~isfinite(startSec) && ismember("drug_elapsed_hms",vars)
                startSec = parse_elapsed_hms( ...
                    string(DrugTimes.drug_elapsed_hms(idx)));
            end

            if isfinite(startSec)
                drugStart = startSec;
                status = "recorded:" + source;
            end
        end
    end
end

if drugStart <= 0
    error('Invalid donepezil injection time for %s: %.3f s.',mouseID,drugStart);
end

preStart = max(0,drugStart-drugWindowSec);

if drugStart < drugWindowSec
    status = status + "/pre-shorter-than-3h";
end

windowStart = [preStart, drugStart, drugStart+drugWindowSec];
windowEnd   = [drugStart, drugStart+drugWindowSec, drugStart+2*drugWindowSec];
end

function x = scalar_to_double(v)

x = NaN;

if isnumeric(v) || islogical(v)
    q = double(v);
    if isscalar(q) && isfinite(q)
        x = q;
    end
    return;
end

s = strtrim(string(v));

if ismissing(s) || strlength(s) == 0
    return;
end

q = str2double(s);

if isfinite(q)
    x = q;
end
end


function sec = parse_elapsed_hms(x)

sec = NaN;

if ismissing(x) || strlength(strtrim(x)) == 0
    return;
end

x = strtrim(string(x));

parts = split(x,":");

if numel(parts) ~= 3
    error('Timing value "%s" must use HH:MM:SS.',x);
end

h = str2double(parts(1));
m = str2double(parts(2));
s = str2double(parts(3));

if any(~isfinite([h m s])) || h < 0 || m < 0 || m >= 60 || s < 0 || s >= 60
    error('Invalid elapsed HH:MM:SS value: %s',x);
end

sec = h*3600 + m*60 + s;
end


function row = metric_row( ...
    sid,mouseID,genotype,sex,cohort,condition, ...
    timeBin,wStart,wEnd,timingStatus,domain,state,metric,value, ...
    validCoverage,recordingCoverage,includePrimary)

row = table( ...
    sid,mouseID,genotype,sex,cohort,condition, ...
    timeBin,wStart,wEnd,timingStatus,domain,state,metric,value, ...
    validCoverage,recordingCoverage,includePrimary, ...
    'VariableNames',{ ...
        'session_id','mouse_id','genotype','sex','cohort','condition', ...
        'time_bin','window_start_s','window_end_s','timing_status', ...
        'domain','state','metric','value', ...
        'valid_coverage_pct','recording_coverage_pct','include_primary'});
end


function scoreOut = reclassify_short_wake(scoreIn,maxDur)

scoreOut = double(scoreIn(:));

isWake = scoreOut == 0;
d = diff([false; isWake; false]);

starts = find(d == 1);
ends_ = find(d == -1)-1;

dur = ends_-starts+1;

short = find(dur <= maxDur);

for i = reshape(short,1,[])
    scoreOut(starts(i):ends_(i)) = 15;
end
end


function B = build_bout_table(score,t)

score = double(score(:));
t = double(t(:));

validCodes = [0 1 2 15];

invalid = ~ismember(score,validCodes);

x = score;
x(invalid) = -999;

change = [true; x(2:end) ~= x(1:end-1)];

starts = find(change);
ends_ = [starts(2:end)-1; numel(x)];

stateCode = x(starts);

keep = ismember(stateCode,validCodes);

starts = starts(keep);
ends_ = ends_(keep);
stateCode = stateCode(keep);

n = numel(starts);

state = strings(n,1);
onset = nan(n,1);
endExclusive = nan(n,1);
duration = nan(n,1);
leftCensored = false(n,1);
rightCensored = false(n,1);
stateBefore = nan(n,1);

for i = 1:n

    s = starts(i);
    e = ends_(i);

    switch stateCode(i)
        case 0
            state(i) = "WK";
        case 1
            state(i) = "NREM";
        case 2
            state(i) = "REM";
        case 15
            state(i) = "MA";
    end

    onset(i) = t(s);
    duration(i) = e-s+1;
    endExclusive(i) = onset(i)+duration(i);

    if s == 1
        leftCensored(i) = true;
    elseif invalid(s-1)
        leftCensored(i) = true;
    else
        stateBefore(i) = score(s-1);
    end

    if e == numel(score)
        rightCensored(i) = true;
    elseif invalid(e+1)
        rightCensored(i) = true;
    end
end

B = table( ...
    stateCode,state,onset,endExclusive,duration, ...
    leftCensored,rightCensored,stateBefore, ...
    'VariableNames',{ ...
        'state_code','state','onset_s','end_exclusive_s', ...
        'duration_s','left_censored','right_censored','state_before'});
end


function ids = paired_ids_for_metric(D)

ids0 = unique(D.mouse_id);

ids = strings(0,1);

for i = 1:numel(ids0)
    m = ids0(i);
    conds = unique(D.condition(D.mouse_id == m));

    if all(ismember(["baseline","donepezil"],conds))
        ids(end+1,1) = m; %#ok<AGROW>
    end
end
end


function ids = paired_ids_phase_survival(D,phaseA,phaseB)

ids = strings(0,1);

for m = unique(D.mouse_id)'

    q = D(D.mouse_id == m,:);
    bins = unique(q.time_bin);

    if ismember(phaseA,bins) && ismember(phaseB,bins)
        ids(end+1,1) = m; %#ok<AGROW>
    end
end
end


function ids = paired_ids_survival(D)

ids0 = unique(D.mouse_id);
ids = strings(0,1);

for i = 1:numel(ids0)
    m = ids0(i);
    conds = unique(D.condition(D.mouse_id == m));

    if all(ismember(["baseline","donepezil"],conds))
        ids(end+1,1) = m; %#ok<AGROW>
    end
end
end


function [p,tstat,df,d] = welch_stats(WT,APP)

WT = double(WT(:));
APP = double(APP(:));

WT = WT(isfinite(WT));
APP = APP(isfinite(APP));

p = NaN;
tstat = NaN;
df = NaN;
d = NaN;

if numel(WT) < 2 || numel(APP) < 2
    return;
end

[~,p,~,stats] = ttest2(APP,WT,'Vartype','unequal');

tstat = stats.tstat;
df = stats.df;

n1 = numel(WT);
n2 = numel(APP);

s1 = std(WT);
s2 = std(APP);

sp = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/max(1,n1+n2-2));

if sp > 0
    d = (mean(APP)-mean(WT))/sp;
end
end


function T = add_sidak_two_phase(T)

T.p_sidak_2phase = nan(height(T),1);

if isempty(T)
    return;
end

groups = unique(T(:,["cohort_view","domain","state","metric"]),'rows');

for i = 1:height(groups)

    g = groups(i,:);

    idx = ...
        T.cohort_view == g.cohort_view & ...
        T.domain == g.domain & ...
        T.state == g.state & ...
        T.metric == g.metric;

    p = double(T.p_welch(idx));

    pAdj = 1-(1-p).^2;
    pAdj(pAdj>1) = 1;

    T.p_sidak_2phase(idx) = pAdj;
end
end


function assert_unique_mouse_metric_rows(T,keyNames,label)

if isempty(T)
    return;
end

K = T(:,cellstr(keyNames));
Ku = unique(K,'rows','stable');

if height(Ku) ~= height(K)
    error('%s contains duplicate mouse/metric/window rows (%d rows, %d unique).', ...
        label,height(K),height(Ku));
end
end


function OUT = safe_lme_anova(Tm,formula)

OUT = table();

% Remove undefined categorical rows defensively.
vars = Tm.Properties.VariableNames;
valid = true(height(Tm),1);

for i = 1:numel(vars)
    x = Tm.(vars{i});
    if isnumeric(x)
        valid = valid & isfinite(x);
    elseif iscategorical(x)
        valid = valid & ~isundefined(x);
    end
end

Tm = Tm(valid,:);

if height(Tm) < 4 || numel(unique(Tm.Mouse)) < 2
    return;
end

try
    lme = fitlme(Tm,formula,'FitMethod','REML');
catch ME
    warning('LME fit failed for "%s": %s',formula,ME.message);
    return;
end

% MATLAB versions differ in both DFMethod support and the type returned by
% anova(lme). Newer releases return a table; older releases may return the
% deprecated Statistics Toolbox "dataset" class.
try
    A = anova(lme,'DFMethod','satterthwaite');
catch
    try
        A = anova(lme);
    catch ME
        warning('LME ANOVA failed for "%s": %s',formula,ME.message);
        return;
    end
end

if isempty(A)
    return;
end

% Convert legacy dataset output to a table when possible. dataset uses
% ObsNames rather than table RowNames, so accessing Properties.RowNames on
% it throws "Unknown dataset property: RowNames" in older MATLAB versions.
legacyObsNames = strings(0,1);

if isa(A,'dataset')
    try
        obs = A.Properties.ObsNames;
        if ~isempty(obs)
            legacyObsNames = string(obs(:));
        end
    catch
        legacyObsNames = strings(0,1);
    end

    try
        A = dataset2table(A);
    catch
        % Keep A as a dataset. The generic helpers below support either
        % dataset or table, so conversion is only a convenience.
    end
end

nA = size(A,1);

% Recover ANOVA term names.
if istable(A)
    if ~isempty(A.Properties.RowNames)
        term = string(A.Properties.RowNames(:));
    elseif ismember('Term',A.Properties.VariableNames)
        term = string(A.Term(:));
    elseif numel(legacyObsNames) == nA
        term = legacyObsNames;
    else
        term = "term_" + string((1:nA)');
    end
elseif isa(A,'dataset')
    try
        obs = A.Properties.ObsNames;
        if ~isempty(obs)
            term = string(obs(:));
        elseif ismember('Term',A.Properties.VarNames)
            term = string(A.Term(:));
        else
            term = "term_" + string((1:nA)');
        end
    catch
        term = "term_" + string((1:nA)');
    end
else
    warning('Unsupported ANOVA output class: %s',class(A));
    return;
end

F   = get_anova_numeric(A,["FStat","F"]);
DF1 = get_anova_numeric(A,["DF1","NumDF"]);
DF2 = get_anova_numeric(A,["DF2","DenDF"]);
p   = get_anova_numeric(A,["pValue","p"]);

OUT = table(term,F,DF1,DF2,p, ...
    'VariableNames',{'term','F','DF1','DF2','p'});
end


function x = get_anova_numeric(A,names)

nA = size(A,1);
x = nan(nA,1);

if istable(A)
    varNames = string(A.Properties.VariableNames);
elseif isa(A,'dataset')
    varNames = string(A.Properties.VarNames);
else
    return;
end

for n = names
    idx = find(varNames == n,1);
    if ~isempty(idx)
        vn = char(varNames(idx));
        try
            x = double(A.(vn));
            x = x(:);
        catch
            x = nan(nA,1);
        end
        return;
    end
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
    nMouse = sum(idx);

    if nMouse > 0
        w(idx) = 1/nMouse;
    end
end

if nMice > 0
    w = w*(n/nMice);
end
end


function [tEvent,S,medianSurvival] = weighted_km(time,event,w)

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

idx = find(S <= 0.5,1,'first');

if isempty(idx)
    medianSurvival = NaN;
else
    medianSurvival = tEvent(idx);
end
end


function FIT = clustered_cox(X,time,censored,cluster,predictorNames)

FIT = struct();
FIT.success = false;
FIT.results = table();
FIT.robustCov = [];
FIT.nClusters = 0;

X = double(X);
time = double(time(:));
censored = logical(censored(:));
cluster = string(cluster(:));
predictorNames = string(predictorNames(:));

if isvector(X)
    X = X(:);
end

valid = ...
    all(isfinite(X),2) & ...
    isfinite(time) & ...
    time>0 & ...
    strlength(cluster)>0;

X = X(valid,:);
time = time(valid);
censored = censored(valid);
cluster = cluster(valid);

clusters = unique(cluster);
G = numel(clusters);
p = size(X,2);

FIT.nClusters = G;

if G <= p+1 || sum(~censored) < p+1
    return;
end

w = mouse_equal_weights(cluster);

try
    [beta,~,~,stats] = coxphfit( ...
        X,time, ...
        'Censoring',censored, ...
        'Frequency',w, ...
        'Ties','breslow', ...
        'Baseline',0);
catch ME
    warning('Cox fit failed: %s',ME.message);
    return;
end

beta = beta(:);

scoreResiduals = double(stats.scores);

clusterScore = zeros(G,p);

for g = 1:G
    idx = cluster == clusters(g);
    clusterScore(g,:) = sum(scoreResiduals(idx,:),1,'omitnan');
end

bread = double(stats.covb);
meat = clusterScore'*clusterScore;

Vrobust = bread*meat*bread;

if G > 1
    Vrobust = Vrobust*(G/(G-1));
end

se = sqrt(max(diag(Vrobust),0));

z = beta./se;

dfCluster = G-1;
pT = 2*(1-tcdf(abs(z),dfCluster));

tcrit = tinv(0.975,dfCluster);

ciLow = exp(beta-tcrit.*se);
ciHigh = exp(beta+tcrit.*se);
HR = exp(beta);

FIT.results = table( ...
    predictorNames,beta,se,HR,ciLow,ciHigh,z,pT, ...
    repmat(G,p,1), ...
    'VariableNames',{ ...
        'predictor','beta','robust_SE','HR','CI_low','CI_high', ...
        'z_robust','p_robust_t','n_clusters'});

FIT.robustCov = Vrobust;
FIT.success = true;
end


function O = linear_combo_hr(FIT,idx,coef)

O = struct('HR',NaN,'CI_low',NaN,'CI_high',NaN,'p',NaN);

if ~FIT.success
    return;
end

idx = idx(:);
coef = double(coef(:));

R = FIT.results;

beta = sum(coef.*double(R.beta(idx)));

V = FIT.robustCov(idx,idx);

variance = coef'*V*coef;
se = sqrt(max(variance,0));

df = FIT.nClusters-1;
tcrit = tinv(0.975,df);

O.HR = exp(beta);
O.CI_low = exp(beta-tcrit*se);
O.CI_high = exp(beta+tcrit*se);

if se > 0
    O.p = 2*(1-tcdf(abs(beta/se),df));
end
end


function O = fit_condition_interaction_cox(D,includeCohort)

O = struct('success',false,'HR',NaN,'p',NaN);

if isempty(D) || ...
        numel(unique(D.genotype))<2 || ...
        numel(unique(D.condition))<2
    return;
end

genotype = double(D.genotype == "APP");
donepezil = double(D.condition == "donepezil");
interaction = genotype.*donepezil;

if includeCohort
    cohort = double(D.cohort == "new");
    X = [genotype donepezil interaction cohort];
    names = ["genotype_APP","donepezil","genotype_x_donepezil","cohort_new"];
else
    X = [genotype donepezil interaction];
    names = ["genotype_APP","donepezil","genotype_x_donepezil"];
end

FIT = clustered_cox( ...
    X,double(D.time_s),logical(D.censored),D.mouse_id,names);

if ~FIT.success
    return;
end

idx = find(FIT.results.predictor == "genotype_x_donepezil");

R = linear_combo_hr(FIT,idx,1);

O.success = true;
O.HR = R.HR;
O.p = R.p;
end


function plot_survival_panel(ax,D,viewName,windowName,stateName)

cla(ax);
hold(ax,'on');

if isempty(D)
    title(ax,stateName + " " + strrep(windowName,"_","-") + " - no data");
    ylim(ax,[0 1]);
    return;
end

WTcolor = [0.30 0.30 0.30];
APPcolor = [0.20 0.40 0.80];

WT = D(D.genotype == "WT",:);
APP = D(D.genotype == "APP",:);

if ~isempty(WT)
    w = mouse_equal_weights(WT.mouse_id);
    [t,s,~] = weighted_km(double(WT.time_s),logical(WT.event),w);

    stairs(ax,[0;t],[1;s], ...
        'LineWidth',1.8,'Color',WTcolor, ...
        'DisplayName',sprintf('WT (%d mice)',numel(unique(WT.mouse_id))));
end

if ~isempty(APP)
    w = mouse_equal_weights(APP.mouse_id);
    [t,s,~] = weighted_km(double(APP.time_s),logical(APP.event),w);

    stairs(ax,[0;t],[1;s], ...
        'LineWidth',1.8,'Color',APPcolor, ...
        'DisplayName',sprintf('APP (%d mice)',numel(unique(APP.mouse_id))));
end

genotype = double(D.genotype == "APP");

if viewName == "all"
    cohort = double(D.cohort == "new");
    FIT = clustered_cox( ...
        [genotype cohort],double(D.time_s),logical(D.censored), ...
        D.mouse_id,["genotype_APP","cohort_new"]);
else
    FIT = clustered_cox( ...
        genotype,double(D.time_s),logical(D.censored), ...
        D.mouse_id,"genotype_APP");
end

xlabel(ax,'Bout duration (s)');
ylabel(ax,'Survival probability');
title(ax,stateName + ", " + strrep(windowName,"_","-"));
ylim(ax,[0 1]);
xlim(ax,[0 survival_xlim(D.time_s)]);
box(ax,'off');
legend(ax,'Location','southwest','Box','off');

if FIT.success
    R = FIT.results(FIT.results.predictor == "genotype_APP",:);

    txt = sprintf( ...
        'APP vs WT HR = %.2f\n95%% CI %.2f-%.2f\nclustered p = %.3g', ...
        R.HR,R.CI_low,R.CI_high,R.p_robust_t);

    text(ax,0.97,0.95,txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',9);
end
end


function xmax = survival_xlim(time)

time = double(time(:));
time = time(isfinite(time) & time>0);

if isempty(time)
    xmax = 1;
    return;
end

xmax = prctile(time,97.5);

if ~isfinite(xmax) || xmax <= 0
    xmax = max(time);
end

if xmax <= 0
    xmax = 1;
end
end


function plot_delta_metric( ...
    DELTA,views,domainName,metricName,states,yLabel,outDir,fileStem)

WTcolor = [0.35 0.35 0.35];
APPcolor = [0.20 0.40 0.80];

bins = ["PRE","DRUG_0_3H","DRUG_3_6H"];
labels = {'Pre','Drug 0-3 h','Drug 3-6 h'};

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        V = DELTA;
        titlePrefix = "Combined cohorts";
    elseif viewName == "thesis"
        V = DELTA(DELTA.cohort == "thesis",:);
        titlePrefix = "Original cohort";
    else
        V = DELTA(DELTA.cohort == "new",:);
        titlePrefix = "New cohort";
    end

    V = V(V.domain == domainName & V.metric == metricName,:);

    if isempty(V)
        continue;
    end

    nStates = numel(states);

    fig = figure('Color','w','Position',[80 80 400*nStates 430]);

    tl = tiledlayout(1,nStates,'TileSpacing','compact','Padding','compact');
    title(tl,titlePrefix + " - " + strrep(metricName,"_"," "));

    for s = 1:nStates

        ax = nexttile;
        hold(ax,'on');

        for b = 1:numel(bins)

            D = V(V.state == states(s) & V.time_bin == bins(b),:);

            WT = double(D.delta_from_pre(D.genotype == "WT"));
            APP = double(D.delta_from_pre(D.genotype == "APP"));

            plot_group(ax,b-0.13,WT,WTcolor);
            plot_group(ax,b+0.13,APP,APPcolor);
        end

        yline(ax,0,'k:');

        xlim(ax,[0.5 3.5]);
        xticks(ax,1:3);
        xticklabels(ax,labels);
        xtickangle(ax,20);

        ylabel(ax,yLabel);
        title(ax,states(s));
        box(ax,'off');
    end

    exportgraphics( ...
        fig, ...
        fullfile(outDir,sprintf('%s_%s.png',char(fileStem),char(viewName))), ...
        'Resolution',300);
end
end


function plot_group(ax,xPos,values,color)

values = double(values(:));
values = values(isfinite(values));

n = numel(values);

if n == 0
    return;
end

if n == 1
    jitter = 0;
else
    jitter = linspace(-0.035,0.035,n)';
end

scatter(ax,xPos+jitter,values,34, ...
    'MarkerEdgeColor',color,'LineWidth',1.1);

m = mean(values);
sem = std(values)/sqrt(n);

errorbar(ax,xPos,m,sem,'o', ...
    'Color',color, ...
    'MarkerFaceColor',color, ...
    'MarkerEdgeColor',color, ...
    'LineWidth',1.6, ...
    'CapSize',8);
end
