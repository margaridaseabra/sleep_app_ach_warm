%% RUN_BASELINE_COMPLETE
% Unified baseline analysis from standardized 1-s scores.
% Windows: WHOLE, ZT3-9, ZT3-6, ZT6-9, ZT9-12.
% WHOLE keeps bouts continuous across bin boundaries.
% Fixed bins require >=95% valid coverage.
% Wake bouts <=15 s are reclassified as MA.
% Outputs state %, bout metrics, fragmentation, survival, cohort models,
% figures, and leave-one-mouse-out survival sensitivity.

clear; clc;

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);
addpath(fullfile(repoDir,'config'));
P = paths_config();

Sessions = readtable(P.sessions,'TextType','string','VariableNamingRule','preserve');
Mice = readtable(P.mice,'TextType','string','VariableNamingRule','preserve');
S = Sessions(Sessions.condition == "baseline",:);

standardDir = fullfile(P.sleep_output,'all','standardized_scores');
outDir = fullfile(P.sleep_output,'baseline_complete');
figDir = fullfile(outDir,'figures');
survDir = fullfile(outDir,'survival');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(survDir), mkdir(survDir); end

windowNames = ["WHOLE","ZT3_9","ZT3_6","ZT6_9","ZT9_12"];
fixedStart = [NaN,0,0,3*3600,6*3600];
fixedEnd   = [NaN,6*3600,3*3600,6*3600,9*3600];
stateNames = ["WK","NREM","REM","MA"];
stateCodes = [0,1,2,15];
carrierStates = ["NREM","REM"];
survivalStates = ["NREM","REM"];
views = ["thesis","new","all"];
minimumValidCoverage = 95;
maThresholdSec = 15;

STATE = table();
METRICS = table();
FRAG = table();
ALLBOUTS = table();

fprintf('\n============================================================\n');
fprintf('UNIFIED BASELINE ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Baseline sessions: %d\n\n',height(S));

for i = 1:height(S)
    sid = string(S.session_id(i));
    mouseID = string(S.mouse_id(i));
    mouseIdx = find(Mice.mouse_id == mouseID,1);
    if isempty(mouseIdx)
        warning('Mouse %s not found in mice.csv; skipping.',mouseID);
        continue;
    end

    genotype = string(Mice.genotype(mouseIdx));
    sex = string(Mice.sex(mouseIdx));
    cohort = string(S.cohort(i));
    scoreFile = fullfile(standardDir,char(sid + ".csv"));
    if ~isfile(scoreFile)
        warning('Missing standardized score file: %s',scoreFile);
        continue;
    end

    T = readtable(scoreFile);
    t = double(T.time_s(:));
    score = double(T.score(:));
    if isempty(t), continue; end
    t = t - min(t);
    scoreMA = reclassify_short_wake(score,maThresholdSec);
    B = build_bout_table(scoreMA,t);

    if ~isempty(B)
        B.session_id = repmat(sid,height(B),1);
        B.mouse_id = repmat(mouseID,height(B),1);
        B.genotype = repmat(genotype,height(B),1);
        B.sex = repmat(sex,height(B),1);
        B.cohort = repmat(cohort,height(B),1);
        ALLBOUTS = [ALLBOUTS; B];
    end

    recordingEnd = floor(max(t)) + 1;
    fprintf('[%02d/%02d] %-8s %-4s  %.2f h\n',i,height(S),mouseID,genotype,recordingEnd/3600);

    for w = 1:numel(windowNames)
        windowName = windowNames(w);
        if windowName == "WHOLE"
            wStart = 0; wEnd = recordingEnd;
        else
            wStart = fixedStart(w); wEnd = fixedEnd(w);
        end

        inWindow = t >= wStart & t < wEnd;
        expectedSec = max(wEnd-wStart,1);
        nAvailable = sum(inWindow);
        validMask = inWindow & ismember(scoreMA,[0 1 2 15]);
        validSec = sum(validMask);
        recordingCoverage = 100*nAvailable/expectedSec;
        validCoverage = 100*validSec/expectedSec;
        includePrimary = validCoverage >= minimumValidCoverage;
        validHours = validSec/3600;

        for s = 1:numel(stateNames)
            stateName = stateNames(s); stateCode = stateCodes(s);
            stateSec = sum(inWindow & scoreMA == stateCode);
            if validSec > 0, percentTime = 100*stateSec/validSec; else, percentTime = NaN; end
            row = table(sid,mouseID,genotype,sex,cohort,windowName,wStart,wEnd,stateName, ...
                stateSec,validSec,recordingCoverage,validCoverage,includePrimary,percentTime, ...
                'VariableNames',{'session_id','mouse_id','genotype','sex','cohort','time_bin', ...
                'window_start_s','window_end_s','state','state_seconds','valid_seconds', ...
                'recording_coverage_pct','valid_coverage_pct','include_primary','percent_time'});
            STATE = [STATE; row];
        end

        for s = 1:numel(stateNames)
            stateName = stateNames(s); stateCode = stateCodes(s);
            stateSec = sum(inWindow & scoreMA == stateCode);
            stateHours = stateSec/3600;
            if isempty(B)
                BR = B;
            else
                BR = B(B.state_code == stateCode & B.onset_s >= wStart & B.onset_s < wEnd,:);
            end

            if isempty(BR)
                nBouts = 0; nDurationBouts = 0; meanDuration = NaN; medianDuration = NaN;
            else
                nBouts = sum(~BR.left_censored);
                if windowName == "WHOLE"
                    useDuration = ~BR.left_censored & ~BR.right_censored;
                else
                    useDuration = ~BR.left_censored & ~BR.right_censored & BR.end_exclusive_s <= wEnd;
                end
                durations = double(BR.duration_s(useDuration));
                nDurationBouts = numel(durations);
                if isempty(durations)
                    meanDuration = NaN; medianDuration = NaN;
                else
                    meanDuration = mean(durations); medianDuration = median(durations);
                end
            end
            if validHours > 0, boutsPerHour = nBouts/validHours; else, boutsPerHour = NaN; end

            row = table(sid,mouseID,genotype,sex,cohort,windowName,wStart,wEnd,stateName, ...
                stateSec,stateHours,validSec,validCoverage,includePrimary,nBouts,boutsPerHour, ...
                nDurationBouts,meanDuration,medianDuration, ...
                'VariableNames',{'session_id','mouse_id','genotype','sex','cohort','time_bin', ...
                'window_start_s','window_end_s','state','state_seconds','state_hours','valid_seconds', ...
                'valid_coverage_pct','include_primary','n_bouts','bouts_per_h','n_duration_bouts', ...
                'mean_bout_duration_s','median_bout_duration_s'});
            METRICS = [METRICS; row];
        end

        if isempty(B)
            maNREM = 0; maREM = 0;
        else
            MA = B(B.state_code == 15 & B.onset_s >= wStart & B.onset_s < wEnd & ~B.left_censored,:);
            maNREM = sum(MA.state_before == 1);
            maREM = sum(MA.state_before == 2);
        end
        nremSec = sum(inWindow & scoreMA == 1); remSec = sum(inWindow & scoreMA == 2);
        nremHours = nremSec/3600; remHours = remSec/3600;
        if nremHours > 0, nremFrag = maNREM/nremHours; else, nremFrag = NaN; end
        if remHours > 0, remFrag = maREM/remHours; else, remFrag = NaN; end

        rowNREM = table(sid,mouseID,genotype,sex,cohort,windowName,wStart,wEnd,string("NREM"), ...
            maNREM,nremSec,nremHours,nremFrag,validCoverage,includePrimary, ...
            'VariableNames',{'session_id','mouse_id','genotype','sex','cohort','time_bin', ...
            'window_start_s','window_end_s','carrier_state','n_MA','state_seconds','state_hours', ...
            'MA_per_h_state','valid_coverage_pct','include_primary'});
        rowREM = table(sid,mouseID,genotype,sex,cohort,windowName,wStart,wEnd,string("REM"), ...
            maREM,remSec,remHours,remFrag,validCoverage,includePrimary, ...
            'VariableNames',{'session_id','mouse_id','genotype','sex','cohort','time_bin', ...
            'window_start_s','window_end_s','carrier_state','n_MA','state_seconds','state_hours', ...
            'MA_per_h_state','valid_coverage_pct','include_primary'});
        FRAG = [FRAG; rowNREM; rowREM];
    end
end

writetable(STATE,fullfile(outDir,'baseline_state_percent_per_mouse.csv'));
writetable(METRICS,fullfile(outDir,'baseline_bout_metrics_per_mouse.csv'));
writetable(FRAG,fullfile(outDir,'baseline_fragmentation_per_mouse.csv'));
writetable(ALLBOUTS,fullfile(outDir,'baseline_bouts.csv'));

SP = STATE(logical(STATE.include_primary),:);
BM = METRICS(logical(METRICS.include_primary),:);
FG = FRAG(logical(FRAG.include_primary),:);

fprintf('\nUsable sessions by window:\n');
for w = 1:numel(windowNames)
    D = SP(SP.time_bin == windowNames(w) & SP.state == "REM",:);
    Nall = sum(STATE.time_bin == windowNames(w) & STATE.state == "REM");
    fprintf('%-8s %d/%d\n',windowNames(w),height(D),Nall);
end

%% Welch APP vs WT statistics
WELCH = table();

for v = 1:numel(views)
    Dview = subset_view(SP,views(v));
    for w = 1:numel(windowNames)
        for s = 1:numel(stateNames)
            D = Dview(Dview.time_bin == windowNames(w) & Dview.state == stateNames(s),:);
            WELCH = [WELCH; make_welch_row(views(v),windowNames(w),"state_percent",stateNames(s),"percent_time",D)];
        end
    end
end

boutMetricNames = ["bouts_per_h","mean_bout_duration_s","median_bout_duration_s"];
for v = 1:numel(views)
    Dview = subset_view(BM,views(v));
    for w = 1:numel(windowNames)
        for s = 1:numel(stateNames)
            D = Dview(Dview.time_bin == windowNames(w) & Dview.state == stateNames(s),:);
            for m = 1:numel(boutMetricNames)
                WELCH = [WELCH; make_welch_row(views(v),windowNames(w),"bout_metric",stateNames(s),boutMetricNames(m),D)];
            end
        end
    end
end

for v = 1:numel(views)
    Dview = subset_view(FG,views(v));
    for w = 1:numel(windowNames)
        for s = 1:numel(carrierStates)
            D = Dview(Dview.time_bin == windowNames(w) & Dview.carrier_state == carrierStates(s),:);
            WELCH = [WELCH; make_welch_row(views(v),windowNames(w),"fragmentation",carrierStates(s),"MA_per_h_state",D)];
        end
    end
end
writetable(WELCH,fullfile(outDir,'baseline_welch_stats.csv'));

%% Combined Genotype x Cohort models
MODEL = table();
for w = 1:numel(windowNames)
    for s = 1:numel(stateNames)
        D = SP(SP.time_bin == windowNames(w) & SP.state == stateNames(s),:);
        MODEL = [MODEL; make_model_row(windowNames(w),"state_percent",stateNames(s),"percent_time",D)];
    end
end

modelBoutMetrics = ["bouts_per_h","mean_bout_duration_s","median_bout_duration_s"];
for w = 1:numel(windowNames)
    for s = 1:numel(stateNames)
        D = BM(BM.time_bin == windowNames(w) & BM.state == stateNames(s),:);
        for m = 1:numel(modelBoutMetrics)
            MODEL = [MODEL; make_model_row(windowNames(w),"bout_metric",stateNames(s),modelBoutMetrics(m),D)];
        end
    end
end

for w = 1:numel(windowNames)
    for s = 1:numel(carrierStates)
        D = FG(FG.time_bin == windowNames(w) & FG.carrier_state == carrierStates(s),:);
        MODEL = [MODEL; make_model_row(windowNames(w),"fragmentation",carrierStates(s),"MA_per_h_state",D)];
    end
end
writetable(MODEL,fullfile(outDir,'baseline_genotype_cohort_models.csv'));

%% Mouse-level figures: WHOLE + all bins
plot_metric_windows(SP,views,windowNames,stateNames,"percent_time","Time in state (%)",figDir,"baseline_state_percent");
plot_metric_windows(BM,views,windowNames,stateNames,"mean_bout_duration_s","Mean bout duration (s)",figDir,"baseline_mean_bout_duration");
plot_metric_windows(BM,views,windowNames,stateNames,"median_bout_duration_s","Median bout duration (s)",figDir,"baseline_median_bout_duration");
plot_metric_windows(BM,views,windowNames,stateNames,"bouts_per_h","Bout onsets per valid recording hour",figDir,"baseline_bout_rate");
plot_fragmentation_windows(FG,views,windowNames,carrierStates,figDir);

%% Build survival dataset
SURV = table();
for w = 1:numel(windowNames)
    windowName = windowNames(w);
    QC = BM(BM.time_bin == windowName & BM.state == "NREM",:);
    primarySessions = unique(QC.session_id(logical(QC.include_primary)));

    for s = 1:numel(survivalStates)
        stateName = survivalStates(s);
        for q = 1:numel(primarySessions)
            sid = primarySessions(q);
            D = ALLBOUTS(ALLBOUTS.session_id == sid & ALLBOUTS.state == stateName,:);
            if isempty(D), continue; end
            D = D(~logical(D.left_censored),:);
            if isempty(D), continue; end

            if windowName == "WHOLE"
                observedTime = double(D.duration_s);
                event = ~logical(D.right_censored);
            else
                wStart = fixedStart(w); wEnd = fixedEnd(w);
                D = D(D.onset_s >= wStart & D.onset_s < wEnd,:);
                if isempty(D), continue; end
                remainingToBoundary = wEnd - double(D.onset_s);
                observedTime = min(double(D.duration_s),remainingToBoundary);
                event = ~logical(D.right_censored) & double(D.end_exclusive_s) <= wEnd;
            end

            censored = ~event;
            valid = isfinite(observedTime) & observedTime > 0;
            D = D(valid,:); observedTime = observedTime(valid); event = event(valid); censored = censored(valid);
            if isempty(D), continue; end

            R = table(D.session_id,D.mouse_id,D.genotype,D.cohort, ...
                repmat(windowName,height(D),1),repmat(stateName,height(D),1), ...
                observedTime,event,censored, ...
                'VariableNames',{'session_id','mouse_id','genotype','cohort','time_bin','state','time_s','event','censored'});
            SURV = [SURV; R];
        end
    end
end
writetable(SURV,fullfile(survDir,'baseline_survival_dataset.csv'));

fprintf('\n============================================================\n');
fprintf('SURVIVAL DATASET\n');
fprintf('============================================================\n');
fprintf('Rows: %d\n',height(SURV));
fprintf('Mice: %d\n\n',numel(unique(SURV.mouse_id)));

%% KM + clustered Cox
COX = table();
KMSUMMARY = table();
for v = 1:numel(views)
    viewName = views(v);
    V = subset_view(SURV,viewName);
    for w = 1:numel(windowNames)
        windowName = windowNames(w);
        for s = 1:numel(survivalStates)
            stateName = survivalStates(s);
            D = V(V.time_bin == windowName & V.state == stateName,:);
            if isempty(D), continue; end

            for genotypeName = ["WT","APP"]
                G = D(D.genotype == genotypeName,:);
                if isempty(G), continue; end
                wg = mouse_equal_weights(G.mouse_id);
                [~,~,medianSurvival] = weighted_km(double(G.time_s),logical(G.event),wg);
                KMSUMMARY = [KMSUMMARY; km_summary_row(viewName,windowName,stateName,genotypeName,G,medianSurvival)];
            end

            if viewName == "all"
                genotype = double(D.genotype == "APP");
                cohort = double(D.cohort == "new");
                FIT = clustered_cox([genotype cohort],double(D.time_s),logical(D.censored),D.mouse_id,["genotype_APP","cohort_new"]);
                if FIT.success
                    R = FIT.results(FIT.results.predictor == "genotype_APP",:);
                    COX = [COX; make_cox_row(viewName,windowName,stateName,D,R)];
                end
            else
                genotype = double(D.genotype == "APP");
                FIT = clustered_cox(genotype,double(D.time_s),logical(D.censored),D.mouse_id,"genotype_APP");
                if FIT.success
                    COX = [COX; make_cox_row(viewName,windowName,stateName,D,FIT.results(1,:))];
                end
            end
        end
    end
end
writetable(KMSUMMARY,fullfile(survDir,'baseline_survival_km_summary.csv'));
writetable(COX,fullfile(survDir,'baseline_survival_cox.csv'));

%% Genotype x Cohort survival models
INTERACTION = table();
for w = 1:numel(windowNames)
    windowName = windowNames(w);
    for s = 1:numel(survivalStates)
        stateName = survivalStates(s);
        D = SURV(SURV.time_bin == windowName & SURV.state == stateName,:);
        if isempty(D), continue; end

        genotype = double(D.genotype == "APP");
        cohort = double(D.cohort == "new");
        interaction = genotype .* cohort;
        FIT = clustered_cox([genotype cohort interaction],double(D.time_s),logical(D.censored),D.mouse_id, ...
            ["genotype_APP","cohort_new","genotype_x_cohort"]);
        if ~FIT.success, continue; end

        R = FIT.results;
        g = R(R.predictor == "genotype_APP",:);
        irow = R(R.predictor == "genotype_x_cohort",:);
        idxG = find(R.predictor == "genotype_APP");
        idxI = find(R.predictor == "genotype_x_cohort");
        betaNew = g.beta + irow.beta;
        varNew = FIT.robustCov(idxG,idxG) + FIT.robustCov(idxI,idxI) + 2*FIT.robustCov(idxG,idxI);
        seNew = sqrt(max(varNew,0));
        dfCluster = FIT.nClusters - 1;
        tcrit = tinv(0.975,dfCluster);
        hrNew = exp(betaNew);
        ciNewLow = exp(betaNew-tcrit*seNew);
        ciNewHigh = exp(betaNew+tcrit*seNew);
        if seNew > 0
            zNew = betaNew/seNew;
            pNew = 2*(1-tcdf(abs(zNew),dfCluster));
        else
            pNew = NaN;
        end

        row = table(windowName,stateName,FIT.nClusters,height(D),sum(D.event), ...
            g.HR,g.CI_low,g.CI_high,g.p_robust_t,hrNew,ciNewLow,ciNewHigh,pNew,irow.HR,irow.p_robust_t, ...
            'VariableNames',{'time_bin','state','n_mice','n_bouts','n_events', ...
            'HR_APP_original','CI_original_low','CI_original_high','p_APP_original', ...
            'HR_APP_new','CI_new_low','CI_new_high','p_APP_new','HR_interaction','p_interaction'});
        INTERACTION = [INTERACTION; row];
    end
end
writetable(INTERACTION,fullfile(survDir,'baseline_survival_cox_interaction.csv'));

%% Leave-one-mouse-out Cox sensitivity
LOO = table();
for v = 1:numel(views)
    viewName = views(v);
    V = subset_view(SURV,viewName);
    for w = 1:numel(windowNames)
        windowName = windowNames(w);
        for s = 1:numel(survivalStates)
            stateName = survivalStates(s);
            D0 = V(V.time_bin == windowName & V.state == stateName,:);
            if isempty(D0), continue; end
            fullFit = fit_genotype_cox_for_view(D0,viewName);
            if ~fullFit.success, continue; end
            miceHere = unique(D0.mouse_id);

            for m = 1:numel(miceHere)
                removedMouse = miceHere(m);
                D = D0(D0.mouse_id ~= removedMouse,:);
                fitOne = fit_genotype_cox_for_view(D,viewName);
                removedGenotype = unique(D0.genotype(D0.mouse_id == removedMouse));
                removedCohort = unique(D0.cohort(D0.mouse_id == removedMouse));
                if fitOne.success
                    hr = fitOne.HR; lo = fitOne.CI_low; hi = fitOne.CI_high; p = fitOne.p;
                else
                    hr = NaN; lo = NaN; hi = NaN; p = NaN;
                end
                row = table(viewName,windowName,stateName,removedMouse,removedGenotype(1),removedCohort(1), ...
                    fullFit.HR,fullFit.p,hr,lo,hi,p,hr-fullFit.HR,abs(hr-fullFit.HR), ...
                    'VariableNames',{'cohort_view','time_bin','state','removed_mouse','removed_genotype', ...
                    'removed_cohort','full_HR','full_p','HR_without_mouse','CI_low','CI_high', ...
                    'p_without_mouse','change_in_HR','abs_change_HR'});
                LOO = [LOO; row];
            end
        end
    end
end
writetable(LOO,fullfile(survDir,'baseline_survival_leave_one_mouse_out.csv'));

%% Survival figures
for v = 1:numel(views)
    viewName = views(v);
    V = subset_view(SURV,viewName);
    if viewName == "all", displayName = "Combined cohorts";
    elseif viewName == "thesis", displayName = "Original cohort";
    else, displayName = "New cohort"; end

    fig = figure('Color','w','Position',[100 100 950 430]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    title(tl,displayName + " - whole baseline recording");
    for s = 1:numel(survivalStates)
        ax = nexttile;
        D = V(V.time_bin == "WHOLE" & V.state == survivalStates(s),:);
        plot_survival_panel(ax,D,viewName,"WHOLE",survivalStates(s));
    end
    exportgraphics(fig,fullfile(figDir,sprintf('baseline_survival_whole_%s.png',char(viewName))),'Resolution',300);

    fixedWindows = ["ZT3_9","ZT3_6","ZT6_9","ZT9_12"];
    fig = figure('Color','w','Position',[60 80 1500 720]);
    tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    title(tl,displayName + " - baseline survival by time bin");
    for s = 1:numel(survivalStates)
        for w = 1:numel(fixedWindows)
            ax = nexttile;
            D = V(V.time_bin == fixedWindows(w) & V.state == survivalStates(s),:);
            plot_survival_panel(ax,D,viewName,fixedWindows(w),survivalStates(s));
        end
    end
    exportgraphics(fig,fullfile(figDir,sprintf('baseline_survival_by_bin_%s.png',char(viewName))),'Resolution',300);
end

%% Console summary
fprintf('\n============================================================\n');
fprintf('CLUSTERED COX RESULTS\n');
fprintf('============================================================\n\n');
disp(COX(:,["cohort_view","time_bin","state","n_mice_WT","n_mice_APP", ...
    "n_bouts_WT","n_bouts_APP","HR_APP_vs_WT","CI_low","CI_high","p_robust_t"]));

fprintf('\n============================================================\n');
fprintf('GENOTYPE x COHORT SURVIVAL MODEL\n');
fprintf('============================================================\n\n');
disp(INTERACTION);

fprintf('\n============================================================\n');
fprintf('MOUSE-BALANCED KM MEDIANS\n');
fprintf('============================================================\n\n');
disp(KMSUMMARY);

fprintf('\n============================================================\n');
fprintf('MOST INFLUENTIAL COMBINED-COHORT SURVIVAL MICE\n');
fprintf('============================================================\n');
if ~isempty(LOO)
    for stateName = ["NREM","REM"]
        for windowName = windowNames
            D = LOO(LOO.state == stateName & LOO.time_bin == windowName & LOO.cohort_view == "all",:);
            if isempty(D), continue; end
            D = sortrows(D,'abs_change_HR','descend');
            fprintf('\n%s %s - top 5 by |change in HR|:\n',stateName,windowName);
            disp(D(1:min(5,height(D)),["removed_mouse","removed_genotype","removed_cohort", ...
                "full_HR","HR_without_mouse","p_without_mouse","abs_change_HR"]));
        end
    end
end

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');
fprintf('Outputs:  %s\n',outDir);
fprintf('Figures:  %s\n',figDir);
fprintf('Survival: %s\n\n',survDir);

%% Local functions
function scoreOut = reclassify_short_wake(scoreIn,maxDur)
scoreOut = double(scoreIn(:));
isWake = scoreOut == 0;
d = diff([false;isWake;false]);
starts = find(d == 1); ends_ = find(d == -1)-1;
dur = ends_ - starts + 1;
short = find(dur <= maxDur);
for i = reshape(short,1,[])
    scoreOut(starts(i):ends_(i)) = 15;
end
end

function B = build_bout_table(score,t)
score = double(score(:)); t = double(t(:));
validCodes = [0 1 2 15];
invalid = ~ismember(score,validCodes);
x = score; x(invalid) = -999;
change = [true;x(2:end) ~= x(1:end-1)];
starts = find(change); ends_ = [starts(2:end)-1;numel(x)];
stateCode = x(starts); keep = ismember(stateCode,validCodes);
starts = starts(keep); ends_ = ends_(keep); stateCode = stateCode(keep);
n = numel(starts);
state = strings(n,1); onset = nan(n,1); endExclusive = nan(n,1); duration = nan(n,1);
leftCensored = false(n,1); rightCensored = false(n,1); stateBefore = nan(n,1);
for i = 1:n
    s = starts(i); e = ends_(i);
    switch stateCode(i)
        case 0, state(i) = "WK";
        case 1, state(i) = "NREM";
        case 2, state(i) = "REM";
        case 15, state(i) = "MA";
    end
    onset(i) = t(s); duration(i) = e-s+1; endExclusive(i) = onset(i)+duration(i);
    if s == 1, leftCensored(i) = true;
    elseif invalid(s-1), leftCensored(i) = true;
    else, stateBefore(i) = score(s-1); end
    if e == numel(score), rightCensored(i) = true;
    elseif invalid(e+1), rightCensored(i) = true; end
end
B = table(stateCode,state,onset,endExclusive,duration,leftCensored,rightCensored,stateBefore, ...
    'VariableNames',{'state_code','state','onset_s','end_exclusive_s','duration_s', ...
    'left_censored','right_censored','state_before'});
end

function D = subset_view(T,viewName)
if viewName == "all", D = T; else, D = T(T.cohort == viewName,:); end
end

function row = make_welch_row(viewName,windowName,domain,stateName,metricName,D)
if isempty(D), WT = []; APP = [];
else
    y = double(D.(metricName)); WT = y(D.genotype == "WT"); APP = y(D.genotype == "APP");
end
WT = WT(isfinite(WT)); APP = APP(isfinite(APP));
[p,tstat,df,d] = welch_stats(WT,APP);
row = table(viewName,windowName,domain,stateName,metricName,numel(WT),numel(APP), ...
    mean(WT,'omitnan'),mean(APP,'omitnan'),mean(APP,'omitnan')-mean(WT,'omitnan'),p,tstat,df,d, ...
    'VariableNames',{'cohort_view','time_bin','domain','state','metric','n_WT','n_APP', ...
    'WT_mean','APP_mean','APP_minus_WT','p_welch','t_stat','df','Cohen_d'});
end

function row = make_model_row(windowName,domain,stateName,metricName,D)
if isempty(D), pg = NaN; pc = NaN; pi = NaN;
else, [pg,pc,pi] = genotype_cohort_model(double(D.(metricName)),D.genotype,D.cohort); end
row = table(windowName,domain,stateName,metricName,pg,pc,pi, ...
    'VariableNames',{'time_bin','domain','state','metric','p_genotype','p_cohort','p_genotype_x_cohort'});
end

function [p,tstat,df,d] = welch_stats(WT,APP)
p = NaN; tstat = NaN; df = NaN; d = NaN;
WT = WT(isfinite(WT)); APP = APP(isfinite(APP));
if numel(WT) < 2 || numel(APP) < 2, return; end
[~,p,~,stats] = ttest2(APP,WT,'Vartype','unequal');
tstat = stats.tstat; df = stats.df;
n1 = numel(WT); n2 = numel(APP); s1 = std(WT); s2 = std(APP);
sp = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/max(1,n1+n2-2));
if sp > 0, d = (mean(APP)-mean(WT))/sp; end
end

function [pg,pc,pi] = genotype_cohort_model(y,genotype,cohort)
pg = NaN; pc = NaN; pi = NaN;
valid = isfinite(y) & strlength(genotype)>0 & strlength(cohort)>0;
y = y(valid); genotype = genotype(valid); cohort = cohort(valid);
if numel(y) < 4, return; end
try
    p = anovan(y,{cellstr(genotype),cellstr(cohort)},'model','interaction','sstype',3, ...
        'display','off','varnames',{'Genotype','Cohort'});
    if numel(p)>=3, pg=p(1); pc=p(2); pi=p(3); end
catch
end
end

function plot_metric_windows(M,views,windows,states,metricName,yLabel,outDir,fileStem)
WTcolor = [0.35 0.35 0.35]; APPcolor = [0.20 0.40 0.80];
windowLabels = {'Whole','ZT3-9','ZT3-6','ZT6-9','ZT9-12'};
for v = 1:numel(views)
    viewName = views(v); Dview = subset_view(M,viewName);
    if viewName == "all", label = "Combined cohorts";
    elseif viewName == "thesis", label = "Original cohort";
    else, label = "New cohort"; end
    fig = figure('Color','w','Position',[80 80 1200 800]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl,label + " - " + strrep(metricName,"_"," "));
    for s = 1:numel(states)
        ax = nexttile; hold(ax,'on');
        for w = 1:numel(windows)
            D = Dview(Dview.time_bin == windows(w) & Dview.state == states(s),:);
            WT = double(D.(metricName)(D.genotype == "WT"));
            APP = double(D.(metricName)(D.genotype == "APP"));
            plot_group(ax,w-0.13,WT,WTcolor); plot_group(ax,w+0.13,APP,APPcolor);
        end
        xlim(ax,[0.5 numel(windows)+0.5]); xticks(ax,1:numel(windows)); xticklabels(ax,windowLabels); xtickangle(ax,20);
        ylabel(ax,yLabel); title(ax,states(s)); box(ax,'off');
        if s == 1
            hWT = plot(ax,nan,nan,'o','MarkerFaceColor',WTcolor,'MarkerEdgeColor',WTcolor);
            hAPP = plot(ax,nan,nan,'o','MarkerFaceColor',APPcolor,'MarkerEdgeColor',APPcolor);
            legend(ax,[hWT hAPP],{'WT','APP/PS1'},'Location','best','Box','off');
        end
    end
    exportgraphics(fig,fullfile(outDir,sprintf('%s_%s.png',char(fileStem),char(viewName))),'Resolution',300);
end
end

function plot_fragmentation_windows(F,views,windows,states,outDir)
WTcolor = [0.35 0.35 0.35]; APPcolor = [0.20 0.40 0.80];
windowLabels = {'Whole','ZT3-9','ZT3-6','ZT6-9','ZT9-12'};
for v = 1:numel(views)
    viewName = views(v); Dview = subset_view(F,viewName);
    if viewName == "all", label = "Combined cohorts";
    elseif viewName == "thesis", label = "Original cohort";
    else, label = "New cohort"; end
    fig = figure('Color','w','Position',[100 100 950 430]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); title(tl,label + " - fragmentation");
    for s = 1:numel(states)
        ax = nexttile; hold(ax,'on');
        for w = 1:numel(windows)
            D = Dview(Dview.time_bin == windows(w) & Dview.carrier_state == states(s),:);
            WT = double(D.MA_per_h_state(D.genotype == "WT")); APP = double(D.MA_per_h_state(D.genotype == "APP"));
            plot_group(ax,w-0.13,WT,WTcolor); plot_group(ax,w+0.13,APP,APPcolor);
        end
        xlim(ax,[0.5 numel(windows)+0.5]); xticks(ax,1:numel(windows)); xticklabels(ax,windowLabels); xtickangle(ax,20);
        ylabel(ax,'MA bouts per hour of state'); title(ax,states(s)); box(ax,'off');
    end
    exportgraphics(fig,fullfile(outDir,sprintf('baseline_fragmentation_%s.png',char(viewName))),'Resolution',300);
end
end

function plot_group(ax,xPos,values,color)
values = values(isfinite(values)); n = numel(values); if n == 0, return; end
if n == 1, jitter = 0; else, jitter = linspace(-0.035,0.035,n)'; end
scatter(ax,xPos+jitter,values,34,'MarkerEdgeColor',color,'LineWidth',1.1);
m = mean(values); sem = std(values)/sqrt(n);
errorbar(ax,xPos,m,sem,'o','Color',color,'MarkerFaceColor',color,'MarkerEdgeColor',color,'LineWidth',1.6,'CapSize',8);
end

function w = mouse_equal_weights(mouseID)
mouseID = string(mouseID); mice = unique(mouseID); n = numel(mouseID); nMice = numel(mice); w = zeros(n,1);
for i = 1:nMice
    idx = mouseID == mice(i); nBoutMouse = sum(idx); if nBoutMouse > 0, w(idx) = 1/nBoutMouse; end
end
if nMice > 0, w = w*(n/nMice); end
end

function [tEvent,S,medianSurvival] = weighted_km(time,event,w)
time = double(time(:)); event = logical(event(:)); w = double(w(:));
valid = isfinite(time) & time>0 & isfinite(w) & w>0;
time = time(valid); event = event(valid); w = w(valid);
eventTimes = sort(unique(time(event))); S = nan(numel(eventTimes),1); surv = 1;
for i = 1:numel(eventTimes)
    tt = eventTimes(i); atRisk = time >= tt; eventsNow = time == tt & event;
    riskWeight = sum(w(atRisk)); eventWeight = sum(w(eventsNow));
    if riskWeight > 0, surv = surv*(1-eventWeight/riskWeight); end
    S(i) = surv;
end
tEvent = eventTimes; idx = find(S <= 0.5,1,'first');
if isempty(idx), medianSurvival = NaN; else, medianSurvival = tEvent(idx); end
end

function FIT = clustered_cox(X,time,censored,cluster,predictorNames)
FIT = struct('success',false,'results',table(),'robustCov',[],'nClusters',0);
X = double(X); time = double(time(:)); censored = logical(censored(:)); cluster = string(cluster(:)); predictorNames = string(predictorNames(:));
if isvector(X), X = X(:); end
valid = all(isfinite(X),2) & isfinite(time) & time>0 & strlength(cluster)>0;
X = X(valid,:); time = time(valid); censored = censored(valid); cluster = cluster(valid);
clusters = unique(cluster); G = numel(clusters); p = size(X,2); FIT.nClusters = G;
if G <= p+1 || sum(~censored) < p+1, return; end
w = mouse_equal_weights(cluster);
try
    [beta,~,~,stats] = coxphfit(X,time,'Censoring',censored,'Frequency',w,'Ties','breslow','Baseline',0);
catch ME
    warning('Cox fit failed: %s',ME.message); return;
end
beta = beta(:); scoreResiduals = double(stats.scores); clusterScore = zeros(G,p);
for g = 1:G
    idx = cluster == clusters(g); clusterScore(g,:) = sum(scoreResiduals(idx,:),1,'omitnan');
end
bread = double(stats.covb); meat = clusterScore'*clusterScore; Vrobust = bread*meat*bread;
if G > 1, Vrobust = Vrobust*(G/(G-1)); end
se = sqrt(max(diag(Vrobust),0)); z = beta./se; pZ = 2*(1-normcdf(abs(z))); dfCluster = G-1;
pT = 2*(1-tcdf(abs(z),dfCluster)); tcrit = tinv(0.975,dfCluster); ciLow = exp(beta-tcrit.*se); ciHigh = exp(beta+tcrit.*se); HR = exp(beta);
R = table(predictorNames,beta,se,HR,ciLow,ciHigh,z,pZ,pT,repmat(G,p,1), ...
    'VariableNames',{'predictor','beta','robust_SE','HR','CI_low','CI_high','z_robust','p_robust_z','p_robust_t','n_clusters'});
FIT.success = true; FIT.results = R; FIT.robustCov = Vrobust;
end

function row = make_cox_row(viewName,windowName,stateName,D,R)
WT = D(D.genotype == "WT",:); APP = D(D.genotype == "APP",:);
row = table(viewName,windowName,stateName,numel(unique(WT.mouse_id)),numel(unique(APP.mouse_id)), ...
    height(WT),height(APP),sum(WT.event),sum(APP.event),R.HR,R.CI_low,R.CI_high,R.p_robust_t,R.p_robust_z,R.beta,R.robust_SE, ...
    'VariableNames',{'cohort_view','time_bin','state','n_mice_WT','n_mice_APP','n_bouts_WT','n_bouts_APP', ...
    'n_events_WT','n_events_APP','HR_APP_vs_WT','CI_low','CI_high','p_robust_t','p_robust_z','beta','robust_SE'});
end

function row = km_summary_row(viewName,windowName,stateName,genotype,D,medianSurvival)
row = table(viewName,windowName,stateName,genotype,numel(unique(D.mouse_id)),height(D),sum(D.event),sum(D.censored),medianSurvival, ...
    'VariableNames',{'cohort_view','time_bin','state','genotype','n_mice','n_bouts','n_events','n_censored','KM_median_s'});
end

function FITOUT = fit_genotype_cox_for_view(D,viewName)
FITOUT = struct('success',false,'HR',NaN,'CI_low',NaN,'CI_high',NaN,'p',NaN);
if isempty(D) || numel(unique(D.genotype)) < 2, return; end
if viewName == "all"
    genotype = double(D.genotype == "APP"); cohort = double(D.cohort == "new");
    FIT = clustered_cox([genotype cohort],double(D.time_s),logical(D.censored),D.mouse_id,["genotype_APP","cohort_new"]);
    if ~FIT.success, return; end
    R = FIT.results(FIT.results.predictor == "genotype_APP",:);
else
    genotype = double(D.genotype == "APP");
    FIT = clustered_cox(genotype,double(D.time_s),logical(D.censored),D.mouse_id,"genotype_APP");
    if ~FIT.success, return; end
    R = FIT.results(1,:);
end
FITOUT.success = true; FITOUT.HR = R.HR; FITOUT.CI_low = R.CI_low; FITOUT.CI_high = R.CI_high; FITOUT.p = R.p_robust_t;
end

function plot_survival_panel(ax,D,viewName,windowName,stateName)
cla(ax); hold(ax,'on');
if isempty(D)
    title(ax,sprintf('%s, %s - no usable data',char(stateName),char(windowName))); ylim(ax,[0 1]); box(ax,'off'); return;
end
WTcolor = [0.30 0.30 0.30]; APPcolor = [0.20 0.40 0.80];
WT = D(D.genotype == "WT",:); APP = D(D.genotype == "APP",:);
if ~isempty(WT)
    wWT = mouse_equal_weights(WT.mouse_id); [tWT,sWT,~] = weighted_km(double(WT.time_s),logical(WT.event),wWT);
    stairs(ax,[0;tWT],[1;sWT],'LineWidth',1.8,'Color',WTcolor,'DisplayName',sprintf('WT (%d mice, %d bouts)',numel(unique(WT.mouse_id)),height(WT)));
end
if ~isempty(APP)
    wAPP = mouse_equal_weights(APP.mouse_id); [tAPP,sAPP,~] = weighted_km(double(APP.time_s),logical(APP.event),wAPP);
    stairs(ax,[0;tAPP],[1;sAPP],'LineWidth',1.8,'Color',APPcolor,'DisplayName',sprintf('APP (%d mice, %d bouts)',numel(unique(APP.mouse_id)),height(APP)));
end
fitResult = fit_genotype_cox_for_view(D,viewName);
xlabel(ax,'Bout duration (s)'); ylabel(ax,'Survival probability');
if windowName == "WHOLE", windowLabel = "Whole recording"; else, windowLabel = strrep(windowName,"_","-"); end
title(ax,stateName + ", " + windowLabel); ylim(ax,[0 1]); xlim(ax,[0 survival_xlim(D.time_s)]); box(ax,'off'); legend(ax,'Location','southwest','Box','off');
if fitResult.success
    txt = sprintf('APP vs WT HR = %.2f\n95%% CI %.2f-%.2f\nclustered p = %.3g',fitResult.HR,fitResult.CI_low,fitResult.CI_high,fitResult.p);
    text(ax,0.97,0.95,txt,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',9);
end
end

function xmax = survival_xlim(time)
time = double(time(:)); time = time(isfinite(time) & time>0);
if isempty(time), xmax = 1; return; end
xmax = prctile(time,97.5); if ~isfinite(xmax) || xmax <= 0, xmax = max(time); end; if xmax <= 0, xmax = 1; end
end
