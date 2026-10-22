%% RUN_04_BASELINE_BOUT_METRICS
%
% Baseline sleep microarchitecture using the current rescored data.
%
% Correct analysis windows:
%   ZT3-6 = 0-3 h from recording start
%   ZT6-9 = 3-6 h from recording start
%
% Outputs:
%   baseline_bout_metrics_per_mouse.csv
%   baseline_fragmentation_per_mouse.csv
%   baseline_bouts.csv
%   baseline_bout_metrics_welch.csv
%   baseline_fragmentation_welch.csv
%   baseline_bout_metrics_genotype_cohort_model.csv
%   figures
%
% Bout-number rule:
%   Count observed bout onsets within each bin.
%   Bouts beginning immediately after undefined data are not counted.
%
% Bout-duration rule:
%   Only fully observed bouts contained within the bin are used.
%   Bouts touching undefined data, recording boundaries, or the end
%   of the analysis bin are excluded from duration estimates.
%
% Fragmentation:
%   Number of MA onsets arising from NREM or REM, normalized by
%   hours spent in the corresponding state.

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

standardDir = fullfile( ...
    P.sleep_output, ...
    'all', ...
    'standardized_scores');

outDir = fullfile( ...
    P.sleep_output, ...
    'baseline');

if ~isfolder(outDir)
    mkdir(outDir);
end

%% Baseline sessions

S = Sessions(Sessions.condition == "baseline", :);

%% Definitions

binNames = ["ZT3_6","ZT6_9"];

binStart = [0, 3*3600];
binEnd   = [3*3600, 6*3600];

stateNames = ["WK","NREM","REM","MA"];
stateCodes = [0,1,2,15];

minimumValidCoverage = 95;
maThresholdSec = 15;

%% Collectors

METRICS = table();
FRAG = table();
ALLBOUTS = table();

fprintf('\n========================================\n');
fprintf('BASELINE BOUT METRICS\n');
fprintf('========================================\n');
fprintf('Sessions: %d\n\n', height(S));

%% Process sessions

for i = 1:height(S)

    sid = S.session_id(i);
    mouseID = S.mouse_id(i);

    mouseIdx = find(Mice.mouse_id == mouseID, 1);

    genotype = Mice.genotype(mouseIdx);
    sex = Mice.sex(mouseIdx);
    cohort = S.cohort(i);

    f = fullfile( ...
        standardDir, ...
        char(sid + ".csv"));

    T = readtable(f);

    t = double(T.time_s(:));
    score = double(T.score(:));

    t = t - min(t);

    %% Reclassify short Wake bouts as MA

    scoreMA = reclassify_short_wake( ...
        score, ...
        maThresholdSec);

    %% Build bout table

    B = build_bout_table(scoreMA, t);

    if ~isempty(B)

        B.session_id = repmat(sid,height(B),1);
        B.mouse_id = repmat(mouseID,height(B),1);
        B.genotype = repmat(genotype,height(B),1);
        B.sex = repmat(sex,height(B),1);
        B.cohort = repmat(cohort,height(B),1);

        ALLBOUTS = [ALLBOUTS; B];
    end

    fprintf('[%02d/%02d] %-8s %-4s\n', ...
        i,height(S),mouseID,genotype);

    %% Analysis bins

    for b = 1:numel(binNames)

        bName = binNames(b);
        bStart = binStart(b);
        bEnd = binEnd(b);

        inBin = t >= bStart & t < bEnd;

        expectedSec = bEnd - bStart;

        validMask = ...
            inBin & ismember(scoreMA,[0 1 2 15]);

        validSec = sum(validMask);

        validCoverage = ...
            100 * validSec / expectedSec;

        includePrimary = ...
            validCoverage >= minimumValidCoverage;

        validHours = validSec / 3600;

        %% State bout metrics

        for s = 1:numel(stateNames)

            stateName = stateNames(s);
            stateCode = stateCodes(s);

            stateSec = sum( ...
                inBin & scoreMA == stateCode);

            stateHours = stateSec / 3600;

            if isempty(B)

                BR = B;

            else

                BR = B( ...
                    B.state_code == stateCode & ...
                    B.onset_s >= bStart & ...
                    B.onset_s < bEnd, :);
            end

            % Reliable observed onsets:
            % exclude bouts beginning at recording start or after
            % undefined data.
            if isempty(BR)

                nBouts = 0;
                nDurationBouts = 0;

                meanDuration = NaN;
                medianDuration = NaN;

            else

                reliableStart = ~BR.left_censored;

                nBouts = sum(reliableStart);

                % Duration must be fully observed and entirely
                % contained within this 3-hour bin.
                useDuration = ...
                    ~BR.left_censored & ...
                    ~BR.right_censored & ...
                    BR.end_exclusive_s <= bEnd;

                durations = ...
                    BR.duration_s(useDuration);

                nDurationBouts = numel(durations);

                if isempty(durations)

                    meanDuration = NaN;
                    medianDuration = NaN;

                else

                    meanDuration = mean(durations);
                    medianDuration = median(durations);
                end
            end

            if validHours > 0
                boutsPerHour = nBouts / validHours;
            else
                boutsPerHour = NaN;
            end

            row = table( ...
                sid, ...
                mouseID, ...
                genotype, ...
                sex, ...
                cohort, ...
                bName, ...
                stateName, ...
                stateSec, ...
                stateHours, ...
                validSec, ...
                validCoverage, ...
                includePrimary, ...
                nBouts, ...
                boutsPerHour, ...
                nDurationBouts, ...
                meanDuration, ...
                medianDuration, ...
                'VariableNames', { ...
                    'session_id', ...
                    'mouse_id', ...
                    'genotype', ...
                    'sex', ...
                    'cohort', ...
                    'time_bin', ...
                    'state', ...
                    'state_seconds', ...
                    'state_hours', ...
                    'valid_seconds', ...
                    'valid_coverage_pct', ...
                    'include_primary', ...
                    'n_bouts', ...
                    'bouts_per_h', ...
                    'n_duration_bouts', ...
                    'mean_bout_duration_s', ...
                    'median_bout_duration_s'});

            METRICS = [METRICS; row];

        end

        %% Fragmentation index

        if isempty(B)

            maNREM = 0;
            maREM = 0;

        else

            MA = B( ...
                B.state_code == 15 & ...
                B.onset_s >= bStart & ...
                B.onset_s < bEnd & ...
                ~B.left_censored, :);

            maNREM = sum(MA.state_before == 1);
            maREM  = sum(MA.state_before == 2);
        end

        nremSec = sum( ...
            inBin & scoreMA == 1);

        remSec = sum( ...
            inBin & scoreMA == 2);

        nremHours = nremSec / 3600;
        remHours = remSec / 3600;

        if nremHours > 0
            nremFrag = maNREM / nremHours;
        else
            nremFrag = NaN;
        end

        if remHours > 0
            remFrag = maREM / remHours;
        else
            remFrag = NaN;
        end

        rowNREM = table( ...
            sid,mouseID,genotype,sex,cohort,bName, ...
            string("NREM"), ...
            maNREM,nremSec,nremHours,nremFrag, ...
            validCoverage,includePrimary, ...
            'VariableNames', { ...
                'session_id','mouse_id','genotype','sex', ...
                'cohort','time_bin','carrier_state', ...
                'n_MA','state_seconds','state_hours', ...
                'MA_per_h_state','valid_coverage_pct', ...
                'include_primary'});

        rowREM = table( ...
            sid,mouseID,genotype,sex,cohort,bName, ...
            string("REM"), ...
            maREM,remSec,remHours,remFrag, ...
            validCoverage,includePrimary, ...
            'VariableNames', { ...
                'session_id','mouse_id','genotype','sex', ...
                'cohort','time_bin','carrier_state', ...
                'n_MA','state_seconds','state_hours', ...
                'MA_per_h_state','valid_coverage_pct', ...
                'include_primary'});

        FRAG = [FRAG; rowNREM; rowREM];

    end
end

%% Save raw outputs

metricsFile = fullfile( ...
    outDir, ...
    'baseline_bout_metrics_per_mouse.csv');

fragFile = fullfile( ...
    outDir, ...
    'baseline_fragmentation_per_mouse.csv');

boutsFile = fullfile( ...
    outDir, ...
    'baseline_bouts.csv');

writetable(METRICS,metricsFile);
writetable(FRAG,fragFile);
writetable(ALLBOUTS,boutsFile);

%% Primary QC subset

M = METRICS(METRICS.include_primary,:);
F = FRAG(FRAG.include_primary,:);

%% Welch statistics

views = ["thesis","new","all"];

metricNames = [ ...
    "bouts_per_h", ...
    "mean_bout_duration_s", ...
    "median_bout_duration_s"];

STATS = table();

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        Dview = M;
    else
        Dview = M(M.cohort == viewName,:);
    end

    for b = 1:numel(binNames)

        for s = 1:numel(stateNames)

            D = Dview( ...
                Dview.time_bin == binNames(b) & ...
                Dview.state == stateNames(s), :);

            for m = 1:numel(metricNames)

                metricName = metricNames(m);

                y = double(D.(metricName));

                WT = y(D.genotype == "WT");
                APP = y(D.genotype == "APP");

                WT = WT(~isnan(WT));
                APP = APP(~isnan(APP));

                [p,tstat,df,d] = ...
                    welch_stats(WT,APP);

                row = table( ...
                    viewName, ...
                    binNames(b), ...
                    stateNames(s), ...
                    metricName, ...
                    numel(WT), ...
                    numel(APP), ...
                    mean(WT,'omitnan'), ...
                    mean(APP,'omitnan'), ...
                    mean(APP,'omitnan') - mean(WT,'omitnan'), ...
                    p,tstat,df,d, ...
                    'VariableNames', { ...
                        'cohort_view','time_bin','state','metric', ...
                        'n_WT','n_APP','WT_mean','APP_mean', ...
                        'APP_minus_WT','p_welch','t_stat','df', ...
                        'Cohen_d'});

                STATS = [STATS; row];

            end
        end
    end
end

statsFile = fullfile( ...
    outDir, ...
    'baseline_bout_metrics_welch.csv');

writetable(STATS,statsFile);

%% Fragmentation statistics

FRAGSTATS = table();

carrierStates = ["NREM","REM"];

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        Dview = F;
    else
        Dview = F(F.cohort == viewName,:);
    end

    for b = 1:numel(binNames)

        for s = 1:numel(carrierStates)

            D = Dview( ...
                Dview.time_bin == binNames(b) & ...
                Dview.carrier_state == carrierStates(s), :);

            WT = double( ...
                D.MA_per_h_state(D.genotype == "WT"));

            APP = double( ...
                D.MA_per_h_state(D.genotype == "APP"));

            WT = WT(~isnan(WT));
            APP = APP(~isnan(APP));

            [p,tstat,df,d] = ...
                welch_stats(WT,APP);

            row = table( ...
                viewName, ...
                binNames(b), ...
                carrierStates(s), ...
                numel(WT), ...
                numel(APP), ...
                mean(WT,'omitnan'), ...
                mean(APP,'omitnan'), ...
                mean(APP,'omitnan') - mean(WT,'omitnan'), ...
                p,tstat,df,d, ...
                'VariableNames', { ...
                    'cohort_view','time_bin','carrier_state', ...
                    'n_WT','n_APP','WT_mean','APP_mean', ...
                    'APP_minus_WT','p_welch','t_stat','df', ...
                    'Cohen_d'});

            FRAGSTATS = [FRAGSTATS; row];

        end
    end
end

fragStatsFile = fullfile( ...
    outDir, ...
    'baseline_fragmentation_welch.csv');

writetable(FRAGSTATS,fragStatsFile);

%% Combined genotype x cohort models

MODEL = table();

modelMetrics = [ ...
    "bouts_per_h", ...
    "mean_bout_duration_s"];

for b = 1:numel(binNames)

    for s = 1:numel(stateNames)

        D = M( ...
            M.time_bin == binNames(b) & ...
            M.state == stateNames(s), :);

        for m = 1:numel(modelMetrics)

            metricName = modelMetrics(m);

            y = double(D.(metricName));

            [pg,pc,pi] = ...
                genotype_cohort_model( ...
                    y, ...
                    D.genotype, ...
                    D.cohort);

            row = table( ...
                binNames(b), ...
                stateNames(s), ...
                metricName, ...
                pg,pc,pi, ...
                'VariableNames', { ...
                    'time_bin','state','metric', ...
                    'p_genotype','p_cohort', ...
                    'p_genotype_x_cohort'});

            MODEL = [MODEL; row];

        end
    end
end

for b = 1:numel(binNames)

    for s = 1:numel(carrierStates)

        D = F( ...
            F.time_bin == binNames(b) & ...
            F.carrier_state == carrierStates(s), :);

        y = double(D.MA_per_h_state);

        [pg,pc,pi] = ...
            genotype_cohort_model( ...
                y,D.genotype,D.cohort);

        row = table( ...
            binNames(b), ...
            carrierStates(s), ...
            string("fragmentation_MA_per_h"), ...
            pg,pc,pi, ...
            'VariableNames', { ...
                'time_bin','state','metric', ...
                'p_genotype','p_cohort', ...
                'p_genotype_x_cohort'});

        MODEL = [MODEL; row];

    end
end

modelFile = fullfile( ...
    outDir, ...
    'baseline_bout_metrics_genotype_cohort_model.csv');

writetable(MODEL,modelFile);

%% Figures

plot_metric_by_view( ...
    M,views,binNames,stateNames, ...
    "mean_bout_duration_s", ...
    "Mean bout duration (s)", ...
    outDir, ...
    "baseline_mean_bout_duration");

plot_metric_by_view( ...
    M,views,binNames,stateNames, ...
    "bouts_per_h", ...
    "Bout onsets per recording hour", ...
    outDir, ...
    "baseline_bout_rate");

plot_fragmentation( ...
    F,views,binNames,carrierStates,outDir);

%% Console summary

fprintf('\n========================================\n');
fprintf('BASELINE BOUT ANALYSIS COMPLETE\n');
fprintf('========================================\n\n');

fprintf('Mean bout duration - NREM and REM:\n\n');

disp(STATS( ...
    STATS.metric == "mean_bout_duration_s" & ...
    ismember(STATS.state,["NREM","REM"]), ...
    ["cohort_view","time_bin","state", ...
     "n_WT","n_APP","WT_mean","APP_mean", ...
     "APP_minus_WT","p_welch","Cohen_d"]));

fprintf('\nFragmentation:\n\n');

disp(FRAGSTATS(:, ...
    ["cohort_view","time_bin","carrier_state", ...
     "n_WT","n_APP","WT_mean","APP_mean", ...
     "APP_minus_WT","p_welch","Cohen_d"]));

fprintf('\nCombined genotype x cohort models:\n\n');

disp(MODEL);

fprintf('\nOutputs saved to:\n%s\n',outDir);

%% Local functions

function scoreOut = reclassify_short_wake(scoreIn,maxDur)

scoreOut = double(scoreIn(:));

isWake = scoreOut == 0;

d = diff([false; isWake; false]);

starts = find(d == 1);
ends_ = find(d == -1) - 1;

dur = ends_ - starts + 1;

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

    duration(i) = e - s + 1;

    endExclusive(i) = ...
        onset(i) + duration(i);

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
    stateCode, ...
    state, ...
    onset, ...
    endExclusive, ...
    duration, ...
    leftCensored, ...
    rightCensored, ...
    stateBefore, ...
    'VariableNames', { ...
        'state_code','state','onset_s', ...
        'end_exclusive_s','duration_s', ...
        'left_censored','right_censored', ...
        'state_before'});

end


function [p,tstat,df,d] = welch_stats(WT,APP)

p = NaN;
tstat = NaN;
df = NaN;
d = NaN;

WT = WT(~isnan(WT));
APP = APP(~isnan(APP));

if numel(WT) < 2 || numel(APP) < 2
    return;
end

[~,p,~,stats] = ...
    ttest2(APP,WT,'Vartype','unequal');

tstat = stats.tstat;
df = stats.df;

n1 = numel(WT);
n2 = numel(APP);

s1 = std(WT);
s2 = std(APP);

sp = sqrt( ...
    ((n1-1)*s1^2 + (n2-1)*s2^2) / ...
    max(1,n1+n2-2));

if sp > 0
    d = (mean(APP)-mean(WT))/sp;
end

end


function [pg,pc,pi] = genotype_cohort_model(y,genotype,cohort)

pg = NaN;
pc = NaN;
pi = NaN;

valid = ...
    ~isnan(y) & ...
    strlength(genotype) > 0 & ...
    strlength(cohort) > 0;

y = y(valid);
genotype = genotype(valid);
cohort = cohort(valid);

try

    p = anovan( ...
        y, ...
        {cellstr(genotype),cellstr(cohort)}, ...
        'model','interaction', ...
        'sstype',3, ...
        'display','off', ...
        'varnames',{'Genotype','Cohort'});

    if numel(p) >= 3
        pg = p(1);
        pc = p(2);
        pi = p(3);
    end

catch
end

end


function plot_metric_by_view( ...
    M,views,bins,states,metricName,yLabel,outDir,fileStem)

WTcolor = [0.35 0.35 0.35];
APPcolor = [0.20 0.40 0.80];

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        Dview = M;
        label = "Combined cohorts";
    elseif viewName == "thesis"
        Dview = M(M.cohort == "thesis",:);
        label = "Original cohort";
    else
        Dview = M(M.cohort == "new",:);
        label = "New cohort";
    end

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 1100 800]);

    tl = tiledlayout(2,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    title(tl,label);

    for s = 1:numel(states)

        ax = nexttile;
        hold(ax,'on');

        for b = 1:numel(bins)

            D = Dview( ...
                Dview.time_bin == bins(b) & ...
                Dview.state == states(s), :);

            WT = double( ...
                D.(metricName)(D.genotype == "WT"));

            APP = double( ...
                D.(metricName)(D.genotype == "APP"));

            plot_group(ax,b-0.13,WT,WTcolor);
            plot_group(ax,b+0.13,APP,APPcolor);
        end

        xlim(ax,[0.5 2.5]);

        xticks(ax,[1 2]);
        xticklabels(ax,{'ZT3-6','ZT6-9'});

        ylabel(ax,yLabel);
        title(ax,states(s));

        box(ax,'off');

    end

    pngFile = fullfile( ...
        outDir, ...
        sprintf('%s_%s.png', ...
        char(fileStem),char(viewName)));

    exportgraphics(fig,pngFile,'Resolution',300);

end

end


function plot_fragmentation(F,views,bins,states,outDir)

WTcolor = [0.35 0.35 0.35];
APPcolor = [0.20 0.40 0.80];

for v = 1:numel(views)

    viewName = views(v);

    if viewName == "all"
        Dview = F;
        label = "Combined cohorts";
    elseif viewName == "thesis"
        Dview = F(F.cohort == "thesis",:);
        label = "Original cohort";
    else
        Dview = F(F.cohort == "new",:);
        label = "New cohort";
    end

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 900 420]);

    tl = tiledlayout(1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    title(tl,label + " - fragmentation");

    for s = 1:numel(states)

        ax = nexttile;
        hold(ax,'on');

        for b = 1:numel(bins)

            D = Dview( ...
                Dview.time_bin == bins(b) & ...
                Dview.carrier_state == states(s), :);

            WT = double( ...
                D.MA_per_h_state(D.genotype == "WT"));

            APP = double( ...
                D.MA_per_h_state(D.genotype == "APP"));

            plot_group(ax,b-0.13,WT,WTcolor);
            plot_group(ax,b+0.13,APP,APPcolor);
        end

        xlim(ax,[0.5 2.5]);

        xticks(ax,[1 2]);
        xticklabels(ax,{'ZT3-6','ZT6-9'});

        ylabel(ax,'MA bouts per hour of state');
        title(ax,states(s));

        box(ax,'off');

    end

    pngFile = fullfile( ...
        outDir, ...
        sprintf('baseline_fragmentation_%s.png', ...
        char(viewName)));

    exportgraphics(fig,pngFile,'Resolution',300);

end

end


function plot_group(ax,xPos,values,color)

values = values(~isnan(values));

n = numel(values);

if n == 0
    return;
end

if n == 1
    jitter = 0;
else
    jitter = linspace(-0.035,0.035,n)';
end

scatter( ...
    ax, ...
    xPos+jitter, ...
    values, ...
    34, ...
    'MarkerEdgeColor',color, ...
    'LineWidth',1.1);

m = mean(values);
sem = std(values)/sqrt(n);

errorbar( ...
    ax,xPos,m,sem,'o', ...
    'Color',color, ...
    'MarkerFaceColor',color, ...
    'MarkerEdgeColor',color, ...
    'LineWidth',1.6, ...
    'CapSize',8);

end
