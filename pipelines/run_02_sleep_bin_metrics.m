%% RUN_02_SLEEP_BIN_METRICS
%
% Compute sleep-state percentages within 3-hour bins.
%
% Bins:
%   ZT3-6   = first 3 h
%   ZT6-9   = second 3 h
%   ZT9-12  = third 3 h
%
% Output is one row per:
% mouse x condition x time bin x state
%
% Bins with <95% valid coverage are flagged and excluded from
% the primary analysis, but are retained in the output table.

clear;
clc;

%% Setup

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir, 'config'));

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

outDir = fullfile(P.sleep_output, 'all');

%% Definitions

binNames = ["ZT3_6","ZT6_9","ZT9_12"];

binStart = [
    0
    3*3600
    6*3600
];

binEnd = [
    3*3600
    6*3600
    9*3600
];

stateNames = ["WK","NREM","REM","MA"];
stateCodes = [0,1,2,15];

minimumValidCoverage = 95;

RESULTS = table();

%% Sessions

for i = 1:height(Sessions)

    sessionID = Sessions.session_id(i);
    mouseID = Sessions.mouse_id(i);

    mouseIdx = find(Mice.mouse_id == mouseID, 1);

    genotype = Mice.genotype(mouseIdx);
    sex = Mice.sex(mouseIdx);

    scoreFile = fullfile( ...
        standardDir, ...
        char(sessionID + ".csv"));

    T = readtable(scoreFile);

    t = double(T.time_s(:));
    score = double(T.score(:));

    % Set beginning of each recording to zero
    t = t - min(t);

    %% Reclassify short Wake bouts as micro-arousals

    scoreMA = score;

    validIndices = find(~isnan(scoreMA));

    if ~isempty(validIndices)

        changePoint = [true; diff(scoreMA) ~= 0];

        runStart = find(changePoint);
        runEnd = [runStart(2:end)-1; numel(scoreMA)];

        for r = 1:numel(runStart)

            idx1 = runStart(r);
            idx2 = runEnd(r);

            runCode = scoreMA(idx1);

            if runCode == 0

                runDuration = ...
                    t(idx2) - t(idx1) + 1;

                if runDuration <= 15
                    scoreMA(idx1:idx2) = 15;
                end

            end
        end
    end

    %% Three-hour bins

    for b = 1:numel(binNames)

        inBin = ...
            t >= binStart(b) & ...
            t < binEnd(b);

        expectedSeconds = ...
            binEnd(b) - binStart(b);

        nAvailable = sum(inBin);

        validMask = ...
            inBin & ismember(scoreMA,[0 1 2 15]);

        nValid = sum(validMask);

        recordingCoverage = ...
            100 * nAvailable / expectedSeconds;

        validCoverage = ...
            100 * nValid / expectedSeconds;

        includePrimary = ...
            validCoverage >= minimumValidCoverage;

        %% State percentages

        for s = 1:numel(stateNames)

            nState = sum( ...
                inBin & ...
                scoreMA == stateCodes(s));

            if nValid > 0
                percentTime = ...
                    100 * nState / nValid;
            else
                percentTime = NaN;
            end

            row = table( ...
                sessionID, ...
                mouseID, ...
                genotype, ...
                sex, ...
                Sessions.cohort(i), ...
                Sessions.condition(i), ...
                binNames(b), ...
                stateNames(s), ...
                nState, ...
                nValid, ...
                recordingCoverage, ...
                validCoverage, ...
                includePrimary, ...
                percentTime, ...
                'VariableNames', { ...
                    'session_id', ...
                    'mouse_id', ...
                    'genotype', ...
                    'sex', ...
                    'cohort', ...
                    'condition', ...
                    'time_bin', ...
                    'state', ...
                    'state_seconds', ...
                    'valid_seconds', ...
                    'recording_coverage_pct', ...
                    'valid_coverage_pct', ...
                    'include_primary', ...
                    'percent_time'});

            RESULTS = [RESULTS; row];

        end
    end
end

%% Save

outFile = fullfile( ...
    outDir, ...
    'sleep_state_3h_bins.csv');

writetable(RESULTS, outFile);

%% Summary

fprintf('\n========================================\n');
fprintf('3-HOUR SLEEP METRICS COMPLETE\n');
fprintf('========================================\n');

fprintf('Rows: %d\n', height(RESULTS));
fprintf('Output:\n%s\n\n', outFile);

fprintf('Usable sessions by time bin:\n');

for b = 1:numel(binNames)

    tmp = RESULTS( ...
        RESULTS.time_bin == binNames(b) & ...
        RESULTS.state == "REM", :);

    fprintf('%-8s : %d/%d\n', ...
        binNames(b), ...
        sum(tmp.include_primary), ...
        height(tmp));
end
