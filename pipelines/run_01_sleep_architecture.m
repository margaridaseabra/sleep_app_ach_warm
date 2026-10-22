%% RUN_01_SLEEP_ARCHITECTURE
%
% Sleep architecture analysis for thesis + new cohorts.
%
% Supported inputs:
%
% legacy_mat
%   old scored MAT files containing sleep_scores
%
% dash_csv
%   new scoring CSV files containing t0_s and final_code
%
% State codes:
%   -1 = undefined
%    0 = Wake
%    1 = NREM
%    2 = REM
%   15 = micro-arousal
%
% Short Wake bouts <= 15 s are reclassified as MA using the
% existing sleep_architecture_from_scores.m implementation.

clear;
clc;

%% Analysis selection

cohortToRun = "all";
conditionToRun = "all";

%% Setup

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir, 'config'));
addpath(fullfile(repoDir, 'Sleep Architecture Analysis'));

P = paths_config();

Mice = readtable(P.mice, ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

Sessions = readtable(P.sessions, ...
    'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

%% Select sessions

includeSleep = to_logical(Sessions.include_sleep);
Sessions = Sessions(includeSleep, :);

if cohortToRun ~= "all"
    Sessions = Sessions(Sessions.cohort == cohortToRun, :);
end

if conditionToRun ~= "all"
    Sessions = Sessions(Sessions.condition == conditionToRun, :);
end

if isempty(Sessions)
    error('No sessions match the requested analysis selection.');
end

%% Output directory

analysisName = cohortToRun;

if conditionToRun ~= "all"
    analysisName = analysisName + "_" + conditionToRun;
end

outDir = fullfile(P.sleep_output, char(analysisName));

if ~isfolder(outDir)
    mkdir(outDir);
end

standardDir = fullfile(outDir, 'standardized_scores');

if ~isfolder(standardDir)
    mkdir(standardDir);
end

%% Sleep definitions

codes = struct( ...
    'WK', 0, ...
    'NREM', 1, ...
    'REM', 2, ...
    'MA', 15);

maThresholdSec = 15;

%% Collectors

OverallAll = table();
PerHourAll = table();
BoutsAll = table();
QC = table();

fprintf('\n========================================\n');
fprintf('SLEEP ARCHITECTURE ANALYSIS\n');
fprintf('========================================\n');
fprintf('Cohort: %s\n', cohortToRun);
fprintf('Condition: %s\n', conditionToRun);
fprintf('Sessions: %d\n\n', height(Sessions));

%% Process recordings

for i = 1:height(Sessions)

    sessionID = Sessions.session_id(i);
    mouseID = Sessions.mouse_id(i);
    condition = Sessions.condition(i);
    cohort = Sessions.cohort(i);
    scoreFormat = Sessions.score_format(i);

    mouseRow = find(Mice.mouse_id == mouseID, 1);

    if isempty(mouseRow)
        error('No metadata found for %s', mouseID);
    end

    genotype = Mice.genotype(mouseRow);
    sex = Mice.sex(mouseRow);

    sourceFile = fullfile( ...
        P.data_root, ...
        char(Sessions.score_file(i)));

    fprintf('[%02d/%02d] %-8s %-4s %-10s ', ...
        i, height(Sessions), mouseID, genotype, condition);

    try

        %% Load and standardize scoring

        [time_s, score] = load_scores( ...
            sourceFile, ...
            scoreFormat);

        nEpochs = numel(score);

        if numel(time_s) ~= nEpochs
            error('Time and score vectors have different lengths.');
        end

        undefinedMask = ...
            isnan(score) | score == -1;

        nUndefined = sum(undefinedMask);

        knownMask = ismember(score, [-1 0 1 2]);

        nUnexpected = sum(~knownMask & ~isnan(score));

        %% Save standardized score file

        Standard = table( ...
            double(time_s(:)), ...
            double(score(:)), ...
            'VariableNames', {'time_s','score'});

        standardFile = fullfile( ...
            standardDir, ...
            char(sessionID + ".csv"));

        writetable(Standard, standardFile);

        %% Run existing scientific core

        R = sleep_architecture_from_scores( ...
            standardFile, ...
            'codes', codes, ...
            'includeMA', true, ...
            'ma_thresh_sec', maThresholdSec, ...
            'reclassify_short_wake_to_MA', true, ...
            'out_prefix', '', ...
            'write_log', false, ...
            'verbose', false);

        %% Overall metrics

        O = R.overall;

        validSeconds = sum(O.total_dur_s);

        if validSeconds > 0
            O.percent_time = ...
                100 .* O.total_dur_s ./ validSeconds;
        else
            O.percent_time = ...
                NaN(height(O),1);
        end

        O.session_id = repmat(sessionID, height(O), 1);
        O.mouse_id = repmat(mouseID, height(O), 1);
        O.genotype = repmat(genotype, height(O), 1);
        O.sex = repmat(sex, height(O), 1);
        O.cohort = repmat(cohort, height(O), 1);
        O.condition = repmat(condition, height(O), 1);
        O.date = repmat(Sessions.date(i), height(O), 1);
        O.score_format = repmat(scoreFormat, height(O), 1);

        OverallAll = [OverallAll; O];

        %% Per-hour metrics

        H = R.per_hour;

        H.session_id = repmat(sessionID, height(H), 1);
        H.mouse_id = repmat(mouseID, height(H), 1);
        H.genotype = repmat(genotype, height(H), 1);
        H.sex = repmat(sex, height(H), 1);
        H.cohort = repmat(cohort, height(H), 1);
        H.condition = repmat(condition, height(H), 1);
        H.date = repmat(Sessions.date(i), height(H), 1);

        PerHourAll = [PerHourAll; H];

        %% Individual bouts

        B = R.runs;

        B.state = strings(height(B),1);

        B.state(B.state_code == 0) = "WK";
        B.state(B.state_code == 1) = "NREM";
        B.state(B.state_code == 2) = "REM";
        B.state(B.state_code == 15) = "MA";
        B.state(B.state_code == -1) = "UNDEFINED";

        unknownState = B.state == "";

        if any(unknownState)
            B.state(unknownState) = "OTHER";
        end

        B.session_id = repmat(sessionID, height(B), 1);
        B.mouse_id = repmat(mouseID, height(B), 1);
        B.genotype = repmat(genotype, height(B), 1);
        B.sex = repmat(sex, height(B), 1);
        B.cohort = repmat(cohort, height(B), 1);
        B.condition = repmat(condition, height(B), 1);
        B.date = repmat(Sessions.date(i), height(B), 1);

        BoutsAll = [BoutsAll; B];

        %% QC summary

        validPercent = ...
            100 * (nEpochs - nUndefined) / nEpochs;

        q = table( ...
            sessionID, ...
            mouseID, ...
            genotype, ...
            sex, ...
            cohort, ...
            condition, ...
            scoreFormat, ...
            nEpochs, ...
            nUndefined, ...
            nUnexpected, ...
            validPercent, ...
            string("OK"), ...
            string(""), ...
            'VariableNames', { ...
                'session_id', ...
                'mouse_id', ...
                'genotype', ...
                'sex', ...
                'cohort', ...
                'condition', ...
                'score_format', ...
                'n_epochs', ...
                'n_undefined', ...
                'n_unexpected_codes', ...
                'valid_percent', ...
                'status', ...
                'error_message'});

        QC = [QC; q];

        fprintf('OK\n');

    catch ME

        fprintf('FAILED\n');

        q = table( ...
            sessionID, ...
            mouseID, ...
            genotype, ...
            sex, ...
            cohort, ...
            condition, ...
            scoreFormat, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            string("FAILED"), ...
            string(ME.message), ...
            'VariableNames', { ...
                'session_id', ...
                'mouse_id', ...
                'genotype', ...
                'sex', ...
                'cohort', ...
                'condition', ...
                'score_format', ...
                'n_epochs', ...
                'n_undefined', ...
                'n_unexpected_codes', ...
                'valid_percent', ...
                'status', ...
                'error_message'});

        QC = [QC; q];

        warning('%s failed: %s', sessionID, ME.message);

    end
end

%% Save outputs

overallFile = fullfile(outDir, 'session_overall.csv');
perHourFile = fullfile(outDir, 'session_per_hour.csv');
boutsFile = fullfile(outDir, 'bouts.csv');
qcFile = fullfile(outDir, 'qc_summary.csv');

if ~isempty(OverallAll)
    writetable(OverallAll, overallFile);
end

if ~isempty(PerHourAll)
    writetable(PerHourAll, perHourFile);
end

if ~isempty(BoutsAll)
    writetable(BoutsAll, boutsFile);
end

writetable(QC, qcFile);

%% Final summary

nSuccess = sum(QC.status == "OK");
nFailed = sum(QC.status == "FAILED");

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Successful sessions: %d\n', nSuccess);
fprintf('Failed sessions: %d\n', nFailed);

fprintf('\nOutput directory:\n%s\n', outDir);

fprintf('\nFiles:\n');
fprintf('%s\n', overallFile);
fprintf('%s\n', perHourFile);
fprintf('%s\n', boutsFile);
fprintf('%s\n', qcFile);

if nFailed > 0
    fprintf('\nFailed sessions:\n');
    disp(QC(QC.status == "FAILED", ...
        ["session_id","error_message"]));
end

%% Local functions

function [time_s, score] = load_scores(filePath, scoreFormat)

scoreFormat = string(scoreFormat);

switch scoreFormat

    case "legacy_mat"

        D = load(filePath, 'sleep_scores');

        if ~isfield(D, 'sleep_scores')
            error('MAT file does not contain sleep_scores.');
        end

        score = double(D.sleep_scores(:));

        time_s = (0:numel(score)-1)';

    case "dash_csv"

        T = readtable(filePath, ...
            'TextType', 'string', ...
            'VariableNamingRule', 'preserve');

        required = ["t0_s","final_code"];

        if ~all(ismember(required, ...
                string(T.Properties.VariableNames)))
            error('Dash CSV requires t0_s and final_code.');
        end

        time_s = double(T.t0_s(:));
        score = double(T.final_code(:));

    otherwise

        error('Unknown score format: %s', scoreFormat);

end

end


function tf = to_logical(x)

if islogical(x)

    tf = x;

elseif isnumeric(x)

    tf = x ~= 0;

else

    x = lower(strtrim(string(x)));

    tf = ismember(x, ...
        ["1","true","yes","y"]);

end

end
