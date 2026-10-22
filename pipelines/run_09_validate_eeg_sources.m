%% RUN_09_VALIDATE_EEG_SOURCES
%
% Validate that every sleep session used for spindle/PAC analysis has a raw
% EEG recording that can be opened by the future microarchitecture pipeline.
%
% This script DOES NOT analyze EEG. It resolves raw recording paths and
% reports which sessions expose EEG + sampling-rate variables.
%
% It never copies raw data into the repository.
%
% Output:
%   outputs/eeg/eeg_source_validation.csv
%
% Run:
%   cd('/Users/margaridaseabra/sleep_app_ach_warm')
%   run('pipelines/run_09_validate_eeg_sources.m')

clear;
clc;

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);
addpath(fullfile(repoDir,'config'));

P = paths_config();

S = readtable(P.sessions, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

S = S(logical(S.include_sleep),:);

outDir = P.eeg_output;
if ~isfolder(outDir), mkdir(outDir); end

OUT = table();

fprintf('\n============================================================\n');
fprintf('EEG SOURCE VALIDATION FOR SPINDLES / PAC\n');
fprintf('============================================================\n\n');

for i = 1:height(S)

    sid = string(S.session_id(i));
    mouseID = string(S.mouse_id(i));
    condition = string(S.condition(i));
    cohort = string(S.cohort(i));

    [relativePath,resolveStatus] = resolve_recording_path(S(i,:),P);

    if strlength(relativePath) == 0
        fullPath = "";
        existsFlag = false;
    else
        fullPath = string(fullfile(P.data_root,char(relativePath)));
        existsFlag = isfile(fullPath);
    end

    hasEEG = false;
    hasFs = false;
    eegVar = "";
    fsVar = "";
    loaderStatus = "missing-file";
    varList = "";

    if existsFlag

        info = whos('-file',char(fullPath));
        names = string({info.name});
        varList = strjoin(names,";");

        if ismember("eeg",names)
            hasEEG = true;
            eegVar = "eeg";
        end

        if ismember("eeg_frequency",names)
            hasFs = true;
            fsVar = "eeg_frequency";
        end

        if hasEEG && hasFs
            loaderStatus = "ready-top-level";
        elseif ismember("recording",names)
            loaderStatus = "recording-struct-needs-inspection";

            % Inspect only the small struct metadata if possible.
            try
                Q = load(char(fullPath),'recording');
                if isstruct(Q.recording)
                    f = string(fieldnames(Q.recording));

                    eegCandidates = ["eeg","EEG","eeg_data","EEG_data"];
                    fsCandidates = ["eeg_frequency","fs_eeg","eeg_fs","Fs","fs"];

                    q = eegCandidates(ismember(eegCandidates,f));
                    if ~isempty(q)
                        hasEEG = true;
                        eegVar = "recording." + q(1);
                    end

                    q = fsCandidates(ismember(fsCandidates,f));
                    if ~isempty(q)
                        hasFs = true;
                        fsVar = "recording." + q(1);
                    end

                    if hasEEG && hasFs
                        loaderStatus = "ready-recording-struct";
                    end
                end
            catch
                % Keep the inspection status. The full pipeline will not
                % silently guess if the structure is unsupported.
            end
        else
            loaderStatus = "no-recognized-eeg-fields";
        end
    end

    ready = existsFlag && hasEEG && hasFs;

    OUT = [OUT; table( ...
        sid,mouseID,condition,cohort, ...
        relativePath,fullPath,string(resolveStatus), ...
        existsFlag,hasEEG,hasFs,eegVar,fsVar,string(loaderStatus),ready,varList, ...
        'VariableNames',{ ...
            'session_id','mouse_id','condition','cohort', ...
            'recording_file_relative','recording_file_full','resolve_status', ...
            'file_exists','has_eeg','has_eeg_frequency', ...
            'eeg_variable','fs_variable','loader_status','ready_for_eeg_analysis', ...
            'top_level_variables'})]; %#ok<AGROW>

    fprintf('[%02d/%02d] %-28s  %-10s  %s\n', ...
        i,height(S),sid,condition,loaderStatus);
end

outFile = fullfile(outDir,'eeg_source_validation.csv');
writetable(OUT,outFile);

fprintf('\n------------------------------------------------------------\n');
fprintf('Ready:   %d/%d sessions\n',sum(OUT.ready_for_eeg_analysis),height(OUT));
fprintf('Missing: %d\n',sum(~OUT.file_exists));
fprintf('Needs loader inspection: %d\n', ...
    sum(OUT.file_exists & ~OUT.ready_for_eeg_analysis));
fprintf('\nOutput:\n%s\n',outFile);
fprintf('============================================================\n\n');

%% ========================================================================
% Local functions
% =========================================================================

function [relativePath,status] = resolve_recording_path(row,P)

relativePath = "";
status = "unresolved";

vars = string(row.Properties.VariableNames);

% 1) Explicit manifest recording path is always preferred.
if ismember("recording_file",vars)
    q = string(row.recording_file(1));

    if ~ismissing(q) && strlength(strtrim(q)) > 0
        relativePath = strtrim(q);
        status = "manifest-recording_file";
        return;
    end
end

% 2) For the June cohort, the companion MAT commonly has exactly the same
%    basename as the score CSV, e.g.
%
%       290526-m59-bas-1_final_scoring_dash.csv
%       290526-m59-bas-1_final_scoring_dash.mat
%
%    Try that exact companion first.
if ~ismember("score_file",vars)
    return;
end

scoreRel = string(row.score_file(1));

if ismissing(scoreRel) || strlength(strtrim(scoreRel)) == 0
    return;
end

[folder,name,ext] = fileparts(char(scoreRel)); %#ok<ASGLU>

sameBaseCandidate = string(fullfile(folder,[name '.mat']));

if isfile(fullfile(P.data_root,char(sameBaseCandidate)))
    relativePath = sameBaseCandidate;
    status = "same-basename-mat";
    return;
end

% 3) Legacy fallback: older exports sometimes used a raw MAT filename
%    without the "_final_scoring_dash" suffix.
base = regexprep(string(name),'_final_scoring_dash$','');

if endsWith(base,".mat")
    rawName = base;
else
    rawName = base + ".mat";
end

legacyCandidate = string(fullfile(folder,char(rawName)));

if isfile(fullfile(P.data_root,char(legacyCandidate)))
    relativePath = legacyCandidate;
    status = "derived-legacy-raw-name";
    return;
end

% 4) Keep the most likely candidate in the output even if missing so the
%    expected location is obvious during QC.
relativePath = sameBaseCandidate;
status = "derived-candidate-missing";
end
