function P = paths_config()
%PATHS_CONFIG Repository and local data paths for the sleep/ACh project.
%
% The repository itself is detected automatically from this file.
%
% Raw/processed experimental data are intentionally kept outside Git.
% To point the pipeline to the data folder, either:
%
%   1) set the environment variable AMBTEMP_DATA_ROOT, or
%   2) place the data at ~/Desktop/AmbTemp_proj (convenient default).
%
% Example from MATLAB for a non-default location:
%   setenv('AMBTEMP_DATA_ROOT','/path/to/AmbTemp_proj')
%
% No user-specific absolute paths should be committed to the repository.

configDir = fileparts(mfilename('fullpath'));
P.repo = fileparts(configDir);

dataRoot = string(getenv('AMBTEMP_DATA_ROOT'));

if strlength(dataRoot) == 0
    homeDir = string(getenv('HOME'));
    defaultRoot = fullfile(homeDir,'Desktop','AmbTemp_proj');

    if isfolder(defaultRoot)
        dataRoot = string(defaultRoot);
    else
        error([ ...
            'Could not locate the experimental data folder. ' ...
            'Set AMBTEMP_DATA_ROOT to the AmbTemp_proj directory, e.g. ' ...
            'setenv(''AMBTEMP_DATA_ROOT'',''/path/to/AmbTemp_proj'').']);
    end
end

P.data_root = char(dataRoot);

P.thesis_scores = fullfile( ...
    P.data_root, ...
    'scored_files_oct25batch');

P.new_scores = fullfile( ...
    P.data_root, ...
    'scored_files_june26batch');

P.mice = fullfile(P.repo,'config','mice.csv');
P.sessions = fullfile(P.repo,'config','sessions.csv');

% Experimental manipulation timing registries.
P.warming_times = fullfile(P.repo,'config','warming_times.csv');
P.drug_times = fullfile(P.repo,'config','drug_times.csv');

P.outputs = fullfile(P.repo,'outputs');
P.sleep_output = fullfile(P.outputs,'sleep_architecture');
P.ach_output = fullfile(P.outputs,'ach');
P.eeg_output = fullfile(P.outputs,'eeg');
P.rem_output = fullfile(P.outputs,'rem_features');
P.survival_output = fullfile(P.outputs,'survival');
end
