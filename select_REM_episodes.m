function select_REM_episodes(varargin)
% select_REM_episodes: find REM bouts in a scored .mat and export selections
%
% INPUTS (name/value, all optional):
%   'file'        : path to .mat (default: UI prompt)
%   'VARS'        : struct with field names in the .mat:
%                   eeg, emg, ach (optional), scores, fs_eeg, fs_ach (opt)
%                   default = struct('eeg','eeg','emg','emg','ach','ne',...
%                                    'scores','sleep_scores',...
%                                    'fs_eeg','eeg_frequency','fs_ach','ne_frequency')
%   'EPOCH_SEC'   : epoch length in seconds used during scoring (default 1)
%   'REM_CODES'   : numeric code(s) that represent REM in sleep_scores (default [2])
%   'OUTDIR'      : folder to save outputs (default: alongside input)
%
% OUTPUT FILES:
%   <input_basename>_REMselection.mat
%       .rem_mask_samples   (logical, length = EEG samples)
%       .rem_bouts_table    (table of bout start/end epochs & times)
%       .bouts              (struct array with EEG/EMG/ACh segments per bout)
%       .params             (struct of run parameters)
%       .eeg_frequency, .fs_eeg  (EEG sampling rate, duplicated for compatibility)
%       .ne_frequency,  .fs_ach  (FP/ACh sampling rate, if present)
%   <input_basename>_REMbouts.csv

%% -------------------- Parse args --------------------
p = inputParser;
addParameter(p,'file','',@(s)ischar(s)||isstring(s));
addParameter(p,'VARS',struct('eeg','eeg','emg','emg','ach','ne', ...
                             'scores','sleep_scores','fs_eeg','eeg_frequency','fs_ach','ne_frequency'));
addParameter(p,'EPOCH_SEC',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'REM_CODES',[2],@(v)isnumeric(v)&&~isempty(v));
addParameter(p,'OUTDIR','',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
FILE       = string(p.Results.file);
VARS       = p.Results.VARS;
EPOCH_SEC  = p.Results.EPOCH_SEC;
REM_CODES  = p.Results.REM_CODES;
OUTDIR     = string(p.Results.OUTDIR);

if FILE=="" || ~isfile(FILE)
    [f,fp] = uigetfile('*.mat','Select scored .mat');
    assert(f~=0,'No file selected.');
    FILE = string(fullfile(fp,f));
end
if OUTDIR=="", OUTDIR = string(fileparts(FILE)); end

%% -------------------- Load & sanity checks --------------------
S = load(FILE);
assert(isfield(S,VARS.scores), 'Missing scores variable "%s".', VARS.scores);
scores = S.(VARS.scores); scores = scores(:);

hasEEG = isfield(S,VARS.eeg);
hasEMG = isfield(S,VARS.emg);
hasACH = isfield(S,VARS.ach);
assert(hasEEG, 'Missing EEG variable "%s".', VARS.eeg);

fs_eeg = []; fs_ach = [];
if isfield(S,VARS.fs_eeg), fs_eeg = S.(VARS.fs_eeg); end
assert(~isempty(fs_eeg) && isnumeric(fs_eeg) && isscalar(fs_eeg), ...
    'Missing or invalid fs_eeg variable "%s".', VARS.fs_eeg);
if hasACH && isfield(S,VARS.fs_ach)
    fs_ach = S.(VARS.fs_ach);
end

EEG = S.(VARS.eeg); EEG = EEG(:);
if hasEMG, EMG = S.(VARS.emg); EMG = EMG(:); else, EMG = []; end
if hasACH, ACH = S.(VARS.ach); ACH = ACH(:); else, ACH = []; end

% Basic length check (1 epoch label per EPOCH_SEC seconds)
total_dur_s   = numel(EEG)/fs_eeg;
expected_ep   = floor(total_dur_s / EPOCH_SEC);
if abs(numel(scores)-expected_ep) > max(5, 0.01*expected_ep)
    warning('scores length (%d) differs from expected (%d). Check EPOCH_SEC/fs_eeg.', numel(scores), expected_ep);
end

%% -------------------- Find REM bouts --------------------
isREM_ep = ismember(scores, REM_CODES);        % logical per epoch
% Convert epoch-wise logical to sample-wise mask
nSamples  = numel(EEG);
samples_per_epoch = round(EPOCH_SEC * fs_eeg);
rem_mask_samples = false(nSamples,1);

nEpochs = numel(scores);
for ep = 1:nEpochs
    if isREM_ep(ep)
        s0 = (ep-1)*samples_per_epoch + 1;
        s1 = min(ep*samples_per_epoch, nSamples);
        rem_mask_samples(s0:s1) = true;
    end
end

% Identify contiguous REM bouts in epoch space
d = diff([false; isREM_ep; false]);  % rising/falling edges
bout_starts = find(d==1);
bout_ends   = find(d==-1)-1;
nBouts = numel(bout_starts);

% Build bout table
start_epoch   = bout_starts(:);
end_epoch     = bout_ends(:);
dur_epochs    = end_epoch - start_epoch + 1;
start_time_s  = (start_epoch-1)*EPOCH_SEC;
end_time_s    = end_epoch*EPOCH_SEC;
dur_s         = dur_epochs*EPOCH_SEC;

rem_bouts_table = table(start_epoch,end_epoch,dur_epochs,start_time_s,end_time_s,dur_s);

%% -------------------- Extract per-bout signals --------------------
bouts = struct('idx_samples',[],'t_s',[],'EEG',[],'EMG',[],'ACH',[]);
t_eeg = (0:nSamples-1)'/fs_eeg;

for k = 1:nBouts
    s0 = max(1, (start_epoch(k)-1)*samples_per_epoch + 1);
    s1 = min(nSamples, end_epoch(k)*samples_per_epoch);
    idx = (s0:s1)';

    bouts(k).idx_samples = idx;
    bouts(k).t_s         = t_eeg(idx);
    bouts(k).EEG         = EEG(idx);
    if ~isempty(EMG), bouts(k).EMG = EMG(idx); end
    if ~isempty(ACH), bouts(k).ACH = ACH_timecrop(ACH, fs_ach, fs_eeg, idx); end
end

%% -------------------- Save outputs --------------------
[~, base, ~] = fileparts(FILE);
new_output_dir = 'N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\REM analysis'; % your desired folder
if ~exist(new_output_dir,'dir'), mkdir(new_output_dir); end

mat_out = fullfile(new_output_dir, base + "_REMselection.mat");
csv_out = fullfile(new_output_dir, base + "_REMbouts.csv");

% --- ADDED: create top-level aliases for sampling rates (no shadowing) ---
eeg_frequency = fs_eeg;
fs_eeg_top    = fs_eeg;        % use a different var name to avoid confusion
save_vars = {'rem_mask_samples','rem_bouts_table','bouts','scores'}; % save scores too

if ~isempty(fs_ach)
    ne_frequency  = fs_ach;
    fs_ach_top    = fs_ach;
end

% --- ADDED: package params for provenance ---
params = struct('file',FILE,'VARS',VARS,'EPOCH_SEC',EPOCH_SEC, ...
                'REM_CODES',REM_CODES,'fs_eeg',eeg_frequency,'fs_ach',[]);
if exist('ne_frequency','var'); params.fs_ach = ne_frequency; end

% --- ADDED: ensure top-level signals exist with standard names ---
eeg = EEG;                          % always save as lower-case 'eeg'
if ~isempty(EMG), emg = EMG; end    % 'emg' if present
if ~isempty(ACH), ne  = ACH;  end   % 'ne'  if present


% Write CSV summary (unchanged)
try
    writetable(rem_bouts_table, csv_out);
catch
    warning('Could not write CSV to %s. Check permissions.', csv_out);
end

% --- ADDED: final save including signals and sampling rates at top level ---
if exist('ne_frequency','var')
    save(mat_out, save_vars{:}, 'params', ...
        'eeg','emg','ne', ...                   % signals
        'eeg_frequency','fs_eeg_top', ...       % EEG sampling rates
        'ne_frequency','fs_ach_top', ...        % FP/ACh sampling rates
        '-v7.3');
else
    % no ACh/FP present
    if exist('emg','var')
        save(mat_out, save_vars{:}, 'params', ...
            'eeg','emg', ...
            'eeg_frequency','fs_eeg_top', ...
            '-v7.3');
    else
        save(mat_out, save_vars{:}, 'params', ...
            'eeg', ...
            'eeg_frequency','fs_eeg_top', ...
            '-v7.3');
    end
end

fprintf('[select_REM_episodes] Saved:\n  %s\n  %s\n', mat_out, csv_out);
end % function


