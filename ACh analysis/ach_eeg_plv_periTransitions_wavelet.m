function PLV = ach_eeg_plv_periTransitions_wavelet(mat_file, scores_csv, ...
                                                   from_state, to_state, varargin)
% ACH_EEG_PLV_PERITRANSITIONS_WAVELET
% -------------------------------------------------------------
% Compute EEG–ACh PLV around state transitions for ONE recording,
% using CWT (wavelet) to extract band-limited complex signals.
%
% Example:
%   PLV = ach_eeg_plv_periTransitions_wavelet( ...
%            '20251003_ambtemp_mouse1_APP.mat', ...
%            '20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv', ...
%            'NREM','REM');
%
% INPUTS
%   mat_file   : .mat with EEG + ACh + their sampling rates
%   scores_csv : 1-Hz scores (time_s, score) or just score column
%   from_state : 'Wake','NREM','REM','MA'
%   to_state   : 'Wake','NREM','REM','MA'
%
% Optional name/value pairs:
%   'codes'    : struct with fields .WK .NREM .REM .MA (default 0,1,2,15)
%   'pre_sec'  : seconds before transition (default 20)
%   'post_sec' : seconds after transition (default 40)
%   'eeg_band' : [f1 f2] Hz for EEG phase (default [5 10], theta)
%   'ach_band' : [f1 f2] Hz for ACh band (default [0.05 0.5])
%
% OUTPUT (struct PLV):
%   .t_rel       : 1×T time axis (s), -pre_sec:dt:post_sec
%   .plv_time    : 1×T PLV(t) across all transitions
%   .plv_events  : nEvents×1 PLV value per event (whole window)
%   .nEvents     : number of transitions used
%   .from_state, .to_state
%   .eeg_band, .ach_band
%
% NOTE: Uses CWT (Wavelet Toolbox). If you don't have it, we can
% switch to a custom Morlet implementation.

% -------------------- defaults --------------------
p = inputParser;
p.addParameter('codes', struct('WK',0,'NREM',1,'REM',2,'MA',15));
p.addParameter('pre_sec', 20);
p.addParameter('post_sec',40);
p.addParameter('eeg_band',[5 10]);       % theta
p.addParameter('ach_band',[0.05 0.5]);   % slow ACh modulation
p.parse(varargin{:});
CODES    = p.Results.codes;
pre_sec  = p.Results.pre_sec;
post_sec = p.Results.post_sec;
eeg_band = p.Results.eeg_band;
ach_band = p.Results.ach_band;

% -------------------- load EEG + ACh --------------------
info  = whos('-file', mat_file);
names = {info.name};

pick = @(cands) cands{find(ismember(cands,names),1,'first')};

% EEG
eeg_name    = pick({'eeg','EEG','eeg1'});
fs_eeg_name = pick({'eeg_frequency','Fs_eeg','fs_eeg','EEG_frequency'});
if isempty(eeg_name) || isempty(fs_eeg_name)
    error('Could not find EEG or its sampling rate in %s', mat_file);
end
S      = load(mat_file, eeg_name, fs_eeg_name);
eeg    = double(S.(eeg_name)(:));
fs_eeg = double(S.(fs_eeg_name));

% ACh
ach_name    = pick({'ach','ACh','ne','dff','dFF','dff_ach'});
fs_ach_name = pick({'ach_frequency','ne_frequency','fs_ach','Fs_ach'});
if isempty(ach_name) || isempty(fs_ach_name)
    error('Could not find ACh or its sampling rate in %s', mat_file);
end
S2     = load(mat_file, ach_name, fs_ach_name);
ach    = double(S2.(ach_name)(:));
fs_ach = double(S2.(fs_ach_name));

% -------------------- load scores --------------------
Msc = readmatrix(scores_csv);
if size(Msc,2) == 1
    score     = Msc(:,1);
    epoch_sec = 1;
    t_scores  = (0:numel(score)-1)' * epoch_sec;
else
    t_scores  = Msc(:,1);
    score     = Msc(:,2);
    dt_s      = diff(t_scores);
    epoch_sec = mode(dt_s(~isnan(dt_s)));
end
score    = score(:);
t_scores = t_scores(:);

% -------------------- state codes --------------------
state2code = containers.Map( ...
    {'Wake','NREM','REM','MA'}, ...
    [CODES.WK, CODES.NREM, CODES.REM, CODES.MA]);

if ~isKey(state2code, from_state) || ~isKey(state2code, to_state)
    error('from_state/to_state must be one of: Wake, NREM, REM, MA');
end

code_from = state2code(from_state);
code_to   = state2code(to_state);

% -------------------- find transitions --------------------
lab = score;
idx_trans = find(lab(1:end-1) == code_from & lab(2:end) == code_to);
if isempty(idx_trans)
    warning('No %s→%s transitions found.', from_state, to_state);
    PLV = struct('t_rel',[],'plv_time',[],'plv_events',[], ...
                 'nEvents',0,'from_state',from_state,'to_state',to_state, ...
                 'eeg_band',eeg_band,'ach_band',ach_band);
    return;
end
t_onsets = t_scores(idx_trans+1);   % first second of new state

% -------------------- common time axis --------------------
dt   = 1/fs_eeg;                % use EEG sampling for time grid
t_rel = -pre_sec : dt : post_sec;
nT   = numel(t_rel);

% -------------------- CWT (wavelet) on full signals --------------------
% EEG: complex time–freq representation
[cfs_eeg, f_eeg] = cwt(eeg, fs_eeg);   % cfs_eeg: nF × nTime
% ACh:
[cfs_ach, f_ach] = cwt(ach, fs_ach);   % cfs_ach: nF × nTime

% Band-averaged complex analytic signals for EEG + ACh
eeg_idx = f_eeg >= eeg_band(1) & f_eeg <= eeg_band(2);
if ~any(eeg_idx)
    error('No EEG wavelet freqs within band %.2f–%.2f Hz', eeg_band(1), eeg_band(2));
end
sig_eeg_band = squeeze(mean(cfs_eeg(eeg_idx, :), 1));   % 1 × nSamples_eeg, complex

ach_idx = f_ach >= ach_band(1) & f_ach <= ach_band(2);
if ~any(ach_idx)
    error('No ACh wavelet freqs within band %.3f–%.3f Hz', ach_band(1), ach_band(2));
end
sig_ach_band = squeeze(mean(cfs_ach(ach_idx, :), 1));   % 1 × nSamples_ach, complex

% -------------------- loop over transitions --------------------
Z          = [];   % complex phase diffs: nEvents × nT
plv_events = [];

for k = 1:numel(t_onsets)
    t0 = t_onsets(k) - pre_sec;
    t1 = t_onsets(k) + post_sec;

    % EEG window indices on analytic band signal
    i0_e = round(t0 * fs_eeg) + 1;
    i1_e = i0_e + nT - 1;
    if i0_e < 1 || i1_e > numel(sig_eeg_band)
        continue;  % too close to edges
    end
    seg_eeg = sig_eeg_band(i0_e:i1_e);

    % ACh indices using same absolute time grid
    t_win = t0 + t_rel;        % absolute times of window
    iA    = round(t_win * fs_ach) + 1;
    if min(iA) < 1 || max(iA) > numel(sig_ach_band)
        continue;
    end
    seg_ach = sig_ach_band(iA);

    % phases directly from complex wavelet coefficients
    phi_eeg = angle(seg_eeg);
    phi_ach = angle(seg_ach);

    dphi = phi_eeg - phi_ach;
    z    = exp(1i * dphi);   % unit complex vectors

    Z(end+1, :)       = z;              %#ok<AGROW>
    plv_events(end+1,1) = abs(mean(z)); %#ok<AGROW>
end

nEvents = size(Z,1);
if nEvents == 0
    warning('All %s→%s transitions were too close to edges.', ...
            from_state, to_state);
    PLV = struct('t_rel',[],'plv_time',[],'plv_events',[], ...
                 'nEvents',0,'from_state',from_state,'to_state',to_state, ...
                 'eeg_band',eeg_band,'ach_band',ach_band);
    return;
end

% PLV as a function of time
plv_time = abs(mean(Z,1));

PLV = struct();
PLV.t_rel      = t_rel;
PLV.plv_time   = plv_time;
PLV.plv_events = plv_events;
PLV.nEvents    = nEvents;
PLV.from_state = from_state;
PLV.to_state   = to_state;
PLV.eeg_band   = eeg_band;
PLV.ach_band   = ach_band;
end
