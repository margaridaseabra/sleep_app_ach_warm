function PLV = ach_eeg_plv_periTransitions_segment(mat_file, scores_csv, ...
    from_state, to_state, t_seg_start, t_seg_end, varargin)
% ACH_EEG_PLV_PERITRANSITIONS_SEGMENT
% -------------------------------------------------------------
% Compute EEG–ACh PLV around state transitions for ONE recording,
% but ONLY using transitions whose onset falls inside a selected
% time window [t_seg_start, t_seg_end] (in seconds).
%
% This matches your "segment" setup used in plot_ach_eeg_segment.
%
% Now also:
%   -> applies 50 Hz notch to EEG for ALL files
%      EXCEPT the following:
%        20251001_baseline_mouse1_APP.mat
%        20251002_baseline_mouse2_WT.mat
%        20251003_ambtemp_mouse1_APP.mat
%        20251005_baseline_mouse8_WT.mat
%        20251006_baseline_mouse4_WT.mat
%
% Example:
%   PLV = ach_eeg_plv_periTransitions_segment( ...
%            mat_file, scores_csv, ...
%            'NREM','REM', ...
%            0, 1800, ...             % segment: first 30 min
%            'codes', CODES);
%
% INPUTS
%   mat_file     : .mat with EEG + ACh + their sampling rates
%   scores_csv   : 1-Hz scores (time_s, score) or just score column
%   from_state   : 'Wake','NREM','REM','MA'
%   to_state     : 'Wake','NREM','REM','MA'
%   t_seg_start  : start of segment (s)
%   t_seg_end    : end of segment (s)
%
% Optional name/value pairs:
%   'codes'      : struct with fields .WK .NREM .REM .MA
%                  (default 0,1,2,15)
%   'pre_sec'    : seconds before transition (default 20)
%   'post_sec'   : seconds after transition (default 40)
%   'eeg_band'   : [f1 f2] Hz for EEG phase (default [5 10], theta)
%   'ach_band'   : [f1 f2] Hz for ACh band (default [0.05 0.5])
%
% OUTPUT (struct PLV):
%   .t_rel       : 1×T time axis (s), -pre_sec:dt:post_sec
%   .plv_time    : 1×T PLV(t) across all transitions in segment
%   .plv_events  : nEvents×1 PLV value per event (whole window)
%   .nEvents     : number of transitions used
%   .from_state, .to_state
%   .eeg_band, .ach_band
%   .t_seg_start, .t_seg_end   : the segment bounds you used

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

% ----- apply 50 Hz notch unless file is on the skip list -----
eeg = apply_50Hz_notch_if_needed(eeg, fs_eeg, mat_file);

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

% -------------------- find all transitions in entire recording --------------------
lab = score;
idx_all = find(lab(1:end-1) == code_from & lab(2:end) == code_to);
if isempty(idx_all)
    warning('No %s→%s transitions found in recording.', from_state, to_state);
    PLV = emptyPLV(from_state, to_state, eeg_band, ach_band, t_seg_start, t_seg_end);
    return;
end
t_onsets_all = t_scores(idx_all+1);   % first second of new state

% -------------------- restrict transitions to segment --------------------
use_idx = t_onsets_all >= t_seg_start & t_onsets_all <= t_seg_end;
t_onsets = t_onsets_all(use_idx);

if isempty(t_onsets)
    warning('No %s→%s transitions inside segment [%.1f, %.1f] s.', ...
        from_state, to_state, t_seg_start, t_seg_end);
    PLV = emptyPLV(from_state, to_state, eeg_band, ach_band, t_seg_start, t_seg_end);
    return;
end

% -------------------- design band-pass filters --------------------
% EEG band
eeg_Wn = eeg_band / (fs_eeg/2);
if eeg_Wn(2) >= 1
    error('EEG band %.2f–%.2f Hz invalid for fs_eeg = %.1f Hz', ...
          eeg_band(1), eeg_band(2), fs_eeg);
end
[beeg, aeeg] = butter(3, eeg_Wn, 'bandpass');

% ACh band
ach_Wn = ach_band / (fs_ach/2);
if ach_Wn(2) >= 1
    error('ACh band %.3f–%.3f Hz invalid for fs_ach = %.1f Hz', ...
          ach_band(1), ach_band(2), fs_ach);
end
[bach, aach] = butter(2, ach_Wn, 'bandpass');

% -------------------- filter full signals once --------------------
eeg_f = filtfilt(beeg, aeeg, eeg);
ach_f = filtfilt(bach, aach, ach);

% -------------------- common time axis for peri-transition window --------------------
dt    = 1/fs_eeg;                        % use EEG sampling
t_rel = -pre_sec : dt : post_sec;        % relative time axis
nT    = numel(t_rel);

% -------------------- loop over transitions in this segment --------------------
Z          = [];   % nEvents × nT complex phase differences
plv_events = [];

for k = 1:numel(t_onsets)
    t0 = t_onsets(k);          % transition time in seconds

    % EEG window indices
    t_start_win = t0 - pre_sec;
    i0_e        = round(t_start_win * fs_eeg) + 1;
    i1_e        = i0_e + nT - 1;

    if i0_e < 1 || i1_e > numel(eeg_f)
        continue;   % too close to edges of recording
    end

    % Enforce full window inside segment (optional but safer)
    if t_start_win < t_seg_start || (t0 + post_sec) > t_seg_end
        continue;
    end

    seg_eeg = eeg_f(i0_e:i1_e);

    % ACh: align to same absolute times as EEG window
    t_win = t_start_win + (0:nT-1)*dt;    % absolute times of each sample
    iA    = round(t_win * fs_ach) + 1;
    if min(iA) < 1 || max(iA) > numel(ach_f)
        continue;
    end
    seg_ach = ach_f(iA);

    % analytic signals via Hilbert
    eeg_h = hilbert(seg_eeg);
    ach_h = hilbert(seg_ach);

    phi_eeg = angle(eeg_h);
    phi_ach = angle(ach_h);

    dphi = phi_eeg - phi_ach;
    z    = exp(1i * dphi);   % unit vectors

    Z(end+1, :)         = z;              %#ok<AGROW>
    plv_events(end+1,1) = abs(mean(z));   %#ok<AGROW>
end

nEvents = size(Z,1);
if nEvents == 0
    warning('No usable %s→%s transitions with full windows inside segment.', ...
            from_state, to_state);
    PLV = emptyPLV(from_state, to_state, eeg_band, ach_band, t_seg_start, t_seg_end);
    return;
end

% PLV as a function of time (average over events)
plv_time = abs(mean(Z,1));

PLV = struct();
PLV.t_rel       = t_rel;
PLV.plv_time    = plv_time;
PLV.plv_events  = plv_events;
PLV.nEvents     = nEvents;
PLV.from_state  = from_state;
PLV.to_state    = to_state;
PLV.eeg_band    = eeg_band;
PLV.ach_band    = ach_band;
PLV.t_seg_start = t_seg_start;
PLV.t_seg_end   = t_seg_end;
end

% -------------------- helper for empty output --------------------
function PLV = emptyPLV(from_state, to_state, eeg_band, ach_band, t_seg_start, t_seg_end)
PLV = struct('t_rel',[],'plv_time',[],'plv_events',[], ...
    'nEvents',0,'from_state',from_state,'to_state',to_state, ...
    'eeg_band',eeg_band,'ach_band',ach_band, ...
    't_seg_start',t_seg_start,'t_seg_end',t_seg_end);
end

% -------------------- helper: 50 Hz notch with skip list --------------------
function eeg_out = apply_50Hz_notch_if_needed(eeg_in, fs_eeg, mat_file)
eeg_out = eeg_in;

% if sampling rate too low, skip
if fs_eeg <= 120
    return;
end

[~, base, ext] = fileparts(mat_file);
fname = [base ext];

skip_list = { ...
    '20251001_baseline_mouse1_APP.mat', ...
    '20251002_baseline_mouse2_WT.mat', ...
    '20251003_ambtemp_mouse1_APP.mat', ...
    '20251005_baseline_mouse8_WT.mat', ...
    '20251006_baseline_mouse4_WT.mat' ...
    };

if ismember(fname, skip_list)
    % do NOT notch these
    return;
end

% design 50 Hz notch
wo = 50 / (fs_eeg/2);        % normalized freq
bw = wo / 35;                % quality factor ~35
[bN, aN] = iirnotch(wo, bw);

eeg_out = filtfilt(bN, aN, eeg_in);
end
