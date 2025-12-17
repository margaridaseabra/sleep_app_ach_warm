function PLV = ach_eeg_plv_state_bouts(mat_file, scores_csv, stateName, ...
                                       t_start, t_end, varargin)
% ACH_EEG_PLV_STATE_BOUTS
% -------------------------------------------------------------
% Compute PLV between EEG band phase and slow ACh phase
% for all bouts of a given STATE (NREM or REM) within [t_start, t_end].
%
% OUTPUT struct PLV (fields ALWAYS present):
%   .stateName    - 'NREM' or 'REM'
%   .band_eeg     - [f1 f2] Hz (e.g. [1 4] for delta)
%   .band_ach     - [f1 f2] Hz (e.g. [0.1 1.0] for slow ACh)
%   .plv_all      - [nBouts×1] PLV per bout  (can be empty)
%   .dur_all      - [nBouts×1] bout duration (s)
%   .t0_all       - [nBouts×1] bout start times (s, global time)

p = inputParser;
p.addParameter('codes', struct('WK',0,'NREM',1,'REM',2,'MA',15));
p.addParameter('eeg_band',[1 4]);      % default: delta
p.addParameter('ach_band',[0.1 1.0]);  % default: slow ACh
p.addParameter('min_bout_sec', 10);    % ignore very short bouts
p.parse(varargin{:});
opt   = p.Results;
CODES = opt.codes;

% ---------- initialise output with EMPTY arrays (important!) ----------
PLV = struct();
PLV.stateName = stateName;
PLV.band_eeg  = opt.eeg_band;
PLV.band_ach  = opt.ach_band;
PLV.plv_all   = [];
PLV.dur_all   = [];
PLV.t0_all    = [];

% ---------- load signals from MAT file ----------
info  = whos('-file', mat_file);
names = {info.name};
pick  = @(cands) cands{find(ismember(cands,names),1,'first')};

eeg_name    = pick({'eeg','EEG','eeg1','eeg_filt'});
fs_eeg_name = pick({'eeg_frequency','fs_eeg','Fs_eeg'});
ach_name    = pick({'ach','ACh','ne','dff','dFF','dff_ach'});
fs_ach_name = pick({'ach_frequency','ne_frequency','fs_ach','Fs_ach'});

S      = load(mat_file, eeg_name, fs_eeg_name, ach_name, fs_ach_name);
eeg    = S.(eeg_name)(:);
fs_eeg = S.(fs_eeg_name);
ach    = S.(ach_name)(:);
fs_ach = S.(fs_ach_name);

% ---------- load 1-Hz scores ----------
M = readmatrix(scores_csv);
if size(M,2) == 1
    score     = M(:,1);
    epoch_sec = 1;
    t_scores  = (0:numel(score)-1)' * epoch_sec;
else
    t_scores  = M(:,1);
    score     = M(:,2);
    dt        = diff(t_scores);
    epoch_sec = mode(dt(~isnan(dt)));
end
score    = score(:);
t_scores = t_scores(:);

% ---------- which state? ----------
switch lower(stateName)
    case 'nrem'
        state_code = CODES.NREM;
    case 'rem'
        state_code = CODES.REM;
    otherwise
        error('stateName must be ''NREM'' or ''REM''.');
end

% ---------- restrict to analysis window & find bouts ----------
inSeg = t_scores >= t_start & t_scores < t_end;
s     = score;
s(~inSeg) = NaN;

isState = (s == state_code);
if ~any(isState)
    warning('No %s bouts in this window.', stateName);
    return;    % PLV struct already has correct fields but empty arrays
end

d  = diff([0; isState; 0]);
i1 = find(d ==  1);          % start indices
i2 = find(d == -1) - 1;      % end indices

t0_all = [];
t1_all = [];

for k = 1:numel(i1)
    t0 = t_scores(i1(k));
    t1 = t_scores(i2(k)) + epoch_sec;

    if t1 <= t_start || t0 >= t_end
        continue;
    end
    t0 = max(t0, t_start);
    t1 = min(t1, t_end);

    if t1 - t0 < opt.min_bout_sec
        continue;
    end

    t0_all(end+1,1) = t0; %#ok<AGROW>
    t1_all(end+1,1) = t1; %#ok<AGROW>
end

nB = numel(t0_all);
if nB == 0
    warning('No %s bouts longer than %.1f s.', stateName, opt.min_bout_sec);
    return;    % again: fields stay, arrays remain empty
end

plv_all = nan(nB,1);
dur_all = nan(nB,1);

% ---------- design filters ----------
% Make sure frequencies are valid
eeg_band = opt.eeg_band;
ach_band = opt.ach_band;

% Adjust bands if needed
if eeg_band(1) < 0.5, eeg_band(1) = 0.5; end
if eeg_band(2) >= fs_eeg/2, eeg_band(2) = fs_eeg/2 - 1; end
if ach_band(1) < 0.05, ach_band(1) = 0.05; end
if ach_band(2) >= fs_ach/2, ach_band(2) = fs_ach/2 - 1; end

% Design filters with error checking
try
    [bE,aE] = butter(4, eeg_band/(fs_eeg/2), 'bandpass');
catch ME
    warning(ME.identifier, 'EEG filter design failed: %s', ME.message);
    fprintf('  fs_eeg=%.1f, band=[%.2f %.2f]\n', fs_eeg, eeg_band);
    return;
end

try
    [bA,aA] = butter(4, ach_band/(fs_ach/2), 'bandpass');
catch ME
    warning(ME.identifier, 'ACh filter design failed: %s', ME.message);
    fprintf('  fs_ach=%.1f, band=[%.2f %.2f]\n', fs_ach, ach_band);
    return;
end

% ---------- compute PLV for each bout ----------
for b = 1:nB
    t0 = t0_all(b);
    t1 = t1_all(b);
    
    fprintf('  Bout %d/%d: %.1f-%.1f s (%.1f s duration)\n', ...
            b, nB, t0, t1, t1-t0);

    % indices in EEG & ACh
    iE0 = max(1, floor(t0*fs_eeg)+1);
    iE1 = min(numel(eeg), floor(t1*fs_eeg));
    iA0 = max(1, floor(t0*fs_ach)+1);
    iA1 = min(numel(ach), floor(t1*fs_ach));

    eeg_seg = eeg(iE0:iE1);
    ach_seg = ach(iA0:iA1);

    % Check segment lengths
    if numel(eeg_seg) < fs_eeg*2 || numel(ach_seg) < fs_ach*2
        fprintf('    Skipping: segments too short\n');
        continue;
    end
    
    % Check for NaN/Inf
    if any(~isfinite(eeg_seg)) || any(~isfinite(ach_seg))
        fprintf('    Skipping: NaN/Inf in segments\n');
        continue;
    end

    % Detrend and normalize
    eeg_seg = detrend(eeg_seg - mean(eeg_seg));
    ach_seg = detrend(ach_seg - mean(ach_seg));
    
    % Filter with error catching
    try
        eeg_f = filtfilt(bE, aE, double(eeg_seg));
        ach_f = filtfilt(bA, aA, double(ach_seg));
    catch ME
        fprintf('    Skipping bout %d: filtering failed (%s)\n', b, ME.message);
        continue;
    end
    
    % Check filtered signals
    if all(abs(eeg_f) < eps) || all(abs(ach_f) < eps)
        fprintf('    Skipping: filtered signal is flat\n');
        continue;
    end

    % Hilbert transform
    try
        phi_eeg = angle(hilbert(eeg_f));
        phi_ach = angle(hilbert(ach_f));
    catch ME
        fprintf('    Skipping bout %d: Hilbert failed (%s)\n', b, ME.message);
        continue;
    end

    % Interpolate ACh phase to EEG time base
    t_eeg = (0:numel(phi_eeg)-1)'/fs_eeg;
    t_ach = (0:numel(phi_ach)-1)'/fs_ach;

    try
        phi_ach_i = interp1(t_ach, unwrap(phi_ach), t_eeg, 'linear', 'extrap');
        phi_ach_i = mod(phi_ach_i, 2*pi);
    catch ME
        fprintf('    Skipping bout %d: interpolation failed (%s)\n', b, ME.message);
        continue;
    end
    
    % Check for valid phases
    if any(~isfinite(phi_eeg)) || any(~isfinite(phi_ach_i))
        fprintf('    Skipping: NaN/Inf in phases\n');
        continue;
    end

    % Compute PLV
    dphi = phi_eeg - phi_ach_i;
    plv_val = abs(mean(exp(1i*dphi)));
    
    fprintf('    PLV = %.3f\n', plv_val);
    
    plv_all(b) = plv_val;
    dur_all(b) = t1 - t0;
end
% ---------- store results ----------
PLV.plv_all = plv_all;
PLV.dur_all = dur_all;
PLV.t0_all  = t0_all;
end
