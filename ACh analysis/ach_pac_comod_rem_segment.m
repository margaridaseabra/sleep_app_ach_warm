function COM = ach_pac_comod_rem_segment(mat_file, scores_csv, t_start, t_end, CODES, varargin)
% ACH_PAC_COMOD_REM_SEGMENT  Robust REM theta–gamma PAC comodulogram
% --------------------------------------------------------------------
% Computes PAC (Modulation Index, Tort et al.) for REM epochs between
% t_start and t_end, using EEG theta phase and gamma amplitude.
%
% This version is robust to zeros in the phase–amplitude histogram and
% converts 0*log(0) situations into MI ~ 0 instead of NaN.
%
% INPUTS
%   mat_file   : .mat file with EEG and sampling rate
%   scores_csv : scoring CSV (time_s, score) or (score only, 1 Hz)
%   t_start    : start time (s)
%   t_end      : end time (s)
%   CODES      : struct with fields .REM (and others, not used here)
%
% NAME/VALUE PAIRS
%   'phase_freqs' : vector of phase frequencies (Hz), e.g. 6:0.5:10
%   'amp_freqs'   : vector of amp frequencies (Hz),   e.g. 30:2:80
%   'nbins'       : number of phase bins, default 18
%   'doPlot'      : logical, plot comodulogram if true
%   'label'       : optional title for plot
%
% OUTPUT
%   COM struct with fields:
%       .MI          : [nPhase x nAmp] modulation index
%       .phase_freqs : phase frequencies (Hz)
%       .amp_freqs   : amplitude frequencies (Hz)
%       .nSamples    : # of REM samples used
%       .fs          : sampling rate (Hz)

% -------------------------------------------------------------
% Parse options
% -------------------------------------------------------------
p = inputParser;
p.addParameter('phase_freqs', 6:0.5:10);
p.addParameter('amp_freqs',   30:2:80);
p.addParameter('nbins',       18);
p.addParameter('doPlot',      false);
p.addParameter('label',       '');
p.parse(varargin{:});
opt = p.Results;

phase_freqs = opt.phase_freqs(:)';
amp_freqs   = opt.amp_freqs(:)';
nbins       = opt.nbins;

% -------------------------------------------------------------
% Load EEG from MAT
% -------------------------------------------------------------
info  = whos('-file', mat_file);
names = {info.name};
pick  = @(cands) cands{find(ismember(cands,names),1,'first')};

eeg_name = pick_if_present(pick, {'eeg','EEG','eeg1','Eeg','eeg_filt'});
fs_name  = pick_if_present(pick, {'eeg_frequency','fs_eeg','Fs_eeg','EEG_frequency'});

if isempty(eeg_name) || isempty(fs_name)
    warning('ach_pac_comod_rem_segment: no EEG or fs in %s', mat_file);
    COM = empty_COM(phase_freqs, amp_freqs);
    return;
end

S   = load(mat_file, eeg_name, fs_name);
eeg = double(S.(eeg_name)(:));
fs  = S.(fs_name);

% Remove NaNs in EEG just in case
eeg(~isfinite(eeg)) = 0;

% -------------------------------------------------------------
% Load scores and build REM mask within [t_start, t_end]
% -------------------------------------------------------------
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

CODE_REM = CODES.REM;

rem_idx = false(size(eeg));  % sample-level mask

for i = 1:numel(score)
    if score(i) ~= CODE_REM
        continue;
    end

    t0 = t_scores(i);
    t1 = t0 + epoch_sec;

    % intersect with analysis window
    if t1 <= t_start || t0 >= t_end
        continue;
    end
    t0c = max(t0, t_start);
    t1c = min(t1, t_end);

    i0 = max(1, floor(t0c * fs) + 1);
    i1 = min(numel(eeg), floor(t1c * fs));
    if i1 > i0
        rem_idx(i0:i1) = true;
    end
end

eeg_rem = eeg(rem_idx);
eeg_rem(~isfinite(eeg_rem)) = 0;

nSamples = numel(eeg_rem);
if nSamples < fs * 5   % require at least ~5 s of REM
    warning('ach_pac_comod_rem_segment: too little REM data in window.');
    COM = empty_COM(phase_freqs, amp_freqs);
    COM.nSamples = nSamples;
    COM.fs       = fs;
    return;
end

% -------------------------------------------------------------
% Compute PAC comodulogram with robust MI
% -------------------------------------------------------------
MI = compute_pac_comod_matrix(eeg_rem, fs, phase_freqs, amp_freqs, nbins);

COM = struct();
COM.MI          = MI;
COM.phase_freqs = phase_freqs;
COM.amp_freqs   = amp_freqs;
COM.nSamples    = nSamples;
COM.fs          = fs;

% -------------------------------------------------------------
% Optional plotting
% -------------------------------------------------------------
if opt.doPlot
    figure('Color','w','Position',[200 200 500 400]);
    imagesc(amp_freqs, phase_freqs, MI);
    axis xy;
    colormap(turbo);
    c = colorbar;
    ylabel(c, 'MI');
    xlabel('Amplitude (Hz)');
    ylabel('Phase (Hz)');
    if ~isempty(opt.label)
        title(opt.label, 'Interpreter','none');
    else
        title('REM PAC (theta phase × gamma amplitude)');
    end
end

end

% =============================================================
% Helper: safe PAC comodulogram
% =============================================================
function MI = compute_pac_comod_matrix(sig, fs, phase_freqs, amp_freqs, nbins)

sig = double(sig(:));
sig(~isfinite(sig)) = 0;

nPhase = numel(phase_freqs);
nAmp   = numel(amp_freqs);
MI     = nan(nPhase, nAmp);

% bandwidths (can be tuned)
phase_bw = 1.0;   % +/- 1 Hz
amp_bw   = 10.0;  % +/- 10 Hz

phase_edges = linspace(-pi, pi, nbins+1);

for ip = 1:nPhase
    fp = phase_freqs(ip);
    f1p = max(0.5, fp - phase_bw);
    f2p = min(fs/2 - 0.5, fp + phase_bw);
    if f2p <= f1p + 0.1
        continue;
    end

    WnP = [f1p f2p] / (fs/2);
    [bP,aP] = butter(4, WnP, 'bandpass');
    try
        sigP = filtfilt(bP, aP, sig);
    catch
        continue;
    end
    phase = angle(hilbert(sigP));
    phase(~isfinite(phase)) = [];

    if numel(phase) < 100
        continue;
    end

    for ia = 1:nAmp
        fa = amp_freqs(ia);
        f1a = max(1, fa - amp_bw);
        f2a = min(fs/2 - 1, fa + amp_bw);
        if f2a <= f1a + 1
            continue;
        end

        WnA = [f1a f2a] / (fs/2);
        [bA,aA] = butter(4, WnA, 'bandpass');
        try
            sigA = filtfilt(bA, aA, sig);
        catch
            continue;
        end
        amp = abs(hilbert(sigA));

        % Keep matching length and finite samples
        n = min(numel(phase), numel(amp));
        if n < 100
            continue;
        end
        ph = phase(1:n);
        am = amp(1:n);
        good = isfinite(ph) & isfinite(am);
        ph = ph(good);
        am = am(good);

        if numel(ph) < 100
            continue;
        end

        % Phase binning
        [~,~,binIdx] = histcounts(ph, phase_edges);
        binMeans = zeros(1, nbins);
        for k = 1:nbins
            thisAmp = am(binIdx == k);
            if isempty(thisAmp)
                binMeans(k) = 0;
            else
                binMeans(k) = mean(thisAmp);
            end
        end

        if all(binMeans == 0)
            MI(ip,ia) = 0;   % no modulation, but not NaN
            continue;
        end

        p = binMeans / sum(binMeans);
        p = p(:);
        p = p / sum(p);     % renormalise robustly
        p(p <= 0) = eps;    % avoid log(0)

        H    = -sum(p .* log(p));      % entropy
        Hmax = log(numel(p));
        MI(ip,ia) = (Hmax - H) / Hmax; % Tort MI, in [0,1]
    end
end

% Replace residual NaNs/Infs with 0 (means "no reliable modulation")
MI(~isfinite(MI)) = 0;

end

% =============================================================
% Small helpers
% =============================================================
function name = pick_if_present(pick, cands)
    try
        name = pick(cands);
    catch
        name = [];
    end
end

function COM = empty_COM(phase_freqs, amp_freqs)
    MI = nan(numel(phase_freqs), numel(amp_freqs));
    COM = struct('MI', MI, ...
                 'phase_freqs', phase_freqs, ...
                 'amp_freqs',   amp_freqs, ...
                 'nSamples',    0, ...
                 'fs',          NaN);
end
