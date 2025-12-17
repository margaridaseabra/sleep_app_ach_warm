function COM = ach_pac_comod_nrem_segment(mat_file, scores_csv, ...
                                           t_start, t_end, CODES, varargin)
% ACH_PAC_COMOD_NREM_SEGMENT
% -------------------------------------------------------------
% Compute a NREM-only PAC comodulogram (phase × amplitude)
% for ONE recording segment and ONE condition.
%
%   - phase frequencies: default 0.5:0.5:4 Hz (slow/delta)
%   - amp   frequencies: default 7:1:25 Hz  (sigma-centered)
%   - state: NREM (CODES.NREM)
%
% INPUTS
%   mat_file   : path to .mat with EEG (+ fs_eeg)
%   scores_csv : 1-Hz scoring CSV (time_s,score) or just score
%   t_start    : segment start (s)
%   t_end      : segment end (s)
%   CODES      : struct with fields .WK .NREM .REM .MA
%
% OPTIONAL name/value:
%   'phase_freqs' : vector of phase freqs (Hz), default 0.5:0.5:4
%   'amp_freqs'   : vector of amp freqs (Hz),   default 7:1:25
%   'nbins'       : #phase bins for PAC,       default 18
%   'doPlot'      : true/false, make figure    default true
%   'label'       : string for figure title    default ''
%
% OUTPUT
%   COM.phase_freqs
%   COM.amp_freqs
%   COM.MI         : nP × nA modulation index
%   COM.nSamples   : # of EEG samples used
%   COM.fs         : sampling rate used
%   COM.state      : 'NREM'
%   COM.t_start, COM.t_end
%
%   (if ROI later: you can compute mean MI in e.g. 0.5–1.5 Hz × 11–16 Hz)

% -------------------------------------------------------------
% parse options
p = inputParser;
p.addParameter('phase_freqs', 0.5:0.5:4);   % slow / delta
p.addParameter('amp_freqs',   7:1:25);      % includes sigma band
p.addParameter('nbins',       18);
p.addParameter('doPlot',      true);
p.addParameter('label',       '');
p.parse(varargin{:});
opt = p.Results;

COM = struct('phase_freqs',opt.phase_freqs, ...
             'amp_freqs',  opt.amp_freqs, ...
             'MI',         [], ...
             'nSamples',   0, ...
             'fs',         [], ...
             'state',      'NREM', ...
             't_start',    t_start, ...
             't_end',      t_end);

if t_end <= t_start
    warning('t_end must be > t_start'); 
    return;
end

% -------------------------------------------------------------
% 1) Load EEG + fs from MAT
% -------------------------------------------------------------
info  = whos('-file', mat_file);
names = {info.name};
pick  = @(cands) cands{find(ismember(cands,names),1,'first')};

eeg_name    = pick({'eeg','EEG','eeg1','Eeg','eeg_filt'});
fs_eeg_name = pick({'eeg_frequency','fs_eeg','Fs_eeg','EEG_frequency'});

if isempty(eeg_name) || isempty(fs_eeg_name)
    warning('No EEG or fs_eeg found in %s', mat_file);
    return;
end

S      = load(mat_file, eeg_name, fs_eeg_name);
eeg    = S.(eeg_name)(:);
fs_eeg = S.(fs_eeg_name);

% apply 50 Hz notch if needed (using your helper)
if exist('apply_50Hz_notch_if_needed','file')
    eeg = apply_50Hz_notch_if_needed(eeg, fs_eeg, mat_file);
end

COM.fs = fs_eeg;

% -------------------------------------------------------------
% 2) Load scores (1 Hz)
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

% -------------------------------------------------------------
% 3) Build a mask of NREM samples within [t_start, t_end]
% -------------------------------------------------------------
codeNREM = CODES.NREM;

idx_epoch = (score == codeNREM) & ...
            (t_scores >= t_start) & (t_scores < t_end);

if ~any(idx_epoch)
    warning('No NREM epochs in this window for %s', mat_file);
    return;
end

mask = false(size(eeg));

bins = find(idx_epoch);
for b = bins'
    t0 = t_scores(b);
    t1 = t0 + epoch_sec;
    i1 = max(1, floor(t0 * fs_eeg) + 1);
    i2 = min(numel(eeg), floor(t1 * fs_eeg));
    mask(i1:i2) = true;
end

eeg_nrem = eeg(mask);
COM.nSamples = numel(eeg_nrem);

if numel(eeg_nrem) < fs_eeg * 5
    warning('Too few NREM samples (<5 s) in this segment, skipping PAC.');
    return;
end

% normalize
eeg_nrem = eeg_nrem(:);
eeg_nrem = detrend(eeg_nrem);
eeg_nrem = eeg_nrem / std(eeg_nrem(~isnan(eeg_nrem)));

% -------------------------------------------------------------
% 4) PAC comodulogram (Tort-like MI) - WITH DIAGNOSTICS
% -------------------------------------------------------------
fP = opt.phase_freqs(:)';   % row
fA = opt.amp_freqs(:)';     % row
nP = numel(fP);
nA = numel(fA);
nbins = opt.nbins;

MI = nan(nP, nA);

fprintf('Computing PAC: %d phase × %d amp frequencies\n', nP, nA);
fprintf('EEG segment: %d samples at %.1f Hz\n', numel(eeg_nrem), fs_eeg);

% helper for bandpass
bp = @(x,fl,fh) bandpass_butter(x, fs_eeg, fl, fh);

n_success = 0;
n_fail = 0;

for ip = 1:nP
    fp = fP(ip);
    % small band around phase freq
    fl_p = max(0.1, fp - 0.5);
    fh_p = min(fs_eeg/2 - 1, fp + 0.5);
    
    % Check band validity
    if fh_p <= fl_p
        fprintf('  Phase freq %.2f Hz: invalid band [%.2f %.2f]\n', fp, fl_p, fh_p);
        continue;
    end

    try
        sig_p = bp(eeg_nrem, fl_p, fh_p);
        ph    = angle(hilbert(sig_p));   % phase
    catch ME
        fprintf('  Phase freq %.2f Hz: filtering failed (%s)\n', fp, ME.message);
        n_fail = n_fail + nA;
        continue;
    end
    
    % Check phase validity
    if any(~isfinite(ph))
        fprintf('  Phase freq %.2f Hz: NaN/Inf in phase\n', fp);
        n_fail = n_fail + nA;
        continue;
    end

    for ia = 1:nA
        fa   = fA(ia);
        fl_a = max(0.5, fa - 2);
        fh_a = min(fs_eeg/2 - 1, fa + 2);
        
        if fh_a <= fl_a
            n_fail = n_fail + 1;
            continue;
        end

        try
            sig_a = bp(eeg_nrem, fl_a, fh_a);
            amp   = abs(hilbert(sig_a)); % amplitude envelope
        catch ME
            n_fail = n_fail + 1;
            continue;
        end
        
        % Check amplitude validity
        if any(~isfinite(amp)) || all(amp < eps)
            n_fail = n_fail + 1;
            continue;
        end

        % Compute MI
        try
            mi_val = tort_MI(ph, amp, nbins);
            if isfinite(mi_val)
                MI(ip, ia) = mi_val;
                n_success = n_success + 1;
            else
                n_fail = n_fail + 1;
            end
        catch ME
            n_fail = n_fail + 1;
        end
    end
end

fprintf('PAC computation: %d successful, %d failed\n', n_success, n_fail);

if n_success == 0
    warning('All PAC computations failed! Check EEG signal quality and frequency bands.');
end

COM.MI = MI;
COM.phase_freqs = fP;
COM.amp_freqs   = fA;

% -------------------------------------------------------------
% 5) Optional figure
% -------------------------------------------------------------
if opt.doPlot
    figure('Color','w','Position',[200 200 600 450]);
    imagesc(fA, fP, MI); axis xy;
    xlabel('Amplitude frequency (Hz)');
    ylabel('Phase frequency (Hz)');
    titleStr = 'NREM PAC comodulogram (EEG)';
    if ~isempty(opt.label)
        titleStr = sprintf('%s – %s', titleStr, opt.label);
    end
    title(titleStr, 'Interpreter','none');
    colormap(turbo);
    c = colorbar;
    ylabel(c,'Modulation Index (Tort)');
end
end

% ==== helper: bandpass with adaptive filtering =====================
function y = bandpass_butter(x, fs, f1, f2)
    % For very low frequencies, use FIR; otherwise use IIR Butterworth
    
    nyq = fs / 2;
    
    % Adjust frequencies if needed
    if f2 >= nyq - 0.5
        f2 = nyq - 0.5;
    end
    if f1 <= 0.1
        f1 = 0.1;
    end
    if f2 <= f1
        error('Invalid frequency band: f2 must be > f1');
    end
    
    x = double(x(:));
    x = detrend(x);  % Remove DC
    
    % Decide: use FIR for very low frequencies to avoid numerical issues
    % Use FIR if f2 < 2 Hz (slow oscillations)
    use_fir = (f2 < 2.0);
    
    if use_fir
        % === FIR filter (more stable for low frequencies) ===
        
        % Design FIR filter using window method
        % Filter order: make it long enough for low frequencies
        % Rule of thumb: order ~ 3*fs/f1
        fir_order = round(3 * fs / f1);
        fir_order = min(fir_order, 10000);  % cap at reasonable value
        fir_order = max(fir_order, 100);     % minimum order
        
        % Make order even
        if mod(fir_order, 2) == 1
            fir_order = fir_order + 1;
        end
        
        try
            % Normalized frequencies
            Wn = [f1 f2] / nyq;
            
            % Design FIR filter
            b_fir = fir1(fir_order, Wn, 'bandpass', hamming(fir_order+1));
            
            % Apply filter (use filtfilt for zero-phase)
            y = filtfilt(b_fir, 1, x);
            
        catch ME
            warning('FIR filter failed for [%.2f %.2f] Hz: %s', f1, f2, ME.message);
            y = nan(size(x));
            return;
        end
        
    else
        % === IIR Butterworth (efficient for higher frequencies) ===
        
        Wn = [f1 f2] / nyq;
        
        % Check Wn validity
        if any(Wn <= 0) || any(Wn >= 1)
            warning('Normalized frequencies out of range: [%.4f %.4f]', Wn(1), Wn(2));
            y = nan(size(x));
            return;
        end
        
        try
            [b, a] = butter(4, Wn, 'bandpass');
            y = filtfilt(b, a, x);
        catch ME
            warning('IIR filter failed for [%.2f %.2f] Hz: %s', f1, f2, ME.message);
            y = nan(size(x));
            return;
        end
    end
    
    % Check output
    if any(~isfinite(y))
        warning('Filtered signal contains %d NaN/Inf values', sum(~isfinite(y)));
        y = nan(size(x));
    end
end
% ==== helper: Tort MI (ROBUST VERSION) =============================
function mi = tort_MI(phase, amp, nbins)
    if nargin < 3
        nbins = 18;
    end
    
    phase = phase(:);
    amp   = amp(:);
    
    % Remove NaN/Inf
    valid = isfinite(phase) & isfinite(amp);
    if sum(valid) < 100  % Need enough samples
        mi = NaN;
        return;
    end
    
    phase = phase(valid);
    amp = amp(valid);
    
    % Create phase bins
    edges = linspace(-pi, pi, nbins+1);
    [~, bin] = histc(phase, edges);
    
    % Fix edge case: values exactly at pi
    bin(bin == nbins+1) = nbins;
    bin(bin == 0) = 1;
    
    % Calculate mean amplitude per bin
    P = zeros(1, nbins);
    for k = 1:nbins
        idx = (bin == k);
        if any(idx)
            P(k) = mean(amp(idx));
        else
            P(k) = 0;
        end
    end
    
    % Check if all bins are empty
    if sum(P) == 0 || all(P == 0)
        mi = NaN;
        return;
    end
    
    % Normalize to probability distribution
    P = P / sum(P);
    
    % Remove zeros to avoid log(0)
    P(P == 0) = eps;
    
    % Calculate entropy
    H = -sum(P .* log(P));
    Hmax = log(nbins);
    
    % Modulation Index
    mi = (Hmax - H) / Hmax;
    
    % Sanity check
    if ~isfinite(mi) || mi < 0 || mi > 1
        mi = NaN;
    end
end
