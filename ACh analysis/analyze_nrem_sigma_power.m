function [SIGMA_RESULTS] = analyze_nrem_sigma_power(mat_file, scores_csv, ...
                                                     t_start, t_end, CODES, varargin)
% ANALYZE_NREM_SIGMA_POWER
% Compute sigma band (spindle) power during NREM sleep including microarousals
%
% INPUTS:
%   mat_file   : path to .mat file with EEG signal
%   scores_csv : path to 1-Hz scoring CSV
%   t_start    : segment start time (sec)
%   t_end      : segment end time (sec)
%   CODES      : struct with .WK .NREM .REM .MA
%
% OPTIONAL name/value:
%   'sigma_band'       : [11 16] Hz (spindle band for mice)
%   'delta_band'       : [1 4] Hz (for delta/sigma ratio)
%   'include_MA'       : true/false, include microarousals (default true)
%   'min_nrem_sec'     : minimum NREM duration to include (default 5 sec)
%   'nfft'             : FFT length for PSD (default 4096)
%   'window_sec'       : window length for Welch (default 4 sec)
%   'overlap_frac'     : overlap fraction (default 0.5)
%
% OUTPUT:
%   SIGMA_RESULTS : struct with fields
%       .sigma_power        : mean sigma power (µV²)
%       .delta_power        : mean delta power (µV²)
%       .delta_sigma_ratio  : delta/sigma ratio
%       .psd_freqs          : frequency vector
%       .psd_power          : full PSD
%       .total_nrem_sec     : total NREM time analyzed
%       .n_ma_epochs        : number of microarousal epochs included
%       .nrem_only_sigma    : sigma power excluding MA epochs
%       .nrem_with_ma_sigma : sigma power including MA epochs

% Parse inputs
p = inputParser;
p.addParameter('sigma_band', [11 16]);
p.addParameter('delta_band', [1 4]);
p.addParameter('include_MA', true);
p.addParameter('min_nrem_sec', 5);
p.addParameter('nfft', 4096);
p.addParameter('window_sec', 4);
p.addParameter('overlap_frac', 0.5);
p.parse(varargin{:});
opt = p.Results;

% Initialize output
SIGMA_RESULTS = struct();
SIGMA_RESULTS.sigma_power = NaN;
SIGMA_RESULTS.delta_power = NaN;
SIGMA_RESULTS.delta_sigma_ratio = NaN;
SIGMA_RESULTS.psd_freqs = [];
SIGMA_RESULTS.psd_power = [];
SIGMA_RESULTS.total_nrem_sec = 0;
SIGMA_RESULTS.n_ma_epochs = 0;
SIGMA_RESULTS.nrem_only_sigma = NaN;
SIGMA_RESULTS.nrem_with_ma_sigma = NaN;

fprintf('\n=== NREM Sigma Band Analysis ===\n');

% ----------------------------------------------------------------
% 1) Load EEG signal
% ----------------------------------------------------------------
info = whos('-file', mat_file);
names = {info.name};
pick = @(cands) cands{find(ismember(cands,names),1,'first')};

eeg_name = pick({'eeg','EEG','eeg1','Eeg','eeg_filt'});
fs_name = pick({'eeg_frequency','fs_eeg','Fs_eeg','EEG_frequency'});

if isempty(eeg_name) || isempty(fs_name)
    warning('No EEG or sampling rate found in %s', mat_file);
    return;
end

S = load(mat_file, eeg_name, fs_name);
eeg = S.(eeg_name)(:);
fs = S.(fs_name);

fprintf('Loaded EEG: %d samples at %.1f Hz\n', numel(eeg), fs);

% ----------------------------------------------------------------
% 2) Load sleep scores
% ----------------------------------------------------------------
M = readmatrix(scores_csv);
if size(M,2) == 1
    scores = M(:,1);
    epoch_sec = 1;
    t_scores = (0:numel(scores)-1)' * epoch_sec;
else
    t_scores = M(:,1);
    scores = M(:,2);
    dt = diff(t_scores);
    epoch_sec = mode(dt(~isnan(dt)));
end

fprintf('Loaded scores: %d epochs at %d-sec resolution\n', ...
        numel(scores), epoch_sec);

% ----------------------------------------------------------------
% 3) Find NREM epochs (with or without MA)
% ----------------------------------------------------------------
nrem_code = CODES.NREM;
ma_code = CODES.MA;

% NREM epochs
idx_nrem = (scores == nrem_code) & ...
           (t_scores >= t_start) & (t_scores < t_end);

% MA epochs
idx_ma = (scores == ma_code) & ...
         (t_scores >= t_start) & (t_scores < t_end);

if opt.include_MA
    idx_combined = idx_nrem | idx_ma;
    fprintf('Including microarousals: %d NREM + %d MA = %d total epochs\n', ...
            sum(idx_nrem), sum(idx_ma), sum(idx_combined));
else
    idx_combined = idx_nrem;
    fprintf('Excluding microarousals: %d NREM epochs only\n', sum(idx_nrem));
end

SIGMA_RESULTS.n_ma_epochs = sum(idx_ma);

if ~any(idx_combined)
    warning('No NREM epochs found in window [%.1f, %.1f]', t_start, t_end);
    return;
end

% ----------------------------------------------------------------
% 4) Extract EEG for NREM+MA epochs
% ----------------------------------------------------------------
mask = false(size(eeg));
bins = find(idx_combined);

for b = bins'
    t0 = t_scores(b);
    t1 = t0 + epoch_sec;
    i1 = max(1, floor(t0 * fs) + 1);
    i2 = min(numel(eeg), floor(t1 * fs));
    mask(i1:i2) = true;
end

eeg_nrem = eeg(mask);
total_sec = numel(eeg_nrem) / fs;

fprintf('Extracted %.1f sec of NREM (%.1f min)\n', total_sec, total_sec/60);

SIGMA_RESULTS.total_nrem_sec = total_sec;

if total_sec < opt.min_nrem_sec
    warning('Insufficient NREM data (%.1f sec < %.1f sec minimum)', ...
            total_sec, opt.min_nrem_sec);
    return;
end

% Normalize EEG
eeg_nrem = eeg_nrem(:);
eeg_nrem = detrend(eeg_nrem);
eeg_nrem = eeg_nrem / std(eeg_nrem(~isnan(eeg_nrem)));

% ----------------------------------------------------------------
% 5) Compute Power Spectral Density (Welch's method)
% ----------------------------------------------------------------
window_samples = round(opt.window_sec * fs);
overlap_samples = round(window_samples * opt.overlap_frac);
nfft = opt.nfft;

fprintf('Computing PSD: window=%.1fs, overlap=%.1f%%, nfft=%d\n', ...
        opt.window_sec, opt.overlap_frac*100, nfft);

[psd, freqs] = pwelch(eeg_nrem, window_samples, overlap_samples, nfft, fs);

SIGMA_RESULTS.psd_freqs = freqs;
SIGMA_RESULTS.psd_power = psd;

fprintf('PSD computed: %d frequency bins from %.2f to %.2f Hz\n', ...
        numel(freqs), freqs(1), freqs(end));

% ----------------------------------------------------------------
% 6) Extract band powers
% ----------------------------------------------------------------
% Sigma band (spindles)
sigma_idx = freqs >= opt.sigma_band(1) & freqs <= opt.sigma_band(2);
sigma_power = mean(psd(sigma_idx));

% Delta band
delta_idx = freqs >= opt.delta_band(1) & freqs <= opt.delta_band(2);
delta_power = mean(psd(delta_idx));

% Ratio
delta_sigma_ratio = delta_power / sigma_power;

SIGMA_RESULTS.sigma_power = sigma_power;
SIGMA_RESULTS.delta_power = delta_power;
SIGMA_RESULTS.delta_sigma_ratio = delta_sigma_ratio;

fprintf('\n--- Band Powers ---\n');
fprintf('Sigma (%.1f-%.1f Hz): %.4e µV²\n', ...
        opt.sigma_band(1), opt.sigma_band(2), sigma_power);
fprintf('Delta (%.1f-%.1f Hz): %.4e µV²\n', ...
        opt.delta_band(1), opt.delta_band(2), delta_power);
fprintf('Delta/Sigma ratio: %.2f\n', delta_sigma_ratio);

% ----------------------------------------------------------------
% 7) Optional: Separate analysis for NREM-only vs NREM+MA
% ----------------------------------------------------------------
if opt.include_MA && sum(idx_ma) > 0
    % NREM only
    mask_nrem_only = false(size(eeg));
    bins_nrem = find(idx_nrem);
    for b = bins_nrem'
        t0 = t_scores(b);
        t1 = t0 + epoch_sec;
        i1 = max(1, floor(t0 * fs) + 1);
        i2 = min(numel(eeg), floor(t1 * fs));
        mask_nrem_only(i1:i2) = true;
    end
    
    eeg_nrem_only = eeg(mask_nrem_only);
    if numel(eeg_nrem_only) > fs * opt.min_nrem_sec
        eeg_nrem_only = detrend(eeg_nrem_only(:));
        eeg_nrem_only = eeg_nrem_only / std(eeg_nrem_only);
        
        [psd_nrem, ~] = pwelch(eeg_nrem_only, window_samples, ...
                               overlap_samples, nfft, fs);
        SIGMA_RESULTS.nrem_only_sigma = mean(psd_nrem(sigma_idx));
    end
    
    % Store the combined result
    SIGMA_RESULTS.nrem_with_ma_sigma = sigma_power;
    
    fprintf('\nComparison:\n');
    fprintf('  NREM only sigma: %.4e µV²\n', SIGMA_RESULTS.nrem_only_sigma);
    fprintf('  NREM+MA sigma:   %.4e µV²\n', SIGMA_RESULTS.nrem_with_ma_sigma);
end

fprintf('=== Analysis Complete ===\n\n');

end