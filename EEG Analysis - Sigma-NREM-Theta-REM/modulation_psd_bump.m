function [f_mod, psd_whiten, psd_plot, total_power, bump_power, peak_f, peak_amp] = ...
    modulation_psd_bump(eeg_state, fs, carrier_band, mod_f_min, mod_f_max, ...
                        env_target_fs, env_win_sec, verbose, label)
% MODULATION_PSD_BUMP
%   - Band-pass EEG in carrier band (sigma/theta)
%   - Hilbert envelope
%   - Downsample envelope
%   - Polynomial detrend (remove slow drift)
%   - Welch PSD of envelope (slow modulations)
%   - Whiten PSD in log–log space (remove 1/f slope)
%   - Return:
%       f_mod         : modulation frequency axis
%       psd_whiten    : whitened PSD (for metrics)
%       psd_plot      : normalized PSD (for plotting)
%       total_power   : area under whitened PSD
%       bump_power    : area around the main bump
%       peak_f        : bump frequency
%       peak_amp      : bump height (whitened PSD)

% 1) Band-pass filter EEG in carrier band
Wn = carrier_band / (fs/2);
[b,a] = butter(4, Wn);                        % 4th-order Butterworth
eeg_bp = filtfilt(b,a,eeg_state);

% 2) Envelope (instantaneous amplitude)
env = abs(hilbert(eeg_bp));

% 3) Downsample envelope
decim = max(1, floor(fs / env_target_fs));
env_ds = decimate(env, decim);
fs_env = fs / decim;

% 4) Remove slow polynomial trend (similar to lab template)
x = (1:numel(env_ds))';
[p,~,mu] = polyfit(x, env_ds, 5);             % 5th-order trend
trend    = polyval(p, x, [], mu);
env_ds_detr = env_ds - trend;

if verbose
    fprintf('%s: env_fs = %.2f Hz, length = %.1f s\n', ...
        label, fs_env, numel(env_ds_detr)/fs_env);
end

% 5) Welch PSD of detrended envelope
L       = numel(env_ds_detr);
win_len = min(round(env_win_sec*fs_env), L);
if win_len < 10
    win_len = L;
end
win   = hamming(win_len);
nover = floor(win_len/2);
nfft  = max(512, 2^nextpow2(win_len));

[psd_all, f_all] = pwelch(env_ds_detr, win, nover, nfft, fs_env);

% 6) Frequency range of interest
idx = (f_all >= mod_f_min) & (f_all <= mod_f_max);
f_mod    = f_all(idx);
psd_raw  = psd_all(idx);

if isempty(f_mod)
    error('%s: no frequency bins between %.3f and %.3f Hz.', ...
          label, mod_f_min, mod_f_max);
end

% 7) Whiten PSD in log–log space (remove 1/f-like slope)
logf = log10(f_mod(:));
logp = log10(psd_raw(:));

% fit straight line logP ≈ a*logf + b
coef   = polyfit(logf, logp, 1);
trend  = polyval(coef, logf);
logp_d = logp - trend;                % detrended (whitened)
psd_whiten = 10.^logp_d;              % back to linear scale

% 8) Summary metrics on whitened PSD
total_power = trapz(f_mod, psd_whiten);

[peak_amp, k] = max(psd_whiten);
peak_f = f_mod(k);

% bump power in ±0.01 Hz window around peak, clipped to [mod_f_min, mod_f_max]
bump_half_width = 0.01;  % adjust if needed
idx_bump = (f_mod >= max(mod_f_min, peak_f - bump_half_width)) & ...
           (f_mod <= min(mod_f_max, peak_f + bump_half_width));

bump_power = trapz(f_mod(idx_bump), psd_whiten(idx_bump));

if verbose
    fprintf('%s: total=%.4g, bump=%.4g, peak_f=%.4f Hz, peak_amp=%.4g\n', ...
        label, total_power, bump_power, peak_f, peak_amp);
end

% 9) Version for plotting: normalize per mouse
psd_plot = psd_whiten / max(psd_whiten);   % or / trapz(f_mod, psd_whiten)

end
