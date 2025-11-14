function [f, psd_raw, band_power, peak_f, peak_amp] = ...
    ach_psd_state(sig, fs, fmin, fmax, win_sec, verbose, label)
% Compute PSD of slow ACh signal in [fmin, fmax] Hz during one state.
% Steps:
%   1) Remove constant offset (detrend)
%   2) Welch PSD
%   3) Keep only f in [fmin, fmax] (e.g. 0.01–0.15 Hz)
%   4) Summarise total power, peak freq, peak height.

if isempty(sig)
    f = []; psd_raw = []; band_power = NaN; peak_f = NaN; peak_amp = NaN;
    return;
end

% Remove DC offset so 0 Hz doesn't dominate
sig = detrend(sig, 'constant');

% Welch parameters
L       = numel(sig);
win_len = min(round(win_sec*fs), L);
if win_len < 10
    win_len = L;
end
win   = hamming(win_len);
nover = floor(win_len/2);
nfft  = max(512, 2^nextpow2(win_len));

[psd_all, f_all] = pwelch(sig, win, nover, nfft, fs);

% Restrict to [fmin, fmax]
idx = (f_all >= fmin) & (f_all <= fmax);
f       = f_all(idx);
psd_raw = psd_all(idx);

if isempty(f)
    % No bins in requested range
    band_power = NaN; peak_f = NaN; peak_amp = NaN;
    return;
end

% Total power in band
band_power = trapz(f, psd_raw);

% Peak in that band
[peak_amp, k] = max(psd_raw);
peak_f = f(k);

if verbose
    fprintf('%s: power=%.4g, peak_f=%.4f Hz, peak_amp=%.4g\n', ...
        label, band_power, peak_f, peak_amp);
end
end
