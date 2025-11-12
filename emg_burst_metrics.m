function [rate_per_min, pct] = emg_burst_metrics(emg, fs)
% Bandpass (20–100 Hz), rectify via Hilbert, smooth 100 ms, robust threshold
emg_bp = bandpass_zero(emg, fs, [20 100]);
env = abs(hilbert(emg_bp));
w = max(1, round(0.10*fs));                  % 100 ms smoothing
env_s = movmean(env, w);

% --- Toolbox-free MAD (unscaled): same as mad(x,1)
med_env = median(env_s);
mad_env = median(abs(env_s - med_env));      % unscaled MAD
thr = med_env + 3 * mad_env;

bin = env_s > thr;

% Clean up bursts: remove <50 ms, merge gaps <50 ms
min_dur = round(0.050*fs);
min_gap = round(0.050*fs);
bin = logical(runlength_prune(bin, min_dur, min_gap));

% Burst metrics
[pct, n_bursts] = burst_stats(bin, fs);
dur_min = numel(emg)/fs/60;
if dur_min>0
    rate_per_min = n_bursts / dur_min;
else
    rate_per_min = NaN;
end
end
