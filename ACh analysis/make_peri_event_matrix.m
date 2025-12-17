function [M, t_rel, event_times_valid] = make_peri_event_matrix(sig, fs, event_times, t_pre, t_post)
% MAKE_PERI_EVENT_MATRIX
%   Build event x time matrix for peri-event analysis.
%
%   IMPORTANT: this version KEEPS ALL EVENTS and pads missing
%   parts of the window with NaNs instead of discarding events
%   near the start/end of the recording.
%
% INPUTS
%   sig         : 1D signal (vector)
%   fs          : sampling rate (Hz)
%   event_times : vector of event times in seconds
%   t_pre       : seconds before event (positive)
%   t_post      : seconds after event (positive)
%
% OUTPUTS
%   M                : [nEvents x nTime] matrix (NaN-padded)
%   t_rel            : relative time axis (seconds)
%   event_times_valid: event times actually used (all events)

sig = sig(:);                 % force column
T_total = numel(sig) / fs;    %#ok<NASGU>

n_samp_pre  = round(t_pre * fs);
n_samp_post = round(t_post * fs);
n_time      = n_samp_pre + n_samp_post + 1;

t_rel = (-n_samp_pre:n_samp_post) / fs;

event_times_valid = event_times(:);
n_events = numel(event_times_valid);

M = nan(n_events, n_time);

for k = 1:n_events
    t0   = event_times_valid(k);
    idx0 = round(t0 * fs) + 1;     % centre sample

    % desired indices in the signal
    i1 = idx0 - n_samp_pre;
    i2 = idx0 + n_samp_post;

    % clamp to signal boundaries
    j1 = max(1, i1);
    j2 = min(numel(sig), i2);

    % corresponding indices in the row
    r1 = 1 + (j1 - i1);
    r2 = n_time - (i2 - j2);

    row = nan(1, n_time);
    row(r1:r2) = sig(j1:j2);

    M(k,:) = row;
end
end
