function cycleHz = compute_cycle_freq(OUT, state)
% Compute ACh cycle frequency from peak-to-peak intervals in raw signal
% 
% INPUTS:
%   OUT   - output structure from ach_analysis
%   state - 'Wake', 'NREM', or 'REM'
%
% OUTPUT:
%   cycleHz - average cycle frequency in Hz (NaN if cannot compute)

    cycleHz = NaN;
    
    % Check if state signal exists
    if ~isfield(OUT, 'sig_state') || ~isfield(OUT.sig_state, state)
        return;
    end
    
    sig = OUT.sig_state.(state);
    if isempty(sig) || numel(sig) < 20
        return;
    end
    
    % Sampling rate
    if isfield(OUT, 'fs')
        fs = OUT.fs;
    else
        warning('OUT.fs missing, assuming 500 Hz');
        fs = 500;
    end
    
    % Detrend + normalize
    sig = detrend(sig(:));
    sig = sig / std(sig);

    % Detect peaks
    [pks, locs] = findpeaks(sig, ...
        'MinPeakProminence', 0.5, ...
        'MinPeakDistance', round(0.25 * fs));  % max 4 Hz
    
    if numel(locs) < 3
        return;
    end
    
    % Peak-to-peak intervals
    intervals = diff(locs) / fs;

    % Remove impossible intervals (noise)
    intervals = intervals(intervals > 0.1 & intervals < 60); % 0.1–60 sec
    
    if isempty(intervals)
        return;
    end

    avg_period = mean(intervals);
    cycleHz = 1 / avg_period;

    % sanity range: ACh rises are ultradian/infraslow
    if cycleHz > 5 || cycleHz < 0.001
        cycleHz = NaN;
    end
end
