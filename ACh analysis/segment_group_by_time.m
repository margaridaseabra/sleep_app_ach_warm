function GROUP_win = segment_group_by_time(GROUP, start_sec, end_sec)
% Segment GROUP data to only include a specific time window
%
% INPUTS:
%   GROUP     - GROUP structure from run_ach_batch_auto
%   start_sec - start time in seconds
%   end_sec   - end time in seconds (can be Inf for rest of recording)
%
% OUTPUT:
%   GROUP_win - GROUP structure with data only from the time window

GROUP_win = GROUP;
nSess = numel(GROUP.sessions);

for k = 1:nSess
    OUT = GROUP.out{k};
    
    % Get sampling frequency
    if isfield(OUT, 'fs')
        fs = OUT.fs;
    elseif isfield(OUT, 'Fs')
        fs = OUT.Fs;
    else
        fs = 1; % Assume 1 Hz if not found
    end
    
    % Convert time to samples
    start_samp = round(start_sec * fs) + 1;
    if isinf(end_sec)
        % Get total length from signal or scores
        if isfield(OUT, 'sig') && ~isempty(OUT.sig)
            end_samp = numel(OUT.sig);
        elseif isfield(OUT, 'scores') && ~isempty(OUT.scores)
            end_samp = numel(OUT.scores);
        else
            warning('Cannot determine recording length for session %d', k);
            continue;
        end
    else
        end_samp = round(end_sec * fs);
    end
    
    % Segment signal
    if isfield(OUT, 'sig') && ~isempty(OUT.sig)
        n_total = numel(OUT.sig);
        end_samp = min(end_samp, n_total);
        if start_samp <= n_total && start_samp < end_samp
            OUT.sig = OUT.sig(start_samp:end_samp);
        else
            OUT.sig = [];
        end
    end
    
    % Segment scores (typically 1 Hz)
    if isfield(OUT, 'scores') && ~isempty(OUT.scores)
        score_start = floor(start_sec) + 1;
        score_end = floor(end_sec);
        if isinf(end_sec)
            score_end = numel(OUT.scores);
        end
        score_end = min(score_end, numel(OUT.scores));
        
        if score_start <= numel(OUT.scores) && score_start < score_end
            OUT.scores = OUT.scores(score_start:score_end);
        else
            OUT.scores = [];
        end
    end
    
    % Segment time vector if present
    if isfield(OUT, 't') && ~isempty(OUT.t)
        if start_samp <= numel(OUT.t) && end_samp <= numel(OUT.t)
            OUT.t = OUT.t(start_samp:end_samp);
            % Adjust to start from 0
            OUT.t = OUT.t - OUT.t(1);
        end
    end
    
    % Segment state-specific signals
    if isfield(OUT, 'sig_state')
        states = fieldnames(OUT.sig_state);
        for s = 1:numel(states)
            state = states{s};
            if ~isempty(OUT.sig_state.(state))
                % For state signals, we need to be more careful
                % These are already filtered by state, so just keep them
                % A more sophisticated approach would track state timing
                % For now, keep them as-is with a warning
                fprintf('  Note: State-specific signals for %s not time-windowed\n', state);
            end
        end
    end
    
    % Segment transitions
    if isfield(OUT, 'transitions') && ~isempty(OUT.transitions)
        for t = 1:numel(OUT.transitions)
            trans = OUT.transitions(t);
            
            % Filter transition times that fall within window
            if isfield(trans, 'times') && ~isempty(trans.times)
                valid_idx = trans.times >= start_sec & trans.times < end_sec;
                
                % Update all transition fields
                if isfield(trans, 'traces') && ~isempty(trans.traces)
                    OUT.transitions(t).traces = trans.traces(:, valid_idx);
                end
                if isfield(trans, 'peaks') && ~isempty(trans.peaks)
                    OUT.transitions(t).peaks = trans.peaks(valid_idx);
                end
                if isfield(trans, 'slopes') && ~isempty(trans.slopes)
                    OUT.transitions(t).slopes = trans.slopes(valid_idx);
                end
                if isfield(trans, 'times')
                    OUT.transitions(t).times = trans.times(valid_idx) - start_sec;
                end
                
                % Recompute mean if traces exist
                if isfield(trans, 'traces') && ~isempty(OUT.transitions(t).traces)
                    OUT.transitions(t).mean = mean(OUT.transitions(t).traces, 2, 'omitnan');
                    if isfield(OUT.transitions(t), 'peaks')
                        OUT.transitions(t).peak_mean = mean(OUT.transitions(t).peaks, 'omitnan');
                    end
                end
            end
        end
    end
    
    % Note: PSD analysis would need to be recomputed for the windowed data
    % For now, we keep the original PSD but add a warning
    if isfield(OUT, 'psd')
        fprintf('  Warning: PSD data not recomputed for time window (session %d)\n', k);
    end
    
    % Update the GROUP structure
    GROUP_win.out{k} = OUT;
end

% Recompute metrics for the windowed data
fprintf('Recomputing metrics for time-windowed data...\n');
GROUP_win = recompute_windowed_metrics(GROUP_win);

end

%% Helper function to recompute metrics
function GROUP = recompute_windowed_metrics(GROUP)
% Recompute transition metrics from windowed data

nSess = numel(GROUP.sessions);

for k = 1:nSess
    OUT = GROUP.out{k};
    
    % Update transition metrics
    % Wake onset
    idxW = [];
    if isfield(OUT, 'transitions')
        idxW = find(strcmp({OUT.transitions.name}, 'Wake_onset'));
    end
    if ~isempty(idxW)
        Tw = OUT.transitions(idxW);
        if isfield(Tw, 'peaks') && ~isempty(Tw.peaks)
            GROUP.metrics(k).WakeOn_peak_mean = mean(Tw.peaks, 'omitnan');
        else
            GROUP.metrics(k).WakeOn_peak_mean = NaN;
        end
        if isfield(Tw, 'slopes') && ~isempty(Tw.slopes)
            GROUP.metrics(k).WakeOn_slope_mean = mean(Tw.slopes, 'omitnan');
        else
            GROUP.metrics(k).WakeOn_slope_mean = NaN;
        end
    else
        GROUP.metrics(k).WakeOn_peak_mean = NaN;
        GROUP.metrics(k).WakeOn_slope_mean = NaN;
    end
    
    % NREM onset
    idxN = [];
    if isfield(OUT, 'transitions')
        idxN = find(strcmp({OUT.transitions.name}, 'NREM_onset'));
    end
    if ~isempty(idxN)
        Tn = OUT.transitions(idxN);
        if isfield(Tn, 'peaks') && ~isempty(Tn.peaks)
            GROUP.metrics(k).NREMOn_peak_mean = mean(Tn.peaks, 'omitnan');
        else
            GROUP.metrics(k).NREMOn_peak_mean = NaN;
        end
        if isfield(Tn, 'slopes') && ~isempty(Tn.slopes)
            GROUP.metrics(k).NREMOn_slope_mean = mean(Tn.slopes, 'omitnan');
        else
            GROUP.metrics(k).NREMOn_slope_mean = NaN;
        end
    else
        GROUP.metrics(k).NREMOn_peak_mean = NaN;
        GROUP.metrics(k).NREMOn_slope_mean = NaN;
    end
    
    % REM onset
    idxR = [];
    if isfield(OUT, 'transitions')
        idxR = find(strcmp({OUT.transitions.name}, 'REM_onset'));
    end
    if ~isempty(idxR)
        Tr = OUT.transitions(idxR);
        if isfield(Tr, 'peaks') && ~isempty(Tr.peaks)
            GROUP.metrics(k).REMOn_peak_mean = mean(Tr.peaks, 'omitnan');
        else
            GROUP.metrics(k).REMOn_peak_mean = NaN;
        end
        if isfield(Tr, 'slopes') && ~isempty(Tr.slopes)
            GROUP.metrics(k).REMOn_slope_mean = mean(Tr.slopes, 'omitnan');
        else
            GROUP.metrics(k).REMOn_slope_mean = NaN;
        end
    else
        GROUP.metrics(k).REMOn_peak_mean = NaN;
        GROUP.metrics(k).REMOn_slope_mean = NaN;
    end
end

% Update metrics table
GROUP.metrics_tbl = struct2table(GROUP.metrics);

end