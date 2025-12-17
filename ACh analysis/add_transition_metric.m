function m = add_transition_metric(m, OUT, trans_name, prefix)
% Add mean peak and slope for a given transition type to metrics struct m.
%
% Example:
%   m = add_transition_metric(m, OUT, 'Wake_onset', 'WakeOn');
%
% This will add fields to m:
%   m.WakeOn_peak_mean
%   m.WakeOn_slope_mean
%
% It expects OUT.transitions(i) to have:
%   .name   = 'Wake_onset' / 'NREM_onset' / 'REM_onset'
%   .peaks  = per-bout peak dF/F values (vector)
%   .slopes = per-bout slopes (optional)

idx = [];
if isfield(OUT,'transitions')
    idx = find(strcmp({OUT.transitions.name}, trans_name));
end

if ~isempty(idx)
    T = OUT.transitions(idx);

    % Peak: average across bouts for this session
    if isfield(T,'peaks') && ~isempty(T.peaks)
        m.([prefix '_peak_mean']) = mean(T.peaks, 'omitnan');
    else
        m.([prefix '_peak_mean']) = NaN;
    end

    % Slope: average across bouts for this session
    if isfield(T,'slopes') && ~isempty(T.slopes)
        m.([prefix '_slope_mean']) = mean(T.slopes, 'omitnan');
    else
        m.([prefix '_slope_mean']) = NaN;
    end
else
    % If this transition type doesn't exist, fill with NaNs
    m.([prefix '_peak_mean'])  = NaN;
    m.([prefix '_slope_mean']) = NaN;
end
end
