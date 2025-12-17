function val = get_state_slope(OUT, stateName)
% Gracefully get mean slope for a state from OUT.state_slopes
%
% OUT.state_slopes.(stateName).mean should exist if slopes were computed
% for that state. If not, return NaN so group stats don't crash.

if isfield(OUT,'state_slopes') && isfield(OUT.state_slopes, stateName) ...
        && isfield(OUT.state_slopes.(stateName), 'mean')

    val = OUT.state_slopes.(stateName).mean;
else
    val = NaN;
end
end
