function [t_events, idx_events] = find_transitions(t_sec, code, from_codes, to_code)
% Find times of transitions where the 1-Hz score changes from
% any of from_codes to to_code.
from_codes = from_codes(:)';
code = code(:);
d = diff(code);
idx = find(ismember(code(1:end-1), from_codes) & code(2:end) == to_code) + 1;
t_events   = t_sec(idx);   % time (s) of first second in new state
idx_events = idx;
end

