function slopes = compute_slopes(traces, t_rel, win)
% Compute per-event linear slope of baseline-corrected traces
% within a specified time window WIN = [t_start t_end] (s)
%
% traces: [nTime x nEvents]
% t_rel : [nTime x 1] time vector (s)

if isempty(traces)
    slopes = [];
    return;
end

t = t_rel(:);
mask = (t >= win(1)) & (t <= win(2));
t_win = t(mask);

nEv = size(traces,2);
slopes = nan(nEv,1);

for j = 1:nEv
    y = traces(mask,j);
    if all(isnan(y)) || numel(y) < 3
        slopes(j) = NaN;
        continue;
    end
    p  = polyfit(t_win, y, 1);   % first-order fit: y = p(1)*t + p(2)
    slopes(j) = p(1);           % slope (dF/F per second)
end
end
