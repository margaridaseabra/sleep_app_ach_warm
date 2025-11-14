function [t_rel, traces] = extract_peri_traces(sig, fs, t_events, win)
% Extract peri-event segments from a continuous signal.
% sig   : continuous ACh trace
% fs    : sampling rate (Hz)
% t_events : vector of event times (s)
% win      : [pre post] (s), e.g. [-100 100]
pre  = win(1);
post = win(2);

nT = numel(t_events);
% Use the first event to define the common time axis
L    = round((post - pre)*fs) + 1;
t_rel= (0:L-1)'/fs + pre;   % relative time axis (s)

traces = nan(L, nT);
for i = 1:nT
    t0 = t_events(i);
    idx0 = round((t0 + pre)*fs) + 1;
    idx1 = idx0 + L - 1;
    if idx0 < 1 || idx1 > numel(sig)
        continue; % should be filtered out earlier by "keep", but just in case
    end
    seg = sig(idx0:idx1);
    traces(:,i) = seg;
end
end

