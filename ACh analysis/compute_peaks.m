function peaks = compute_peaks(traces, t_rel, resp_win)
% Compute peak value in response window for each event.
mask_resp = (t_rel >= resp_win(1)) & (t_rel <= resp_win(2));
peaks = nan(1,size(traces,2));
for i = 1:size(traces,2)
    x = traces(:,i);
    if all(isnan(x)), continue; end
    peaks(i) = max(x(mask_resp));
end
peaks = peaks(:);
end

