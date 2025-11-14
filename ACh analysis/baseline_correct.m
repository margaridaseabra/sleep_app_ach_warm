function traces_bc = baseline_correct(traces, t_rel, base_win)
% Baseline-correct each peri-event trace: subtract mean in base_win.
L = size(traces,1);
traces_bc = traces;
mask_base = (t_rel >= base_win(1)) & (t_rel <= base_win(2));
for i = 1:size(traces,2)
    x = traces(:,i);
    if all(isnan(x)), continue; end
    base_val = mean(x(mask_base),'omitnan');
    traces_bc(:,i) = x - base_val;
end
end

