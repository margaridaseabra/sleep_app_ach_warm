function [pkHz, pkPow] = local_band_peak(F, Pavg, ZorP, fband)
mask = F>=fband(1) & F<=fband(2);
if ~any(mask)
    pkHz = NaN; pkPow = NaN; return;
end
[~, relIdx] = max(Pavg(mask));
idx0 = find(mask,1,'first') + relIdx - 1;
pkHz = F(idx0);
pkPow = mean(ZorP(idx0,:), 'omitnan');
end
