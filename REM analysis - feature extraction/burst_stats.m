function [pct, n_bursts] = burst_stats(bin, fs)
bin = bin(:);
pct = 100 * sum(bin)/numel(bin);
% count rising edges
d = diff([false; bin; false]);
starts = find(d==1);  %#ok<NASGU>
ends   = find(d==-1)-1; %#ok<NASGU>
n_bursts = sum(d==1);
end