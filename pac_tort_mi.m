function mi = pac_tort_mi(x, fs, f_phase, f_amp)
% Tort MI (KL divergence of amp-by-phase histogram), toolbox-free
% f_phase: [low high] (e.g., [5 9]); f_amp: [30 60] or [60 100]

x = double(x(:));

% Phase (theta)
xp  = bandpass_zero(x, fs, f_phase);
phi = angle(hilbert(xp));

% Amplitude (gamma)
xa  = bandpass_zero(x, fs, f_amp);
amp = abs(hilbert(xa));

% Bin amplitude by phase
nbins = 18;
edges = linspace(-pi, pi, nbins+1);
[~,~,bin] = histcounts(phi, edges);

A = zeros(1, nbins);
for b = 1:nbins
    ab = amp(bin==b);
    if ~isempty(ab)
        A(b) = mean(ab);
    end
end

if ~any(A)                               % all zeros -> undefined MI
    mi = NaN; 
    return;
end

P = A / sum(A);                          % probability over phase bins

% ---- toolbox-free entropy (ignore NaNs/zeros) ----
mask = isfinite(P) & (P > 0);
H = -sum(P(mask) .* log(P(mask)));       % natural log
Hmax = log(nbins);

mi = (Hmax - H) / Hmax;                  % Tort's normalized MI in [0,1]
end
