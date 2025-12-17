function EP = compute_epoch_spectra(eeg, fs, states, epochs_t, varargin)
% compute_epoch_spectra
% -------------------------------------------------------------
% Compute band powers per epoch and return a table EP.
%
% Required:
%   eeg      : [nSamples x 1] notched EEG
%   fs       : sampling rate (Hz)
%   states   : [nEpochs x 1] state labels for each epoch
%   epochs_t : [nEpochs x 1] epoch start times (s)
%
% Optional name-value:
%   'mouse','geno','cond' : metadata added as columns

ip = inputParser;
addParameter(ip,'mouse',"",@(x)ischar(x)||isstring(x));
addParameter(ip,'geno',"",@(x)ischar(x)||isstring(x));
addParameter(ip,'cond',"",@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
mouse = string(ip.Results.mouse);
geno  = string(ip.Results.geno);
cond  = string(ip.Results.cond);

nEpochs = numel(states);

% --- infer epoch length in samples (assumes constant) --------------------
if nEpochs > 1
    epoch_len_s = epochs_t(2) - epochs_t(1);
else
    % fallback: guess from length of signal
    epoch_len_s = numel(eeg)/fs;
end
epoch_len_samples = round(epoch_len_s * fs);

% frequency bands (adjust to your exact definitions)
bands = struct( ...
    'delta',[0.5 4], ...
    'theta',[6 10], ...
    'sigma',[11 16], ...
    'beta', [15 30]);

% preallocate
delta_pow = nan(nEpochs,1);
theta_pow = nan(nEpochs,1);
sigma_pow = nan(nEpochs,1);
beta_pow  = nan(nEpochs,1);

% --- loop epochs --------------------------------------------------------
for iE = 1:nEpochs
    idx_start = (iE-1)*epoch_len_samples + 1;
    idx_end   = min(iE*epoch_len_samples, numel(eeg));
    x = eeg(idx_start:idx_end);
    
    if numel(x) < fs   % too short? skip
        continue;
    end
    
    % Welch PSD
    window    = hamming(round(fs));     % 1-second window
    noverlap  = round(fs * 0.5);        % 50% overlap
    nfft      = [];                     % let MATLAB choose
    [ Pxx, f ] = pwelch(x, window, noverlap, nfft, fs);


    % Tell bandpower that Pxx,f are a PSD:
    delta_pow(iE) = bandpower(Pxx, f, bands.delta, 'psd');
    theta_pow(iE) = bandpower(Pxx, f, bands.theta, 'psd');
    sigma_pow(iE) = bandpower(Pxx, f, bands.sigma, 'psd');
    beta_pow(iE)  = bandpower(Pxx, f, bands.beta, 'psd');

end

% --- build table --------------------------------------------------------
EP = table;
EP.mouse    = repmat(mouse, nEpochs, 1);
EP.geno     = repmat(geno,  nEpochs, 1);
EP.cond     = repmat(cond,  nEpochs, 1);
EP.t_start  = epochs_t(:);
EP.state    = string(states(:));

EP.delta_pow = delta_pow;
EP.theta_pow = theta_pow;
EP.sigma_pow = sigma_pow;
EP.beta_pow  = beta_pow;
