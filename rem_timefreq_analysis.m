function rem_timefreq_analysis(varargin)
% rem_timefreq_analysis: high-res time-frequency spectra for REM bouts
%
% INPUTS (name/value, all optional):
%   'file'        : path to *_REMselection.mat or *_REMselection_FULL.mat
%   'method'      : 'spectrogram' (default) or 'cwt'
%   'fRange'      : [1 100] (Hz) default
%   'norm'        : 'db' (default), 'zrem', or 'none'
%   'thetaBand'   : [5 9] Hz (default)
%   'bands'       : struct with fields: delta=[1 4], theta=[5 9], sigma=[10 15],
%                   beta=[15 30], lgamma=[30 60], hgamma=[60 100]
%   'outdir'      : output folder (default: alongside file in _REM_TFR/)
%   'saveMatrices': true/false (save TFR matrices per bout as .mat)
%
% OUTPUTS:
%   - Per-bout PNG spectrograms + CSV with band means and peak features.

% ---------- Parse inputs ----------
p = inputParser;
addParameter(p,'file','',@(s)ischar(s)||isstring(s));
addParameter(p,'method','spectrogram',@(s)any(strcmpi(s,{'spectrogram','cwt'})));
addParameter(p,'fRange',[1 100],@(v)isnumeric(v)&&numel(v)==2);
addParameter(p,'norm','db',@(s)any(strcmpi(s,{'db','zrem','none'})));
addParameter(p,'thetaBand',[5 9],@(v)isnumeric(v)&&numel(v)==2);
defaultBands = struct('delta',[1 4], 'theta',[5 9], 'sigma',[10 15], ...
                      'beta',[15 30], 'lgamma',[30 60], 'hgamma',[60 100]);
addParameter(p,'bands',defaultBands,@isstruct);
addParameter(p,'outdir','',@(s)ischar(s)||isstring(s));
addParameter(p,'saveMatrices',false,@islogical);
parse(p,varargin{:});

FILE   = string(p.Results.file);
METHOD = lower(p.Results.method);
FR     = p.Results.fRange;
NORM   = lower(p.Results.norm);
TH     = p.Results.thetaBand;
BANDS  = p.Results.bands;
OUTDIR = string(p.Results.outdir);
SAVE_M = p.Results.saveMatrices; %#ok<NASGU>

if FILE=="" || ~isfile(FILE)
    [f,fp]=uigetfile('*REMselection*.mat','Select REM selection .mat');
    assert(f~=0,'No file selected.'); FILE=string(fullfile(fp,f));
end

S = load(FILE);

% ---------- Locate required fields ----------
fs_eeg = tryget(S,{'eeg_frequency','fs_eeg','Fs','fs'});
assert(~isempty(fs_eeg),'Could not find sampling rate (eeg_frequency/fs_eeg).');
fs_eeg = double(fs_eeg);

% EEG at top level is optional; we'll fall back to bout segments if missing
EEG_name = fieldnames(S);
EEG_name = EEG_name(ismember(lower(EEG_name), {'eeg', 'EEG'}));
EEG_name = EEG_name{1}; 
EEG = []; % Initialize EEG to avoid potential errors if no field is found
if ~isempty(EEG_name)
    EEG = S.(EEG_name);
end
bouts  = S.(findfield(S,{'bouts'}));
remtbl = S.(findfield(S,{'rem_bouts_table'}));

[folder,base] = fileparts(FILE);
if OUTDIR=="", OUTDIR = fullfile(folder, base + "_REM_TFR"); end
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% ---------- CWT fallback if toolbox missing ----------
if strcmpi(METHOD,'cwt') && ~(exist('cwt','file')==2)
    warning('Wavelet Toolbox not found. Falling back to ''spectrogram''.');
    METHOD = 'spectrogram';
end

% ---------- Baseline for 'zrem' normalization ----------
if strcmpi(NORM,'zrem')
    if ~isempty(EEG)
        allIdx = vertcat(bouts.idx_samples);
        baselineSig = double(EEG(allIdx));
    else
        baselineSig = cell2mat(arrayfun(@(b) double(b.EEG(:)), bouts, 'UniformOutput', false));
    end
end

% ---------- Peak header suffix ----------
if strcmpi(NORM,'db')
    pow_suffix = '_db';
elseif strcmpi(NORM,'zrem')
    pow_suffix = '_z';
else
    pow_suffix = '_raw';
end

% ---------- Feature prealloc ----------
nB = numel(bouts);
bands_list = fieldnames(BANDS)';           
nBands     = numel(bands_list);
nPeakCols  = 7;
nNewCols   = 4;                 % EMG + PAC (if present)
nThetaCols = 2;                 % theta_peak_sd, theta_bw_Hz
nFragCols  = 2;                 % ibi_prev_s, ibi_next_s
feat_cols  = 5 + nBands + nPeakCols + nNewCols + nThetaCols + nFragCols;
feat_rows  = cell(nB, feat_cols);

% --- Inter-REM intervals (fragmentation proxies) ---
% Use the bout table to compute the time from previous REM end to this REM start, and to next REM start.
start_times = remtbl.start_time_s(:);
end_times   = remtbl.end_time_s(:);

ibi_prev_s = nan(nB,1);
ibi_next_s = nan(nB,1);
if nB >= 2
    ibi_prev_s(2:end) = start_times(2:end) - end_times(1:end-1);       % gap from previous REM end
    ibi_next_s(1:end-1) = start_times(2:end) - end_times(1:end-1);     % same gap but assigned to current row for "next"
end


% ---------- Loop bouts ----------
for k = 1:nB
    sig = double(bouts(k).EEG(:));
    % t   = bouts(k).t_s(:);  % not needed for calc; plotting uses T returned by spectrogram/cwt

    % --- TFR ---
    switch METHOD
        case 'spectrogram'
            winSec = 0.5;
            wlen   = max(64, round(winSec*fs_eeg));
            nover  = round(0.9*wlen);
            nfft   = 2^nextpow2(max(wlen, 4*wlen));
            [Sxx,F,T] = spectrogram(sig, wlen, nover, nfft, fs_eeg, 'yaxis');
            P = abs(Sxx).^2;
        case 'cwt'
            [cfs,F,T] = cwt(sig, fs_eeg, 'FrequencyLimits', FR);
            P = abs(cfs).^2;
    end

    % Restrict frequency range
    fmask = F>=FR(1) & F<=FR(2);
    F = F(fmask); P = P(fmask,:);

    % --- Normalization ---
    switch NORM
        case 'db'
            ref = median(P,2); ref(ref==0)=eps;
            Z = 10*log10(P./ref);
        case 'zrem'
            [Pb,Fb] = quick_ps(baselineSig, fs_eeg, METHOD, FR);
            mu = mean(Pb,2); sd = std(Pb,0,2); sd(sd==0)=eps;
            mu = interp1(Fb,mu,F,'linear','extrap');
            sd = interp1(Fb,sd,F,'linear','extrap');
            Z  = (P - mu)./sd;
        case 'none'
            Z = P;
    end

    % Pick matrix for "power" outputs (raw when none, else normalized)
    if strcmpi(NORM,'none')
        PforPow = P;
    else
        PforPow = Z;
    end

    % --- Theta peak (frequency via raw power average over time) ---
    Pavg = mean(P,2);                   % raw spectrum avg over time
    thmask = F>=TH(1) & F<=TH(2);
    if any(thmask)
        [~,imax] = max(Pavg(thmask));
        theta_peak_Hz = F(find(thmask,1,'first') + imax - 1);
    else
        theta_peak_Hz = NaN;
    end

    % ---------- Theta peak stability (SD over time) ----------
thmask = F>=TH(1) & F<=TH(2);
Fth = F(thmask);
if any(thmask)
    % frame-wise theta peak: find max in theta band for each time bin (use raw P)
    [~, ix_time] = max(P(thmask,:), [], 1);      % 1 x nFrames (indices within Fth)
    theta_peak_traj = Fth(ix_time);              % Hz per frame
    theta_peak_sd = std(theta_peak_traj, 0, 'omitnan');
else
    theta_peak_sd = NaN;
end

% ---------- Theta -3 dB bandwidth around bout theta peak ----------
% Use time-averaged RAW spectrum Pavg to find bandwidth around the bout's theta_peak_Hz
if ~isnan(theta_peak_Hz)
    % index of closest freq to the bout theta peak
    [~, i0] = min(abs(F - theta_peak_Hz));
    S0 = Pavg(i0);                % reference power at the peak
    if isfinite(S0) && S0>0
        halfPow = S0/2;           % -3 dB point in power units
        % search within theta band on each side
        inTheta = (F>=TH(1) & F<=TH(2));
        % left edge
        i_left = i0;
        while i_left>1 && inTheta(i_left-1) && Pavg(i_left-1) >= halfPow
            i_left = i_left-1;
        end
        % right edge
        i_right = i0;
        while i_right<numel(F) && inTheta(i_right+1) && Pavg(i_right+1) >= halfPow
            i_right = i_right+1;
        end
        theta_bw_Hz = F(i_right) - F(i_left);
    else
        theta_bw_Hz = NaN;
    end
else
    theta_bw_Hz = NaN;
end


    % --- Band averages (match normalization space in outputs) ---
    bandvals = nan(nBands,1);
    for b = 1:nBands
        fr = BANDS.(bands_list{b});
        bm = F>=fr(1) & F<=fr(2);
        if any(bm)
            if strcmpi(NORM,'none')
                bandvals(b) = mean(P(bm,:),'all','omitnan');
            else
                bandvals(b) = mean(Z(bm,:),'all','omitnan');
            end
        end
    end

    % --- Peak features ---
    % Global peak across full FR
    [peak_global_freq_Hz, peak_global_power] = local_band_peak(F, Pavg, PforPow, FR);

    % Theta peak power (at theta_peak_Hz)
    if ~isnan(theta_peak_Hz)
        [~, idx_theta] = min(abs(F - theta_peak_Hz));
        theta_peak_power = mean(PforPow(idx_theta,:), 'omitnan');
    else
        theta_peak_power = NaN;
    end

    % Low- and high-gamma peaks
    [lgamma_peak_freq_Hz, lgamma_peak_power] = local_band_peak(F, Pavg, PforPow, [30 60]);
    [hgamma_peak_freq_Hz, hgamma_peak_power] = local_band_peak(F, Pavg, PforPow, [60 100]);
    
    % ---------- EMG burst metrics (robust, per-bout) ----------
    % If EMG present in bouts, compute; else NaN
    if isfield(bouts(k),'EMG') && ~isempty(bouts(k).EMG)
        emg = double(bouts(k).EMG(:));
        [emg_burst_rate_per_min, emg_burst_pct] = emg_burst_metrics(emg, fs_eeg);
    else
        emg_burst_rate_per_min = NaN;
        emg_burst_pct = NaN;
    end

    % ---------- θ–γ PAC (Tort MI; θ phase 5–9 Hz vs γ amplitude) ----------
    % Use your TH for theta; low/high gamma fixed at [30 60] and [60 100]
    pac_theta_lgamma_MI = pac_tort_mi(double(bouts(k).EEG(:)), fs_eeg, TH, [30 60]);
    pac_theta_hgamma_MI = pac_tort_mi(double(bouts(k).EEG(:)), fs_eeg, TH, [60 100]);

    % --- Store features row ---
    band_cells = num2cell(bandvals(:).');   % 1×nBands
    peak_cells = { ...
        peak_global_freq_Hz, peak_global_power, ...
        theta_peak_power, ...
        lgamma_peak_freq_Hz, lgamma_peak_power, ...
        hgamma_peak_freq_Hz, hgamma_peak_power ...
    };

    theta_cells = {theta_peak_sd, theta_bw_Hz};
    frag_cells  = {ibi_prev_s(k), ibi_next_s(k)};   % fragmentation proxies

    feat_rows(k, :) = [ { ...
        k, remtbl.start_time_s(k), remtbl.end_time_s(k), remtbl.dur_s(k), theta_peak_Hz ...
    }, band_cells, peak_cells, { emg_burst_rate_per_min, emg_burst_pct, pac_theta_lgamma_MI, pac_theta_hgamma_MI }, ...
    theta_cells, frag_cells ];


    % --- Plot ---
    hf = figure('Color','w','Position',[100 100 900 420]);
    imagesc(T,F,Z); axis xy;
    ylim([1 30]);                      % focus on delta–beta
    if strcmpi(NORM,'none')
        lo = prctile(Z(:), 5);         % robust autoscale
        hi = prctile(Z(:), 95);
    clim([lo hi]);
    end
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(sprintf('REM Bout %d  |  %s  |  %s', k, upper(METHOD), upper(NORM)));
    cb = colorbar;
    if strcmpi(NORM,'db'), cb_label='Power (dB)';
    elseif strcmpi(NORM,'zrem'), cb_label='Z-score';
    else, cb_label='Power'; end
    cb.Label.String = cb_label;
    hold on;
    if ~isnan(theta_peak_Hz)
        plot([T(1) T(end)],[theta_peak_Hz theta_peak_Hz],'--','LineWidth',1);
    end
    hold off;

    png_out = fullfile(OUTDIR, sprintf('%s_REMbout_%03d_%s_%s.png', base, k, METHOD, NORM));
    saveas(hf, png_out); close(hf);
end

% ---------- Save feature table ----------
vars = {'bout','start_time_s','end_time_s','dur_s','theta_peak_Hz'};
band_headers = strcat(bands_list, ['_' NORM]);
peak_headers = { ...
    'peak_global_freq_Hz', ['peak_global_power' pow_suffix], ...
    ['theta_peak_power' pow_suffix], ...
    'lgamma_peak_freq_Hz', ['lgamma_peak_power' pow_suffix], ...
    'hgamma_peak_freq_Hz', ['hgamma_peak_power' pow_suffix] ...
};

% If you have EMG/PAC:
extra_headers1 = { 'emg_burst_rate_per_min','emg_burst_pct','pac_theta_lgamma_MI','pac_theta_hgamma_MI' };

% New theta/fragmentation headers:
extra_headers2 = { 'theta_peak_sd','theta_bw_Hz','ibi_prev_s','ibi_next_s' };

vars = [vars, band_headers, peak_headers, extra_headers1, extra_headers2];

assert( size(feat_rows,2) == numel(vars), ...
    'Header/feature width mismatch: %d vs %d', size(feat_rows,2), numel(vars));


Tfeat = cell2table(feat_rows, 'VariableNames', vars);

if OUTDIR=="", OUTDIR = folder; end
writetable(Tfeat, fullfile(OUTDIR, sprintf('%s_REM_TFR_features.csv', base)));

fprintf('[OK] Saved TFR figures and features to: %s\n', OUTDIR);
end

