function OUT = psd_sigma_theta_analysis(eeg_mat_file, scores_csv, varargin)
% psd_sigma_theta_analysis
% -------------------------------------------------------------
% For each mouse/session, this function:
%   1) Takes NREM EEG, filters it in the sigma band (e.g. 10–15 Hz),
%      computes the sigma power ENVELOPE, and then the PSD of that
%      envelope in the slow range (0–0.15 Hz).
%      -> "How sigma power fluctuates over time" during NREM.
%
%   2) Takes REM EEG, filters it in the theta band (e.g. 5–9 Hz),
%      computes the theta power ENVELOPE, and then the PSD of that
%      envelope in 0–0.15 Hz.
%      -> "How theta power fluctuates over time" during REM.
%
%   3) Extracts three summary metrics from those modulation PSDs:
%         - Total modulation power (area under PSD curve)
%         - Peak modulation frequency (Hz)
%         - Peak modulation amplitude (height at that peak)
%
%   4) Plots:
%         Left column:   normalized modulation PSDs (sigma NREM,
%                        theta REM) vs 0–0.15 Hz (A.U.).
%         Right column:  three separate bar plots:
%                        (a) sigma vs theta modulation power
%                        (b) sigma vs theta peak modulation frequency
%                        (c) sigma vs theta peak modulation amplitude
%
%   5) Saves a CSV with the six metrics and a PNG of the figure.
%
% This is conceptually similar to Fig. 3.9 in the paper you showed,
% where the x-axis 0–0.15 Hz reflects the FREQUENCY OF POWER
% MODULATIONS (not the carrier 10–15 Hz sigma itself).
%
% -------------------------------------------------------------
% Example call:
%
% OUT = psd_sigma_theta_analysis('mouse1_baseline.mat', ...
%       'mouse1_base_scores_1Hz.csv', ...
%       'codes', struct('WK',1,'NREM',4,'REM',9,'MA',15), ...
%       'ma_thresh_sec', 15, ...
%       'reclassify_short_wake_to_MA', true, ...
%       'mouse_id', 'M1', ...
%       'session', 'baseline', ...
%       'out_prefix', 'M1_base', ...
%       'out_dir', pwd, ...
%       'verbose', true);
%
% -------------------------------------------------------------

%% --------------- Parse inputs (similar style to other script) ---------
p = inputParser;
addRequired(p,'eeg_mat_file',@ischar);
addRequired(p,'scores_csv',@ischar);

% Sleep codes
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);

% 1 Hz scoring epoch length
addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);

% Micro-arousal (MA) handling: reclassify short Wake to MA
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);
addParameter(p,'reclassify_short_wake_to_MA',true,@islogical);

% Carrier bands for the envelope (Hz)
bands.sigma = [10 15];   % NREM sigma carrier
bands.theta = [5  9];    % REM theta carrier
addParameter(p,'bands',bands,@isstruct);

% Modulation frequency range (Hz) for the envelope PSD
addParameter(p,'mod_f_max',0.15,@(x)isscalar(x)&&x>0);

% Envelope downsampling target (Hz): we downsample the envelope to
% something like 20–50 Hz because its modulations are very slow.
addParameter(p,'env_target_fs',20,@(x)isscalar(x)&&x>0);

% Window length (sec) for Welch on the envelope
addParameter(p,'env_win_sec',300,@(x)isscalar(x)&&x>0);  % e.g. 5 min

% Book-keeping
addParameter(p,'mouse_id','',@ischar);
addParameter(p,'session','',@ischar);
addParameter(p,'out_prefix','',@ischar);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'verbose',true,@islogical);

parse(p, eeg_mat_file, scores_csv, varargin{:});
S = p.Results;
C = S.codes;
verbose = S.verbose;

if isempty(S.out_dir)
    S.out_dir = fileparts(scores_csv);
    if isempty(S.out_dir), S.out_dir = pwd; end
end
if ~exist(S.out_dir,'dir'); mkdir(S.out_dir); end

if S.reclassify_short_wake_to_MA && ~isfield(C,'MA')
    warning('reclassify_short_wake_to_MA=true but codes.MA missing; disabling.');
    S.reclassify_short_wake_to_MA = false;
end

OUT = struct();
OUT.success = false;
OUT.files   = struct();

if verbose
    fprintf('\n=== PSD sigma/theta modulation analysis ===\n');
    fprintf('EEG file   : %s\n', eeg_mat_file);
    fprintf('Scores CSV : %s\n', scores_csv);
    fprintf('Output dir : %s\n', S.out_dir);
end

%% --------------- Load 1 Hz scores and reclassify short Wake ---------
T = readtable(scores_csv);
assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
    'Expected columns time_s and score in %s', scores_csv);

t    = double(T.time_s(:));   % seconds (1 Hz grid)
code = double(T.score(:));    % integer state codes

n_epoch      = numel(code);
dur_scores_s = n_epoch * S.epoch_sec;

if verbose
    fprintf('Score vector: %d samples (%.1f min)\n', n_epoch, dur_scores_s/60);
    fprintf('Unique codes in scores: %s\n', mat2str(unique(code)'));
end

% ---- MA reclassification: short Wake bouts -> MA ---------------------
% We operate directly on the 1-Hz code vector.
if S.reclassify_short_wake_to_MA && isfield(C,'MA')
    if verbose
        fprintf('Reclassifying Wake runs <= %d s to MA...\n', S.ma_thresh_sec);
    end
    [st0,en0] = runs_from_codes(code);  % contiguous runs on the 1-Hz code vector
    rc   = code(st0);                   % run code
    rdur = t(en0) - t(st0) + 1;         % run duration (s), inclusive
    is_short_wake = (rc == C.WK) & (rdur <= S.ma_thresh_sec);
    for k = find(is_short_wake).'
        code(st0(k):en0(k)) = C.MA;
    end
end

%% --------------- Load EEG + sampling rate ---------------------------
EEG_struct = load(eeg_mat_file);

% Try common EEG field names
cand_eeg = {'eeg','EEG','EEG_rawtrace','EEG_raw'};
eeg = [];
eeg_name = '';
for k = 1:numel(cand_eeg)
    if isfield(EEG_struct, cand_eeg{k})
        eeg = EEG_struct.(cand_eeg{k});
        eeg_name = cand_eeg{k};
        break;
    end
end
if isempty(eeg)
    error('Could not find EEG variable in %s', eeg_mat_file);
end
eeg = double(eeg(:));

% Try common sampling-rate field names
cand_fs = {'eeg_frequency','fs_eeg','fs','Fs'};
fs = [];
fs_name = '';
for k = 1:numel(cand_fs)
    if isfield(EEG_struct, cand_fs{k})
        fs = EEG_struct.(cand_fs{k});
        fs_name = cand_fs{k};
        break;
    end
end
if isempty(fs)
    error('Could not find sampling rate in %s', eeg_mat_file);
end
fs = double(fs);

if verbose
    fprintf('Using EEG variable: %s (n = %d samples)\n', eeg_name, numel(eeg));
    fprintf('Using sampling rate: %s = %.3f Hz\n', fs_name, fs);
end

%% --------------- Align continuous EEG to 1 Hz labels -----------------
% Assume EEG starts at same time as scoring.
n_eeg     = numel(eeg);
dur_eeg_s = n_eeg / fs;
t_eeg     = (0:n_eeg-1)'/fs;                 % seconds from start
sec_idx   = floor(t_eeg / S.epoch_sec) + 1;  % which 1-Hz epoch each sample belongs to

valid  = sec_idx >= 1 & sec_idx <= n_epoch;
eeg    = eeg(valid);
sec_idx = sec_idx(valid);
state_vec = code(sec_idx);

if verbose
    fprintf('EEG duration total : %.1f min\n', dur_eeg_s/60);
    fprintf('Overlap with scores: %.1f min (used)\n', max(t_eeg(valid))/60);
end

%% --------------- Extract NREM and REM EEG segments -------------------
mask_nrem = (state_vec == C.NREM);
mask_rem  = (state_vec == C.REM);

eeg_nrem = eeg(mask_nrem);
eeg_rem  = eeg(mask_rem);

L_nrem = numel(eeg_nrem);
L_rem  = numel(eeg_rem);

if verbose
    fprintf('NREM samples: %d (%.1f s, %.2f min)\n', L_nrem, L_nrem/fs, L_nrem/fs/60);
    fprintf('REM  samples: %d (%.1f s, %.2f min)\n', L_rem,  L_rem/fs,  L_rem/fs/60);
end

if L_nrem == 0
    error('No NREM samples found (code %d).', C.NREM);
end
if L_rem == 0
    error('No REM samples found (code %d).', C.REM);
end
if verbose
    if L_nrem/fs < 60
        fprintf('[NOTE] <60 s of NREM: sigma modulation PSD may be noisy.\n');
    end
    if L_rem/fs < 60
        fprintf('[NOTE] <60 s of REM: theta modulation PSD may be noisy.\n');
    end
end

%% --------------- Sigma envelope modulation PSD (NREM) ----------------
[sig_f_mod, sig_psd_mod, sig_power, sig_peak_f, sig_peak_amp] = ...
    modulation_psd(eeg_nrem, fs, S.bands.sigma, S.mod_f_max, ...
                   S.env_target_fs, S.env_win_sec, verbose, 'NREM sigma');

%% --------------- Theta envelope modulation PSD (REM) -----------------
[th_f_mod, th_psd_mod, th_power, th_peak_f, th_peak_amp] = ...
    modulation_psd(eeg_rem, fs, S.bands.theta, S.mod_f_max, ...
                   S.env_target_fs, S.env_win_sec, verbose, 'REM theta');

%% --------------- Pack metrics into OUT struct ------------------------
OUT.sigma.mod_f       = sig_f_mod;
OUT.sigma.mod_psd     = sig_psd_mod;       % raw PSD of sigma envelope
OUT.sigma.mod_power   = sig_power;         % area under curve
OUT.sigma.mod_peak_f  = sig_peak_f;        % peak modulation frequency
OUT.sigma.mod_peak_amp= sig_peak_amp;      % peak amplitude

OUT.theta.mod_f       = th_f_mod;
OUT.theta.mod_psd     = th_psd_mod;
OUT.theta.mod_power   = th_power;
OUT.theta.mod_peak_f  = th_peak_f;
OUT.theta.mod_peak_amp= th_peak_amp;

OUT.params = S;
OUT.mouse_id = S.mouse_id;
OUT.session  = S.session;
OUT.eeg_file = eeg_mat_file;
OUT.scores_file = scores_csv;

%% --------------- Save metrics table for group stats ------------------
if isempty(S.out_prefix)
    base = sprintf('PSDmod_%s_%s', safe_str(S.mouse_id), safe_str(S.session));
else
    base = sprintf('%s_PSDmod', safe_str(S.out_prefix));
end
csv_name = fullfile(S.out_dir, [base '.csv']);

metrics_tbl = table( ...
    {S.mouse_id}, {S.session}, ...
    sig_power, sig_peak_f, sig_peak_amp, ...
    th_power,  th_peak_f,  th_peak_amp, ...
    'VariableNames', {'mouse','session', ...
                      'sigma_mod_power','sigma_mod_peak_f','sigma_mod_peak_amp', ...
                      'theta_mod_power','theta_mod_peak_f','theta_mod_peak_amp'});
writetable(metrics_tbl, csv_name);
OUT.files.metrics_csv = csv_name;

%% --------------- Plot modulation PSDs + 3 summary bar plots ----------
fig = figure('Name','Sigma/Theta modulation PSD', ...
             'Color','w','Visible','on');

tstr = sprintf('%s – %s', S.mouse_id, S.session);
sgtitle(tstr);

% --- Normalize PSDs for plotting (A.U.), but keep raw values for metrics
sig_psd_norm = sig_psd_mod / max(sig_psd_mod);
th_psd_norm  = th_psd_mod  / max(th_psd_mod);

% Layout: 2 rows x 3 columns
% (1) NREM sigma modulation PSD (0–0.15 Hz)
subplot(2,3,1);
plot(sig_f_mod, sig_psd_norm, 'LineWidth',1.5);
xlabel('Frequency (Hz)');                 % modulation frequency
ylabel('Sigma power (A.U.)');
xlim([0 S.mod_f_max]);
ylim([0 1.05]);
title('NREM – sigma power modulation');
box off;

% (4) REM theta modulation PSD (0–0.15 Hz)
subplot(2,3,4);
plot(th_f_mod, th_psd_norm, 'LineWidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Theta power (A.U.)');
xlim([0 S.mod_f_max]);
ylim([0 1.05]);
title('REM – theta power modulation');
box off;

% --- Three separate bar plots comparing sigma vs theta ---

% (2) Modulation power
subplot(2,3,2);
bar_data_power = [sig_power; th_power];
bar(bar_data_power);
set(gca,'XTick',1:2,'XTickLabel',{'Sigma (NREM)','Theta (REM)'});
ylabel('Modulation power (raw units)');
title('Total modulation power');
xtickangle(20);
box off;

% (3) Peak modulation frequency
subplot(2,3,3);
bar_data_freq = [sig_peak_f; th_peak_f];
bar(bar_data_freq);
set(gca,'XTick',1:2,'XTickLabel',{'Sigma (NREM)','Theta (REM)'});
ylabel('Peak modulation freq (Hz)');
title('Peak modulation frequency');
xtickangle(20);
box off;

% (5) Peak modulation amplitude
subplot(2,3,5);
bar_data_amp = [sig_peak_amp; th_peak_amp];
bar(bar_data_amp);
set(gca,'XTick',1:2,'XTickLabel',{'Sigma (NREM)','Theta (REM)'});
ylabel('Peak modulation amplitude');
title('Peak modulation amplitude');
xtickangle(20);
box off;

% Leave subplot(2,3,6) empty or use for notes/legend if you want:
subplot(2,3,6);
axis off;
text(0,0.5, ...
    sprintf('Carrier bands:\nSigma: %.1f–%.1f Hz\nTheta: %.1f–%.1f Hz', ...
    S.bands.sigma(1), S.bands.sigma(2), ...
    S.bands.theta(1), S.bands.theta(2)), ...
    'FontSize',9);

% Save figure
fig_name = fullfile(S.out_dir, [base '_figure.png']);
saveas(fig, fig_name);
OUT.files.fig_png = fig_name;

if verbose
    fprintf('\nFiles saved:\n');
    fprintf('  Metrics CSV : %s\n', csv_name);
    fprintf('  Figure PNG  : %s\n', fig_name);
    fprintf('Analysis complete.\n\n');
end

OUT.success = true;
end  % main function


%% ---------------------- Helper: safe_str for filenames ---------------
function s = safe_str(s)
if isempty(s), s = 'NA'; return; end
s = regexprep(s,'\s+','_');
s = regexprep(s,'[^\w-]','');
end

%% ---------------------- Helper: contiguous runs from codes ----------
function [st,en] = runs_from_codes(code)
% Given a vector of integer codes, returns start/end indices of each run
code = code(:);
if isempty(code)
    st = []; en = [];
    return;
end
d = diff(code);
change_idx = [1; find(d ~= 0) + 1];
st = change_idx;
en = [change_idx(2:end)-1; numel(code)];
end

%% ---------------------- Helper: modulation PSD of envelope ----------
function [f_mod, psd_mod, total_power, peak_f, peak_amp] = ...
    modulation_psd(eeg_state, fs, carrier_band, mod_f_max, ...
                   env_target_fs, env_win_sec, verbose, label)
% 1) Band-pass filter EEG in the carrier band (sigma/theta)
Wn = carrier_band / (fs/2);
[b,a] = butter(4, Wn);                        % 4th-order Butterworth
eeg_bp = filtfilt(b,a,eeg_state);

% 2) Get the instantaneous amplitude (envelope) => "power" over time
env = abs(hilbert(eeg_bp));

% 3) Downsample the envelope to a manageable sampling rate, e.g. ~20 Hz
decim = max(1, floor(fs / env_target_fs));
env_ds = decimate(env, decim);
fs_env = fs / decim;

% Remove DC / slow offset so we don't get a massive spike at ~0 Hz
env_ds = detrend(env_ds, 'constant');   % subtract mean


if verbose
    fprintf('%s: env_fs = %.2f Hz, length = %.1f s\n', ...
        label, fs_env, numel(env_ds)/fs_env);
end

% 4) Compute PSD of the envelope (slow modulations) using Welch
L = numel(env_ds);
win_len = min(round(env_win_sec*fs_env), L);
if win_len < 10
    % if very little data, just use the whole segment as one window
    win_len = L;
end
win     = hamming(win_len);
nover   = floor(win_len/2);
nfft    = max(512, 2^nextpow2(win_len));

[psd_all, f_all] = pwelch(env_ds, win, nover, nfft, fs_env);

% 5) Keep only modulation frequencies between mod_f_min and mod_f_max
mod_f_min = 0.01;  % Hz, i.e. ignore fluctuations slower than 1/0.01 = 100 s
idx = (f_all >= mod_f_min) & (f_all <= mod_f_max);

f_mod   = f_all(idx);
psd_mod = psd_all(idx);


if isempty(f_mod)
    error('%s: no frequency bins <= %.3f Hz in envelope PSD.', label, mod_f_max);
end

% 6) Summary metrics on this modulation PSD
total_power = trapz(f_mod, psd_mod);     % area under curve
[peak_amp, k] = max(psd_mod);
peak_f = f_mod(k);

if verbose
    fprintf('%s: total_power=%.4g, peak_f=%.4f Hz, peak_amp=%.4g\n', ...
        label, total_power, peak_f, peak_amp);
end
end
