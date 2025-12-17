function OUT = plot_ach_eeg_segment(mat_file, scores_csv, t_start, t_end, varargin)
% PLOT_ACH_EEG_SEGMENT
% -------------------------------------------------------------
% Science-style overview for ONE recording and time window:
%   - EEG trace with state-coloured line
%   - EEG spectrogram (0–30 Hz)
%   - EMG trace with state-coloured line
%   - ACh trace with state-coloured line
%   - ACh power spectrogram (0–1 Hz)
%
% Also computes ACh onset features (canonical stable transitions):
%   - Wake onset: NREM/REM -> Wake
%   - NREM onset: Wake -> NREM
%   - REM  onset: NREM -> REM
%
% Per state (Wake/NREM/REM):
%   n               : # onsets in window
%   deltaF_all      : ΔF/F per onset
%   slope_all       : slope per onset (ΔF/F per s)
%   deltaF_mean     : mean ΔF/F
%   slope_mean      : mean slope
%   cycleHz         : ACh cycle frequency estimate in that state+window
%
% For each state a SEPARATE FIGURE with feature histograms is created.
%
% USAGE:
%   OUT = plot_ach_eeg_segment(mat_file, scores_csv, 300, 900);
%
% OPTIONAL NAME/VALUE:
%   'codes'   : struct with fields .WK .NREM .REM .MA (default 0,1,2,15)
%   't_pre'   : secs before onset for features (default 5)
%   't_post'  : secs after onset for features (default 5)

% -------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------
p = inputParser;
p.addParameter('codes', struct('WK',0,'NREM',1,'REM',2,'MA',15));
p.addParameter('t_pre',  5);
p.addParameter('t_post', 5);
p.parse(varargin{:});
opt = p.Results;
CODES = opt.codes;

if t_end <= t_start
    error('t_end must be greater than t_start.');
end
dur = t_end - t_start;

% -------------------------------------------------------------
% Load signals from MAT file
% -------------------------------------------------------------
info = whos('-file', mat_file);
names = {info.name};

pick = @(cands) cands{find(ismember(cands, names),1,'first')};

eeg_name    = pick({'eeg','EEG','eeg1','Eeg','eeg_filt'});
fs_eeg_name = pick({'eeg_frequency','fs_eeg','Fs_eeg','EEG_frequency'});

emg_name    = pick({'emg','EMG','emg1','EMG1','emg_filt'});
fs_emg_name = pick({'eeg_frequency','fs_emg','Fs_emg','EMG_frequency'});

ach_name    = pick({'ach','ACh','ne','dff','dFF','dff_ach'});
fs_ach_name = pick({'ach_frequency','ne_frequency','fs_ach','Fs_ach'});

S = load(mat_file);

if isempty(eeg_name) || isempty(fs_eeg_name)
    error('Could not find EEG + sampling rate in %s', mat_file);
end
if isempty(emg_name) || isempty(fs_emg_name)
    error('Could not find EMG + sampling rate in %s', mat_file);
end
if isempty(ach_name) || isempty(fs_ach_name)
    error('Could not find ACh + sampling rate in %s', mat_file);
end

eeg    = S.(eeg_name)(:);
fs_eeg = S.(fs_eeg_name);
emg    = S.(emg_name)(:);
fs_emg = S.(fs_emg_name);
ach    = S.(ach_name)(:);
fs_ach = S.(fs_ach_name);

% 50 Hz notch for EEG, skipping specific files
eeg = apply_50Hz_notch_if_needed(eeg, fs_eeg, mat_file);

% -------------------------------------------------------------
% Load 1-Hz scores
% -------------------------------------------------------------
M = readmatrix(scores_csv);
if size(M,2) == 1
    score     = M(:,1);
    epoch_sec = 1;
    t_scores  = (0:numel(score)-1)' * epoch_sec;
else
    t_scores  = M(:,1);
    score     = M(:,2);
    dt = diff(t_scores);
    epoch_sec = mode(dt(~isnan(dt)));
end
score    = score(:);
t_scores = t_scores(:);

% -------------------------------------------------------------
% Clip signals to [t_start, t_end] and build time axes (relative)
% -------------------------------------------------------------
% EEG
i1_eeg = max(1, floor(t_start*fs_eeg)+1);
i2_eeg = min(numel(eeg), floor(t_end*fs_eeg));
eeg_win = eeg(i1_eeg:i2_eeg);
t_eeg = ((i1_eeg:i2_eeg)-1)/fs_eeg - t_start;

% EMG
i1_emg = max(1, floor(t_start*fs_emg)+1);
i2_emg = min(numel(emg), floor(t_end*fs_emg));
emg_win = emg(i1_emg:i2_emg);
t_emg = ((i1_emg:i2_emg)-1)/fs_emg - t_start;

% ACh
i1_ach = max(1, floor(t_start*fs_ach)+1);
i2_ach = min(numel(ach), floor(t_end*fs_ach));
ach_win = ach(i1_ach:i2_ach);
t_ach = ((i1_ach:i2_ach)-1)/fs_ach - t_start;

% Scores window
idx_scores   = t_scores >= t_start & t_scores < t_end;
t_scores_win = t_scores(idx_scores) - t_start;
score_win    = score(idx_scores);

% State runs (for colouring)
if ~isempty(score_win)
    s = score_win;
    run_starts = [1; find(diff(s)~=0)+1];
    run_ends   = [run_starts(2:end)-1; numel(s)];
    run_codes  = s(run_starts);
else
    run_starts = [];
    run_ends   = [];
    run_codes  = [];
end

% -------------------------------------------------------------
% State colours
% -------------------------------------------------------------
col_wake = [0.8 0.3 0.8];
col_nrem = [0.97 0.6 0.6];
col_rem  = [0.6 0.9 0.6];
col_ma   = [1.0 0.95 0.6];

stateColors = containers.Map( ...
   {CODES.WK,      CODES.NREM,     CODES.REM,     CODES.MA}, ...
   {col_wake,      col_nrem,       col_rem,       col_ma});

% -------------------------------------------------------------
% FIGURE LAYOUT  (no separate "state-only" axis)
% -------------------------------------------------------------
figure('Color','w','Position',[50 50 1100 750]);

nRows  = 5;
ax_eeg   = subplot(nRows,1,1);
ax_specE = subplot(nRows,1,2);
ax_emg   = subplot(nRows,1,3);
ax_ach   = subplot(nRows,1,4);
ax_specA = subplot(nRows,1,5);

% -------------------------------------------------------------
% 2) EEG trace with coloured segments
% -------------------------------------------------------------
axes(ax_eeg); cla; hold on;
plot_state_coloured_trace(t_eeg, eeg_win, ...
    t_scores_win, epoch_sec, run_starts, run_ends, run_codes, stateColors);
xlim([0 dur]);
ylabel('EEG (a.u.)');
box on;

% -------------------------------------------------------------
% 3) EEG spectrogram
% -------------------------------------------------------------
axes(ax_specE); cla; hold on;
win_e = round(fs_eeg * 2);       % 2-s window
ov_e  = round(win_e * 0.75);
f_e   = 0:0.5:30;                % 0–30 Hz
if numel(eeg_win) > win_e
    [S,F,T] = spectrogram(eeg_win, win_e, ov_e, f_e, fs_eeg);
    P = 10*log10(abs(S).^2 + eps);
    clim = [prctile(P(:),5) prctile(P(:),95)];
    imagesc(T, F, P, clim);
    axis xy;
    ylabel('Freq (Hz)');
    xlim([0 dur]);
    colormap(ax_specE, 'turbo');
    colorbar;
else
    text(0.5,0.5,'EEG too short for spectrogram','Units','normalized');
    axis off;
end
box on;

% -------------------------------------------------------------
% 4) EMG trace with coloured segments
% -------------------------------------------------------------
axes(ax_emg); cla; hold on;
plot_state_coloured_trace(t_emg, emg_win, ...
    t_scores_win, epoch_sec, run_starts, run_ends, run_codes, stateColors);
xlim([0 dur]);
ylabel('EMG (a.u.)');
box on;

% -------------------------------------------------------------
% 5) ACh trace with coloured segments
% -------------------------------------------------------------
axes(ax_ach); cla; hold on;
plot_state_coloured_trace(t_ach, ach_win, ...
    t_scores_win, epoch_sec, run_starts, run_ends, run_codes, stateColors);
xlim([0 dur]);
ylabel('\DeltaF/F ACh');
box on;

% -------------------------------------------------------------
% 6) ACh power spectrogram (0–1 Hz)
% -------------------------------------------------------------
axes(ax_specA); cla; hold on;
win_a = max(round(fs_ach*5), 10);   % 5-s window
ov_a  = round(win_a*0.75);
f_a   = 0:0.05:1.0;                 % 0–1 Hz
if numel(ach_win) > win_a
    [S2,F2,T2] = spectrogram(ach_win, win_a, ov_a, f_a, fs_ach);
    P2 = 10*log10(abs(S2).^2 + eps);
    clim2 = [prctile(P2(:),5) prctile(P2(:),95)];
    imagesc(T2, F2, P2, clim2);
    axis xy;
    xlabel('Time (s, window relative)');
    ylabel('Freq (Hz)');
    title('ACh power spectrogram (0–1 Hz)');
    xlim([0 dur]);
    colormap(ax_specA, 'hot');
    colorbar;
else
    text(0.5,0.5,'ACh too short for spectrogram','Units','normalized');
    axis off;
end
box on;

% Link x-axes (no ax_state now)
linkaxes([ax_eeg ax_specE ax_emg ax_ach ax_specA],'x');
set(ax_specA,'XLim',[0 dur]);

% -------------------------------------------------------------
% 7) ACh onset features inside this window
% -------------------------------------------------------------
OUT = struct();
OUT.mat_file   = mat_file;
OUT.scores_csv = scores_csv;
OUT.t_start    = t_start;
OUT.t_end      = t_end;
OUT.codes      = CODES;
OUT.epoch_sec  = epoch_sec;

FEATURES = struct();

CODE_WAKE = CODES.WK;
CODE_NREM = CODES.NREM;
CODE_REM  = CODES.REM;

minPreBins  = round(10 / epoch_sec);   % 10 s stable prev
minPostBins = round(20 / epoch_sec);   % 20 s stable next

idx_W_fromN = find_stable_transitions(score, CODE_NREM, CODE_WAKE, ...
    minPreBins, minPostBins);
idx_W_fromR = find_stable_transitions(score, CODE_REM,  CODE_WAKE, ...
    minPreBins, minPostBins);
idx_Wake    = sort([idx_W_fromN; idx_W_fromR]);

idx_NREM = find_stable_transitions(score, CODE_WAKE, CODE_NREM, ...
    minPreBins, minPostBins);
idx_REM  = find_stable_transitions(score, CODE_NREM, CODE_REM, ...
    minPreBins, minPostBins);

t_Wake = (idx_Wake-1)*epoch_sec;
t_NREM = (idx_NREM-1)*epoch_sec;
t_REM  = (idx_REM -1)*epoch_sec;

t_Wake = t_Wake(t_Wake >= t_start & t_Wake <= t_end);
t_NREM = t_NREM(t_NREM >= t_start & t_NREM <= t_end);
t_REM  = t_REM(t_REM  >= t_start & t_REM  <= t_end);

FEATURES.Wake = compute_onset_features(ach, fs_ach, t_Wake,  opt.t_pre, opt.t_post);
FEATURES.NREM = compute_onset_features(ach, fs_ach, t_NREM,  opt.t_pre, opt.t_post);
FEATURES.REM  = compute_onset_features(ach, fs_ach, t_REM,   opt.t_pre, opt.t_post);

FEATURES.Wake.cycleHz = estimate_cycle_freq_state(ach, fs_ach, score, t_scores, CODES.WK,  t_start, t_end);
FEATURES.NREM.cycleHz = estimate_cycle_freq_state(ach, fs_ach, score, t_scores, CODES.NREM,t_start, t_end);
FEATURES.REM.cycleHz  = estimate_cycle_freq_state(ach, fs_ach, score, t_scores, CODES.REM, t_start, t_end);

OUT.features = FEATURES;

disp('--- ACh onset features in this window ---');
disp(FEATURES);

% short text box in main figure (anchor on EEG axis now)
ax_txt = ax_eeg;
axes(ax_txt);
txt = sprintf(['Wake: n=%d, \\DeltaF=%.2f, slope=%.3f, f=%.3f Hz\n' ...
               'NREM: n=%d, \\DeltaF=%.2f, slope=%.3f, f=%.3f Hz\n' ...
               'REM:  n=%d, \\DeltaF=%.2f, slope=%.3f, f=%.3f Hz'], ...
    FEATURES.Wake.n, FEATURES.Wake.deltaF_mean, FEATURES.Wake.slope_mean, FEATURES.Wake.cycleHz, ...
    FEATURES.NREM.n, FEATURES.NREM.deltaF_mean, FEATURES.NREM.slope_mean, FEATURES.NREM.cycleHz, ...
    FEATURES.REM.n,  FEATURES.REM.deltaF_mean,  FEATURES.REM.slope_mean,  FEATURES.REM.cycleHz);
annotation('textbox',[0.72 0.75 0.27 0.2], 'String', txt, ...
    'FitBoxToText','on','EdgeColor','none','FontSize',9);

% -------------------------------------------------------------
% 8) Separate feature figures per state
% -------------------------------------------------------------
plot_state_features(FEATURES.Wake,'Wake');
plot_state_features(FEATURES.NREM,'NREM');
plot_state_features(FEATURES.REM,'REM');

% -------------------------------------------------------------
% 9) Context-labelled bouts (NREM_preREM vs NREM_preWake, etc.)
% -------------------------------------------------------------
bouts = build_context_bouts(score, t_scores, CODES);

idxNREM_preREM  = find(strcmp({bouts.state},'NREM') & strcmp({bouts.next_state},'REM'));
idxNREM_preWake = find(strcmp({bouts.state},'NREM') & strcmp({bouts.next_state},'Wake'));

idxWake_postREM  = find(strcmp({bouts.state},'Wake') & strcmp({bouts.prev_state},'REM'));
idxWake_postNREM = find(strcmp({bouts.state},'Wake') & strcmp({bouts.prev_state},'NREM'));

OUT.context.bouts          = bouts;
OUT.context.idxNREM_preREM = idxNREM_preREM;
OUT.context.idxNREM_preWake = idxNREM_preWake;
OUT.context.idxWake_postREM  = idxWake_postREM;
OUT.context.idxWake_postNREM = idxWake_postNREM;

end

% ==== helper: plot signal with state-coloured line segments ============
function plot_state_coloured_trace(t_sig, sig, ...
    t_scores_win, epoch_sec, run_starts, run_ends, run_codes, stateColors)

hold on;
if isempty(run_starts)
    plot(t_sig, sig, 'k');
    return;
end

for r = 1:numel(run_starts)
    k1 = run_starts(r);
    k2 = run_ends(r);
    code = run_codes(r);
    if ~isKey(stateColors, code), continue; end
    c = stateColors(code);
    t0 = t_scores_win(k1);
    t1 = t_scores_win(k2) + epoch_sec;
    idx = t_sig >= t0 & t_sig <= t1;
    if nnz(idx) >= 2
        plot(t_sig(idx), sig(idx), 'Color', c, 'LineWidth', 1);
    end
end
end

% ==== helper: find stable transitions ===================================
function idx_keep = find_stable_transitions(score, prev_code, next_code, ...
                                            minPreBins, minPostBins)
n = numel(score);
idx_candidates = find(score(2:end) == next_code & ...
                      score(1:end-1) == prev_code) + 1;
idx_keep = [];
for k = 1:numel(idx_candidates)
    i = idx_candidates(k);
    i_pre_start  = max(1, i - minPreBins);
    pre_segment  = score(i_pre_start:i-1);
    i_post_end   = min(n, i + minPostBins - 1);
    post_segment = score(i:i_post_end);
    if numel(pre_segment)  < minPreBins,  continue; end
    if numel(post_segment) < minPostBins, continue; end
    if all(pre_segment == prev_code) && all(post_segment == next_code)
        idx_keep(end+1,1) = i; %#ok<AGROW>
    end
end
end

% ==== helper: onset features per list of times ==========================
function F = compute_onset_features(ach, fs, t_events, t_pre, t_post)
F = struct('n',0,'deltaF_mean',NaN,'slope_mean',NaN,...
           'deltaF_all',[],'slope_all',[]);
if isempty(t_events)
    return;
end

deltaF = nan(numel(t_events),1);
slope  = nan(numel(t_events),1);

for k = 1:numel(t_events)
    t0 = t_events(k);
    idx0 = round(t0*fs) + 1;
    i1 = max(1, idx0 - round(t_pre*fs));
    i2 = min(numel(ach), idx0 + round(t_post*fs));
    t_win = ((i1:i2)-1)/fs - t0;
    y = ach(i1:i2);

    base_idx = t_win < 0;
    post_idx = t_win >= 0 & t_win <= t_post;
    if ~any(base_idx) || ~any(post_idx)
        continue;
    end
    base = mean(y(base_idx),'omitnan');
    peak = max(y(post_idx));
    deltaF(k) = peak - base;

    X = [t_win(:) ones(numel(t_win),1)];
    b = X\y(:);
    slope(k) = b(1);
end

valid = ~isnan(deltaF) & ~isnan(slope);
F.n            = sum(valid);
F.deltaF_all   = deltaF(valid);
F.slope_all    = slope(valid);
F.deltaF_mean  = mean(F.deltaF_all,'omitnan');
F.slope_mean   = mean(F.slope_all,'omitnan');
end

% ==== helper: estimate cycle frequency for one state ====================
function cycleHz = estimate_cycle_freq_state(ach, fs, score, t_scores, stateCode, t_start, t_end)
cycleHz = NaN;

idx = (t_scores >= t_start) & (t_scores <= t_end) & (score == stateCode);
if ~any(idx), return; end

mask = false(size(ach));
epoch_sec = mode(diff(t_scores));
bins = find(idx);
for b = bins'
    t0 = t_scores(b);
    t1 = t0 + epoch_sec;
    i1 = max(1, floor(t0*fs)+1);
    i2 = min(numel(ach), floor(t1*fs));
    mask(i1:i2) = true;
end

ach_state = ach(mask);
if numel(ach_state) < fs*5
    return;
end

sig = detrend(ach_state(:));
sig = sig / std(sig);

[pks,locs] = findpeaks(sig, ...
    'MinPeakHeight', 0.5, ...
    'MinPeakDistance', round(fs*0.2), ...
    'MinPeakProminence', 0.3);

if numel(locs) < 3, return; end

ipi = diff(locs)/fs;
med_i = median(ipi);
sd_i  = std(ipi);
valid = abs(ipi - med_i) < 3*sd_i;
if sum(valid) < 2, return; end

period = mean(ipi(valid));
cycleHz = 1/period;
if cycleHz < 0.01 || cycleHz > 5
    cycleHz = NaN;
end
end

% ==== helper: feature plots per state ===================================
function plot_state_features(F, stateName)
if F.n == 0
    return;
end

figure('Color','w','Position',[200 200 700 300]);
subplot(1,2,1);
histogram(F.deltaF_all, 'FaceColor',[0.3 0.3 0.8]);
xlabel('\DeltaF/F at onset');
ylabel('Count');
title(sprintf('%s: \\DeltaF/F (n=%d)', stateName, F.n));

subplot(1,2,2);
histogram(F.slope_all, 'FaceColor',[0.8 0.3 0.3]);
xlabel('Slope (\DeltaF/F per s)');
ylabel('Count');
title(sprintf('%s: slopes', stateName));

sgtitle(sprintf('%s ACh onset features (cycle f = %.3f Hz)', ...
    stateName, F.cycleHz), 'FontWeight','bold');
end

% ==== helper: 50 Hz notch with skip list ================================
function eeg_out = apply_50Hz_notch_if_needed(eeg_in, fs_eeg, mat_file)
% APPLY_50HZ_NOTCH_IF_NEEDED
% -------------------------------------------------------------
% Notch-filter EEG at 50 Hz for all files EXCEPT a skip list.

    eeg_out = eeg_in;

    % if sampling rate is too low, skip
    if fs_eeg <= 120
        return;
    end

    [~, base, ext] = fileparts(mat_file);
    fname = [base ext];

    skip_list = { ...
        '20251001_baseline_mouse1_APP.mat', ...
        '20251002_baseline_mouse2_WT.mat', ...
        '20251003_ambtemp_mouse1_APP.mat', ...
        '20251005_baseline_mouse8_WT.mat', ...
        '20251006_baseline_mouse4_WT.mat' ...
        };

    % Do NOT notch these files
    if ismember(fname, skip_list)
        return;
    end

    % 50 Hz notch
    wo = 50 / (fs_eeg/2);   % normalized frequency
    bw = wo / 35;          % Q ~ 35
    [bN,aN] = iirnotch(wo, bw);

    eeg_out = filtfilt(bN, aN, eeg_in);
end
