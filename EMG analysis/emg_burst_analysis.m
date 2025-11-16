function OUT = emg_burst_analysis(mat_file, scores_csv, varargin)
% EMG_BURST_ANALYSIS
% -------------------------------------------------------------------------
% Detect EMG bursts and quantify them per sleep bout and per state.
%
% INPUTS
%   mat_file   : .mat file with EMG trace and its sampling rate.
%                Expected fields (any of these):
%                    EMG trace : 'emg','EMG'
%                    fs (Hz)   : 'emg_frequency','eeg_frequency','Fs_emg'
%
%   scores_csv : CSV with 1-Hz sleep scores, with columns:
%                   time_s  (0,1,2,...)
%                   score   (integer code per second)
%
% OPTIONAL NAME/VALUE PAIRS
%   'codes'        : struct with fields .WK .NREM .REM .MA (defaults 0,1,2,15)
%   'epoch_sec'    : epoch length of scoring (default 1 s)
%   'hp_cutoff'    : EMG high-pass cutoff (Hz), default 10
%   'lp_env'       : envelope low-pass cutoff (Hz), default 10
%   'burst_z_thr'  : z-threshold on envelope, default 2.5
%   'burst_min_dur': minimum burst duration (s), default 0.05
%   'burst_merge_gap': max gap to merge bursts (s), default 0.05
%
%   'mouse_id'     : string label
%   'session'      : string label
%   'out_prefix'   : base name for export files
%   'out_dir'      : directory to save outputs (default: alongside CSV)
%   'verbose'      : true/false (default true)
%
% OUTPUT (struct OUT)
%   OUT.bursts : table of individual bursts (start_s,end_s,dur_s,peak_amp,
%                peak_z, state)
%   OUT.bouts  : table of sleep bouts with EMG burst metrics:
%                   state, state_name, start_s, end_s, dur_s,
%                   n_bursts, bursts_per_min, mean_peak_amp
%   OUT.state_summary : table summarising metrics per state (Wake/NREM/REM)
%   OUT.files.fig     : path to per-session EMG burst figure
%
% Requires helper: runs_from_codes.m

%% -------------------- Parse inputs ---------------------------
p = inputParser;
addRequired(p,'mat_file',@ischar);
addRequired(p,'scores_csv',@ischar);

addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);

addParameter(p,'hp_cutoff',10,@(x)isscalar(x)&&x>0);
addParameter(p,'lp_env',10,@(x)isscalar(x)&&x>0);
addParameter(p,'burst_z_thr',2.5,@(x)isscalar(x)&&x>0);
addParameter(p,'burst_min_dur',0.05,@(x)isscalar(x)&&x>0);
addParameter(p,'burst_merge_gap',0.05,@(x)isscalar(x)&&x>=0);

addParameter(p,'mouse_id','',@ischar);
addParameter(p,'session','',@ischar);
addParameter(p,'out_prefix','',@ischar);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'notch_freqs',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'notch_Q',30,@(x)isscalar(x)&&x>0);

parse(p, mat_file, scores_csv, varargin{:});
S       = p.Results;
C       = S.codes;
verbose = S.verbose;

if isempty(S.out_dir)
    S.out_dir = fileparts(scores_csv);
    if isempty(S.out_dir), S.out_dir = pwd; end
end
if ~exist(S.out_dir,'dir'); mkdir(S.out_dir); end

OUT        = struct();
OUT.files  = struct();
OUT.success = false;

if verbose
    fprintf('\n=== EMG burst analysis ===\n');
    fprintf('MAT file    : %s\n', mat_file);
    fprintf('Scores CSV  : %s\n', scores_csv);
    fprintf('Output dir  : %s\n', S.out_dir);
end

%% -------------------- Load scores (1 Hz) ---------------------
T = readtable(scores_csv);
assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
       'Expected columns time_s and score in %s', scores_csv);

t_sec = double(T.time_s(:));   % 0,1,2,...
code  = double(T.score(:));

n_epoch      = numel(code);
dur_scores_s = n_epoch * S.epoch_sec;

if verbose
    fprintf('Scores: %d epochs (%.1f min)\n', n_epoch, dur_scores_s/60);
    fprintf('Unique codes: %s\n', mat2str(unique(code)'));
end

%% -------------------- Load EMG trace -------------------------
Smat = load(mat_file);

cand_emg   = {'emg','EMG'};
cand_fs    = {'emg_frequency','eeg_frequency','Fs_emg','fs_emg'};

emg      = [];
emg_name = '';
for k = 1:numel(cand_emg)
    if isfield(Smat,cand_emg{k})
        emg      = Smat.(cand_emg{k});
        emg_name = cand_emg{k};
        break;
    end
end
if isempty(emg)
    error('Could not find EMG variable in %s', mat_file);
end
emg = double(emg(:));

fs_emg   = [];
fs_name  = '';
for k = 1:numel(cand_fs)
    if isfield(Smat,cand_fs{k})
        fs_emg  = Smat.(cand_fs{k});
        fs_name = cand_fs{k};
        break;
    end
end
if isempty(fs_emg)
    error('Could not find EMG sampling rate in %s', mat_file);
end
fs_emg = double(fs_emg);

if verbose
    fprintf('Using EMG variable: %s (n = %d samples)\n', emg_name, numel(emg));
    fprintf('EMG sampling rate : %s = %.3f Hz\n', fs_name, fs_emg);
end

%% -------------------- Align EMG with scores ------------------
n_emg    = numel(emg);
t_emg    = (0:n_emg-1)' / fs_emg;
sec_idx  = floor(t_emg / S.epoch_sec) + 1;
valid    = sec_idx >= 1 & sec_idx <= n_epoch;

emg       = emg(valid);
t_emg     = t_emg(valid);
sec_idx   = sec_idx(valid);
state_vec = code(sec_idx);

if verbose
    fprintf('EMG duration total : %.1f min\n', n_emg/fs_emg/60);
    fprintf('Overlap with scores: %.1f min\n', max(t_emg)/60);
end

%% -------------------- Preprocess EMG (envelope) ---------------
hp = S.hp_cutoff;
lp = S.lp_env;
% -------------------- Optional notch filtering -------------------------
emg_filt = emg;

if ~isempty(S.notch_freqs)
    nf = S.notch_freqs(:)';
    for f0 = nf
        if f0 <= 0 || f0 >= fs_emg/2
            continue; % skip impossible frequencies
        end
        wo = f0/(fs_emg/2);
        bw = wo/S.notch_Q;
        [bn,an] = iirnotch(wo, bw);
        emg_filt = filtfilt(bn, an, emg_filt);
    end
end

% -------------------- Band-pass & envelope -----------------------------
[b,a] = butter(2, [hp 200]/(fs_emg/2), 'bandpass');
emg_bp = filtfilt(b,a, emg_filt);

env_raw = abs(emg_bp);

win_samp = round(0.05 * fs_emg);  % 50 ms smoothing
if win_samp < 1, win_samp = 1; end
env_smooth = movmean(env_raw, win_samp);

% baseline from NREM if available, otherwise all data
base_mask = (state_vec == C.NREM);
if nnz(base_mask) < fs_emg * 10   % less than 10 s of NREM → use all
    base_mask = true(size(env_smooth));
end

mu_env  = mean(env_smooth(base_mask));
sd_env  = std(env_smooth(base_mask));
env_z   = (env_smooth - mu_env) / max(sd_env, eps);

%% -------------------- Detect bursts ---------------------------
thr     = S.burst_z_thr;
min_dur = S.burst_min_dur;
merge_gap = S.burst_merge_gap;

is_high = env_z > thr;
d       = diff([0; is_high; 0]);
idx_start = find(d==1);
idx_end   = find(d==-1) - 1;

% durations
dur = (idx_end - idx_start + 1) / fs_emg;
keep = dur >= min_dur;
idx_start = idx_start(keep);
idx_end   = idx_end(keep);

% merge bursts that are too close
if merge_gap > 0 && ~isempty(idx_start)
    new_start = idx_start(1);
    new_end   = idx_end(1);
    merged_st = [];
    merged_en = [];
    for i = 2:numel(idx_start)
        gap_s = (idx_start(i) - new_end - 1) / fs_emg;
        if gap_s <= merge_gap
            new_end = idx_end(i);  % extend current burst
        else
            merged_st(end+1,1) = new_start; %#ok<AGROW>
            merged_en(end+1,1) = new_end;   %#ok<AGROW>
            new_start = idx_start(i);
            new_end   = idx_end(i);
        end
    end
    merged_st(end+1,1) = new_start;
    merged_en(end+1,1) = new_end;
    idx_start = merged_st;
    idx_end   = merged_en;
end

n_bursts = numel(idx_start);
if verbose
    fprintf('Detected %d EMG bursts.\n', n_bursts);
end

% build burst table
if n_bursts == 0
    bursts_tbl = table();
else
    start_s = t_emg(idx_start);
    end_s   = t_emg(idx_end);
    dur_s   = end_s - start_s;

    peak_amp = nan(n_bursts,1);
    peak_z   = nan(n_bursts,1);
    state_mid= nan(n_bursts,1);

    for i = 1:n_bursts
        seg = idx_start(i):idx_end(i);
        [peak_amp(i), ~] = max(env_smooth(seg));
        peak_z(i)        = max(env_z(seg));
        mid_idx          = round(mean(seg));
        state_mid(i)     = state_vec(mid_idx);
    end

    bursts_tbl = table(start_s, end_s, dur_s, peak_amp, peak_z, state_mid, ...
        'VariableNames',{'start_s','end_s','dur_s','peak_amp','peak_z','state'});
end
OUT.bursts = bursts_tbl;

%% -------------------- Bout-level metrics ----------------------
[st_b, en_b] = runs_from_codes(code);
n_bouts = numel(st_b);

state_name = cell(n_bouts,1);
start_b_s  = zeros(n_bouts,1);
end_b_s    = zeros(n_bouts,1);
dur_b_s    = zeros(n_bouts,1);
n_burst_b  = zeros(n_bouts,1);
rate_b     = zeros(n_bouts,1);
mean_amp_b = nan(n_bouts,1);

for i = 1:n_bouts
    s_code = code(st_b(i));
    switch s_code
        case C.WK,   state_name{i} = 'Wake';
        case C.NREM, state_name{i} = 'NREM';
        case C.REM,  state_name{i} = 'REM';
        case C.MA,   state_name{i} = 'MA';
        otherwise,   state_name{i} = sprintf('code%d',s_code);
    end

    start_b_s(i) = t_sec(st_b(i));
    end_b_s(i)   = t_sec(en_b(i));
    dur_b_s(i)   = end_b_s(i) - start_b_s(i) + 1;

    if n_bursts > 0
        mid_burst_s = (bursts_tbl.start_s + bursts_tbl.end_s)/2;
        in_bout = mid_burst_s >= start_b_s(i) & mid_burst_s <= end_b_s(i);
        n_burst_b(i) = nnz(in_bout);
        rate_b(i)    = n_burst_b(i) / (dur_b_s(i)/60);  % bursts/min
        mean_amp_b(i)= mean(bursts_tbl.peak_amp(in_bout),'omitnan');
    else
        n_burst_b(i) = 0;
        rate_b(i)    = 0;
        mean_amp_b(i)= NaN;
    end
end

bouts_tbl = table(code(st_b), state_name, start_b_s, end_b_s, dur_b_s, ...
                  n_burst_b, rate_b, mean_amp_b, ...
       'VariableNames', {'state_code','state_name','start_s','end_s', ...
                         'dur_s','n_bursts','bursts_per_min','mean_peak_amp'});
OUT.bouts = bouts_tbl;

%% -------------------- State-wise summary ----------------------
state_labels = {'Wake','NREM','REM'};
state_codes  = [C.WK,  C.NREM,  C.REM];

nStates = numel(state_labels);
state_col     = cell(nStates,1);
n_bouts_col   = zeros(nStates,1);
rate_mean_col = nan(nStates,1);
rate_sem_col  = nan(nStates,1);
amp_mean_col  = nan(nStates,1);
amp_sem_col   = nan(nStates,1);

for i = 1:nStates
    sc = state_codes(i);
    mask = bouts_tbl.state_code == sc;
    state_col{i}   = state_labels{i};
    n_bouts_col(i) = nnz(mask);

    rates = bouts_tbl.bursts_per_min(mask);
    amps  = bouts_tbl.mean_peak_amp(mask);

    rate_mean_col(i) = mean(rates,'omitnan');
    amp_mean_col(i)  = mean(amps,'omitnan');

    rate_sem_col(i) = std(rates,'omitnan') / max(1,sqrt(sum(~isnan(rates))));
    amp_sem_col(i)  = std(amps,'omitnan') / max(1,sqrt(sum(~isnan(amps))));
end

OUT.state_summary = table(state_col, n_bouts_col, ...
    rate_mean_col, rate_sem_col, amp_mean_col, amp_sem_col, ...
    'VariableNames', {'state','n_bouts', ...
                      'bursts_per_min_mean','bursts_per_min_sem', ...
                      'peak_amp_mean','peak_amp_sem'});

%% -------------------- Simple per-session figure ----------------
fig = figure('Name','EMG bursts','Color','w','Visible','on');

% Panel 1: short EMG snippet with bursts
subplot(1,2,1); hold on;
if n_emg > fs_emg*60
    % last 60 s
    t_plot = t_emg(t_emg >= max(t_emg)-60);
    emg_plot = emg(end-numel(t_plot)+1:end);
else
    t_plot = t_emg;
    emg_plot = emg;
end
plot(t_plot, emg_plot);
xlabel('Time (s)');
ylabel('EMG (raw)');
title('Example EMG segment');
box off;

% mark bursts within this segment
if n_bursts > 0
    for i = 1:n_bursts
        if bursts_tbl.start_s(i) >= t_plot(1) && bursts_tbl.start_s(i) <= t_plot(end)
            xline(bursts_tbl.start_s(i),'r:');
        end
    end
end

% Panel 2: bursts per minute per state
subplot(1,2,2); hold on;
sm = OUT.state_summary;
cats = categorical(sm.state);
bar(cats, sm.bursts_per_min_mean, 'FaceColor',[0.6 0.6 0.6]);
errorbar(cats, sm.bursts_per_min_mean, sm.bursts_per_min_sem, ...
         'k','LineStyle','none','CapSize',8);
ylabel('Bursts per min');
title('EMG bursts per state');
box off;

sgtitle(sprintf('%s – %s : EMG bursts', S.mouse_id, S.session));

if isempty(S.out_prefix)
    base_fig = sprintf('EMGBursts_%s_%s', S.mouse_id, S.session);
else
    base_fig = sprintf('%s_EMGBursts', S.out_prefix);
end
fig_file = fullfile(S.out_dir, [base_fig '.png']);
saveas(fig, fig_file);
OUT.files.fig = fig_file;

OUT.success = true;
end
