function run_phase1_full_analysis()
% =================== Phase-1 FULL ANALYSIS (single .m) ===================
% Requires the scored .mat to include:
%   eeg, emg, ne (or ACh ΔF/F), sleep_scores, eeg_frequency, ne_frequency
% You will map score codes (WAKE/NREM/REM) via a quick dialog.
%
% Output:
%   _results_phase1/*.png      _processed_phase1/*.csv
% ========================================================================

%% -------------------- User Settings --------------------
VARS.eeg    = 'eeg';
VARS.emg    = 'emg';
VARS.ach    = 'ne';               % ACh/NE ΔF/F
VARS.scores = 'sleep_scores';     % vector of integer codes
VARS.fs_eeg = 'eeg_frequency';
VARS.fs_ach = 'ne_frequency';     % e.g., 10.1725 Hz

% epoch length (seconds) of your sleep scoring (use 1 for 1 Hz labels)
EPOCH_SEC = 1;

% analysis params
SIGMA_BAND     = [10 15];         % Hz
ACH_PSD_MAXHZ  = 0.15;            % Hz (NREM slow component)
REM_POST_WIN_S = [40 60];         % seconds after REM onset
MT_TAPS        = 3.5;             % pmtm time-bandwidth (auto-capped)

EEG_BANDS = struct('delta',[0.5 4], 'theta',[5 9], 'sigma',[10 15], ...
                   'beta',[16 30], 'lowgamma',[30 55], 'highgamma',[65 100]);

DO_SPINDLES = true;
SP_PARAMS.sigma_band    = SIGMA_BAND;
SP_PARAMS.env_smooth_s  = 0.2;    % s
SP_PARAMS.thr_z         = 2.0;    % z-threshold
SP_PARAMS.min_dur_s     = 0.4;    % s
SP_PARAMS.max_dur_s     = 3.0;    % s

% PSD fallbacks
VERY_SHORT_N = 16;    % <16 samples -> periodogram
SHORT_SEC    = 5;     % <5 s       -> Welch

%% -------------------- Pick one .mat --------------------
[fname, fdir] = uigetfile({'*.mat','MAT-files (*.mat)'}, 'Select one scored .mat file');
if isequal(fname,0), error('Canceled. No file selected.'); end
fpath = fullfile(fdir, fname);
ses_id = erase(fname, ".mat");

% Output folders next to the file
% --- Output folders (session subfolder) ---
safe_id = regexprep(ses_id, '[^\w\-]', '_');

DATA_ROOT = fullfile(fdir, '_processed_phase1');
FIG_ROOT  = fullfile(fdir, '_results_phase1');
DATA_OUT  = fullfile(DATA_ROOT, safe_id);
FIG_OUT   = fullfile(FIG_ROOT,  safe_id);

if ~exist(DATA_OUT,'dir'), mkdir(DATA_OUT); end
if ~exist(FIG_OUT,'dir'),  mkdir(FIG_OUT);  end


%% -------------------- Load & extract --------------------
S = load(fpath);

eeg    = must_get(S, VARS.eeg);
emg    = must_get(S, VARS.emg);
ach    = must_get(S, VARS.ach);
scores = must_get(S, VARS.scores);

fs_eeg = must_get_scalar(S, VARS.fs_eeg);
fs_ach = must_get_scalar(S, VARS.fs_ach);

% Sizes & time vectors
Neeg = numel(eeg); Nemg = numel(emg); Nach = numel(ach);
t_eeg = (0:Neeg-1).' / fs_eeg;                 % seconds
t_emg = (0:Nemg-1).' / fs_eeg;                 % assuming EMG at fs_eeg
t_ach = (0:Nach-1).' / fs_ach;                 % seconds (e.g., 10.1725 Hz)

fprintf('EEG duration: %.1f min @ %.4f Hz | ACh duration: %.1f min @ %.4f Hz\n', ...
    Neeg/fs_eeg/60, fs_eeg, Nach/fs_ach/60, fs_ach);

%% -------------------- Expand epoch labels to EEG length --------------------
scores_ep = scores(:);
Nepochs   = numel(scores_ep);
epoch_len_samples = max(1, round(fs_eeg * EPOCH_SEC));

% repeat each label to sample rate, then trim/pad to EEG length
scores_samp = repelem(scores_ep, epoch_len_samples);
if numel(scores_samp) > Neeg
    scores_samp = scores_samp(1:Neeg);
elseif numel(scores_samp) < Neeg
    scores_samp(end+1:Neeg,1) = scores_samp(end);
end
fprintf('Epoch set to %.3f s (%d samples). Nepochs=%d\n', EPOCH_SEC, epoch_len_samples, Nepochs);

% Map codes explicitly (no guessing) + optional EXCLUDE
[isWake_EEG, isNREM_EEG, isREM_EEG, codeMap, excludeMask] = map_scores_prompt(scores_samp);

% Report unlabeled time (not any of the 3 states, not excluded)
unlabeled = ~(isWake_EEG | isNREM_EEG | isREM_EEG) & ~excludeMask;
if nnz(unlabeled) > 0
    miss_sec = nnz(unlabeled)/fs_eeg;
    fprintf('Warning: %.2f s (%.2f min) not labeled as Wake/NREM/REM.\n', miss_sec, miss_sec/60);
end

% Keep only not-excluded (but still keep unlabeled for plotting info)
keep_eeg    = ~excludeMask;
isWake_EEG  = isWake_EEG & keep_eeg;
isNREM_EEG  = isNREM_EEG & keep_eeg;
isREM_EEG   = isREM_EEG  & keep_eeg;

% Map masks to ACh time (nearest)
isWake_ACh = logical(interp1(t_eeg, double(isWake_EEG), t_ach, 'nearest', 0));
isNREM_ACh = logical(interp1(t_eeg, double(isNREM_EEG), t_ach, 'nearest', 0));
isREM_ACh  = logical(interp1(t_eeg, double(isREM_EEG),  t_ach, 'nearest', 0));

% QC printouts
tot_labeled_samp = nnz(isWake_EEG)+nnz(isNREM_EEG)+nnz(isREM_EEG);
tot_min = tot_labeled_samp/fs_eeg/60;
pctW = 100*nnz(isWake_EEG)/max(tot_labeled_samp,1);
pctN = 100*nnz(isNREM_EEG)/max(tot_labeled_samp,1);
pctR = 100*nnz(isREM_EEG) /max(tot_labeled_samp,1);
fprintf('Architecture (EEG labeled): %.1f min | %%Wake/%%NREM/%%REM = %.1f / %.1f / %.1f\n', tot_min, pctW, pctN, pctR);

%% -------------------- Full-Session Traces by State (no spectrogram) ----
dash_fig = plot_fullsession_traces_shaded( ...
    t_eeg, double(eeg), double(emg), t_ach, double(ach), ...
    isWake_EEG, isNREM_EEG, isREM_EEG, fs_eeg, ses_id);
saveas(dash_fig, fullfile(FIG_OUT, [ses_id '_traces_by_state.png']));
close(dash_fig);

%% -------------------- Analyses --------------------
% 1) Sleep architecture (EEG timeline)
[arch_tbl, arch_fig] = sleep_architecture(isWake_EEG, isNREM_EEG, isREM_EEG, fs_eeg, ses_id, codeMap);
saveas(arch_fig, fullfile(FIG_OUT, [ses_id '_architecture.png'])); close(arch_fig);

% 2) EEG PSDs per state (robust MT/Welch fallback)
[~, psdEEG_fig] = eeg_psd_statewise(eeg, fs_eeg, ...
    {isWake_EEG, isNREM_EEG, isREM_EEG}, {'Wake','NREM','REM'}, MT_TAPS, SHORT_SEC, VERY_SHORT_N);
saveas(psdEEG_fig, fullfile(FIG_OUT, [ses_id '_EEG_PSD.png'])); close(psdEEG_fig);

% 3) EEG band powers
band_tbl = eeg_bandpowers(eeg, fs_eeg, ...
    {isWake_EEG, isNREM_EEG, isREM_EEG}, {'Wake','NREM','REM'}, EEG_BANDS);

% 4) ACh/NE PSD in NREM (ACh timebase)
[psdACh_tbl, psdACh_fig] = ach_psd_nrem_TIME(ach, fs_ach, isNREM_ACh, MT_TAPS, ACH_PSD_MAXHZ, SHORT_SEC, VERY_SHORT_N);
saveas(psdACh_fig, fullfile(FIG_OUT, [ses_id '_ACh_PSD_NREM.png'])); close(psdACh_fig);

% 5) Sigma & simple spindles (EEG timeline)
[sigma_tbl, sp_fig] = sigma_summary_and_spindles(eeg, fs_eeg, isNREM_EEG, SIGMA_BAND, DO_SPINDLES, SP_PARAMS);
saveas(sp_fig, fullfile(FIG_OUT, [ses_id '_SigmaSpindles.png'])); close(sp_fig);

% 6) REM-locked ACh (ACh timebase)
rem_tbl = rem_locked_ach_TIME(ach, fs_ach, isREM_ACh, REM_POST_WIN_S, ses_id);

%% -------------------- Save tables --------------------
writetable(arch_tbl,   fullfile(DATA_OUT, [ses_id '_architecture.csv']));
writetable(band_tbl,   fullfile(DATA_OUT, [ses_id '_EEG_bandpowers.csv']));
writetable(psdACh_tbl, fullfile(DATA_OUT, [ses_id '_ACh_PSD_NREM.csv']));
writetable(sigma_tbl,  fullfile(DATA_OUT, [ses_id '_sigma_spindles.csv']));
writetable(rem_tbl,    fullfile(DATA_OUT, [ses_id '_REM_locked_ACh.csv']));

% One-row summary
summary_row = summary_pack(ses_id, arch_tbl, band_tbl, psdACh_tbl, sigma_tbl, rem_tbl);
writetable(struct2table(summary_row), fullfile(DATA_OUT, [ses_id '_SUMMARY_phase1.csv']));

disp('Done ✔  Dashboard + CSVs/PNGs saved next to your .mat');
end

% ===================== Helpers =====================
function x = must_get(S, name)
    assert(isfield(S,name), 'Missing variable "%s" in .mat', name);
    x = S.(name); x = x(:);
end
function fs = must_get_scalar(S, name)
    assert(isfield(S,name), 'Missing "%s" in .mat', name);
    fs = S.(name);
    assert(isscalar(fs) && isfinite(fs) && fs>0, 'Invalid "%s"', name);
end

function [isW,isN,isR,codeMap,excludeMask] = map_scores_prompt(scores_samp)
    scores_samp = scores_samp(:);
    u = sort(unique(scores_samp(~isnan(scores_samp))));
    answ = inputdlg({sprintf('Unique codes found: %s\n\nCode for WAKE:', mat2str(u(:)')),...
                     'Code for NREM:','Code for REM:','Code to EXCLUDE (optional):'}, ...
                     'Map sleep score codes (explicit)', 1, {'','','',''});
    if isempty(answ), error('Canceled by user.'); end
    codeMap = struct('WAKE',str2double(answ{1}), 'NREM',str2double(answ{2}), 'REM',str2double(answ{3}));
    excl = NaN; if ~isempty(answ{4}), excl = str2double(answ{4}); end
    want = [codeMap.WAKE, codeMap.NREM, codeMap.REM];
    if any(~ismember(want, u)), error('Provided codes not found in data.'); end
    if numel(unique(want))<3, error('Codes for WAKE/NREM/REM must be distinct.'); end
    excludeMask = false(size(scores_samp));
    if ~isnan(excl), excludeMask = (scores_samp==excl); end
    isW = (scores_samp==codeMap.WAKE) & ~excludeMask;
    isN = (scores_samp==codeMap.NREM) & ~excludeMask;
    isR = (scores_samp==codeMap.REM)  & ~excludeMask;
end

function [T, fig] = sleep_architecture(isW,isN,isR, fs, session_name, codeMap)
    secsW = nnz(isW)/fs; secsN = nnz(isN)/fs; secsR = nnz(isR)/fs; total = secsW+secsN+secsR;
    pct = 100*[secsW,secsN,secsR]/max(total,eps);
    [cntW,durW] = bouts_and_durations(isW, fs);
    [cntN,durN] = bouts_and_durations(isN, fs);
    [cntR,durR] = bouts_and_durations(isR, fs);
    T = table( string(session_name), total/60, secsW/60, secsN/60, secsR/60, ...
        pct(1), pct(2), pct(3), cntW, median(durW), cntN, median(durN), cntR, median(durR), ...
        codeMap.WAKE, codeMap.NREM, codeMap.REM, ...
        'VariableNames', {'session','minutes_total','minutes_wake','minutes_nrem','minutes_rem', ...
        'pct_wake','pct_nrem','pct_rem','bouts_wake','boutDurMed_wake','bouts_nrem','boutDurMed_nrem','bouts_rem','boutDurMed_rem', ...
        'code_WAKE','code_NREM','code_REM'});
    fig = figure('Color','w','Name','Sleep Architecture');
    tiledlayout(1,2,'padding','compact'); nexttile;
    bar([pct(1),pct(2),pct(3)]); set(gca,'XTickLabel',{'Wake','NREM','REM'}); ylabel('% of time'); title(strrep(session_name,'_','\_'));
    nexttile; bar([cntW,cntN,cntR]); set(gca,'XTickLabel',{'Wake','NREM','REM'}); ylabel('Bout count');
end

% Additional helpers for other analyses (ACh PSD, Spindles, etc.)
