%% ==============================================================
%  ACh onset + whole-bout stability + NREM PLV + NREM PAC
%  mouse1 APP, ONE episode per condition:
%     - baseline : BL-late
%     - ambtemp  : AMB-early
%     - drugs    : DR-late
% ==============================================================

clearvars; clear functions; close all; clc;

% --------------------------------------------------------------
% 1) Paths
% --------------------------------------------------------------
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';

CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

% --------------------------------------------------------------
% 2) Define ONE episode per condition
% --------------------------------------------------------------
segs = struct([]);

% === BASELINE: use BL-late ===
segs(1).mat_file   = fullfile(sigDir,  '20251001_baseline_mouse1_APP.mat');
segs(1).scores_csv = fullfile(scoreDir,'20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
segs(1).t_start    = 19704;      % BL-late
segs(1).t_end      = 21000;
segs(1).cond       = 'baseline';
segs(1).label      = 'BL-late';

% === AMBTEMP: use AMB-early ===
segs(2).mat_file   = fullfile(sigDir,  '20251003_ambtemp_mouse1_APP.mat');
segs(2).scores_csv = fullfile(scoreDir,'20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv');
segs(2).t_start    = 10966;      % AMB-early
segs(2).t_end      = 12770;
segs(2).cond       = 'ambtemp';
segs(2).label      = 'AMB-early';

% === DRUGS: use DR-late ===
segs(3).mat_file   = fullfile(sigDir,  '20251022_drugs_mouse1_APP.mat');
segs(3).scores_csv = fullfile(scoreDir,'20251022_drugs_mouse1_APP_scored_scores_1Hz.csv');
segs(3).t_start    = 22422;      % DR-late
segs(3).t_end      = 24576;
segs(3).cond       = 'drugs';
segs(3).label      = 'DR-late';

nSeg = numel(segs);

% --------------------------------------------------------------
% 3) Overview plots + onset features (one OUT per episode)
%    Use CELL ARRAY to avoid struct mismatch issues
% --------------------------------------------------------------
OUT = cell(1, nSeg);

for k = 1:nSeg
    fprintf('\n=== Episode %d: %s (%s) ===\n', ...
            k, segs(k).label, segs(k).cond);

    tmp = plot_ach_eeg_segment( ...
              segs(k).mat_file, ...
              segs(k).scores_csv, ...
              segs(k).t_start, ...
              segs(k).t_end, ...
              'codes', CODES);

    tmp.cond  = segs(k).cond;
    tmp.label = segs(k).label;
    tmp.codes = CODES;

    OUT{k} = tmp;
end

episode_labels = {segs.label};   % {'BL-late','AMB-early','DR-late'}
cond_labels    = {segs.cond};    % {'baseline','ambtemp','drugs'}

% Convert once to struct array for downstream functions
OUTs = [OUT{:}];   % 1×3 struct

% --------------------------------------------------------------
% 4) Onset-based comparison (per episode & per condition)
% --------------------------------------------------------------
t_pre_peri  = 20;       % sec before onset
t_post_peri = 40;       % sec after onset
dt          = 0.2;      % peri-onset time step

% (a) per-episode comparison
STATS_ep = compare_ach_segments_onsets_features( ...
               OUTs, episode_labels, ...
               t_pre_peri, t_post_peri);

% (b) per-condition comparison (baseline vs ambtemp vs drugs)
[COND_onset, STATS_onset] = ach_compare_conditions_onsets( ...
                               OUTs, cond_labels, ...
                               t_pre_peri, t_post_peri, dt);

ach_plot_condition_summary_allStates(COND_onset, STATS_onset);

% --------------------------------------------------------------
% 5) Whole-bout stability features (ACh over whole bouts)
% --------------------------------------------------------------
[BOUT, STATS_BOUT] = ach_compute_bout_features_conditions(OUTs, cond_labels);
ach_plot_bout_features_allStates(BOUT, STATS_BOUT);

% --------------------------------------------------------------
% 6) NREM whole-bout PLV: delta EEG (1–4 Hz) vs slow ACh (0.1–1 Hz)
%    One PLV struct per condition
% --------------------------------------------------------------
PLV_NREM = cell(1, nSeg);   % CELL ARRAY to avoid struct mismatch

for k = 1:nSeg
    PLV_NREM{k} = ach_eeg_plv_state_bouts( ...
        OUTs(k).mat_file, ...
        OUTs(k).scores_csv, ...
        'NREM', ...
        OUTs(k).t_start, ...
        OUTs(k).t_end, ...
        'codes',        OUTs(k).codes, ...
        'eeg_band',     [1 4], ...      % delta phase
        'ach_band',     [0.1 1.0], ...  % slow ACh phase
        'min_bout_sec', 10);
end

% Now convert CELL -> struct array (all have same fields)
PLV_NREM_s = [PLV_NREM{:}];

% Plot summary (PLV per bout per condition + PLV vs duration)
ach_plot_plv_bout_summary(PLV_NREM_s, cond_labels, ...
    'NREM delta–ACh PLV across conditions');

% --------------------------------------------------------------
% 7) NREM PAC comodulogram per condition (slow-phase × sigma-amp)
% --------------------------------------------------------------
phase_freqs = 0.5:0.5:4;    % slow/delta phase
amp_freqs   = 7:1:25;       % amplitude around sigma

COM_cells = cell(1, nSeg);

for k = 1:nSeg
    COM_cells{k} = ach_pac_comod_nrem_segment( ...
                       OUTs(k).mat_file, ...
                       OUTs(k).scores_csv, ...
                       OUTs(k).t_start, ...
                       OUTs(k).t_end, ...
                       OUTs(k).codes, ...
                       'phase_freqs', phase_freqs, ...
                       'amp_freqs',   amp_freqs, ...
                       'nbins',       18, ...
                       'doPlot',      true, ...
                       'label',       sprintf('%s (%s)', ...
                                              OUTs(k).label, ...
                                              OUTs(k).cond));
end

COM = [COM_cells{:}];

% --------------------------------------------------------------
% 8) NREM PAC ROI summary: SO-phase (0.5–1.5 Hz) × sigma-amp (11–16 Hz)
% --------------------------------------------------------------
phase_roi = [0.5 1.5];
amp_roi   = [11 16];

mi_roi = nan(1, nSeg);

for k = 1:nSeg
    if k > numel(COM) || isempty(COM(k).MI)
        continue;
    end
    pf = COM(k).phase_freqs;
    af = COM(k).amp_freqs;
    ip = pf >= phase_roi(1) & pf <= phase_roi(2);
    ia = af >= amp_roi(1)   & af <= amp_roi(2);
    roi_vals = COM(k).MI(ip, ia);
    mi_roi(k) = mean(roi_vals(:), 'omitnan');
end

figure('Color','w','Position',[300 300 450 350]);
bar(mi_roi);
set(gca,'XTick',1:nSeg,'XTickLabel',cond_labels);
ylabel('SO–sigma PAC (MI)');
title('NREM PAC ROI (slow-phase 0.5–1.5 Hz, sigma 11–16 Hz)');
box on; grid on;
