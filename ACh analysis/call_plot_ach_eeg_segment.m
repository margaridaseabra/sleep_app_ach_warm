%% ==============================================================
%  Full ACh + PLV analysis for mouse1 APP
%  Conditions: baseline, ambtemp, drugs
%  Episodes: 2 per condition (early / late)
%
%  PLV transitions (θ–ACh):
%     - Wake → NREM
%     - NREM → REM
%     - REM  → Wake
%
%  Extra PLV band (σ–ACh in NREM):
%     - Wake → NREM, EEG sigma [11–16 Hz]
% ==============================================================

clearvars; clear functions; close all; clc;

% --------------------------------------------------------------
% 1) Paths and common settings
% --------------------------------------------------------------
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';

CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

% --------------------------------------------------------------
% 2) Define all episodes (segments) for this mouse
% --------------------------------------------------------------
segs = struct([]);

% === BASELINE episodes ===
segs(1).mat_file   = fullfile(sigDir,  '20251001_baseline_mouse1_APP.mat');
segs(1).scores_csv = fullfile(scoreDir,'20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
segs(1).t_start    = 19600;
segs(1).t_end      = 20400;
segs(1).cond       = 'baseline';
segs(1).label      = 'BL-early';

segs(2)            = segs(1);
segs(2).t_start    = 19704;
segs(2).t_end      = 21000;
segs(2).label      = 'BL-late';

% === AMBTEMP episodes ===
segs(3).mat_file   = fullfile(sigDir,  '20251003_ambtemp_mouse1_APP.mat');
segs(3).scores_csv = fullfile(scoreDir,'20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv');
segs(3).t_start    = 10966;
segs(3).t_end      = 12770;
segs(3).cond       = 'ambtemp';
segs(3).label      = 'AMB-early';

segs(4)            = segs(3);
segs(4).t_start    = 16428;
segs(4).t_end      = 18163;
segs(4).label      = 'AMB-late';

% === DRUGS episodes ===
segs(5).mat_file   = fullfile(sigDir,  '20251022_drugs_mouse1_APP.mat');
segs(5).scores_csv = fullfile(scoreDir,'20251022_drugs_mouse1_APP_scored_scores_1Hz.csv');
segs(5).t_start    = 18766;
segs(5).t_end      = 20027;
segs(5).cond       = 'drugs';
segs(5).label      = 'DR-early';

segs(6)            = segs(5);
segs(6).t_start    = 22422;
segs(6).t_end      = 24576;
segs(6).label      = 'DR-late';

nSeg = numel(segs);

% --------------------------------------------------------------
% 3) Run ACh/EEG overview + onset features for each episode
%    Use CELL ARRAY to avoid "dissimilar structure" errors.
% --------------------------------------------------------------
OUT = cell(1, nSeg);   % <--- cell, not struct

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

    OUT{k} = tmp;   % store in cell
end

% Handy label arrays
episode_labels = {segs.label};
cond_labels    = {segs.cond};

% Convert ONCE to struct array for functions that want structs
OUTs = [OUT{:}];    % 1×nSeg struct array

% --------------------------------------------------------------
% 4) Per-episode ACh onset comparison
% --------------------------------------------------------------
t_pre_peri  = 20;
t_post_peri = 40;

% STATS_ep = compare_ach_segments_onsets_features( ...
%                 OUTs, episode_labels, ...
%                 t_pre_peri, t_post_peri);

% --------------------------------------------------------------
% 5) Condition-level ACh onset comparison
% --------------------------------------------------------------
dt = 0.2;

[COND_onset, STATS_onset] = ach_compare_conditions_onsets( ...
                              OUTs, cond_labels, ...
                              t_pre_peri, t_post_peri, dt);

ach_plot_condition_summary_allStates(COND_onset, STATS_onset);

% --------------------------------------------------------------
% 6) Bout-wise ACh features across whole bouts
% --------------------------------------------------------------
[BOUT, STATS_BOUT] = ach_compute_bout_features_conditions(OUTs, cond_labels);
ach_plot_bout_features_allStates(BOUT, STATS_BOUT);

% --------------------------------------------------------------
% 7) PLV around multiple transitions, per episode
%    Transitions (θ-band):
%      - Wake → NREM   (theta)
%      - NREM → REM    (theta)
%      - REM  → Wake   (theta)
%    Extra transition (σ-band):
%      - Wake → NREM   (sigma)
% --------------------------------------------------------------
trans_specs = struct([]);

% θ-band transitions
trans_specs(1).from     = 'Wake';
trans_specs(1).to       = 'NREM';
trans_specs(1).field    = 'Wake_NREM_theta';
trans_specs(1).label    = 'Wake→NREM θ–ACh PLV';
trans_specs(1).eeg_band = [5 10];       % theta
trans_specs(1).ach_band = [0.05 0.5];   % slow ACh

trans_specs(2).from     = 'NREM';
trans_specs(2).to       = 'REM';
trans_specs(2).field    = 'NREM_REM_theta';
trans_specs(2).label    = 'NREM→REM θ–ACh PLV';
trans_specs(2).eeg_band = [5 10];
trans_specs(2).ach_band = [0.05 0.5];

trans_specs(3).from     = 'REM';
trans_specs(3).to       = 'Wake';
trans_specs(3).field    = 'REM_Wake_theta';
trans_specs(3).label    = 'REM→Wake θ–ACh PLV';
trans_specs(3).eeg_band = [5 10];
trans_specs(3).ach_band = [0.05 0.5];

% σ-band (NREM-related) – Wake→NREM in sigma band
trans_specs(4).from     = 'Wake';
trans_specs(4).to       = 'NREM';
trans_specs(4).field    = 'Wake_NREM_sigma';
trans_specs(4).label    = 'Wake→NREM σ–ACh PLV';
trans_specs(4).eeg_band = [11 16];      % sigma / spindles-ish
trans_specs(4).ach_band = [0.05 0.5];

% Compute PLV per episode (still using cell OUT for convenience)
for k = 1:nSeg
    fprintf('\n--- PLV for episode %d: %s (%s) ---\n', ...
            k, segs(k).label, segs(k).cond);

    for tt = 1:numel(trans_specs)
        spec = trans_specs(tt);

        PLV_tmp = ach_eeg_plv_periTransitions_segment( ...
                        segs(k).mat_file, ...
                        segs(k).scores_csv, ...
                        spec.from, spec.to, ...
                        segs(k).t_start, segs(k).t_end, ...
                        'codes', CODES, ...
                        'pre_sec', 20, 'post_sec', 40, ...
                        'eeg_band', spec.eeg_band, ...
                        'ach_band', spec.ach_band);

        OUT{k}.plv.(spec.field) = PLV_tmp;  % cell version
        fprintf('  %s: %d transitions used.\n', ...
                spec.field, PLV_tmp.nEvents);
    end
end

% Convert again to struct array now that plv fields exist
OUTs = [OUT{:}];

% --------------------------------------------------------------
% 8) Condition-level PLV comparison for EACH transition
%     (Baseline vs Ambtemp vs Drugs)
% --------------------------------------------------------------
for tt = 1:numel(trans_specs)
    spec = trans_specs(tt);

    % Use cell to avoid "dissimilar structures" problem
    PLV_cells = {};
    cond_plv  = {};
    idx       = 0;

    for k = 1:nSeg
        if isfield(OUTs(k),'plv') && isfield(OUTs(k).plv, spec.field)
            P = OUTs(k).plv.(spec.field);
            if ~isempty(P) && isfield(P,'t_rel') && ~isempty(P.t_rel) ...
                    && isfield(P,'nEvents') && P.nEvents > 0

                idx = idx + 1;
                PLV_cells{idx} = P;             %#ok<SAGROW>
                cond_plv{idx}  = OUTs(k).cond;  %#ok<SAGROW>
            end
        end
    end

    if isempty(PLV_cells)
        fprintf('\n[PLV] No usable data for %s across episodes, skipping.\n', ...
                spec.label);
        continue;
    end

    % Convert cell -> struct array once
    PLV_list = [PLV_cells{:}];   % 1×N struct

    fprintf('\n[PLV] Condition comparison for %s\n', spec.label);
    STATS_PLV = ach_plv_plot_condition_summary( ...
                    PLV_list, cond_plv, spec.label);

    % STATS_PLV.p_anova is your ANOVA p across conditions
end
%%
%% ==============================================================
%  Full ACh + PLV analysis for mouse1 APP
%  Conditions: baseline, ambtemp, drugs
%  Episodes: 2 per condition (early / late)
%
%  PLV transitions (θ–ACh):
%     - Wake → NREM
%     - NREM → REM
%     - REM  → Wake
%
%  Extra PLV band (σ–ACh in NREM):
%     - Wake → NREM, EEG sigma [11–16 Hz]
% ==============================================================

clearvars; clear functions; close all; clc;

%% 1) Paths and common settings
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';

CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

%% 2) Define all episodes (segments) for this mouse

segs = struct([]);

% === BASELINE episodes ===
segs(1).mat_file   = fullfile(sigDir,  '20251001_baseline_mouse1_APP.mat');
segs(1).scores_csv = fullfile(scoreDir,'20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
segs(1).t_start    = 19600;
segs(1).t_end      = 20400;
segs(1).cond       = 'baseline';
segs(1).label      = 'BL-early';

segs(2)            = segs(1);
segs(2).t_start    = 19704;
segs(2).t_end      = 21000;
segs(2).label      = 'BL-late';

% === AMBTEMP episodes ===
segs(3).mat_file   = fullfile(sigDir,  '20251003_ambtemp_mouse1_APP.mat');
segs(3).scores_csv = fullfile(scoreDir,'20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv');
segs(3).t_start    = 10966;
segs(3).t_end      = 12770;
segs(3).cond       = 'ambtemp';
segs(3).label      = 'AMB-early';

segs(4)            = segs(3);
segs(4).t_start    = 16428;
segs(4).t_end      = 18163;
segs(4).label      = 'AMB-late';

% === DRUGS episodes ===
segs(5).mat_file   = fullfile(sigDir,  '20251022_drugs_mouse1_APP.mat');
segs(5).scores_csv = fullfile(scoreDir,'20251022_drugs_mouse1_APP_scored_scores_1Hz.csv');
segs(5).t_start    = 18766;
segs(5).t_end      = 20027;
segs(5).cond       = 'drugs';
segs(5).label      = 'DR-early';

segs(6)            = segs(5);
segs(6).t_start    = 22422;
segs(6).t_end      = 24576;
segs(6).label      = 'DR-late';

nSeg = numel(segs);

%% 3) Run ACh/EEG overview + onset features for each episode

OUT = cell(1,nSeg);

for k = 1:nSeg
    fprintf('\n=== Episode %d: %s (%s) ===\n', ...
            k, segs(k).label, segs(k).cond);

    tmp = plot_ach_eeg_segment( ...
              segs(k).mat_file, ...
              segs(k).scores_csv, ...
              segs(k).t_start, ...
              segs(k).t_end, ...
              'codes', CODES);

    % tag metadata
    tmp.cond  = segs(k).cond;
    tmp.label = segs(k).label;
    tmp.codes = CODES;

    OUT{k} = tmp;
end

% convenience arrays
episode_labels = {segs.label};
cond_labels    = {segs.cond};
OUTs           = [OUT{:}];     % struct array version

%% 4) Condition-level ACh onset comparison (WT vs conditions)

t_pre_peri  = 20;   % s before onset
t_post_peri = 40;   % s after onset
dt          = 0.2;  % common time step

[COND_onset, STATS_onset] = ach_compare_conditions_onsets( ...
                              OUTs, cond_labels, ...
                              t_pre_peri, t_post_peri, dt);

ach_plot_condition_summary_allStates(COND_onset, STATS_onset);
% → gives you Wake / NREM / REM onset plots + peak/slope bars with p values

%% 5) Bout-wise ACh features across whole bouts

[BOUT, STATS_BOUT] = ach_compute_bout_features_conditions(OUTs, cond_labels);
ach_plot_bout_features_allStates(BOUT, STATS_BOUT);
% → one figure per state: mean ACh, slope, amplitude, duration per condition

%% 6) PLV around multiple transitions, per episode

trans_specs = struct([]);

% θ-band transitions (EEG 5–10 Hz, ACh 0.05–0.5 Hz)
trans_specs(1).from     = 'Wake';
trans_specs(1).to       = 'NREM';
trans_specs(1).field    = 'Wake_NREM_theta';
trans_specs(1).label    = 'Wake→NREM θ–ACh PLV';
trans_specs(1).eeg_band = [5 10];
trans_specs(1).ach_band = [0.05 0.5];

trans_specs(2).from     = 'NREM';
trans_specs(2).to       = 'REM';
trans_specs(2).field    = 'NREM_REM_theta';
trans_specs(2).label    = 'NREM→REM θ–ACh PLV';
trans_specs(2).eeg_band = [5 10];
trans_specs(2).ach_band = [0.05 0.5];

trans_specs(3).from     = 'REM';
trans_specs(3).to       = 'Wake';
trans_specs(3).field    = 'REM_Wake_theta';
trans_specs(3).label    = 'REM→Wake θ–ACh PLV';
trans_specs(3).eeg_band = [5 10];
trans_specs(3).ach_band = [0.05 0.5];

% σ-band (NREM spindles-ish): Wake→NREM, EEG 11–16 Hz
trans_specs(4).from     = 'Wake';
trans_specs(4).to       = 'NREM';
trans_specs(4).field    = 'Wake_NREM_sigma';
trans_specs(4).label    = 'Wake→NREM σ–ACh PLV';
trans_specs(4).eeg_band = [11 16];
trans_specs(4).ach_band = [0.05 0.5];

% --- compute PLV per episode ---
for k = 1:nSeg
    fprintf('\n--- PLV for episode %d: %s (%s) ---\n', ...
            k, segs(k).label, segs(k).cond);

    for tt = 1:numel(trans_specs)
        spec = trans_specs(tt);

        PLV_tmp = ach_eeg_plv_periTransitions_segment( ...
                        segs(k).mat_file, ...
                        segs(k).scores_csv, ...
                        spec.from, spec.to, ...
                        segs(k).t_start, segs(k).t_end, ...
                        'codes', CODES, ...
                        'pre_sec', 20, 'post_sec', 40, ...
                        'eeg_band', spec.eeg_band, ...
                        'ach_band', spec.ach_band);

        OUT{k}.plv.(spec.field) = PLV_tmp;
        fprintf('  %s: %d transitions used.\n', ...
                spec.field, PLV_tmp.nEvents);
    end
end

OUTs = [OUT{:}];   % update struct array with .plv

%% 7) PLV condition-level visualisation: heatmaps + mean traces

for tt = 1:numel(trans_specs)
    spec = trans_specs(tt);

    % collect PLV structs per condition
    PLV_cells = {};
    cond_plv  = {};
    idx       = 0;

    for k = 1:nSeg
        if isfield(OUTs(k),'plv') && isfield(OUTs(k).plv, spec.field)
            P = OUTs(k).plv.(spec.field);
            if ~isempty(P) && isfield(P,'t_rel') && ~isempty(P.t_rel) ...
                    && isfield(P,'plv_mat') && ~isempty(P.plv_mat) ...
                    && isfield(P,'nEvents') && P.nEvents > 0

                idx = idx + 1;
                PLV_cells{idx} = P;             %#ok<SAGROW>
                cond_plv{idx}  = OUTs(k).cond;  %#ok<SAGROW>
            end
        end
    end

    if isempty(PLV_cells)
        fprintf('\n[PLV] No usable data for %s across episodes, skipping.\n', ...
                spec.label);
        continue;
    end

    PLV_list = [PLV_cells{:}];   % struct array

    % --- pretty PLV figure: event heatmaps + mean±SEM per condition ---
    ach_plv_heatmap_per_condition(PLV_list, cond_plv, spec.label);

    % If you still want numeric p-values across conditions, you can also call:
    % STATS_PLV = ach_plv_plot_condition_summary(PLV_list, cond_plv, spec.label);
end

%% 8) OPTIONAL: context-labelled bouts (NREM_preREM vs NREM_preWake, etc.)
% Example for baseline only – adapt as you like.

do_context = false;  % set true if you want to run this block

if do_context
    figure('Color','w','Position',[100 100 800 400]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    conds_unique = unique(cond_labels,'stable');

    for c = 1:numel(conds_unique)
        thisCond = conds_unique{c};
        idxSeg   = find(strcmp(cond_labels, thisCond));

        NREM_preREM  = 0;
        NREM_preWake = 0;
        Wake_postREM  = 0;
        Wake_postNREM = 0;

        for kk = idxSeg
            % reload scores for this segment
            Msc = readmatrix(segs(kk).scores_csv);
            if size(Msc,2) == 1
                score     = Msc(:,1);
                t_scores  = (0:numel(score)-1)';
            else
                t_scores  = Msc(:,1);
                score     = Msc(:,2);
            end

            bouts = build_context_bouts(score, t_scores, CODES);

            NREM_preREM  = NREM_preREM  + sum(strcmp({bouts.state},'NREM') & strcmp({bouts.next_state},'REM'));
            NREM_preWake = NREM_preWake + sum(strcmp({bouts.state},'NREM') & strcmp({bouts.next_state},'Wake'));
            Wake_postREM  = Wake_postREM  + sum(strcmp({bouts.state},'Wake') & strcmp({bouts.prev_state},'REM'));
            Wake_postNREM = Wake_postNREM + sum(strcmp({bouts.state},'Wake') & strcmp({bouts.prev_state},'NREM'));
        end

        % simple bar plot summary
        nexttile(1);
        hold on;
        bar((c-1)*4 + (1:4), ...
            [NREM_preREM NREM_preWake Wake_postREM Wake_postNREM]);
        xticks(1:4*numel(conds_unique));
        xticklabels({'NREM→REM','NREM→Wake','REM→Wake','NREM→Wake'});
        xtickangle(45);
        ylabel('# bouts');
        title('Context-labelled bouts (all episodes)');
        legend(conds_unique,'Location','northwest');
    end
end
