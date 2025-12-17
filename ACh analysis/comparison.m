%% ==============================================================
%  GROUP SUMMARY: APP vs WT – ALL KEY METRICS
%
%  Metrics:
%  --------------------------------------------------------------
%  1) Onset ACh metrics (Wake, NREM, REM onsets)
%       - Peak ACh (deltaF_all from ONSET features)
%       - Slope around onset (slope_all)
%       - Mean peri-onset ACh (mean of peri traces)
%
%  2) State ACh metrics (from bouts)
%       - Mean ACh per state (Wake, NREM, REM)
%       - NREM bout slope (mean of slope_all)
%       - NREM bout "variance" (variance of meanACh_all across bouts)
%
%  3) Coupling metrics
%       - NREM PLV (delta EEG phase × slow ACh phase)
%       - REM PAC ROI (theta phase × gamma amplitude)
%       - REM PAC COMODULOGRAMS (2 x 3 grid: genotype x condition)
%
%  Animals:
%      APP : mouse1
%      WT  : mouse4
%
%  Conditions:
%      baseline, ambtemp, drugs
%
%  NOTE:
%    - n = 1 per genotype → descriptive plots only, no p-values
%    - Colours:
%         APP baseline/ambtemp/drugs = light→dark red
%         WT  baseline/ambtemp/drugs = light→dark grey
% ==============================================================

clearvars; clear functions; close all; clc;

%% --------------------------------------------------------------
% 1) Paths + scoring codes
% --------------------------------------------------------------
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';

CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

conditions   = {'baseline','ambtemp','drugs'};
cond_labels  = {'BL','AMB','DR'};   % shorthand, only for text if needed
nCond        = numel(conditions);

state_names  = {'Wake','NREM','REM'};
nStates      = numel(state_names);

%% --------------------------------------------------------------
% 2) Define animals and windows  (APP mouse1 + WT mouse4)
% --------------------------------------------------------------
animals = struct([]);

% ---------- APP mouse ----------
animals(1).id        = 'mouse1_APP';
animals(1).genotype  = 'APP';

animals(1).baseline_mat    = fullfile(sigDir,  '20251001_baseline_mouse1_APP.mat');
animals(1).baseline_score  = fullfile(scoreDir,'20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
animals(1).baseline_window = [19704, 21000];

animals(1).ambtemp_mat     = fullfile(sigDir,  '20251003_ambtemp_mouse1_APP.mat');
animals(1).ambtemp_score   = fullfile(scoreDir,'20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv');
animals(1).ambtemp_window  = [16428, 18163];

animals(1).drugs_mat       = fullfile(sigDir,  '20251022_drugs_mouse1_APP.mat');
animals(1).drugs_score     = fullfile(scoreDir,'20251022_drugs_mouse1_APP_scored_scores_1Hz.csv');
animals(1).drugs_window    = [18766, 20027];

% ---------- WT mouse ----------
animals(2).id        = 'mouse4_WT';
animals(2).genotype  = 'WT';

animals(2).baseline_mat    = fullfile(sigDir,  '20251006_baseline_mouse4_WT.mat');
animals(2).baseline_score  = fullfile(scoreDir,'20251006_baseline_mouse4_WT_scored_scores_1Hz.csv');
animals(2).baseline_window = [14821, 17023];

animals(2).ambtemp_mat     = fullfile(sigDir,  '20251007_ambtemp_mouse4_WT.mat');
animals(2).ambtemp_score   = fullfile(scoreDir,'20251007_ambtemp_mouse4_WT_scored_scores_1Hz.csv');
animals(2).ambtemp_window  = [16015, 17831];

animals(2).drugs_mat       = fullfile(sigDir,  '20251022_drugs_mouse4_WT.mat');
animals(2).drugs_score     = fullfile(scoreDir,'20251022_drugs_mouse4_WT_scored_scores_1Hz.csv');
animals(2).drugs_window    = [15548, 17206];

nAnimals    = numel(animals);
geno_labels = {animals.genotype};
isAPP       = strcmp(geno_labels,'APP');
isWT        = strcmp(geno_labels,'WT');

%% --------------------------------------------------------------
% 3) Colour scheme (consistent for ALL plots)
% --------------------------------------------------------------
% condition order: {'baseline','ambtemp','drugs'}

% APP: light -> darker red
colors.APP = [
    1.00 0.60 0.60;  % baseline
    0.90 0.30 0.30;  % ambtemp
    0.60 0.00 0.00]; % drugs

% WT: light -> darker grey
colors.WT = [
    0.85 0.85 0.85;  % baseline
    0.55 0.55 0.55;  % ambtemp
    0.25 0.25 0.25]; % drugs

%% --------------------------------------------------------------
% 4) Preallocate metric matrices
% --------------------------------------------------------------
% Onset metrics: [animal x state x condition]
onset_peak  = nan(nAnimals, nStates, nCond);  % from deltaF_all
onset_slope = nan(nAnimals, nStates, nCond);  % from slope_all
onset_mean  = nan(nAnimals, nStates, nCond);  % mean peri-onset ACh

% State ACh (mean of bout means) : [animal x state x condition]
state_mean_ach = nan(nAnimals, nStates, nCond);

% NREM bout stability: [animal x condition]
nrem_bout_slope = nan(nAnimals, nCond);   % mean slope across NREM bouts
nrem_bout_var   = nan(nAnimals, nCond);   % variance of meanACh_all across NREM bouts

% Coupling metrics: [animal x condition]
PLV_NREM_mean = nan(nAnimals, nCond);     % NREM PLV (delta–slow ACh)
PAC_REM_roi   = nan(nAnimals, nCond);     % REM PAC ROI (theta–gamma)

% Store REM comodulograms for plotting later
COM_REM(nAnimals, nCond) = struct('MI',[],'phase_freqs',[],'amp_freqs',[]);

%% --------------------------------------------------------------
% 5) Parameters for peri-onset, PLV and PAC
% --------------------------------------------------------------
% Peri-onset window
t_pre_peri  = 20;   % sec before onset
t_post_peri = 40;   % sec after onset
dt_peri     = 0.2;  % step for common time axis

% NREM PLV: delta EEG phase (1–4 Hz) vs slow ACh phase (0.1–1 Hz)
eeg_band_NREM = [1 4];
ach_band_NREM = [0.1 1.0];
min_bout_sec  = 10;

% REM PAC: theta phase (6–10 Hz) × gamma amplitude (30–80 Hz)
phase_freqs_rem = 6:0.5:10;   % theta phase
amp_freqs_rem   = 30:2:80;    % gamma amplitude

% ROI for REM PAC summary
phase_roi_rem = [6 9];        % 6–9 Hz theta phase
amp_roi_rem   = [40 80];      % 40–80 Hz gamma

%% --------------------------------------------------------------
% 6) MAIN LOOP: per animal – build OUTs once, then extract metrics
% --------------------------------------------------------------
fprintf('=== Running onset, state, PLV and REM PAC metrics ===\n');

for a = 1:nAnimals
    fprintf('\n==============================================\n');
    fprintf('Animal %d/%d: %s (%s)\n', ...
        a, nAnimals, animals(a).id, animals(a).genotype);
    fprintf('==============================================\n');

    % --------- build per-animal segs (baseline / ambtemp / drugs) ---------
    segs = struct([]);
    for c = 1:nCond
        cond = conditions{c};
        segs(c).mat_file   = animals(a).([cond '_mat']);
        segs(c).scores_csv = animals(a).([cond '_score']);
        win                = animals(a).([cond '_window']);
        segs(c).t_start    = win(1);
        segs(c).t_end      = win(2);
        segs(c).cond       = cond;
        segs(c).label      = cond_labels{c};
    end

    % --------- make OUTs by calling plot_ach_eeg_segment ---------
    OUT = cell(1, nCond);
    for c = 1:nCond
        fprintf('  Episode %d: %s (%s)\n', c, segs(c).label, segs(c).cond);

        tmp = plot_ach_eeg_segment( ...
                  segs(c).mat_file, ...
                  segs(c).scores_csv, ...
                  segs(c).t_start, ...
                  segs(c).t_end, ...
                  'codes', CODES);

        tmp.cond  = segs(c).cond;
        tmp.label = segs(c).label;
        tmp.codes = CODES;

        OUT{c} = tmp;
    end

    OUTs              = [OUT{:}];          % 1×3 struct
    cond_labels_local = conditions;       % {'baseline','ambtemp','drugs'}

    % ==========================================================
    % 6A) Onset-based metrics (all states, all conditions)
    %      using ach_compare_conditions_onsets
    % ==========================================================
    [COND_onset, ~] = ach_compare_conditions_onsets( ...
                         OUTs, cond_labels_local, ...
                         t_pre_peri, t_post_peri, dt_peri);

    for c = 1:nCond
        for s = 1:nStates
            st_struct = COND_onset(c).state(s);

            % Peak ACh (deltaF_all)
            if isfield(st_struct, 'deltaF_all') && ~isempty(st_struct.deltaF_all)
                onset_peak(a,s,c) = mean(st_struct.deltaF_all, 'omitnan');
            end

            % Slope (slope_all)
            if isfield(st_struct, 'slope_all') && ~isempty(st_struct.slope_all)
                onset_slope(a,s,c) = mean(st_struct.slope_all, 'omitnan');
            end

            % Mean peri-onset ACh (mean over peri traces)
            if isfield(st_struct, 'traces') && ~isempty(st_struct.traces)
                evt_means = mean(st_struct.traces, 2, 'omitnan'); % one mean per event
                onset_mean(a,s,c) = mean(evt_means, 'omitnan');
            end
        end
    end

    % ==========================================================
    % 6B) State/bout metrics using ach_compute_bout_features_conditions
    % ==========================================================
    [BOUT, ~] = ach_compute_bout_features_conditions( ...
                     OUTs, cond_labels_local, ...
                     'plv_eeg_band', [1 4], ...
                     'plv_ach_band', [0.1 1.0], ...
                     'minBoutDur',   5);

    for c = 1:nCond
        % Mean ACh per state
        for s = 1:nStates
            st_struct = BOUT(c).state(s);
            if ~isempty(st_struct.meanACh_all)
                state_mean_ach(a,s,c) = mean(st_struct.meanACh_all, 'omitnan');
            end
        end

        % NREM bout stability metrics (state index 2 = NREM)
        nrem_struct = BOUT(c).state(2);

        if ~isempty(nrem_struct.slope_all)
            nrem_bout_slope(a,c) = mean(nrem_struct.slope_all, 'omitnan');
        end

        if ~isempty(nrem_struct.meanACh_all)
            nrem_bout_var(a,c) = var(nrem_struct.meanACh_all, 'omitnan');
        end
    end

    % ==========================================================
    % 6C) Coupling metrics per condition (NREM PLV + REM PAC)
    % ==========================================================
    for c = 1:nCond
        cond       = conditions{c};
        mat_file   = animals(a).([cond '_mat']);
        scores_csv = animals(a).([cond '_score']);
        win        = animals(a).([cond '_window']);
        t_start    = win(1);
        t_end      = win(2);

        if ~exist(mat_file,'file') || ~exist(scores_csv,'file')
            warning('    Missing files for %s %s, skipping coupling metrics.', animals(a).id, cond);
            continue;
        end

        % ---- NREM PLV (delta EEG – slow ACh) ----
        try
            PLV = ach_eeg_plv_state_bouts( ...
                mat_file, scores_csv, 'NREM', ...
                t_start, t_end, ...
                'codes',        CODES, ...
                'eeg_band',     eeg_band_NREM, ...
                'ach_band',     ach_band_NREM, ...
                'min_bout_sec', min_bout_sec);

            PLV_NREM_mean(a,c) = mean(PLV.plv_all, 'omitnan');
            fprintf('    [%s] NREM PLV mean = %.4f\n', cond, PLV_NREM_mean(a,c));
        catch ME
            warning('    NREM PLV failed for %s %s: %s', ...
                animals(a).id, cond, ME.message);
        end

        % ---- REM PAC (theta–gamma ROI) ----
        try
            COMr = ach_pac_comod_rem_segment( ...
                       mat_file, scores_csv, ...
                       t_start, t_end, CODES, ...
                       'phase_freqs', phase_freqs_rem, ...
                       'amp_freqs',   amp_freqs_rem, ...
                       'nbins',       18, ...
                       'doPlot',      false);

            if ~isempty(COMr) && ~isempty(COMr.MI)
                pf = COMr.phase_freqs;
                af = COMr.amp_freqs;

                ip = pf >= phase_roi_rem(1) & pf <= phase_roi_rem(2);
                ia = af >= amp_roi_rem(1)   & af <= amp_roi_rem(2);

                roi_vals         = COMr.MI(ip, ia);
                PAC_REM_roi(a,c) = mean(roi_vals(:), 'omitnan');

                fprintf('    [%s] REM PAC ROI MI = %.4g\n', cond, PAC_REM_roi(a,c));

                COM_REM(a,c).MI          = COMr.MI;
                COM_REM(a,c).phase_freqs = COMr.phase_freqs;
                COM_REM(a,c).amp_freqs   = COMr.amp_freqs;
            else
                fprintf('    [%s] REM PAC: no MI data.\n', cond);
            end
        catch ME
            warning('    REM PAC failed for %s %s: %s', ...
                animals(a).id, cond, ME.message);
        end
    end
end

fprintf('\n=== All metric matrices built ===\n');

%% --------------------------------------------------------------
% 7) Split by genotype (APP vs WT)
% --------------------------------------------------------------
% Onset metrics
onset_peak_APP  = onset_peak(isAPP,:,:);
onset_peak_WT   = onset_peak(isWT,:,:);
onset_slope_APP = onset_slope(isAPP,:,:);
onset_slope_WT  = onset_slope(isWT,:,:);
onset_mean_APP  = onset_mean(isAPP,:,:);
onset_mean_WT   = onset_mean(isWT,:,:);

% State mean ACh
state_mean_APP = state_mean_ach(isAPP,:,:);
state_mean_WT  = state_mean_ach(isWT,:,:);

% NREM bout stability
nrem_slope_APP = nrem_bout_slope(isAPP,:);
nrem_slope_WT  = nrem_bout_slope(isWT,:);
nrem_var_APP   = nrem_bout_var(isAPP,:);
nrem_var_WT    = nrem_bout_var(isWT,:);

% Coupling
PLV_APP = PLV_NREM_mean(isAPP,:);
PLV_WT  = PLV_NREM_mean(isWT,:);

PAC_APP = PAC_REM_roi(isAPP,:);
PAC_WT  = PAC_REM_roi(isWT,:);

%% --------------------------------------------------------------
% 8) PLOTTING – onset metrics (peak, slope, mean) per state
%     Each figure: x = condition, bars = APP vs WT
% --------------------------------------------------------------
for s = 1:nStates
    st_label = state_names{s};

    % Peak ACh
    M_APP = squeeze(onset_peak_APP(:,s,:));   % [nAPP x nCond]
    M_WT  = squeeze(onset_peak_WT(:,s,:));    % [nWT  x nCond]
    plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ...
        sprintf('ACh peak around %s onset (\\DeltaF/F)', st_label), ...
        sprintf('%s onset: ACh peak (preliminary)', st_label));

    % Slope
    M_APP = squeeze(onset_slope_APP(:,s,:));
    M_WT  = squeeze(onset_slope_WT(:,s,:));
    plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ...
        sprintf('ACh slope around %s onset (\\DeltaF/F per s)', st_label), ...
        sprintf('%s onset: ACh slope (preliminary)', st_label));

    % Mean peri-onset ACh
    M_APP = squeeze(onset_mean_APP(:,s,:));
    M_WT  = squeeze(onset_mean_WT(:,s,:));
    plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ...
        sprintf('Mean ACh around %s onset (\\DeltaF/F)', st_label), ...
        sprintf('%s onset: mean ACh (preliminary)', st_label));
end

%% --------------------------------------------------------------
% 9) PLOTTING – state ACh level & NREM stability
% --------------------------------------------------------------
for s = 1:nStates
    st_label = state_names{s};
    M_APP = squeeze(state_mean_APP(:,s,:));
    M_WT  = squeeze(state_mean_WT(:,s,:));

    plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ...
        sprintf('Mean ACh in %s bouts (\\DeltaF/F)', st_label), ...
        sprintf('%s mean ACh per bout (preliminary)', st_label));
end

% NREM bout slope
plot_APP_WT_metric(nrem_slope_APP, nrem_slope_WT, conditions, colors, ...
    'NREM bout ACh slope (\\DeltaF/F per s)', ...
    'NREM bout ACh slope (stability, preliminary)');

% NREM bout variance (across meanACh of NREM bouts)
plot_APP_WT_metric(nrem_var_APP, nrem_var_WT, conditions, colors, ...
    'NREM bout ACh variance (across bouts)', ...
    'NREM bout ACh variance (stability, preliminary)');

%% --------------------------------------------------------------
% 10) PLOTTING – coupling metrics (NREM PLV + REM PAC)
% --------------------------------------------------------------
% NREM PLV (delta–slow ACh)
plot_APP_WT_metric(PLV_APP, PLV_WT, conditions, colors, ...
    'Mean NREM PLV (delta EEG – slow ACh)', ...
    'NREM PLV: APP vs WT (preliminary)');

% REM PAC ROI (theta–gamma)
plot_APP_WT_metric(PAC_APP, PAC_WT, conditions, colors, ...
    'REM theta–gamma PAC (MI, ROI)', ...
    'REM PAC ROI: APP vs WT (preliminary)');

%% --------------------------------------------------------------
% 11) PLOTTING – REM PAC comodulograms (genotype x condition)
% --------------------------------------------------------------
plot_REM_PAC_comodulograms(COM_REM, animals, conditions);

fprintf('\n=== All plotting complete ===\n');

%% ==============================================================
% 12) Helper plotting function (bars + dots, condition × genotype)
% ==============================================================

function plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ylab, ttl)
% M_APP, M_WT: [nAnimals x nCond] matrices (can be 1 x nCond)
% Plots bars grouped by condition, with APP vs WT per condition.
% Colours:
%   APP: colors.APP (nCond x 3)
%   WT : colors.WT  (nCond x 3)

    if isempty(M_APP) || isempty(M_WT)
        warning('plot_APP_WT_metric: empty APP or WT matrix, skipping plot (%s).', ttl);
        return;
    end

    % Make sure both are 2D (nAnimals x nCond)
    M_APP = squeeze(M_APP);
    M_WT  = squeeze(M_WT);
    if isvector(M_APP), M_APP = reshape(M_APP, 1, []); end
    if isvector(M_WT),  M_WT  = reshape(M_WT,  1, []); end

    % Check same number of conditions
    nCond_app = size(M_APP, 2);
    nCond_wt  = size(M_WT,  2);
    if nCond_app ~= nCond_wt
        error('plot_APP_WT_metric: APP and WT have different number of columns (%d vs %d).', ...
              nCond_app, nCond_wt);
    end
    nCond = nCond_app;

    % Trim condition labels if needed
    cond_plot = conditions(1:nCond);

    % Sanity check for colour matrices
    if size(colors.APP,1) < nCond || size(colors.WT,1) < nCond
        error('plot_APP_WT_metric: not enough colour rows for %d conditions.', nCond);
    end

    % Means and SEM across animals
    mean_APP = mean(M_APP, 1, 'omitnan');
    mean_WT  = mean(M_WT,  1, 'omitnan');

    sem_APP = std(M_APP, 0, 1, 'omitnan') ./ max(1, sqrt(sum(~isnan(M_APP),1)));
    sem_WT  = std(M_WT,  0, 1, 'omitnan') ./ max(1, sqrt(sum(~isnan(M_WT),1)));

    % X positions: group per condition
    x_base = 1:nCond;
    width  = 0.35;
    x_APP  = x_base - width/2;
    x_WT   = x_base + width/2;

    figure('Color','w','Position',[200 200 600 450]); 
    hold on;

    % Dummy handles for legend (we'll overwrite on first iteration)
    hAPP = [];
    hWT  = [];

    % Loop over conditions and draw APP/WT bars manually
    for c = 1:nCond
        % --- APP bar ---
        if c == 1
            hAPP = bar(x_APP(c), mean_APP(c), width, ...
                'FaceColor', colors.APP(c,:), ...
                'EdgeColor', 'k');
        else
            bar(x_APP(c), mean_APP(c), width, ...
                'FaceColor', colors.APP(c,:), ...
                'EdgeColor', 'k');
        end

        % --- WT bar ---
        if c == 1
            hWT = bar(x_WT(c), mean_WT(c), width, ...
                'FaceColor', colors.WT(c,:), ...
                'EdgeColor', 'k');
        else
            bar(x_WT(c), mean_WT(c), width, ...
                'FaceColor', colors.WT(c,:), ...
                'EdgeColor', 'k');
        end

        % --- Error bars ---
        errorbar(x_APP(c), mean_APP(c), sem_APP(c), 'k', ...
                 'LineStyle','none', 'LineWidth',1);
        errorbar(x_WT(c),  mean_WT(c),  sem_WT(c),  'k', ...
                 'LineStyle','none', 'LineWidth',1);

        % --- Scatter individual animals ---
        jitter = 0.05;
        for i = 1:size(M_APP,1)
            if ~isnan(M_APP(i,c))
                plot(x_APP(c) + (rand-0.5)*2*jitter, M_APP(i,c), 'o', ...
                     'MarkerFaceColor', colors.APP(c,:), ...
                     'MarkerEdgeColor', 'k');
            end
        end
        for i = 1:size(M_WT,1)
            if ~isnan(M_WT(i,c))
                plot(x_WT(c) + (rand-0.5)*2*jitter, M_WT(i,c), 'o', ...
                     'MarkerFaceColor', colors.WT(c,:), ...
                     'MarkerEdgeColor', 'k');
            end
        end
    end

    set(gca,'XTick',x_base,'XTickLabel',cond_plot,'FontSize',12);
    xlim([0.5, nCond+0.5]);
    ylabel(ylab,'FontSize',13);
    if ~isempty(hAPP) && ~isempty(hWT)
        legend([hAPP hWT], {'APP','WT'},'Location','best');
    end
    title(ttl,'FontSize',14);
    grid on; box on;
end

%% ==============================================================
% 13) Helper: REM PAC comodulograms (2 x 3 grid: genotype x condition)
% ==============================================================

function plot_REM_PAC_comodulograms(COM_REM, animals, conditions)
% COM_REM: [nAnimals x nCond] struct with fields .MI, .phase_freqs, .amp_freqs
% Figure: rows = genotype (one row per animal here), cols = conditions

    nAnimals = numel(animals);
    nCond    = numel(conditions);

    % Find global max MI for consistent colour scale, ignoring NaNs
    mi_max = NaN;
    for a = 1:nAnimals
        for c = 1:nCond
            M = COM_REM(a,c).MI;
            if ~isempty(M)
                val = M(:);
                val = val(~isnan(val));
                if ~isempty(val)
                    this_max = max(val);
                    if isnan(mi_max)
                        mi_max = this_max;
                    else
                        mi_max = max(mi_max, this_max);
                    end
                end
            end
        end
    end

    if isnan(mi_max)
        warning('plot_REM_PAC_comodulograms: PAC MI is all NaN for all animals/conditions.');
        % You can still make an empty figure if you want:
        figure('Color','w'); title('REM PAC: all MI values are NaN');
        return;
    end

    figure('Color','w','Position',[100 100 1200 500]);
    tlo = tiledlayout(nAnimals, nCond, 'TileSpacing','compact', 'Padding','compact');
    title(tlo, 'REM PAC Comodulograms (theta phase × gamma amplitude)', ...
          'FontSize', 16, 'FontWeight', 'bold');

    for a = 1:nAnimals
        for c = 1:nCond
            nexttile;

            COM = COM_REM(a,c);
            if isempty(COM.MI)
                text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
                axis off;
                continue;
            end

            M = COM.MI;
            imagesc(COM.amp_freqs, COM.phase_freqs, M);
            set(gca,'YDir','normal');
            colormap(turbo);
            caxis([0 mi_max]);  % NaNs will appear as missing / blank

            if a == nAnimals
                xlabel('Amplitude (Hz)');
            end
            if c == 1
                ylabel(sprintf('%s phase (Hz)', animals(a).genotype));
            end
            title(sprintf('%s', conditions{c}), 'FontSize', 11);
        end
    end

    cb = colorbar('Location','eastoutside');
    cb.Label.String = 'Modulation Index (MI)';
end
