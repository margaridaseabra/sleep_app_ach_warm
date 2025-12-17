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
segs(1).mat_file   = fullfile(sigDir,  '20251006_baseline_mouse4_WT.mat');
segs(1).scores_csv = fullfile(scoreDir,'20251006_baseline_mouse4_WT_scored_scores_1Hz.csv');
segs(1).t_start    = 14821;      % BL-late
segs(1).t_end      = 17023;
segs(1).cond       = 'baseline';
segs(1).label      = 'BL';

% === AMBTEMP: use AMB-early ===
segs(2).mat_file   = fullfile(sigDir,  '20251007_ambtemp_mouse4_WT.mat');
segs(2).scores_csv = fullfile(scoreDir,'20251007_ambtemp_mouse4_WT_scored_scores_1Hz.csv');
segs(2).t_start    = 16015;      % AMB-early
segs(2).t_end      = 17831;
segs(2).cond       = 'ambtemp';
segs(2).label      = 'AMB';
 
% === DRUGS: use DR-late ===
segs(3).mat_file   = fullfile(sigDir,  '20251022_drugs_mouse4_WT.mat');
segs(3).scores_csv = fullfile(scoreDir,'20251022_drugs_mouse4_WT_scored_scores_1Hz.csv');
segs(3).t_start    = 15548;      % DR-late
segs(3).t_end      = 17206;
segs(3).cond       = 'drugs';
segs(3).label      = 'DR';

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

%% --------------------------------------------------------------
% 7) NREM PAC comodulogram per condition (slow-phase × sigma-amp)
%    WITH DIAGNOSTICS
% --------------------------------------------------------------
phase_freqs = 0.5:0.5:4;    % slow/delta phase
amp_freqs   = 7:1:25;       % amplitude around sigma

COM_cells = cell(1, nSeg);

fprintf('\n=== PAC COMODULOGRAM ANALYSIS ===\n');

for k = 1:nSeg
    fprintf('\nComputing PAC for %s (%s)...\n', segs(k).label, segs(k).cond);
    
    COM_cells{k} = ach_pac_comod_nrem_segment( ...
                       segs(k).mat_file, ...
                       segs(k).scores_csv, ...
                       segs(k).t_start, ...
                       segs(k).t_end, ...
                       CODES, ...
                       'phase_freqs', phase_freqs, ...
                       'amp_freqs',   amp_freqs, ...
                       'nbins',       18, ...
                       'doPlot',      true, ...
                       'label',       sprintf('%s (%s)', ...
                                              segs(k).label, ...
                                              segs(k).cond));
    
    % Check results
    if isempty(COM_cells{k}) || isempty(COM_cells{k}.MI)
        warning('No PAC data for %s', segs(k).label);
    else
        fprintf('  PAC matrix: %d × %d\n', size(COM_cells{k}.MI));
        fprintf('  MI range: [%.4f, %.4f]\n', ...
                min(COM_cells{k}.MI(:)), max(COM_cells{k}.MI(:)));
        fprintf('  NREM samples used: %d (%.1f sec)\n', ...
                COM_cells{k}.nSamples, COM_cells{k}.nSamples/COM_cells{k}.fs);
    end
end

COM = [COM_cells{:}];

fprintf('\n=== PAC COMPUTATION COMPLETE ===\n');
fprintf('Generated %d comodulograms\n', numel(COM));
%% --------------------------------------------------------------
% 8) NREM PAC ROI summary: SO-phase (0.5–1.5 Hz) × sigma-amp (11–16 Hz)
%    ENHANCED VERSION
% --------------------------------------------------------------
phase_roi = [0.5 1.5];   % slow oscillation
amp_roi   = [11 16];     % spindles/sigma

mi_roi = nan(1, nSeg);
mi_full = cell(1, nSeg);  % store full MI matrices for visualization

fprintf('\n=== PAC ROI ANALYSIS ===\n');
fprintf('ROI: Phase %.1f-%.1f Hz × Amplitude %d-%d Hz\n', ...
        phase_roi(1), phase_roi(2), amp_roi(1), amp_roi(2));

for k = 1:nSeg
    if k > numel(COM) || isempty(COM(k).MI)
        fprintf('%s: No PAC data\n', cond_labels{k});
        continue;
    end
    
    pf = COM(k).phase_freqs;
    af = COM(k).amp_freqs;
    
    % Find ROI indices
    ip = pf >= phase_roi(1) & pf <= phase_roi(2);
    ia = af >= amp_roi(1)   & af <= amp_roi(2);
    
    % Extract ROI values
    roi_vals = COM(k).MI(ip, ia);
    mi_roi(k) = mean(roi_vals(:), 'omitnan');
    mi_full{k} = COM(k).MI;
    
    fprintf('%s: Mean SO-spindle MI = %.4f\n', cond_labels{k}, mi_roi(k));
    fprintf('  Phase freqs in ROI: %.1f-%.1f Hz (%d bins)\n', ...
            min(pf(ip)), max(pf(ip)), sum(ip));
    fprintf('  Amp freqs in ROI: %d-%d Hz (%d bins)\n', ...
            min(af(ia)), max(af(ia)), sum(ia));
end

% ============ FIGURE 1: ROI Bar Plot ============
figure('Color','w','Position',[100 100 600 450]);

% Bar plot with individual data points
x = 1:nSeg;
bar(x, mi_roi, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','k', 'LineWidth',1.5);
hold on;

% Add error bars or SEM if you have multiple episodes per condition
% (currently you have 1 episode per condition, so no error bars)

% Formatting
set(gca, 'XTick', x, 'XTickLabel', cond_labels, 'FontSize', 12);
ylabel('Modulation Index (MI)', 'FontSize', 14);
title(sprintf('NREM SO-Spindle PAC\n(Phase: %.1f-%.1f Hz, Amp: %d-%d Hz)', ...
              phase_roi(1), phase_roi(2), amp_roi(1), amp_roi(2)), ...
      'FontSize', 14, 'FontWeight', 'bold');
ylim([0 max(mi_roi)*1.2]);
box on; grid on;

% Add value labels on bars
for k = 1:nSeg
    if ~isnan(mi_roi(k))
        text(k, mi_roi(k) + max(mi_roi)*0.03, sprintf('%.3f', mi_roi(k)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
end

% ============ FIGURE 2: All Comodulograms Side-by-Side ============
figure('Color','w','Position',[150 150 1400 400]);
tiledlayout(1, nSeg, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:nSeg
    nexttile; 
    
    if isempty(mi_full{k})
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
        title(cond_labels{k});
        continue;
    end
    
    % Plot comodulogram
    imagesc(amp_freqs, phase_freqs, mi_full{k});
    axis xy;
    
    % Mark ROI with rectangle
    hold on;
    rectangle('Position', [amp_roi(1), phase_roi(1), ...
                          diff(amp_roi), diff(phase_roi)], ...
             'EdgeColor', 'w', 'LineWidth', 2.5, 'LineStyle', '--');
    
    xlabel('Amplitude (Hz)', 'FontSize', 11);
    if k == 1
        ylabel('Phase (Hz)', 'FontSize', 11);
    end
    title(cond_labels{k}, 'FontSize', 13, 'FontWeight', 'bold');
    
    colormap(turbo);
    c = colorbar;
    if k == nSeg
        ylabel(c, 'MI', 'FontSize', 11);
    end
    
    % Same color scale for all
    clim([0 max(cellfun(@(x) max(x(:)), mi_full(~cellfun(@isempty, mi_full))))]);
end

sgtitle('NREM PAC Comodulograms (White box = SO-Spindle ROI)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% ============ FIGURE 3: ROI Comparison with Statistical Test ============
if nSeg > 2
    figure('Color','w','Position',[200 200 500 400]);
    
    % Box plot or violin plot would be better with multiple bouts
    % but with 1 value per condition, just show bars with pattern
    
    bar(1:nSeg, mi_roi, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');
    hold on;
    
    % Add connecting lines to show trend
    plot(1:nSeg, mi_roi, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r');
    
    set(gca, 'XTick', 1:nSeg, 'XTickLabel', cond_labels);
    ylabel('SO-Spindle PAC (MI)', 'FontSize', 13);
    title('PAC Strength Across Conditions', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; box on;
    
    % Note about statistics
    text(0.5, max(mi_roi)*0.95, ...
         'Note: Statistical test requires multiple episodes per condition', ...
         'Units', 'normalized', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
end

fprintf('\n=== PAC ROI SUMMARY COMPLETE ===\n');