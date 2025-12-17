%% ==============================================================
%  COMPREHENSIVE SLEEP ANALYSIS: ALL ANIMALS, NREM + REM
%  Includes: PLV (NREM & REM) + PAC (NREM & REM)
%  Genotypes: APP and WT
% ==============================================================

clearvars; clear functions; close all; clc;

% --------------------------------------------------------------
% 1) Define all animals and their files
% --------------------------------------------------------------
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';

CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

% Define your animals
animals = struct([]);

% ========== APP MICE ==========
% Mouse 1 - APP
animals(1).id = 'mouse1_APP';
animals(1).genotype = 'APP';
animals(1).baseline_mat   = fullfile(sigDir, '20251001_baseline_mouse1_APP.mat');
animals(1).baseline_score = fullfile(scoreDir, '20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
animals(1).baseline_window = [19704, 21000];  % [t_start, t_end]

animals(1).ambtemp_mat   = fullfile(sigDir, '20251003_ambtemp_mouse1_APP.mat');
animals(1).ambtemp_score = fullfile(scoreDir, '20251003_ambtemp_mouse1_APP_scored_scores_1Hz.csv');
animals(1).ambtemp_window = [10966, 12770];

animals(1).drugs_mat   = fullfile(sigDir, '20251022_drugs_mouse1_APP.mat');
animals(1).drugs_score = fullfile(scoreDir, '20251022_drugs_mouse1_APP_scored_scores_1Hz.csv');
animals(1).drugs_window = [22422, 24576];

% ========== WT MICE ==========
% Mouse 2 - WT (add your WT mouse info)
animals(2).id = 'mouse2_WT';
animals(2).genotype = 'WT';
animals(2).baseline_mat   = fullfile(sigDir, '20251006_baseline_mouse4_WT.mat');
animals(2).baseline_score = fullfile(scoreDir, '20251006_baseline_mouse4_WT_scored_scores_1Hz.csv');
animals(2).baseline_window = [14821, 17023];  % Adjust to your data

animals(2).ambtemp_mat   = fullfile(sigDir, '20251007_ambtemp_mouse4_WT.mat');
animals(2).ambtemp_score = fullfile(scoreDir, '20251007_ambtemp_mouse4_WT_scored_scores_1Hz.csv');
animals(2).ambtemp_window = [16015, 17831];

animals(2).drugs_mat   = fullfile(sigDir, '20251022_drugs_mouse4_WT.mat');
animals(2).drugs_score = fullfile(scoreDir, '20251022_drugs_mouse4_WT_scored_scores_1Hz.csv');
animals(2).drugs_window = [15548, 17206];


nAnimals = numel(animals);
conditions = {'baseline', 'ambtemp', 'drugs'};
nCond = numel(conditions);

fprintf('=== ANALYZING %d ANIMALS ===\n', nAnimals);
for i = 1:nAnimals
    fprintf('  %d. %s (%s)\n', i, animals(i).id, animals(i).genotype);
end
fprintf('\n');

% --------------------------------------------------------------
% 2) Initialize storage for all results
% --------------------------------------------------------------
ALL_RESULTS = struct();

for a = 1:nAnimals
    animal_id = animals(a).id;
    ALL_RESULTS.(animal_id).genotype = animals(a).genotype;
    
    for c = 1:nCond
        cond = conditions{c};
        ALL_RESULTS.(animal_id).(cond) = struct();
    end
end

% --------------------------------------------------------------
% 3) MAIN ANALYSIS LOOP: Each animal × each condition
% --------------------------------------------------------------

for a = 1:nAnimals
    animal = animals(a);
    animal_id = animal.id;
    
    fprintf('\n########################################\n');
    fprintf('ANIMAL %d/%d: %s (%s)\n', a, nAnimals, animal_id, animal.genotype);
    fprintf('########################################\n\n');
    
    for c = 1:nCond
        cond = conditions{c};
        
        fprintf('--- Condition: %s ---\n', upper(cond));
        
        % Get file paths and window for this condition
        mat_file   = animal.([cond '_mat']);
        scores_csv = animal.([cond '_score']);
        time_window = animal.([cond '_window']);
        t_start = time_window(1);
        t_end   = time_window(2);
        
        % Check if files exist
        if ~exist(mat_file, 'file')
            warning('File not found: %s', mat_file);
            continue;
        end
        if ~exist(scores_csv, 'file')
            warning('File not found: %s', scores_csv);
            continue;
        end
        
        % ==============================================
        % A) NREM PLV (Delta EEG × Delta ACh)
        % ==============================================
        fprintf('  Computing NREM PLV...\n');
        try
            PLV_NREM = ach_eeg_plv_state_bouts( ...
                mat_file, scores_csv, 'NREM', t_start, t_end, ...
                'codes', CODES, ...
                'eeg_band', [1 4], ...      % Delta
                'ach_band', [1 4], ...      % Delta (matched for stability)
                'min_bout_sec', 5);
            
            % Store results
            ALL_RESULTS.(animal_id).(cond).NREM_PLV = PLV_NREM;
            
            % Summary stats
            n_bouts = numel(PLV_NREM.plv_all);
            valid_plv = sum(~isnan(PLV_NREM.plv_all));
            mean_plv = mean(PLV_NREM.plv_all, 'omitnan');
            
            fprintf('    Found %d NREM bouts, %d valid PLV (%.1f%%)\n', ...
                    n_bouts, valid_plv, 100*valid_plv/n_bouts);
            fprintf('    Mean PLV: %.4f\n', mean_plv);
        catch ME
            warning('%s', ['NREM PLV failed: ' ME.message]);
            ALL_RESULTS.(animal_id).(cond).NREM_PLV = struct('plv_all',[],'dur_all',[],'t0_all',[]);
        end
        
        % ==============================================
        % B) REM PLV (Theta EEG × Delta ACh)
        % ==============================================
        fprintf('  Computing REM PLV...\n');
        try
            PLV_REM = ach_eeg_plv_state_bouts( ...
                mat_file, scores_csv, 'REM', t_start, t_end, ...
                'codes', CODES, ...
                'eeg_band', [6 10], ...     % Theta for REM
                'ach_band', [1 4], ...      % Slow modulation
                'min_bout_sec', 5);
            
            ALL_RESULTS.(animal_id).(cond).REM_PLV = PLV_REM;
            
            n_bouts = numel(PLV_REM.plv_all);
            valid_plv = sum(~isnan(PLV_REM.plv_all));
            mean_plv = mean(PLV_REM.plv_all, 'omitnan');
            
            fprintf('    Found %d REM bouts, %d valid PLV (%.1f%%)\n', ...
                    n_bouts, valid_plv, 100*valid_plv/n_bouts);
            fprintf('    Mean PLV: %.4f\n', mean_plv);
        catch ME
            warning('%s', ['NREM PLV failed: ' ME.message]);
            ALL_RESULTS.(animal_id).(cond).REM_PLV = struct('plv_all',[],'dur_all',[],'t0_all',[]);
        end
        
        % ==============================================
        % C) NREM PAC (Delta phase × Spindle amplitude)
        % ==============================================
        fprintf('  Computing NREM PAC...\n');
        try
            PAC_NREM = ach_pac_comod_nrem_segment( ...
                mat_file, scores_csv, t_start, t_end, CODES, ...
                'phase_freqs', 0.5:0.5:4, ...   % Slow/delta
                'amp_freqs', 7:1:25, ...        % Spindles
                'nbins', 18, ...
                'doPlot', false);               % Don't plot for each animal
            
            ALL_RESULTS.(animal_id).(cond).NREM_PAC = PAC_NREM;
            
            fprintf('    MI range: [%.4f, %.4f]\n', ...
                    min(PAC_NREM.MI(:)), max(PAC_NREM.MI(:)));
        catch ME
                warning(ME.identifier, '%s', ME.message);
            ALL_RESULTS.(animal_id).(cond).NREM_PAC = struct('MI',[],'phase_freqs',[],'amp_freqs',[]);
        end
        
        % ==============================================
        % D) REM PAC (Theta phase × Gamma amplitude)
        % ==============================================
        fprintf('  Computing REM PAC...\n');
        try
            PAC_REM = ach_pac_comod_rem_segment( ...
                mat_file, scores_csv, t_start, t_end, CODES, ...
                'phase_freqs', 6:0.5:10, ...    % Theta
                'amp_freqs', 30:2:80, ...       % Gamma
                'nbins', 18, ...
                'doPlot', false);
            
            ALL_RESULTS.(animal_id).(cond).REM_PAC = PAC_REM;
            
            fprintf('    MI range: [%.4f, %.4f]\n', ...
                    min(PAC_REM.MI(:)), max(PAC_REM.MI(:)));
        catch ME
                warning(ME.identifier, '%s', ME.message);
            ALL_RESULTS.(animal_id).(cond).REM_PAC = struct('MI',[],'phase_freqs',[],'amp_freqs',[]);
        end
        
        fprintf('\n');
    end
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');
%%
% --------------------------------------------------------------
% 4) Save results
% --------------------------------------------------------------
save_file = 'all_animals_plv_pac_results.mat';
save(save_file, 'ALL_RESULTS', 'animals', 'CODES', 'conditions');
fprintf('Results saved to: %s\n', save_file);
% Before the plotting lines, add:
run('plot_summary_functions.m');

% Then call the plots:
plot_summary_nrem_plv(ALL_RESULTS, animals, conditions);
plot_summary_rem_plv(ALL_RESULTS, animals, conditions);
plot_summary_nrem_pac(ALL_RESULTS, animals, conditions);
plot_summary_rem_pac(ALL_RESULTS, animals, conditions);
plot_genotype_comparison(ALL_RESULTS, animals, conditions);
%%
% --------------------------------------------------------------
% 5) Generate summary statistics and plots
% --------------------------------------------------------------
fprintf('\n=== GENERATING SUMMARY PLOTS ===\n');

% Make sure results are loaded
load('all_animals_plv_pac_results.mat')

% Run the simple plotting script
simple_plotting_script