% run_REM_ambtemp_3h_master.m
% -------------------------------------------------------
% 3 h REM + MA analysis for BASELINE vs AMBTEMP only,
% using CROPPED scoring files (*_crop.csv).
%
% Produces:
%   - ALL_REM_3h.mat, ALL_SUMMARY_3h.mat
%   - Bar plots of REM metrics (3 h baseline vs ambtemp, WT vs APP)
%   - Survival-style REM duration plots baseline vs ambtemp

clear; clc;

% ----------- PATHS: ADAPT TO YOUR SETUP -----------
data_dir   = '/Users/margaridaseabra/24.11notchedsignal';      % EEG .mat (notched)
scores_dir = '/Users/margaridaseabra/24.11_cropped_ambtemp';
out_dir    = '/Users/margaridaseabra/24.11_REM_ambtemp_3h';    % output folder

if ~isfolder(out_dir); mkdir(out_dir); end

addpath(fullfile(pwd,'functions'));    % where your helper functions live
% (read_scores_csv, compute_epoch_spectra, build_bout_table,
%  annotate_REM_with_MA, cluster_REM_bouts, add_REM_size_category, ...)

% ----------- 1) Run REM+MA pipeline on cropped data -----------
[ALL_REM_3h, ALL_SUMMARY_3h] = run_REM_MA_pipeline_ambtemp_cropped( ...
                                    data_dir, scores_dir, out_dir);

save(fullfile(out_dir,'ALL_REM_3h.mat'),      'ALL_REM_3h');
save(fullfile(out_dir,'ALL_SUMMARY_3h.mat'),  'ALL_SUMMARY_3h');

%% ----------- 2) Bar plots + genotype×condition stats -----------
OUT_REMbars = plot_REM_3h_baseline_vs_ambtemp_APPvsWT( ...
                    ALL_SUMMARY_3h, out_dir, ...
                    'useFDR', true, ...
                    'minNperGroupForStats', 3);

% ----------- 3) REM duration "survival" plots + stats ----------
OUT_surv = plot_REM_survival_3h_baseline_vs_ambtemp( ...
                    ALL_REM_3h, out_dir);

fprintf('\n✓ REM ambtemp 3 h analysis finished. Outputs in: %s\n', out_dir);
%% If I want the other metrics to test for similarity
%% ----------- 2) Bar plots + genotype×condition stats -----------
OUT_REMbars = plot_REM_3h_baseline_vs_ambtemp_APPvsWT_similarity( ...
                    ALL_SUMMARY_3h, out_dir, ...
                    'useFDR', true, ...
                    'minNperGroupForStats', 3);

% ----------- 3) REM duration "survival" plots + stats ----------
OUT_surv = plot_REM_survival_3h_baseline_vs_ambtemp( ...
                    ALL_REM_3h, out_dir);

fprintf('\n✓ REM ambtemp 3 h analysis finished. Outputs in: %s\n', out_dir);