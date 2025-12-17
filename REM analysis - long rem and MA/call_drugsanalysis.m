% run_REM_drugss_3h_master.m
% -------------------------------------------------------
% 3 h REM + MA analysis for BASELINE vs drugss only,
% using CROPPED scoring files (*_crop.csv).
%
% Produces:
%   - ALL_REM_3h.mat, ALL_SUMMARY_3h.mat
%   - Bar plots of REM metrics (3 h baseline vs drugss, WT vs APP)
%   - Survival-style REM duration plots baseline vs drugss

clear; clc;

% ----------- PATHS: ADAPT TO YOUR SETUP -----------
data_dir   = '/Users/margaridaseabra/24.11notchedsignal';      % EEG .mat (notched)
scores_dir = '/Users/margaridaseabra/24.11_cropped_3h_washout';
out_dir    = '/Users/margaridaseabra/24.11_REM_drugs_6h';      % output folder

if ~isfolder(out_dir); mkdir(out_dir); end

addpath(fullfile(pwd,'functions'));

% ----------- 1) Run REM+MA pipeline on cropped data -----------
[ALL_REM_6h, ALL_SUMMARY_6h] = run_REM_MA_pipeline_drugs_6h( ...
                                    data_dir, scores_dir, out_dir);

save(fullfile(out_dir,'ALL_REM_6h.mat'),     'ALL_REM_6h');
save(fullfile(out_dir,'ALL_SUMMARY_6h.mat'), 'ALL_SUMMARY_6h');

%% ----------- 2) Bar plots + genotype×condition stats -----------
OUT_REMbars = plot_REM_6h_baseline_vs_drugs_APPvsWT( ...
                    ALL_SUMMARY_6h, out_dir, ...
                    'useFDR', true, ...
                    'minNperGroupForStats', 3);

% ----------- 3) REM duration "survival" plots + stats ----------
OUT_surv = plot_REM_survival_6h_baseline_vs_drugs( ...
                    ALL_REM_6h, out_dir);

fprintf('\n✓ REM drugss 6 h analysis finished. Outputs in: %s\n', out_dir);

%% If I want the other metrics to test for similarity
% %% ----------- 2) Bar plots + genotype×condition stats -----------
% OUT_REMbars = plot_REM_6h_baseline_vs_drugss_APPvsWT_similarity( ...
%                     ALL_SUMMARY_6h, out_dir, ...
%                     'useFDR', true, ...
%                     'minNperGroupForStats', 3);
% 
% % ----------- 3) REM duration "survival" plots + stats ----------
% OUT_surv = plot_REM_survival_6h_baseline_vs_drugss( ...
%                     ALL_REM_3h, out_dir);
% 
% fprintf('\n✓ REM drugss 3 h analysis finished. Outputs in: %s\n', out_dir);