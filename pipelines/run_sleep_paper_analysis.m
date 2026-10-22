%% run_sleep_paper_analysis.m
% One-command master pipeline for the sleep / ambient-temperature paper.
%
% Run from the repository root:
%
%   cd('/Users/margaridaseabra/sleep_app_ach_warm')
%   run('pipelines/run_sleep_paper_analysis.m')
%
% This performs a clean publication-style rerun of:
%   1. metadata / source-file validation
%   2. standardized sleep architecture extraction
%   3. complete baseline analysis
%   4. baseline GraphPad Prism export
%   5. complete ambient-warming analysis
%   6. warming GraphPad Prism export
%
% NOTE:
% Existing child scripts use CLEAR, so this master script deliberately
% contains no state that needs to survive between stages.

fprintf('\n============================================================\n');
fprintf('SLEEP PAPER MASTER ANALYSIS\n');
fprintf('============================================================\n\n');

assert(isfile(fullfile('config','paths_config.m')), ...
    ['Repository root not detected. First cd to the sleep_app_ach_warm ' ...
     'repository, then run this script.']);

assert(isfile(fullfile('pipelines','run_00_prepare_data.m')), ...
    'Missing pipelines/run_00_prepare_data.m');
assert(isfile(fullfile('pipelines','run_01_sleep_architecture.m')), ...
    'Missing pipelines/run_01_sleep_architecture.m');
assert(isfile(fullfile('pipelines','run_baseline_complete.m')), ...
    'Missing pipelines/run_baseline_complete.m');
assert(isfile(fullfile('pipelines','run_06_export_baseline_prism.m')), ...
    'Missing pipelines/run_06_export_baseline_prism.m');
assert(isfile(fullfile('pipelines','run_warming_complete.m')), ...
    'Missing pipelines/run_warming_complete.m');
assert(isfile(fullfile('pipelines','run_07_export_warming_prism.m')), ...
    'Missing pipelines/run_07_export_warming_prism.m');

fprintf('1/6  Validate metadata and source files\n');
run('pipelines/run_00_prepare_data.m');

fprintf('\n2/6  Build standardized sleep architecture outputs\n');
run('pipelines/run_01_sleep_architecture.m');

fprintf('\n3/6  Run complete baseline analysis\n');
run('pipelines/run_baseline_complete.m');

fprintf('\n4/6  Export baseline data for GraphPad Prism\n');
run('pipelines/run_06_export_baseline_prism.m');

fprintf('\n5/6  Run complete ambient-warming analysis\n');
run('pipelines/run_warming_complete.m');

fprintf('\n6/6  Export warming data for GraphPad Prism\n');
run('pipelines/run_07_export_warming_prism.m');

fprintf('\n============================================================\n');
fprintf('MASTER PIPELINE COMPLETE\n');
fprintf('============================================================\n');
fprintf('Prism outputs:\n');
fprintf('  results/prism/baseline/\n');
fprintf('  results/prism/warming/\n\n');
fprintf('Not yet included in this master pipeline:\n');
fprintf('  - spindle analysis\n');
fprintf('  - phase-amplitude coupling / coupling analysis\n');
fprintf('============================================================\n');
