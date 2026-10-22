%% run_sleep_paper_analysis_v2.m
% One-command macroarchitecture pipeline: baseline + warming + donepezil.
%
% Spindle/PAC processing is intentionally not launched until
% run_09_validate_eeg_sources.m confirms that raw EEG is available.

fprintf('\n============================================================\n');
fprintf('SLEEP PAPER MASTER ANALYSIS V2\n');
fprintf('============================================================\n\n');

assert(isfile(fullfile('config','paths_config.m')), ...
    'cd to the sleep_app_ach_warm repository before running this script.');

run('pipelines/run_00_prepare_data.m');
run('pipelines/run_01_sleep_architecture.m');

run('pipelines/run_baseline_complete.m');
run('pipelines/run_06_export_baseline_prism.m');

run('pipelines/run_warming_complete.m');
run('pipelines/run_07_export_warming_prism.m');

run('pipelines/run_donepezil_complete.m');
run('pipelines/run_08_export_donepezil_prism.m');

fprintf('\n============================================================\n');
fprintf('MACROARCHITECTURE PIPELINE COMPLETE\n');
fprintf('============================================================\n');
fprintf('Prism:\n');
fprintf('  results/prism/baseline/\n');
fprintf('  results/prism/warming/\n');
fprintf('  results/prism/donepezil/\n\n');
fprintf('Next EEG check:\n');
fprintf('  run(''pipelines/run_09_validate_eeg_sources.m'')\n');
fprintf('============================================================\n');
