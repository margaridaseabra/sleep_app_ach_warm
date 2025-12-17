%% Baseline ACh analysis: WT vs APP under baseline condition only

% --- EDIT THESE PATHS ---
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';
outDir   = fullfile(sigDir, 'ACh_baseline_group');
% -------------------------

if ~exist(outDir, 'dir'); mkdir(outDir); end

% 1) Run the batch analysis (all conditions present in those folders)
GROUP_all = run_ach_batch_auto(sigDir, scoreDir);
%%
% 2) Keep only 'baseline' sessions
GROUP_baseline = subset_group_by_cond(GROUP_all, 'baseline');

% 3) Save for safety
save(fullfile(outDir, 'ACh_GROUP_baseline_only.mat'), ...
     'GROUP_baseline','-v7.3');

% 4) Summary bar plots + t-tests (the figure you already got)
ach_stats_baseline_ttest(GROUP_baseline, outDir);

%% 5) Wake / NREM / REM onset traces + peak bar plots + stats
group_plot_all_transitions(GROUP_baseline, outDir,false);

% 6) NREM ACh PSD traces + "peak oscillation" bars (power, freq, amp)
ach_plot_nrem_psd_baseline(GROUP_baseline, outDir,false);

% 7) Wake ACh PSD traces + "peak oscillation" bars (power, freq, amp)
ach_plot_Wake_psd_baseline(GROUP_baseline, outDir,false);
