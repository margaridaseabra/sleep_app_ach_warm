load(fullfile('REM_MA_out','ALL_SUMMARY.mat'));  % adjust path if needed

STATS_baseline = plot_baseline_APP_vs_WT_REM_stats(ALL_SUMMARY, 'REM_MA_out', ...
    'baseline_cond', 'baseline', ...
    'remove_outliers', true);

%%
SURV_BASE = plot_baseline_REM_survival_APPvsWT(ALL_REM, 'REM_MA_out', ...
    'baseline_cond', 'baseline');   % or whatever your label is
%%

