function call_run_REM_NREM_MA_transition_analysis
% Wrapper to run the REM–NREM–MA transition analysis
% and make all the visualisation plots.

    % ---- 1) Point to your scores folder ----
    input_dir = '/Users/margaridaseabra/24.11scores';

    % ---- 2) Run transition analysis ----
    OUT_TRANS = run_REM_NREM_MA_transition_analysis(input_dir, ...
        'pattern', '*_scores_1Hz.csv', ...
        'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), ...
        'quantiles', [33 66], ...      % Short / Medium / Long
        'ma_window_sec', 120, ...      % window for MA before REM
        'wake_include_MA', true, ...   % count MA as wake after REM
        'out_dir', fullfile(input_dir, 'rem_transitions'));

    % ---- 3) Extract tables and output dir ----
    REM_TRANS  = OUT_TRANS.REM_TRANS;
    NREM_TRANS = OUT_TRANS.NREM_TRANS;
    outDir     = OUT_TRANS.params.out_dir;

    % ---- 4) Make all the plots ----
    plot_prev_nrem_before_REM_by_type(REM_TRANS, outDir);
    plot_MA_before_REM_by_type(REM_TRANS, outDir);
    plot_wake_after_REM_by_type(REM_TRANS, outDir);
    plot_NREM_sequence_patterns(NREM_TRANS, outDir);

    fprintf('✅ All REM–NREM–MA transition plots saved in: %s\n', outDir);
end
