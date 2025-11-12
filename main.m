%% REM selection
select_REM_episodes('file','N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\scored files\20251007_ambtemp_mouse4_notched_scored.mat');
select_REM_for_TFR_context('file', 'N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\scored files\20251007_ambtemp_mouse4_notched_scored.mat', 'PRE_SEC', 5, 'POST_SEC', 5);

%% REM analysis
rem_timefreq_analysis('file','N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\scored files\20251007_ambtemp_mouse4_notched_scored_REMselection_CONTEXT_forTFR.mat', ...
                      'method','spectrogram', 'norm','none');

%% Conclusions
run_rem_conclusions('N:\SUN-IN-Kjaerby\0_Personal folders\Margarida\REM analysis\features');
