% If REM is coded as 2 in sleep_scores:
%plot_REM_bout_full('20251002-baseline-mouse2_newlyscored.mat', 2, 1);

% Try another bout, e.g. the one with a lot of EMG bursts:
plot_REM_bout_full('20251002-baseline-mouse2_newlyscored.mat', 2, 29);

%% Example: REM coded as 2, look at bout 15 (same as your screenshot):
analyze_REM_burst_locked('20251002-baseline-mouse2_newlyscored.mat', 2, 29);
