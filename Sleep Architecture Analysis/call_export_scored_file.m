R = export_sleep_scores_with_MA('/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/20251001_mouse1_baseline_scored.mat','mouse2_base', ...
    'out_dir','/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files', ...
    'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), ...
    'ma_thresh_sec',15);

% Quick check
if R.success
    disp(R.files)                  % paths of the exported files
    fprintf('MA bouts: %d\n', R.n_MA_bouts);
end