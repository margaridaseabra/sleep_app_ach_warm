meta_csv   = '/Users/margaridaseabra/sleep_app_ach_warm/Ambtemp prep-3h-selection/ambtemp-times-all-mice.xlsx';   % your Mouse / Genotype / times
scores_dir = '/Users/margaridaseabra/24.11scores';
mats_dir   = '/Users/margaridaseabra/24.11 signalnotscored';
out_dir    = '/Users/margaridaseabra/24.11_cropped_ambtemp';

OUT = crop_baseline_ambtemp_segments(meta_csv, scores_dir, mats_dir, out_dir, ...
    'eegVar','eeg', ...
    'fsVar','eeg_frequency', ...
    'ambVar','ambtemp', ...      % if this exists in your .mat
    'ambFsVar','ambtemp_fs');    % or '' if same fs as EEG
