% After run_all_mice_eeg_psd_auto(...)
plot_group_psd_curves('EEG_PSD_AllMice');
csvFile = fullfile('EEG_PSD_AllMice', 'EEG_band_power_allmice.csv');
outPng  = fullfile('EEG_PSD_AllMice', 'GroupBandPower_grid.png');

plot_bandpower_grid(csvFile, outPng);
csvFile = fullfile('EEG_PSD_AllMice', 'EEG_band_power_allmice.csv');
outDir  = fullfile('EEG_PSD_AllMice', 'BandPower_ByCondition');

plot_bandpower_per_condition(csvFile, outDir);
