%% Notch the files that are not notched 
batch_notch_50Hz_to_folder('/Users/margaridaseabra/24.11 signalnotscored', ...
                           '/Users/margaridaseabra/24.11notchedsignal');
%% BASELINE
%% PSD Analysis for notched files
run_all_mice_eeg_psd_auto( ...
    '/Users/margaridaseabra/24.11notchedsignal', ...
    '/Users/margaridaseabra/24.11scores');
%%
plot_group_psd_curves('EEG_PSD_AllMice');
%%
plot_bandpower_baseline_stats( ...
    'EEG_PSD_AllMice/EEG_band_power_allmice.csv', []);
%%
bandpower_ratio_baseline('EEG_PSD_AllMice/EEG_band_power_allmice.csv', []);
%%
plot_baseline_bandpower_focus('EEG_PSD_AllMice/EEG_band_power_allmice.csv', []);
%% AMBIENT TEMPERATURE
run_all_mice_eeg_psd_auto_crop( ...
    '/Users/margaridaseabra/24.11_cropped_ambtemp', ...
    '/Users/margaridaseabra/24.11_cropped_ambtemp');
%%
csvFile = 'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv';
%%
plot_group_psd_curves('EEG_PSD_AllMice_AmbTemp');
%%
% Ambtemp-only bar plots with stars:
plot_bandpower_cond_stats_ambtemp( ...
  'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
  'ambtemp', []);
%%
plot_bandpower_cond_stats_ambtemp( ...
    'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
    'APP', [], 'condition');
%%
plot_bandpower_cond_stats_ambtemp( ...
    'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
    'WT', [], 'condition');

%%
csvFile = 'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv';

plot_bandpower_geno_cond_anova( ...
    csvFile, ...
    "NREM", ["Sigma","Beta","lGamma1"], ...
    'EEG_PSD_AllMice_AmbTemp/NREM_baseline_vs_ambtemp_geno_cond.png');

%%
plot_bandpower_geno_cond_anova( ...
    csvFile, ...
    "REM", ["Sigma","Beta","lGamma1","lGamma2"], ...
    'EEG_PSD_AllMice_AmbTemp/REM_baseline_vs_ambtemp_geno_cond.png');

%%
csvFile = 'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv';

% NREM
plot_bandpower_change_by_genotype(csvFile, "NREM", ["Sigma","Beta","lGamma1"], []);

% REM
plot_bandpower_change_by_genotype(csvFile, "REM", ["Sigma","Beta","lGamma1","lGamma2"], []);

%% Compare WT baseline vs WT ambient temperature


%% Compare APP baseline vs APP ambient temperature

