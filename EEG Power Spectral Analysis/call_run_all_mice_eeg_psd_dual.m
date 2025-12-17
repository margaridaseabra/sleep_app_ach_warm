%% EEG BwASELINE ANALYSIS: WT vs APP
eegDir    = '/Users/margaridaseabra/24.11notchedsignal';
scoreDirB = '/Users/margaridaseabra/24.11scores';

% 1) Per-epoch band power for the "baseline dataset"
run_all_mice_eeg_psd_auto_dual(eegDir, scoreDirB, ...
    '/Users/margaridaseabra/EEG_PSD_Baseline');
bandCsvB = '/Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv';



%%
bandCsv_base = '/Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv';
outDir_base  = '/Users/margaridaseabra/EEG_group_baseline_simple';

OUT_EEG_base = run_EEG_bandpower_APPvsWT_simple( ...
    bandCsv_base, outDir_base, ...
    'condition','baseline', ...   % only baseline rows
    'useFDR', true, ...           % FDR ON for stars
    'minN', 3);                   % require >=3 mice per genotype






%% 2) Group stats: BASELINE ONLY, APP vs WT;

bandCsv = '/Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv';
outDir  = '/Users/margaridaseabra/EEG_group_baseline';

OUT_EEG_base = run_EEG_bandpower_group_stats(bandCsv, outDir);

%% 
scoreDir = '/Users/margaridaseabra/24.11scores';
outDirPSD  = '/Users/margaridaseabra/EEG_group_baseline';
OUT_psd_base = run_group_psd_APPvsWT(eegDir, scoreDir, outDirPSD, ...
    'condition','baseline');

%% EEG AMBTEMP ANALYSIS
% 1) Per-epoch band power for 3 h windows
run_all_mice_eeg_psd_auto_dual(eegDir, scoreDir3h, ... 
    '/Users/margaridaseabra/EEG_PSD_Amb3h');

bandCsv3h = '/Users/margaridaseabra/EEG_PSD_Amb3h/EEG_band_power_allmice.csv';

% 2) Group stats: baseline-3h vs ambtemp-3h, APP vs WT
OUT_EEG_3h = run_EEG_bandpower_group_stats_manual( ...
    bandCsv3h, '/Users/margaridaseabra/EEG_group_ambtemp3h');
%%
% --- 1. Compute PSD and band power per mouse ---
eegDir   = '/Users/margaridaseabra/24.11notchedsignal';   % your .mat
scoreDir = '/Users/margaridaseabra/24.11scores';          % all scoring CSVs

run_all_mice_eeg_psd_auto_dual(eegDir, scoreDir);
% -> writes /Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv
%%
% 1) Create per-epoch band power table (already working)
eegDir   = '/Users/margaridaseabra/24.11notchedsignal';   % raw EEG
scoreDir = '/Users/margaridaseabra/24.11scores';          % full-night scores

run_all_mice_eeg_psd_auto(eegDir, scoreDir);  
% -> /Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv

% 2) Group stats + band bars
bandCsv = '/Users/margaridaseabra/EEG_PSD_Baseline/EEG_band_power_allmice.csv';
outDir  = '/Users/margaridaseabra/EEG_group_baseline';

OUT_EEG_base = run_EEG_bandpower_group_stats(bandCsv, outDir);
