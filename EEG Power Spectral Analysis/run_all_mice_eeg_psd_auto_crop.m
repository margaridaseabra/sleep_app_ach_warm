function run_all_mice_eeg_psd_auto_crop(eegDir, scoreDir)
    % Batch analysis of all EEG files in eegDir, using score CSVs in scoreDir
    %
    % EEG example:
    %   20251023-drugs-mouse3_APP.mat
    %
    % Score CSV example:
    %   20251023_drugs_mouse3_APP_scored_scores_1Hz.csv

    baseOutDir = 'EEG_PSD_AllMice_AmbTemp';
    if ~exist(baseOutDir, 'dir')
        mkdir(baseOutDir);
    end

    eegFiles = dir(fullfile(eegDir, '*.mat'));
    allTbl   = table();

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegNameFull, ~] = fileparts(eegFiles(k).name);
        % eegNameFull: e.g. '20251023-drugs-mouse3_APP'

        % 1) Strip any trailing processing suffixes, including _crop
        %    Works for names like:
        %    20251013-ambtemp-mouse2_WT_notched
        %    20251013_ambtemp_mouse2_WT_crop
        eegBase = regexprep(eegNameFull, '(_scored.*|_notched.*|_crop.*)$', '');

        % 2) Replace '-' with '_' to match CSV style
        %    20251013-ambtemp-mouse2_WT -> 20251013_ambtemp_mouse2_WT
        patternBase = strrep(eegBase, '-', '_');

        % 3) Score CSV pattern: allow extra suffix like "_crop"
        %    20251013_ambtemp_mouse2_WT_scored_scores_1Hz.csv
        %    20251013_ambtemp_mouse2_WT_scored_scores_1Hz_crop.csv
        scorePattern = sprintf('%s*scored_scores_1Hz*.csv', patternBase);

        scoreCandidates = dir(fullfile(scoreDir, scorePattern));

        if isempty(scoreCandidates)
            warning('No score CSV found for EEG %s with pattern %s', ...
                    eegNameFull, scorePattern);
            continue;
        elseif numel(scoreCandidates) > 1
            warning('Multiple CSVs found for EEG %s, using the first one.', ...
                    eegNameFull);
        end

        scoreFile = fullfile(scoreDir, scoreCandidates(1).name);

        % 4) Run single-session analysis
        [sessionTbl, ~, meta] = eeg_each_band_from_csv_ambtemp( ...
            eegFile, scoreFile, baseOutDir);

        fprintf(' Done: %s (%s, %s) from %s\n', ...
                meta.mouseID, meta.genotype, meta.condition, eegNameFull);

        allTbl = [allTbl; sessionTbl];
    end

    % 5) Save big summary CSV
    outCsv = fullfile(baseOutDir, 'EEG_band_power_allmice.csv');
    writetable(allTbl, outCsv);
    fprintf(' Saved group summary to %s\n', outCsv);
end
