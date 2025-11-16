function run_all_mice_eeg_psd_auto(eegDir, scoreDir)
    % Batch analysis of all EEG files in eegDir, using score CSVs in scoreDir

    baseOutDir = 'EEG_PSD_AllMice';
    if ~exist(baseOutDir, 'dir'), mkdir(baseOutDir); end

    % Adjust pattern if needed – this matches your example filenames
    eegFiles = dir(fullfile(eegDir, '*notched50Hz_bw3*.mat'));
    allTbl   = table();

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegName, ~] = fileparts(eegFiles(k).name);

        % 1) Remove '_notched...' suffix:
        eegBase = regexprep(eegName, '_notched.*$', '');
        % Ex: '20251023-drugs_mouse3'

        % 2) Replace '_' with '-' to match score filename style:
        % '20251023-drugs_mouse3' -> '20251023-drugs-mouse3'
        patternBase = strrep(eegBase, '_', '-');

        % 3) Score pattern: 20251023-drugs-mouse3-WT_scored_scores_1Hz.csv
        scorePattern = sprintf('%s-*_scored_scores_1Hz*.csv', patternBase);
        scoreCandidates = dir(fullfile(scoreDir, scorePattern));

        if isempty(scoreCandidates)
            warning('No score CSV found for EEG %s with pattern %s', ...
                    eegName, scorePattern);
            continue;
        elseif numel(scoreCandidates) > 1
            warning('Multiple CSVs found for EEG %s, using the first one.', eegName);
        end

        scoreFile = fullfile(scoreDir, scoreCandidates(1).name);

        % Run single-session analysis
        [sessionTbl, ~, meta] = eeg_each_band_from_csv( ...
            eegFile, scoreFile, baseOutDir);

        fprintf('📄 Done: %s (%s, %s) from %s\n', ...
                meta.mouseID, meta.genotype, meta.condition, eegName);

        allTbl = [allTbl; sessionTbl]; %#ok<AGROW>
    end

    % Save big summary CSV
    outCsv = fullfile(baseOutDir, 'EEG_band_power_allmice.csv');
    writetable(allTbl, outCsv);
    fprintf('📁 Saved group summary to %s\n', outCsv);
end
