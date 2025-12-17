function run_all_mice_eeg_psd_auto(eegDir, scoreDir)
    % Batch analysis of all EEG files in eegDir, using score CSVs in scoreDir
    %
    % EEG example:
    %   20251023-drugs-mouse3_APP.mat
    %
    % Score CSV example:
    %   20251023_drugs_mouse3_APP_scored_scores_1Hz.csv

    baseOutDir = 'EEG_PSD_AllMice';
    if ~exist(baseOutDir, 'dir')
        mkdir(baseOutDir);
    end

    eegFiles = dir(fullfile(eegDir, '*.mat'));
    allTbl   = table();

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegNameFull, ~] = fileparts(eegFiles(k).name);
        % eegNameFull: e.g. '20251023-drugs-mouse3_APP'

        % 1) Strip any trailing processing suffix if it ever exists
        %    (safe even if there is none)
        eegBase = regexprep(eegNameFull, '(_scored.*|_notched.*)$', '');

        % 2) Replace '-' with '_' to match score naming
        %    '20251023-drugs-mouse3_APP' -> '20251023_drugs_mouse3_APP'
        patternBase = strrep(eegBase, '-', '_');

        % 3) Score CSV pattern: exact match of base + '_scored_scores_1Hz.csv'
        %    -> '20251023_drugs_mouse3_APP_scored_scores_1Hz.csv'
        scorePattern = sprintf('%s_scored_scores_1Hz.csv', patternBase);

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
        [sessionTbl, ~, meta] = eeg_each_band_from_csv( ...
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
