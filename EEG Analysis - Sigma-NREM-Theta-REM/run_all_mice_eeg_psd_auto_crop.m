function run_all_mice_eeg_psd_auto_crop(eegDir, scoreDir)
    % Batch analysis of all EEG files in eegDir, using score CSVs in scoreDir
    %
    % EEG example (cropped/notched):
    %   20251013_ambtemp_mouse2_WT_crop_notched.mat
    %
    % Score CSV example:
    %   20251013_ambtemp_mouse2_WT_scored_scores_1Hz_crop.csv

    baseOutDir = 'EEG_PSD_AllMice_AmbTemp';
    if ~exist(baseOutDir, 'dir')
        mkdir(baseOutDir);
    end

    eegFiles = dir(fullfile(eegDir, '*.mat'));
    if isempty(eegFiles)
        error('No EEG .mat files found in %s', eegDir);
    end

    allTbl   = table();
    fprintf('Found %d EEG files in %s\n', numel(eegFiles), eegDir);

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegNameFull, ~] = fileparts(eegFiles(k).name);

        fprintf('\n=== Processing EEG file %d/%d: %s ===\n', ...
                k, numel(eegFiles), eegNameFull);

        % 1) Strip any trailing processing suffixes, including _crop / _notched
        eegBase = regexprep(eegNameFull, '(_scored.*|_notched.*|_crop.*)$', '');
        fprintf('  eegBase: %s\n', eegBase);

        % 2) Replace '-' with '_' to match CSV style
        patternBase = strrep(eegBase, '-', '_');
        fprintf('  patternBase: %s\n', patternBase);

        % 3) Score CSV pattern: allow extra suffix like "_crop"
        scorePattern = sprintf('%s*scored_scores_1Hz*.csv', patternBase);
        fprintf('  Looking for score files with pattern: %s\n', scorePattern);

        scoreCandidates = dir(fullfile(scoreDir, scorePattern));
        if isempty(scoreCandidates)
            warning('  No score CSV found for EEG %s with pattern %s', ...
                    eegNameFull, scorePattern);
            continue;
        elseif numel(scoreCandidates) > 1
            warning('  Multiple score CSVs found, using the first one: %s', ...
                    scoreCandidates(1).name);
        end

        scoreFile = fullfile(scoreDir, scoreCandidates(1).name);
        fprintf('  Using score file: %s\n', scoreCandidates(1).name);

        % 4) Run single-session analysis
        [sessionTbl, ~, meta] = eeg_each_band_from_csv_ambtemp( ...
            eegFile, scoreFile, baseOutDir);

        fprintf('  sessionTbl rows: %d\n', height(sessionTbl));

        if isempty(sessionTbl)
            warning('  No bandpower rows produced for %s (check eeg_each_band_from_csv).', ...
                    eegNameFull);
            continue;
        end

        fprintf('  Done: %s (%s, %s) from %s\n', ...
                meta.mouseID, meta.genotype, meta.condition, eegNameFull);

        allTbl = [allTbl; sessionTbl];
    end

    fprintf('\nTotal rows in allTbl: %d\n', height(allTbl));

    % 5) Save big summary CSV (only if non-empty)
    outCsv = fullfile(baseOutDir, 'EEG_band_power_allmice.csv');
    if isempty(allTbl)
        warning('AllTbl is empty – writing an empty CSV would break later steps. NOT saving.');
    else
        writetable(allTbl, outCsv);
        fprintf(' Saved group summary to %s\n', outCsv);
    end
end
