function run_all_mice_rem_onset_wavelets(eegDir, scoreDir, outBaseDir)

    if nargin < 3 || isempty(outBaseDir)
        outBaseDir = 'EEG_REM_wavelets';
    end
    if ~exist(outBaseDir,'dir'); mkdir(outBaseDir); end

    eegFiles = dir(fullfile(eegDir, '*.mat'));

    for k = 1:numel(eegFiles)
        eegFile = fullfile(eegDir, eegFiles(k).name);
        [~, eegNameFull] = fileparts(eegFiles(k).name);

        % same naming logic as your PSD script:
        eegBase     = regexprep(eegNameFull, '(_scored.*|_notched.*)$', '');
        patternBase = strrep(eegBase, '-', '_');
        scorePattern = sprintf('%s_scored_scores_1Hz.csv', patternBase);

        scoreCandidates = dir(fullfile(scoreDir, scorePattern));
        if isempty(scoreCandidates)
            warning('No score CSV found for %s', eegNameFull);
            continue;
        end
        scoresCsv = fullfile(scoreDir, scoreCandidates(1).name);

        % per-session output folder (optional)
        outDir = fullfile(outBaseDir, eegBase);
        if ~exist(outDir,'dir'); mkdir(outDir); end

        fprintf('\n=== REM-onset wavelets for %s ===\n', eegNameFull);
        rem_onset_wavelet_session(eegFile, scoresCsv, outDir);
    end
end
