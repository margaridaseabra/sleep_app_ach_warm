function batch_notch_eeg_50Hz(inputDir, outputDir)
% batch_notch_eeg_50Hz
% -----------------------------------------------------------
% For every .mat file in inputDir:
%   - load variables
%   - apply 50 Hz notch (bw ≈ 3 Hz) to 'eeg' using eeg_frequency
%   - overwrite S.eeg with the notched signal
%   - save to outputDir with suffix '_notched50Hz_bw3.mat'
%
% Original files in inputDir are NOT modified.
%
% EXAMPLE:
%   inputDir  = '/path/to/20.11 scored files';
%   outputDir = '/path/to/20.11 scored files_notched';
%   batch_notch_eeg_50Hz(inputDir, outputDir);

    if nargin < 1 || isempty(inputDir)
        error('You must provide inputDir (folder with .mat files).');
    end

    if nargin < 2 || isempty(outputDir)
        % default: create sibling folder "<inputDir>_notched"
        outputDir = [inputDir '_notched'];
    end

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    % Adjust the pattern if you want (e.g. '*_scored.mat')
    matFiles = dir(fullfile(inputDir, '*.mat'));

    if isempty(matFiles)
        warning('No .mat files found in %s', inputDir);
        return;
    end

    for k = 1:numel(matFiles)
        inFile  = fullfile(inputDir, matFiles(k).name);
        [~, baseName, ~] = fileparts(matFiles(k).name);

        fprintf('Processing %s...\n', matFiles(k).name);

        % Load everything from the .mat file
        S = load(inFile);

        % Check that we have eeg + eeg_frequency
        if ~isfield(S, 'eeg') || ~isfield(S, 'eeg_frequency')
            warning('Skipping %s: missing eeg or eeg_frequency.', matFiles(k).name);
            continue;
        end

        x  = double(S.eeg);
        Fs = double(S.eeg_frequency);

        % -------- 50 Hz notch filter --------
        f0 = 50;    % Hz
        bw = 3;     % Hz, approximate bandwidth

        wo = f0 / (Fs/2);   % normalized notch freq (0–1)
        bwN = bw / (Fs/2);  % normalized bandwidth

        [b, a] = iirnotch(wo, bwN);

        x_notched = filtfilt(b, a, x);

        % Cast back to original class to save space (usually 'single')
        origClass = class(S.eeg);
        S.eeg = cast(x_notched, origClass);

        % Optional: add some metadata so you remember it's notched
        S.eeg_notch_meta = struct( ...
            'method', 'iirnotch + filtfilt', ...
            'f0_Hz',  f0, ...
            'bandwidth_Hz', bw);

        % Output filename
        outFile = fullfile(outputDir, [baseName '_notched50Hz_bw3.mat']);

        % Use -v7.3 for large files (lots of samples)
        save(outFile, '-struct', 'S', '-v7.3');

        fprintf('  -> saved %s\n', outFile);
    end

    fprintf('✅ Done. Notched EEG files saved in: %s\n', outputDir);
end
