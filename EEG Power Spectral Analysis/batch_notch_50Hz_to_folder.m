function batch_notch_50Hz_to_folder(eegDir, outDir)
% Create a new folder with 50 Hz–notched EEG.
% Files in "alreadyNotched" are just copied (no extra filter).
%
% Usage:
%   batch_notch_50Hz_to_folder('/path/to/raw_eeg', '/path/to/eeg_notched');

    if nargin < 1 || isempty(eegDir)
        eegDir = pwd;
    end
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(eegDir, 'eeg_notched_50Hz');
    end
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % --- 5 files that are already notched ---
    alreadyNotched = { ...
        '20251001_baseline_mouse1_APP', ...
        '20251002_baseline_mouse2_WT', ...
        '20251003_ambtemp_mouse1_APP', ...
        '20251005_baseline_mouse8_WT', ...
        '20251006_baseline_mouse4_WT' ...
    };

    files = dir(fullfile(eegDir, '*.mat'));

    for k = 1:numel(files)
        inFile = fullfile(eegDir, files(k).name);
        [~, base, ext] = fileparts(files(k).name);
        outFile = fullfile(outDir, [base ext]);

        % 1) If this file is in the "already notched" list → just copy
        if ismember(base, alreadyNotched)
            copyfile(inFile, outFile);
            fprintf('⏭️  Copied already-notched file: %s\n', base);
            continue;
        end

        % 2) Otherwise: load and notch
        S = load(inFile);
        if ~isfield(S,'eeg') || ~isfield(S,'eeg_frequency')
            warning('File %s lacks eeg or eeg_frequency, skipping.', base);
            continue;
        end

        eeg_raw = double(S.eeg(:));
        fs      = S.eeg_frequency;

        mains = 50;      % Hz
        bw_Hz = 3;       % bandwidth
        harmonics = mains:mains:floor(fs/2);

        eeg_notched = eeg_raw;
        for f0 = harmonics
            wo = f0/(fs/2);
            bw = bw_Hz/(fs/2);
            [b,a] = iirnotch(wo, bw);
            eeg_notched = filtfilt(b,a,eeg_notched);
        end

        % Keep original as eeg_raw, store notched in eeg
        S.eeg_raw   = S.eeg;
        S.eeg       = single(eeg_notched);
        S.notch_info = struct('mains',mains, ...
                              'bw_Hz',bw_Hz, ...
                              'fs',fs, ...
                              'harmonics',harmonics);

        save(outFile, '-struct', 'S', '-v7.3');
        fprintf('✅ Notched and saved: %s\n', base);
    end
end
