function [sessionTbl, result, meta] = eeg_each_band_from_csv_ambtemp( ...
            eegFile, scoreFile, outputRootDir)

    % ---------- 1. Parse metadata from score filename ----------
    meta = parse_meta_from_scorefile_ambtemp(scoreFile);
    mouseID   = meta.mouseID;
    genotype  = meta.genotype;
    condition = meta.condition;

    % Output dir: root/genotype/condition/mouseID
    outputDir = fullfile(outputRootDir, genotype, condition, mouseID);
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    % ---------- 2. Load EEG ----------
    eegData = load(eegFile);
    if ~isfield(eegData, 'eeg') || ~isfield(eegData, 'eeg_frequency')
        error('EEG file %s must contain variables "eeg" and "eeg_frequency".', eegFile);
    end
    EEG = eegData.eeg;
    Fs  = eegData.eeg_frequency;
    % -------- Optional 50 Hz notch filter --------
    doNotch = false;   % set to false if your EEG is already notched

    if doNotch
        % 50 Hz notch with ~3 Hz bandwidth
        wo = 50/(Fs/2);      % normalised notch frequency
        bw = 3/(Fs/2);       % normalised bandwidth (≈ your "bw3")
        [b, a] = iirnotch(wo, bw);

        % zero-phase filtering so you don't smear transitions
        EEG = filtfilt(b, a, double(EEG));
    end


    % ---------- 3. Load scoring from CSV ----------
    T = readtable(scoreFile);

    % Try to find the score column – EDIT list if needed
    candidateCols = {'score','Score','state','State','Stage','SleepState','state_1Hz','score_1Hz'};
    scoreCol = '';
    for i = 1:numel(candidateCols)
        if ismember(candidateCols{i}, T.Properties.VariableNames)
            scoreCol = candidateCols{i};
            break;
        end
    end
    if isempty(scoreCol)
        error('Could not find a score column in %s', scoreFile);
    end

    rawScore = T.(scoreCol);

    % Convert to numeric score vector
    if iscell(rawScore) || isstring(rawScore) || ischar(rawScore)
        % Text labels case (e.g. W/N/R/MA)
        s = string(rawScore);
        score = zeros(numel(s),1);
        score(ismember(lower(s), {'w','wake','wk'}))  = 0;
        score(ismember(lower(s), {'n','nrem','nr'}))  = 1;
        score(ismember(lower(s), {'r','rem'}))        = 2;
        score(ismember(lower(s), {'ma','micro','microarousal'})) = 15;
    else
        % Already numeric
        score = rawScore;
    end

    % 1 Hz scoring
    epochLength     = 1;                  % seconds
    samplesPerEpoch = round(Fs * epochLength);

    % ---------- 4. Define state codes + merge MA → NREM ----------
    % ⚠️ IMPORTANT: If your numeric codes are different, EDIT these lines.
    % Example for simple coding: 1=Wake, 2=NREM, 3=REM, 4=MA
    C.WK   = 0;
    C.NREM = 1;
    C.REM  = 2;
    C.MA   = 15;

    % If your CSV uses e.g. 1=Wake, 4=NREM, 9=REM, 15=MA, change to:
    % C.WK = 1; C.NREM = 4; C.REM = 9; C.MA = 15;

    % Merge MA into NREM
    score(score == C.MA) = C.NREM;

    states = struct('Wake', C.WK, 'NREM', C.NREM, 'REM', C.REM);
    stateNames = fieldnames(states);

    % ---------- 5. Define frequency bands ----------
    bands = struct( ...
        'Delta',   [1,   4], ...
        'Theta',   [5,   9], ...
        'Sigma',   [10, 15], ...
        'Beta',    [15, 30], ...
        'lGamma1', [30, 49], ...
        'lGamma2', [51, 60], ...
        'hGamma',  [60,100]);
    bandNames = fieldnames(bands);

    % ---------- 6. Loop over states ----------
    result = struct();
    rows   = [];

    for s = 1:numel(stateNames)
        name = stateNames{s};
        val  = states.(name);

        epochs = find(score == val);
        if isempty(epochs)
            warning('%s | %s | %s: no epochs for %s', ...
                     mouseID, genotype, condition, name);
            continue;
        end

        allPxx = [];
        F_ref  = [];

        for i = 1:numel(epochs)
            startIdx = (epochs(i)-1)*samplesPerEpoch + 1;
            endIdx   = startIdx + samplesPerEpoch - 1;
            if endIdx > numel(EEG), continue; end

            segment = EEG(startIdx:endIdx);
            [pxx, F] = pwelch(segment, hamming(512), [], [], Fs);

            if isempty(F_ref), F_ref = F; end
            if numel(pxx) == numel(F_ref)
                allPxx(:, end+1) = pxx; %#ok<AGROW>
            end
        end

        if isempty(allPxx)
            warning('%s | %s | %s: no valid PSD for %s', ...
                     mouseID, genotype, condition, name);
            continue;
        end

        meanPxx = mean(allPxx, 2);
        result.(name).F   = F_ref;
        result.(name).PSD = meanPxx;

        % ---------- PSD figure ----------
        fig = figure('Visible','off');
        semilogy(F_ref, meanPxx, 'k', 'LineWidth', 1.5);
        xline(50, '--r', '50 Hz Notch', 'LabelOrientation','horizontal');
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density (µV^2/Hz)');
        title(sprintf('PSD - %s - %s (%s, %s)', ...
              name, mouseID, genotype, condition));
        grid on; xlim([0 100]);
        saveas(fig, fullfile(outputDir, sprintf('PSD_%s.png', name)));
        close(fig);

        % ---------- Band power integration ----------
        bandDb = struct();
        for b = 1:numel(bandNames)
            bname = bandNames{b};
            range = bands.(bname);
            idx   = F_ref >= range(1) & F_ref <= range(2);
            if any(idx)
                powerLinear = trapz(F_ref(idx), meanPxx(idx));
                bandDb.(bname) = 10 * log10(powerLinear);
            else
                bandDb.(bname) = NaN;
            end
        end
        result.(name).bandPowerDb = bandDb;

        % Long-format rows for this state
        for b = 1:numel(bandNames)
            bname = bandNames{b};
            valDb = bandDb.(bname);

            newRow = table( ...
                string(mouseID), ...
                string(genotype), ...
                string(condition), ...
                string(name), ...
                string(bname), ...
                valDb, ...
                'VariableNames', {'MouseID','Genotype','Condition', ...
                                  'State','Band','Power_dB'});
            rows = [rows; newRow]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        sessionTbl = table();
    else
        sessionTbl = rows;
    end

    save(fullfile(outputDir, 'PSD_data.mat'), 'result', 'Fs', ...
         'mouseID', 'genotype', 'condition');

    fprintf('✅ Finished EEG PSD for %s (%s, %s)\n', ...
             mouseID, genotype, condition);
end
