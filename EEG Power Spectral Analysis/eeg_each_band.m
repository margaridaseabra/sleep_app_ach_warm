function eeg_each_band(eegFile, scoreFile, outputDir, mouseID)

    % Create output directory
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    % Load EEG and scoring data
    eegData = load(eegFile);
    scoreData = load(scoreFile);

    EEG = eegData.eeg;
    Fs = eegData.eeg_frequency;
    TT = scoreData.TT;
    C = scoreData.C;
    score = TT{:,1};
    epochLength = scoreData.epochSec;
    samplesPerEpoch = round(Fs * epochLength);

    % Merge MA (15) into NREM (1)
    score(score == C.MA) = C.NREM;

    % Define sleep states
    states = struct('Wake', C.WK, 'NREM', C.NREM, 'REM', C.REM);
    stateNames = fieldnames(states);

    % Define frequency bands (once)
    bands = struct(...
        'Delta',  [1, 4], ...
        'Theta',  [5, 9], ...
        'Sigma',  [10, 15], ...
        'Beta',   [15, 30], ...
        'lGamma1', [30, 49], ...
        'lGamma2', [51, 60], ...
        'hGamma', [60, 100]);

    bandNames = fieldnames(bands);

    % Initialize result and summary table
    result = struct();
    csvTable = table(string(mouseID), 'VariableNames', {'MouseID'});


    % Loop over states
    for s = 1:length(stateNames)
        name = stateNames{s};
        val = states.(name);
        epochs = find(score == val);

        if isempty(epochs)
            warning('%s: No epochs found for %s', mouseID, name);
            continue;
        end

        % Compute PSD per epoch
        allPxx = [];
        F_ref = [];
        for i = 1:length(epochs)
            startIdx = (epochs(i)-1)*samplesPerEpoch + 1;
            endIdx   = startIdx + samplesPerEpoch - 1;
            if endIdx <= length(EEG)
                segment = EEG(startIdx:endIdx);
                [pxx, F] = pwelch(segment, hamming(512), [], [], Fs);
                if isempty(F_ref), F_ref = F; end
                if length(pxx) == length(F_ref)
                    allPxx(:, end+1) = pxx;
                end
            end
        end

        % Average PSD
        meanPxx = mean(allPxx, 2);
        result.(name).F = F_ref;
        result.(name).PSD = meanPxx;

        % --- PSD Plot ---
        fig = figure('Visible','off');
        semilogy(F_ref, meanPxx, 'k', 'LineWidth', 1.5);
        xline(50, '--r', '50 Hz Notch', 'LabelOrientation','horizontal');
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density (uV^2/Hz)');
        title(['PSD - ' name ' - ' mouseID]);
        grid on; xlim([0 100]);
        saveas(fig, fullfile(outputDir, ['PSD_' name '.png']));
        close(fig);

        % --- Band power computation (in dB) ---
        bandDb = struct();
        for b = 1:length(bandNames)
            bname = bandNames{b};
            range = bands.(bname);
            idx = F_ref >= range(1) & F_ref <= range(2);
            if any(idx)
                powerLinear = trapz(F_ref(idx), meanPxx(idx));
                powerDb = 10 * log10(powerLinear);
                bandDb.(bname) = powerDb;
            else
                bandDb.(bname) = NaN;
            end
        end
        result.(name).bandPowerDb = bandDb;

        % --- Bar plot per state ---
        fig = figure('Visible','off');
        barData = cell2mat(struct2cell(bandDb));
        bar(barData);
        set(gca, 'XTickLabel', fieldnames(bandDb), 'XTickLabelRotation', 45);
        ylabel('Power (dB)');
        title(['Band Power - ' name ' - ' mouseID]);
        grid on;
        saveas(fig, fullfile(outputDir, ['BandPower_' name '.png']));
        close(fig);

        % --- Add band data to CSV summary ---
        for b = 1:length(bandNames)
            csvTable.([name '_' bandNames{b} '_dB']) = bandDb.(bandNames{b});
        end
    end

    % Save results
    save(fullfile(outputDir, 'PSD_data.mat'), 'result', 'Fs', 'mouseID');
    writetable(csvTable, fullfile(outputDir, 'PSD_data.csv'));

    fprintf('✅ Analysis complete for %s\n', mouseID);
end
