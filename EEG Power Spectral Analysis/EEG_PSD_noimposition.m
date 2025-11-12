% EEG PSD Analysis per State (Unbiased, Data-Driven)

% ---- Setup ----
eegFile = '/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/20251002-mouse2-baseline_test.mat';
scoreFile = '/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/mouse2_base_scores_tt.mat';

% Load data
eegData = load(eegFile);
scoreData = load(scoreFile);

% Extract EEG signal and sampling rate
EEG = eegData.eeg;
Fs = eegData.eeg_frequency;

% Scoring info
TT = scoreData.TT;       % Timetable with state labels
C = scoreData.C;         % Contains state codes like WK, NREM, REM, MA
score = TT{:,1};         % Numeric score vector
epochLength = scoreData.epochSec;
samplesPerEpoch = round(Fs * epochLength);

% Define states (include MA in NREM)
states = struct('Wake', C.WK, 'NREM', [C.NREM, C.MA], 'REM', C.REM);
stateNames = fieldnames(states);

% Create output folder
outputFolder = 'EEG_PSD_Per_State';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% ---- Main Analysis ----
PSD_all = struct();

for s = 1:length(stateNames)
    stateName = stateNames{s};
    stateVals = states.(stateName);
    
    % Find relevant epochs
    epochs = find(ismember(score, stateVals));
    fprintf('%s: %d epochs\n', stateName, length(epochs));
    if isempty(epochs), continue; end

    % Compute PSD for each epoch and accumulate
    allPxx = [];
    F_ref = [];
    for i = 1:length(epochs)
        startIdx = round((epochs(i)-1)*samplesPerEpoch + 1);
        endIdx   = round(startIdx + samplesPerEpoch - 1);
        if endIdx <= length(EEG)
            segment = EEG(startIdx:endIdx);
            [pxx, F] = pwelch(segment, hamming(512), [], [], Fs);

            if isempty(F_ref)
                F_ref = F;  % store frequency vector
            end

            if length(pxx) ~= length(F_ref)
                continue;
            end
            allPxx(:, end+1) = pxx;
        end
    end

    % Average PSD
    meanPxx = mean(allPxx, 2);
    PSD_all.(stateName) = meanPxx;

    % Plot individual PSD
    fig = figure('Visible','off');
    semilogy(F_ref, meanPxx, 'k', 'LineWidth', 1.5);
    xline(50, '--r', '50 Hz Notch');
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (uV^2/Hz)');
    title(['PSD - ' stateName]);
    xlim([0 100]);
    saveas(fig, fullfile(outputFolder, [stateName '_PSD.png']));
    close(fig);
end

% ---- Overlay PSD Comparison ----
fig = figure;
hold on;
colors = lines(length(stateNames));
for i = 1:length(stateNames)
    stateName = stateNames{i};
    semilogy(F_ref, PSD_all.(stateName), 'LineWidth', 1.5, 'Color', colors(i,:));
end
xline(50, '--r', '50 Hz Notch');
legend(stateNames, 'Location', 'best');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (uV^2/Hz)');
title('PSD Comparison Across Sleep States');
xlim([0 100]);
grid on;
saveas(fig, fullfile(outputFolder, 'PSD_Overlay.png'));
close(fig);

disp('Unbiased EEG PSD analysis complete.');
