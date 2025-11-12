% Load EEG and scoring data
eegFile = '/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/20251002-mouse2-baseline_test.mat';
scoreFile = '/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/mouse2_base_scores_tt.mat';
%%
eegData = load(eegFile);
scoreData = load(scoreFile);
%%
% Extract EEG signal and sampling frequency
EEG = eegData.eeg;
Fs=eegData.eeg_frequency;

TT = scoreData.TT;       % Timetable containing state annotations
C = scoreData.C;         % Contains state codes and parameters

% Convert timetable to state vector (assumes first variable is the score)
score = TT{:,1};         % Get numeric score vector from timetable
epochLength = scoreData.epochSec;
samplesPerEpoch = round(Fs * epochLength);
% Store sleep stage labels (reverse lookup)
states = struct('Wake', C.WK, 'NREM', C.NREM, 'REM', C.REM);

% Frequency bands
bands = struct('Delta', [1, 4], 'Theta', [5, 9], ...
               'Alpha', [10, 15], 'Beta', [15, 30], ...
               'lGamma1', [30, 49],'lGamma2', [51, 60], 'hGamma' ,[60 100] );

states = struct('Wake', 0, 'NREM', 1, 'REM', 2);
outputFolder = 'EEG_Power_Analysis';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

stateName = stateNames{s};
stateVal = states.(stateName);
epochs = find(score == stateVal);
if isempty(epochs); end

% --- PSD averaging across epochs ---
allPxx = [];
F_ref = [];
for i = 1:length(epochs)
    startIdx = round((epochs(i)-1)*samplesPerEpoch + 1);
    endIdx   = round(startIdx + samplesPerEpoch - 1);
    if endIdx <= length(EEG)
        segment = EEG(startIdx:endIdx);
        [pxx, F] = pwelch(segment, hamming(512), [], [], Fs);

        if isempty(F_ref)
            F_ref = F;  % set reference frequencies
        end

        if length(pxx) ~= length(F_ref)
            warning('Skipping epoch %d: mismatched PSD size.', i);
            continue;
        end

        allPxx(:, end+1) = pxx;
    end
end

% Average PSD
meanPxx = mean(allPxx, 2);

% --- Plot PSD ---
fig = figure('Visible','off');
semilogy(F_ref, meanPxx, 'k', 'LineWidth', 1.5);
xline(50, '--r', '50 Hz Notch', 'LabelOrientation','horizontal', ...
      'LabelVerticalAlignment','bottom');
grid on;
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (uV^2/Hz)');
title(['PSD - ' stateName]);
xlim([0 100]);
saveas(fig, fullfile(outputFolder, [stateName '_PSD.png']));
close(fig);

% --- Compute relative power per band ---
totalPower = trapz(F_ref, meanPxx);

bandPowers = struct();
bandNames = fieldnames(bands);
for b = 1:length(bandNames)
    bandName = bandNames{b};
    range = bands.(bandName);
    if strcmp(bandName, 'lGamma')
        idx = (F_ref >= 30 & F_ref < 49.5) | (F_ref > 50.5 & F_ref <= 60);
    else
        idx = F_ref >= range(1) & F_ref <= range(2);
    end
    bandPower = trapz(F_ref(idx), meanPxx(idx));
    bandPowers.(bandName) = 100 * (bandPower / totalPower);
end

% --- Save results ---
txtFile = fullfile(outputFolder, [stateName '_band_power.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, 'Relative band power for %s\n', stateName);
for b = 1:length(bandNames)
    fprintf(fid, '%s: %.2f %%\n', bandNames{b}, bandPowers.(bandNames{b}));
end
fclose(fid);

disp('EEG power and band analysis complete.');

