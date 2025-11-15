%% Load the data
data = load('/Users/margaridaseabra/Downloads/scored files/20251001-baseline-mouse1_scored.mat');

% Extract relevant variables
EEG_rawtrace = data.eeg;  % Continuous EEG data
EMG_rawtrace = data.emg;  % Continuous EMG data
sleep_score = data.sleep_scores;  % Sleep scoring (0 = Wake, 1 = NREM, 2 = REM)
sampling_freq = data.eeg_frequency;  % Sampling frequency
start_time = data.start_time;  % Start time

%% Choosing a specific time window (if needed)
% If you want to crop, enable this section; otherwise, it uses the full data.
optional_crop = false;  % Set to true to enable cropping, false to use full data

if optional_crop
    start_time = 0;  % Define your start time in seconds
    end_time = 1200;  % Define your end time in seconds

    % Convert times to sample indices
    start_idx = round(start_time * sampling_freq);
    end_idx = round(end_time * sampling_freq);

    % Crop the EEG and EMG data
    EEG_cropped = EEG_rawtrace(start_idx:end_idx);
    EMG_cropped = EMG_rawtrace(start_idx:end_idx);

    % Crop the sleep score vector as well
    sleep_score_cropped = sleep_score(start_idx:end_idx);
else
    % Use full data if cropping is not enabled
    EEG_cropped = EEG_rawtrace;
    EMG_cropped = EMG_rawtrace;
    sleep_score_cropped = sleep_score;
end

%% Generate time vector (to match the cropped data)
time_vector = (0:length(EEG_cropped)-1) / sampling_freq;

%% Plot EEG Spectrogram (Top)
% Create the spectrogram for EEG using a short-time Fourier transform (STFT)
figure;

subplot(4,1,1);  % Spectrogram on the top
[S, F, T] = spectrogram(EEG_cropped, 256, 250, 256, sampling_freq);
imagesc(T, F, 10*log10(abs(S)));
axis xy;
title('EEG Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
colormap jet;
caxis([-115 -95]);  % Set color limits for better contrast

%% Plot EEG Signal (Second)
% Plot the EEG signal with highlighted sleep states
subplot(4,1,2);  % EEG signal plot
plot(time_vector, EEG_cropped, 'k');
hold on;

% Highlight the sleep states on EEG signal
wake_idx = find(sleep_score_cropped == 0);
nrem_idx = find(sleep_score_cropped == 1);
rem_idx = find(sleep_score_cropped == 2);

fill([time_vector(wake_idx)'; flipud(time_vector(wake_idx)')], ...
    [ones(length(wake_idx), 1) * max(EEG_cropped)'; ones(length(wake_idx), 1) * min(EEG_cropped)'], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Wake in red
fill([time_vector(nrem_idx)'; flipud(time_vector(nrem_idx)')], ...
    [ones(length(nrem_idx), 1) * max(EEG_cropped)'; ones(length(nrem_idx), 1) * min(EEG_cropped)'], ...
    'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % NREM in green
fill([time_vector(rem_idx)'; flipud(time_vector(rem_idx)')], ...
    [ones(length(rem_idx), 1) * max(EEG_cropped)'; ones(length(rem_idx), 1) * min(EEG_cropped)'], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % REM in blue

title('EEG Signal');
xlabel('Time (s)');
ylabel('EEG Amplitude');

%% Plot EMG Signal (Third)
% Plot the EMG signal with highlighted sleep states
subplot(4,1,3);  % EMG signal plot
plot(time_vector, EMG_cropped, 'k');
hold on;

% Highlight the sleep states on EMG signal
fill([time_vector(wake_idx)'; flipud(time_vector(wake_idx)')], ...
    [ones(length(wake_idx), 1) * max(EMG_cropped)'; ones(length(wake_idx), 1) * min(EMG_cropped)'], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Wake in red
fill([time_vector(nrem_idx)'; flipud(time_vector(nrem_idx)')], ...
    [ones(length(nrem_idx), 1) * max(EMG_cropped)'; ones(length(nrem_idx), 1) * min(EMG_cropped)'], ...
    'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % NREM in green
fill([time_vector(rem_idx)'; flipud(time_vector(rem_idx)')], ...
    [ones(length(rem_idx), 1) * max(EMG_cropped)'; ones(length(rem_idx), 1) * min(EMG_cropped)'], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % REM in blue

title('EMG Signal');
xlabel('Time (s)');
ylabel('EMG Amplitude');

%% Plot Sleep Stage Transitions (Bottom)
subplot(4,1,4);  % Sleep Stage Transitions
plot(time_vector, sleep_score_cropped, 'LineWidth', 2);
title('Sleep Stage Transitions');
xlabel('Time (s)');
ylabel('Sleep State');
ylim([-0.5, 2.5]);  % Adjust the y-axis to fit the scoring range
yticks([0 1 2]);
yticklabels({'Wake', 'NREM', 'REM'});

