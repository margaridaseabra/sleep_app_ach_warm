
%% Loading the data
% Load your .mat file
data = load('/Users/margaridaseabra/Downloads/scored files/20251001-baseline-mouse1_scored.mat');

% Extract relevant variables
EEG_rawtrace = data.eeg;   % Continuous EEG data
EMG_rawtrace = data.emg;   % Continuous EMG data
sleep_score = data.sleep_scores;  % Sleep scoring (e.g., 0 = Wake, 1 = NREM, 2 = REM)
sampling_freq = data.eeg_frequency;  % Sampling frequency (assuming this is the correct variable)
start_time = data.start_time;  % Start time (this may be useful for aligning or cropping)


%% Choosing a specific time window if necessary
% Define the time window (in seconds)
optional_crop = false;  % Set to true to enable cropping, false to use full data
if optional_crop
    start_time = 0;  % e.g., starting at 600 seconds
    end_time = 200000;   % e.g., ending at 1200 seconds

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

%% Plotting the EEG and EMG
% Identify gaps in the sleep scores and replace them with NaN
% Assuming a value of 0 means no scoring (can be adjusted based on your data)
sleep_score_cropped(isnan(sleep_score_cropped)) = NaN;

% Generate the time vector to match the length of the cropped data
time_vector = (0:length(EEG_cropped)-1) / sampling_freq;

% Plot EEG and EMG with sleep states on top
figure;

% Plot EEG Signal
subplot(3,1,1);
plot(time_vector, EEG_cropped);
title('EEG Signal');
xlabel('Time (s)');
ylabel('EEG Amplitude');
hold on;

% Highlight sleep states on EEG signal (use shaded areas)
% Wake = 0 (red), NREM = 1 (green), REM = 2 (blue)
wake_idx = find(sleep_score_cropped == 0);
nrem_idx = find(sleep_score_cropped == 1);
rem_idx = find(sleep_score_cropped == 2);

% Highlight the sleep states on the EEG signal
fill([time_vector(wake_idx)'; flipud(time_vector(wake_idx)')], ...
    [ones(length(wake_idx), 1) * max(EEG_cropped)'; ones(length(wake_idx), 1) * min(EEG_cropped)'], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Wake in red
fill([time_vector(nrem_idx)'; flipud(time_vector(nrem_idx)')], ...
    [ones(length(nrem_idx), 1) * max(EEG_cropped)'; ones(length(nrem_idx), 1) * min(EEG_cropped)'], ...
    'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % NREM in green
fill([time_vector(rem_idx)'; flipud(time_vector(rem_idx)')], ...
    [ones(length(rem_idx), 1) * max(EEG_cropped)'; ones(length(rem_idx), 1) * min(EEG_cropped)'], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % REM in blue

% Plot EMG Signal
subplot(3,1,2);
plot(time_vector, EMG_cropped);
title('EMG Signal');
xlabel('Time (s)');
ylabel('EMG Amplitude');
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

% Plot Sleep Stage Transitions (0 = Wake, 1 = NREM, 2 = REM)
subplot(3,1,3);
plot(time_vector, sleep_score_cropped, 'LineWidth', 2);
title('Sleep Stage Transitions');
xlabel('Time (s)');
ylabel('Sleep State');
ylim([-0.5, 2.5]);  % Adjust the y-axis to fit the scoring range



%% Microarousals Detection

% Define microarousal duration threshold 
MA_maxdur = 15;  % in seconds

% Find wake bouts
wake_bouts = find(sleep_score_cropped == 0);  % Sleep score 1 is wakefulness

% Calculate bout durations
bout_durations = diff(wake_bouts);

% Identify microarousals (short wake bouts)
MA_idx = find(bout_durations < MA_maxdur);  % Short bouts of wakefulness
MA_onset = wake_bouts(MA_idx);
MA_duration = bout_durations(MA_idx);

%% Sleep stages Transitions
% identifies transitions between sleep stages and creates a state
% transition matrix
wake_to_nrem = find(diff(sleep_score_cropped == 0) == 1);  % Wake to NREM
nrem_to_rem = find(diff(sleep_score_cropped == 1) == 2);  % NREM to REM

% Store onset and offset for each transition
wake_to_nrem_onset = wake_to_nrem;
nrem_to_rem_onset = nrem_to_rem;

%% Analysis of Sleep Architecture -- UNDER CONSTRUCTION
% Find the onsets and offsets of wake, NREM, and REM bouts
wake_onsets = find(diff(sleep_score_cropped == 0) == 1);  % Transition to wake (0)
wake_offsets = find(diff(sleep_score_cropped == 0) == -1);  % Transition from wake (0)

nrem_onsets = find(diff(sleep_score_cropped == 1) == 1);  % Transition to NREM (1)
nrem_offsets = find(diff(sleep_score_cropped == 1) == -1);  % Transition from NREM (1)

rem_onsets = find(diff(sleep_score_cropped == 2) == 1);  % Transition to REM (2)
rem_offsets = find(diff(sleep_score_cropped == 2) == -1);  % Transition from REM (2)

% Calculate the duration of each bout
wake_durations = wake_offsets - wake_onsets;  % Durations of wake bouts
nrem_durations = nrem_offsets - nrem_onsets;  % Durations of NREM bouts
rem_durations = rem_offsets - rem_onsets;    % Durations of REM bouts

% Total time spent in each state (in seconds)
total_wake_time = sum(wake_durations);
total_nrem_time = sum(nrem_durations);
total_rem_time = sum(rem_durations);

% Total duration of the analysis window (in seconds)
total_time = length(sleep_score_cropped);

% Proportion of time spent in each state
wake_proportion = total_wake_time / total_time;
nrem_proportion = total_nrem_time / total_time;
rem_proportion = total_rem_time / total_time;



%% Spectral Analysis
% Compute PSD for EEG during different sleep states (e.g., Wake, NREM, REM)
fs = sampling_freq;  % Sampling frequency

% Extract segments of EEG corresponding to each state (wake, NREM, REM)
eeg_wake = EEG_cropped(sleep_score_cropped == 0);
eeg_nrem = EEG_cropped(sleep_score_cropped == 1);
eeg_rem = EEG_cropped(sleep_score_cropped == 2);

% Compute PSD for each segment
[pxx_wake, f_wake] = pwelch(eeg_wake, [], [], [], fs);
[pxx_nrem, f_nrem] = pwelch(eeg_nrem, [], [], [], fs);
[pxx_rem, f_rem] = pwelch(eeg_rem, [], [], [], fs);

% Plot PSDs
figure;
subplot(3,1,1);
plot(f_wake, 10*log10(pxx_wake));
title('Wake EEG Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');

subplot(3,1,2);
plot(f_nrem, 10*log10(pxx_nrem));
title('NREM EEG Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');

subplot(3,1,3);
plot(f_rem, 10*log10(pxx_rem));
title('REM EEG Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');

%% Visualizations
% Plot the duration of each bout for each state
figure;
subplot(3,1,1);
histogram(wake_durations, 20);  % Histogram of wake bout durations
title('Wake Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

subplot(3,1,2);
histogram(nrem_durations, 20);  % Histogram of NREM bout durations
title('NREM Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

subplot(3,1,3);
histogram(rem_durations, 20);  % Histogram of REM bout durations
title('REM Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

% Plot the total proportion of time spent in each state
figure;
bar([wake_proportion, nrem_proportion, rem_proportion]);
set(gca, 'XTickLabel', {'Wake', 'NREM', 'REM'});
title('Proportion of Time Spent in Each State');
ylabel('Proportion');

% Create a figure for plotting EEG and EMG signals aligned with sleep states
% Plot EEG and EMG with sleep states
figure;

% Plot EEG Signal
subplot(3,1,1);
plot(time_vector, EEG_cropped);
title('EEG Signal');
xlabel('Time (s)');
ylabel('EEG Amplitude');

% Plot EMG Signal
subplot(3,1,2);
plot(time_vector, EMG_cropped);
title('EMG Signal');
xlabel('Time (s)');
ylabel('EMG Amplitude');

% Plot Sleep Stage Transitions (0 = Wake, 1 = NREM, 2 = REM)
subplot(3,1,3);
plot(time_vector, sleep_score_cropped, 'LineWidth', 2);
title('Sleep Stage Transitions');
xlabel('Time (s)');
ylabel('Sleep State');


% Plot the duration of each bout for each state
figure;
subplot(3,1,1);
histogram(wake_durations, 20);  % Histogram of wake bout durations
title('Wake Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

subplot(3,1,2);
histogram(nrem_durations, 20);  % Histogram of NREM bout durations
title('NREM Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

subplot(3,1,3);
histogram(rem_durations, 20);  % Histogram of REM bout durations
title('REM Bout Durations');
xlabel('Duration (s)');
ylabel('Frequency');

% Plot the total proportion of time spent in each state
figure;
bar([wake_proportion, nrem_proportion, rem_proportion]);
set(gca, 'XTickLabel', {'Wake', 'NREM', 'REM'});
title('Proportion of Time Spent in Each State');
ylabel('Proportion');
