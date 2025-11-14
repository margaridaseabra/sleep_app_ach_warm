% Load your .mat files
data1 = load('/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/20251001_mouse1_baseline_automaticallynotched.mat');
data2 = load('/Users/margaridaseabra/Library/CloudStorage/OneDrive-UniversityofCopenhagen/scored files/20251007_ambtemp_mouse4_missingnotch.mat');

signal_notched = data1.eeg;
signal_raw = data2.eeg;
Fs = data2.eeg_frequency;


% Design 50 Hz notch filter
f0   = 50;                   % notch frequency (Hz)
Q    = 35;                   % quality factor (adjust as needed)
f0_norm = f0 / (Fs/2);       % normalized (0‑1) frequency relative to Nyquist

% Get coefficients directly
[B, A] = designNotchPeakIIR( ...
    "Response", "notch", ...
    "CenterFrequency", f0_norm, ...
    "QualityFactor", Q );

% Now apply zero‑phase filtering
signal_filtered = filtfilt(B, A, signal_raw);

% Find the minimum length
L = min([length(signal_raw), length(signal_filtered), length(signal_notched)]);

% Trim all to same length
t = (0:L-1)/Fs;
signal_raw      = signal_raw(1:L);
signal_filtered = signal_filtered(1:L);
signal_notched  = signal_notched(1:L);

% Plot
figure;
subplot(3,1,1); plot(t, signal_raw);     title('Raw Signal');
subplot(3,1,2); plot(t, signal_filtered); title('Filtered (50 Hz Notch)');
subplot(3,1,3); plot(t, signal_notched);  title('Auto-Notched Signal');



% Export filtered signal to .mat file
filtered_signal = signal_filtered;
save('20251007_mouse4_filtered_50Hznotch.mat', 'filtered_signal', 'Fs');