function plot_REM_episodes_from_mat(mat_file, REM_CODE)
% plot_REM_episodes_from_mat(mat_file, REM_CODE)
%
% For each REM episode in a scored recording, plot:
%   (1) EEG trace
%   (2) EMG trace
%   (3) EEG PSD for that REM episode
%
% Designed for files like:
%   20251002-baseline-mouse2_newlyscored.mat
%
% Expected variables inside MAT:
%   eeg           : [1 x N] or [N x 1] double
%   emg           : [1 x N] or [N x 1] double
%   eeg_frequency : scalar sampling rate (Hz)
%   sleep_scores  : [1 x T] or [T x 1] double, 1 value per second
%
% Sleep-score codes:
%   We assume REM is coded as REM_CODE (default = 2).
%   Change REM_CODE if your scoring uses a different label.

    % ---------- inputs & defaults ----------
    if nargin < 1 || isempty(mat_file)
        [f,p] = uigetfile('*.mat', 'Select .mat file with EEG/EMG/sleep_scores');
        if isequal(f,0)
            disp('User cancelled.');
            return;
        end
        mat_file = fullfile(p,f);
    end

    if nargin < 2 || isempty(REM_CODE)
        REM_CODE = 2;   % <-- change here if REM is 0 or 1
    end

    % ---------- load data ----------
    S = load(mat_file);

    required_vars = {'eeg','emg','eeg_frequency','sleep_scores'};
    for i = 1:numel(required_vars)
        if ~isfield(S, required_vars{i})
            error('MAT file is missing required variable: %s', required_vars{i});
        end
    end

    EEG          = S.eeg(:)';           % force row
    EMG          = S.emg(:)';           % force row
    Fs           = S.eeg_frequency;     % scalar
    sleep_scores = S.sleep_scores(:)';  % 1 value per second

    if ~isscalar(Fs)
        error('eeg_frequency must be a scalar.');
    end

    Nsamples = numel(EEG);
    if numel(EMG) ~= Nsamples
        error('EEG and EMG must have same number of samples.');
    end

    Nsec = numel(sleep_scores);   % number of scored seconds

    % ---------- sanity check of alignment ----------
    est_N = round(Nsec * Fs);
    fprintf('EEG samples: %d, Nsec*Fs ≈ %d (diff = %d samples)\n', ...
            Nsamples, est_N, Nsamples - est_N);

    % Time vectors
    t_eeg = (0:Nsamples-1) / Fs;        % seconds
    t_scores = 0:(Nsec-1);              % second indices for scoring

    % ---------- find REM episodes in score domain ----------
    isREM_sec = (sleep_scores == REM_CODE);

    d = diff([0, isREM_sec, 0]);    % edges
    rem_start_sec_idx = find(d ==  1);   % indices into sleep_scores
    rem_end_sec_idx   = find(d == -1) - 1;

    nREM = numel(rem_start_sec_idx);
    if nREM == 0
        warning('No REM episodes found with REM_CODE = %d.', REM_CODE);
        return;
    end

    % Print durations
    rem_durs_sec = rem_end_sec_idx - rem_start_sec_idx + 1;
    fprintf('Found %d REM episodes. Median duration = %.1f s (range %.1f–%.1f s)\n', ...
        nREM, median(rem_durs_sec), min(rem_durs_sec), max(rem_durs_sec));

    % ---------- PSD parameters ----------
    % We'll adapt the window length to the segment length later.
    base_win_sec   = 4;          % nominal Welch window length (s)
    base_win_samps = round(base_win_sec * Fs);
    fmax_plot = 50;              % max frequency for PSD plot (Hz)

    % ---------- loop over REM episodes ----------
    for k = 1:nREM
        % ----- REM episode in seconds -----
        sec1 = rem_start_sec_idx(k);   % 1-based second index
        sec2 = rem_end_sec_idx(k);     % inclusive

        % Convert seconds -> sample indices
        % sec 1 covers [0,1) s, sec 2 covers [1,2) s, etc.
        samp1 = floor((sec1-1) * Fs) + 1;
        samp2 = min(round(sec2 * Fs), Nsamples);

        eeg_seg = EEG(samp1:samp2);
        emg_seg = EMG(samp1:samp2);
        t_seg   = t_eeg(samp1:samp2);

        dur_seg = t_seg(end) - t_seg(1);

        % Skip very short REMs if you want (optional)
        % if dur_seg < 3
        %     continue;
        % end

        % ----- PSD for this segment -----
        seg_len = numel(eeg_seg);
        win   = min(seg_len, base_win_samps);
        if win < 4  % too short, skip
            continue;
        end
        nover = floor(0.5 * win);
        nfft  = max(2^nextpow2(win), win);

        [Pxx, F] = pwelch(eeg_seg, win, nover, nfft, Fs, 'onesided');
        idx_f = F <= fmax_plot;

        % ---------- plotting ----------
        figure('Name', sprintf('REM %d', k), 'Color', 'w');

        % (1) EEG trace
        subplot(3,1,1);
        plot(t_seg, eeg_seg);
        xlabel('Time (s)');
        ylabel('EEG (a.u.)');
        title(sprintf('REM episode %d – EEG (%.1f s)', k, dur_seg));
        grid on;

        % (2) EMG trace
        subplot(3,1,2);
        plot(t_seg, emg_seg);
        xlabel('Time (s)');
        ylabel('EMG (a.u.)');
        title('EMG');
        grid on;

        % (3) PSD
        subplot(3,1,3);
        plot(F(idx_f), 10*log10(Pxx(idx_f)));
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        title('EEG PSD during REM episode');
        grid on;

        % If you want to automatically save each figure, uncomment:
        % [p,f,~] = fileparts(mat_file);
        % outname = sprintf('%s_REM%02d.png', f, k);
        % exportgraphics(gcf, fullfile(p, outname), 'Resolution', 300);

        drawnow;
    end
end
