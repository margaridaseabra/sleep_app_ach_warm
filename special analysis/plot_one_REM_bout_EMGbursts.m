function plot_one_REM_bout_EMGbursts(mat_file, REM_CODE, bout_idx)
% plot_one_REM_bout_EMGbursts(mat_file, REM_CODE, bout_idx)
%
% Plot one specific REM bout with:
%   (1) EEG trace
%   (2) EMG trace
%   (3) EMG envelope with threshold + EMG-burst markers
%   (4) EEG PSD for that REM bout
%
% Designed for files like: 20251002-baseline-mouse2_newlyscored.mat
%
% EXPECTED VARIABLES IN MAT:
%   eeg           : [1 x N] or [N x 1], EEG
%   emg           : [1 x N] or [N x 1], EMG
%   eeg_frequency : scalar, Fs in Hz
%   sleep_scores  : [1 x T] or [T x 1], 1 value per second (0/1/2 etc.)
%
% INPUTS
%   mat_file : string, path to .mat
%   REM_CODE : value in sleep_scores that means REM (default = 2)
%   bout_idx : index of the REM bout you want to plot (1-based)
%

    % ---------- inputs ----------
    if nargin < 1 || isempty(mat_file)
        [f,p] = uigetfile('*.mat','Select .mat file');
        if isequal(f,0); return; end
        mat_file = fullfile(p,f);
    end
    if nargin < 2 || isempty(REM_CODE)
        REM_CODE = 2;   % <-- change if REM is coded differently
    end
    if nargin < 3
        error('You must provide bout_idx (which REM bout to plot).');
    end

    % ---------- load data ----------
    S = load(mat_file);
    EEG          = S.eeg(:)';          % row
    EMG          = S.emg(:)';          % row
    Fs           = S.eeg_frequency;
    sleep_scores = S.sleep_scores(:)'; % row, 1 / second

    Nsamples = numel(EEG);
    if numel(EMG) ~= Nsamples
        error('EEG and EMG must have same length.');
    end
    Nsec   = numel(sleep_scores);
    t_eeg  = (0:Nsamples-1)/Fs;

    % ---------- map 1-second scores to per-sample state ----------
    est_N = round(Nsec * Fs);
    fprintf('EEG samples = %d, Nsec*Fs ≈ %d (diff = %d)\n', ...
        Nsamples, est_N, Nsamples-est_N);

    score_per_sample = zeros(1, Nsamples);
    for s = 1:Nsec
        i1 = floor((s-1)*Fs) + 1;
        i2 = min(round(s*Fs), Nsamples);
        score_per_sample(i1:i2) = sleep_scores(s);
    end

    isREM = (score_per_sample == REM_CODE);

    % ---------- find REM bouts (in samples) ----------
    d_rem    = diff([0 isREM 0]);
    rem_start = find(d_rem ==  1);   % sample indices
    rem_end   = find(d_rem == -1) - 1;

    nREM = numel(rem_start);
    if nREM == 0
        error('No REM bouts found (REM_CODE = %d).', REM_CODE);
    end

    rem_durs = (rem_end - rem_start + 1) / Fs;

    if bout_idx < 1 || bout_idx > nREM
        error('bout_idx must be between 1 and %d (you asked for %d).', ...
              nREM, bout_idx);
    end

    fprintf('Total REM bouts: %d\n', nREM);
    fprintf('Bout %d duration = %.1f s\n', bout_idx, rem_durs(bout_idx));

    % ---------- EMG envelope + burst detection (global, based on REM) ----------
    bpFilt = designfilt('bandpassiir', ...
                        'FilterOrder',4, ...
                        'HalfPowerFrequency1',20, ...
                        'HalfPowerFrequency2',300, ...
                        'SampleRate',Fs);
    emg_bp = filtfilt(bpFilt, EMG);

    emg_env = abs(emg_bp);
    win_env = round(0.1 * Fs);         % 100 ms smoothing
    emg_env = movmean(emg_env, win_env);

    rem_env = emg_env(isREM);
    mu_rem  = mean(rem_env);
    sd_rem  = std(rem_env);
    thr     = mu_rem + 2.5*sd_rem;     % tweak factor if needed

    fprintf('EMG env in REM: mean = %.3f, SD = %.3f, thr = %.3f\n', ...
        mu_rem, sd_rem, thr);

    isBurst = (emg_env > thr) & isREM;
    d_b     = diff([0 isBurst 0]);
    burst_start = find(d_b ==  1);
    burst_end   = find(d_b == -1) - 1;

    min_dur = 0.1; % seconds
    good = (burst_end - burst_start + 1) >= min_dur*Fs;
    burst_start = burst_start(good);
    burst_end   = burst_end(good);

    nBursts = numel(burst_start);
    fprintf('Detected %d EMG bursts in REM (>= %.1f s)\n', nBursts, min_dur);

    % ---------- extract THIS REM bout ----------
    i1 = rem_start(bout_idx);
    i2 = rem_end(bout_idx);

    tt      = t_eeg(i1:i2);
    eeg_seg = EEG(i1:i2);
    emg_seg = EMG(i1:i2);
    env_seg = emg_env(i1:i2);

    % bursts within this REM bout
    b_in_seg = burst_start(burst_start >= i1 & burst_start <= i2);

    % ---------- compute PSD of EEG for this REM bout ----------
    seg_len = numel(eeg_seg);
    win   = min(seg_len, round(4*Fs));    % 4 s window or shorter
    nover = floor(0.5 * win);
    nfft  = max(2^nextpow2(win), win);

    [Pxx, F] = pwelch(eeg_seg, win, nover, nfft, Fs, 'onesided');
    fmax_plot = 50;
    idx_f = F <= fmax_plot;

    % ---------- plot ----------
    figure('Color','w', 'Name', sprintf('REM bout %d', bout_idx));

    % 1) EEG
    subplot(4,1,1);
    plot(tt, eeg_seg);
    xlabel('Time (s)');
    ylabel('EEG');
    title(sprintf('REM bout %d (%.1f s) – EEG', bout_idx, rem_durs(bout_idx)));
    grid on;

    % 2) EMG
    subplot(4,1,2);
    plot(tt, emg_seg);
    xlabel('Time (s)');
    ylabel('EMG');
    title('Raw EMG');
    grid on;

    % 3) EMG envelope with threshold + bursts
    subplot(4,1,3);
    plot(tt, env_seg); hold on;
    yline(thr, 'r--', 'Threshold');
    if ~isempty(b_in_seg)
        plot(t_eeg(b_in_seg), emg_env(b_in_seg), ...
             'mo','MarkerFaceColor','m','DisplayName','Burst onset');
        legend('show');
    end
    xlabel('Time (s)');
    ylabel('EMG env');
    title('EMG envelope (REM) with bursts');
    grid on;

    % 4) EEG PSD
    subplot(4,1,4);
    plot(F(idx_f), 10*log10(Pxx(idx_f)));
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('EEG PSD for this REM bout');
    grid on;

end
