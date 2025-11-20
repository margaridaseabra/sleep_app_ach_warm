function explore_REM_EMG_bursts(mat_file, REM_CODE, maxREM)
% explore_REM_EMG_bursts(mat_file, REM_CODE, maxREM)
%
% Inside REM sleep:
%   - detect EMG bursts (using EMG envelope threshold)
%   - for EACH REM bout, plot:
%       1) EEG trace
%       2) EMG trace
%       3) EMG envelope (with threshold + burst markers in that bout)
%   - also compute a mean peri-burst spectrogram (kept from before)
%
% INPUTS
%   mat_file : .mat file with eeg, emg, eeg_frequency, sleep_scores
%   REM_CODE : value in sleep_scores that means REM (default = 2)
%   maxREM   : (optional) max number of REM bouts to plot (default = Inf)
%
% EXPECTED VARIABLES IN MAT:
%   eeg           : [1 x N] or [N x 1], EEG
%   emg           : [1 x N] or [N x 1], EMG
%   eeg_frequency : scalar, Fs in Hz
%   sleep_scores  : [1 x T] or [T x 1], 1 value per second (0/1/2 etc.)

    % ---------- inputs ----------
    if nargin < 1 || isempty(mat_file)
        [f,p] = uigetfile('*.mat','Select .mat file');
        if isequal(f,0); return; end
        mat_file = fullfile(p,f);
    end
    if nargin < 2 || isempty(REM_CODE)
        REM_CODE = 2;   % <-- change if your REM label is different
    end
    if nargin < 3 || isempty(maxREM)
        maxREM = Inf;   % no limit by default
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

    % Map 1-second scores to per-sample state
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
    d_rem  = diff([0 isREM 0]);
    rem_start = find(d_rem ==  1);   % sample indices
    rem_end   = find(d_rem == -1) - 1;

    nREM = numel(rem_start);
    if nREM == 0
        warning('No REM bouts found (REM_CODE = %d).', REM_CODE);
        return;
    end

    rem_durs = (rem_end - rem_start + 1) / Fs;
    fprintf('Found %d REM bouts. Median duration = %.1f s (range %.1f–%.1f s)\n', ...
        nREM, median(rem_durs), min(rem_durs), max(rem_durs));

    % ---------- EMG envelope and burst detection ----------
    % band-pass EMG to 20–300 Hz (muscle band)
    bpFilt = designfilt('bandpassiir', ...
                        'FilterOrder',4, ...
                        'HalfPowerFrequency1',20, ...
                        'HalfPowerFrequency2',300, ...
                        'SampleRate',Fs);
    emg_bp = filtfilt(bpFilt, EMG);

    emg_env = abs(emg_bp);
    win_env = round(0.1 * Fs);           % 100 ms smoothing
    emg_env = movmean(emg_env, win_env);

    % threshold based on REM samples only
    rem_env = emg_env(isREM);
    mu_rem  = mean(rem_env);
    sd_rem  = std(rem_env);
    thr     = mu_rem + 2.5*sd_rem;       % adjust factor if needed

    fprintf('EMG env in REM: mean = %.3f, SD = %.3f, thr = %.3f\n', ...
        mu_rem, sd_rem, thr);

    % burst detection in REM
    isBurst = (emg_env > thr) & isREM;
    d_b     = diff([0 isBurst 0]);
    burst_start = find(d_b ==  1);
    burst_end   = find(d_b == -1) - 1;

    min_dur = 0.1;  % seconds
    good = (burst_end - burst_start + 1) >= min_dur*Fs;
    burst_start = burst_start(good);
    burst_end   = burst_end(good);

    nBursts = numel(burst_start);
    fprintf('Detected %d EMG bursts in REM (>= %.1f s)\n', nBursts, min_dur);
    if nBursts == 0
        warning('No bursts found; try lowering threshold or min_dur.');
    end

    % ---------- 1) VISUALIZE ALL REM BOUTS ----------
    nToPlot = min(nREM, maxREM);
    fprintf('Plotting %d REM bouts (maxREM = %g)\n', nToPlot, maxREM);

    for r = 1:nToPlot
        i1 = rem_start(r);
        i2 = rem_end(r);

        tt      = t_eeg(i1:i2);
        eeg_seg = EEG(i1:i2);
        emg_seg = EMG(i1:i2);
        env_seg = emg_env(i1:i2);

        % bursts within this REM bout
        b_in_seg = burst_start(burst_start >= i1 & burst_start <= i2);

        figure('Color','w', 'Name', sprintf('REM bout %d', r));

        % EEG
        subplot(3,1,1);
        plot(tt, eeg_seg);
        xlabel('Time (s)');
        ylabel('EEG');
        title(sprintf('REM bout %d (%.1f s) – EEG', r, rem_durs(r)));
        grid on;

        % raw EMG
        subplot(3,1,2);
        plot(tt, emg_seg);
        xlabel('Time (s)');
        ylabel('EMG');
        title('Raw EMG');
        grid on;

        % EMG envelope with threshold and burst markers
        subplot(3,1,3);
        plot(tt, env_seg); hold on;
        yline(thr, 'r--', 'Threshold');
        if ~isempty(b_in_seg)
            plot(t_eeg(b_in_seg), emg_env(b_in_seg), ...
                 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Burst onset');
        end
        xlabel('Time (s)');
        ylabel('EMG env');
        title('EMG envelope (REM) with bursts');
        grid on;

        drawnow;
    end

    % ---------- 2) KEEP THE MEAN PERI-BURST SPECTROGRAM (optional) ----------
    if nBursts > 0
        win_sec  = 5;                    % +/- 5 s around burst
        win_samp = round(win_sec * Fs);

        spec_win   = round(0.5*Fs);      % 0.5 s window
        spec_nover = round(0.25*Fs);
        spec_nfft  = max(2^nextpow2(spec_win), spec_win);

        % example burst to get time/freq grid
        first_burst = burst_start(1);
        i1 = max(1, first_burst - win_samp);
        i2 = min(Nsamples, first_burst + win_samp);
        [S_ex,F_spec,T_spec] = spectrogram(EEG(i1:i2), ...
                                           spec_win, spec_nover, spec_nfft, Fs);
        T_spec = T_spec - win_sec;       % roughly center at 0
        valid_t = (T_spec >= -win_sec) & (T_spec <= win_sec);
        T_spec = T_spec(valid_t);

        F_lim = F_spec(F_spec <= 50);
        nF    = numel(F_lim);
        nT    = numel(T_spec);
        S_stack = zeros(nF, nT, nBursts);

        for b = 1:nBursts
            center = burst_start(b);
            i1 = max(1, center - win_samp);
            i2 = min(Nsamples, center + win_samp);
            eeg_win = EEG(i1:i2);

            [S_b,F_b,T_b] = spectrogram(eeg_win, spec_win, spec_nover, spec_nfft, Fs);
            idxF = F_b <= 50;
            S_b  = S_b(idxF,:);
            F_b  = F_b(idxF);

            % re-time relative to burst center
            T_b = T_b + ( (i1-1)/Fs ) - (center-1)/Fs;

            % map to common T_spec grid (nearest)
            [~, idxT] = arrayfun(@(x) min(abs(T_b - x)), T_spec);
            S_b = S_b(:, idxT);

            S_stack(:,:,b) = S_b;
        end

        S_mean = 10*log10(mean(abs(S_stack).^2, 3) + eps);

        figure('Color','w','Name','Mean peri-burst spectrogram (EEG)');
        imagesc(T_spec, F_lim, S_mean);
        axis xy;
        xlabel('Time from EMG burst (s)');
        ylabel('Frequency (Hz)');
        title('Mean EEG spectrogram aligned to EMG bursts during REM');
        colorbar;
    end
end
