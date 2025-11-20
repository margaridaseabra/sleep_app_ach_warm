function analyze_REM_burst_locked(mat_file, REM_CODE, bout_idx)
% analyze_REM_burst_locked(mat_file, REM_CODE, bout_idx)
%
% For ONE REM bout:
%   - detect EMG bursts inside that bout
%   - show EMG envelope + burst markers over the whole bout
%   - compute burst-locked averages:
%       * mean EEG trace around bursts
%       * mean EEG spectrogram around bursts
%   - compare pre vs post burst PSD (within this bout)
%
% INPUTS
%   mat_file : .mat file with fields:
%                eeg, emg, eeg_frequency, sleep_scores
%   REM_CODE : value in sleep_scores that == REM (default = 2)
%   bout_idx: index of REM bout to analyze (1-based)
%
% Assumes sleep_scores is one score per second.

    % ---------------- inputs ----------------
    if nargin < 1 || isempty(mat_file)
        [f,p] = uigetfile('*.mat','Select .mat file'); 
        if isequal(f,0); return; end
        mat_file = fullfile(p,f);
    end
    if nargin < 2 || isempty(REM_CODE)
        REM_CODE = 2;
    end
    if nargin < 3
        error('You must provide bout_idx.');
    end

    % ---------------- load data ----------------
    S = load(mat_file);
    EEG          = S.eeg(:)'; 
    EMG          = S.emg(:)'; 
    Fs           = S.eeg_frequency;
    sleep_scores = S.sleep_scores(:)';

    Nsamples = numel(EEG);
    Nsec     = numel(sleep_scores);
    t_eeg    = (0:Nsamples-1)/Fs;

    % --- map scores (1 per s) to per-sample labels ---
    score_per_sample = zeros(1,Nsamples);
    for s = 1:Nsec
        i1 = floor((s-1)*Fs) + 1;
        i2 = min(round(s*Fs), Nsamples);
        score_per_sample(i1:i2) = sleep_scores(s);
    end

    isREM = (score_per_sample == REM_CODE);

    % --- find REM bouts ---
    d_rem    = diff([0 isREM 0]);
    rem_start = find(d_rem ==  1);
    rem_end   = find(d_rem == -1) - 1;

    nREM = numel(rem_start);
    if nREM == 0
        error('No REM bouts found (REM_CODE = %d).',REM_CODE);
    end
    if bout_idx < 1 || bout_idx > nREM
        error('bout_idx must be between 1 and %d.', nREM);
    end

    i1 = rem_start(bout_idx);
    i2 = rem_end(bout_idx);

    eeg_seg = EEG(i1:i2);
    emg_seg = EMG(i1:i2);
    t_seg   = t_eeg(i1:i2);
    dur_seg = t_seg(end) - t_seg(1);

    fprintf('Analyzing REM bout %d (%.1f s)\n', bout_idx, dur_seg);

    % ==========================================================
    % 1) EMG envelope & burst detection inside THIS REM bout
    % ==========================================================
    bpFilt = designfilt('bandpassiir', ...
                        'FilterOrder',4, ...
                        'HalfPowerFrequency1',20, ...
                        'HalfPowerFrequency2',300, ...
                        'SampleRate',Fs);
    emg_bp = filtfilt(bpFilt, emg_seg);
    emg_env = abs(emg_bp);
    win_env = round(0.1*Fs);
    emg_env = movmean(emg_env, win_env);

        % ---------- EMG stats for this REM bout ----------
    mu_rem = mean(emg_env);
    sd_rem = std(emg_env);
    max_rem = max(emg_env);
    fprintf('EMG env stats (REM bout %d): mean = %.4g, SD = %.4g, max = %.4g\n', ...
            bout_idx, mu_rem, sd_rem, max_rem);

    % ---------- choose threshold ----------
    % Option A: percentile-based (recommended)
    thr_percentile = 97.5;   % try 95–99; smaller = more bursts
    thr = prctile(emg_env, thr_percentile);

    % Option B: mean + k*SD (comment out A and uncomment this if you prefer)
    % thr_k = 1.5;           % was 2.5 before, now more lenient
    % thr = mu_rem + thr_k*sd_rem;

    fprintf('Using threshold thr = %.4g (%.1fth percentile)\n', thr, thr_percentile);

    % ---------- raw burst candidates from threshold crossings ----------
    isBurst_local = (emg_env > thr);
    d_b = diff([0 isBurst_local 0]);
    b_start = find(d_b ==  1);
    b_end   = find(d_b == -1) - 1;

    % compute times within REM bout
    burst_start_s = (b_start-1) / Fs;   % 0 = REM bout start
    burst_end_s   = (b_end-1)   / Fs;
    burst_dur_s   = burst_end_s - burst_start_s;

    fprintf('\nRaw EMG-burst candidates in REM bout %d (before filtering):\n', bout_idx);
    fprintf('  Burst   Start_s   End_s   Dur_s\n');
    for k = 1:numel(b_start)
        fprintf('  %3d    %7.3f  %7.3f  %7.3f\n', ...
                k, burst_start_s(k), burst_end_s(k), burst_dur_s(k));
    end
    if isempty(b_start)
        fprintf('  (none above threshold)\n');
    end
    fprintf('\n');

    % ---------- 1) minimum duration filter ----------
    min_dur_s = 0.1;   % 30 ms (shorter than your ~0.07 s bursts)
    keep_duration = burst_dur_s >= min_dur_s;

    % ---------- 2) edge margin filter ----------
    edge_margin_s = 2;  % ignore bursts within 2 s of start/end of REM bout
    keep_edge = (burst_start_s > edge_margin_s) & ...
                (burst_end_s   < (dur_seg - edge_margin_s));

    % combine filters
    keep = keep_duration & keep_edge;

    b_start       = b_start(keep);
    b_end         = b_end(keep);
    burst_start_s = burst_start_s(keep);
    burst_end_s   = burst_end_s(keep);
    burst_dur_s   = burst_dur_s(keep);

    % ---------- final bursts after filtering ----------
    nBursts = numel(b_start);
    fprintf('Detected %d EMG bursts in this REM bout (after filtering).\n', nBursts);
    if nBursts == 0
        warning('No bursts left after filtering – lower thr_percentile, min_dur_s or edge_margin_s.');
        return;
    end

    fprintf('\nEMG bursts kept in REM bout %d:\n', bout_idx);
    fprintf('  Burst   Start_s   End_s   Dur_s  (within bout)\n');
    for k = 1:nBursts
        fprintf('  %3d    %7.3f  %7.3f  %7.3f\n', ...
                k, burst_start_s(k), burst_end_s(k), burst_dur_s(k));
    end
    fprintf('\n');



    % ==========================================================
    % 2) Burst-locked window extraction
    % ==========================================================
    preT  = 2;   % seconds before burst
    postT = 3;   % seconds after burst
    win_samp = round((preT+postT)*Fs);
    t_win   = linspace(-preT, postT, win_samp);

    % we will stack segments as [time x burst]
    eeg_win_all = [];
    emg_win_all = [];

    for k = 1:nBursts
        c = b_start(k); % center index in rem-bout coordinates
        j1 = c - round(preT*Fs);
        j2 = c + round(postT*Fs) - 1;

        if j1 < 1 || j2 > numel(eeg_seg)
            % skip bursts too close to edges
            continue;
        end

        eeg_win_all(:,end+1) = eeg_seg(j1:j2); %#ok<AGROW>
        emg_win_all(:,end+1) = emg_seg(j1:j2); %#ok<AGROW>
    end

    nUsed = size(eeg_win_all,2);
    fprintf('Using %d bursts for peri-burst average (others too close to edges).\n', nUsed);
    if nUsed == 0
        warning('All bursts were too close to edges – no peri-burst windows.');
        return;
    end

    % mean and SEM (EEG)
    eeg_mean = mean(eeg_win_all,2);
    eeg_sem  = std(eeg_win_all,0,2)/sqrt(nUsed);

    % simple rectified EMG for visualisation
    emg_env_win = abs(emg_win_all);
    emg_mean = mean(emg_env_win,2);
    emg_sem  = std(emg_env_win,0,2)/sqrt(nUsed);

    % ==========================================================
    % 3) Burst-locked spectrogram (average)
    % ==========================================================
    spec_win   = round(0.5*Fs);  % 0.5 s window
    spec_nover = round(0.25*Fs);
    spec_nfft  = max(2^nextpow2(spec_win), spec_win);

    % get time/freq grid from first window
    eeg_first = eeg_win_all(:,1);
    [S0,F0,T0] = spectrogram(eeg_first, spec_win, spec_nover, spec_nfft, Fs);
    T0 = T0 + t_win(1);  % approximate alignment
    T0 = T0 - preT;      % center at 0
    fmax_spec = 50;
    idxF = F0 <= fmax_spec;
    F0 = F0(idxF);
    S_stack = zeros(numel(F0), numel(T0), nUsed);

    for k = 1:nUsed
        eegk = eeg_win_all(:,k);
        [Sk,Fk,Tk] = spectrogram(eegk, spec_win, spec_nover, spec_nfft, Fs);
        Tk = Tk + t_win(1);
        Tk = Tk - preT;

        % restrict freqs
        Sk = Sk(idxF,:);

        % map onto common Tk grid (nearest)
        [~, idxT] = arrayfun(@(x) min(abs(Tk - x)), T0);
        Sk = Sk(:, idxT);

        S_stack(:,:,k) = Sk;
    end

    S_mean = 10*log10(mean(abs(S_stack).^2,3) + eps);

    % ==========================================================
    % 4) Plot EMG envelope + burst-locked averages + spectrogram
    % ==========================================================
    figure('Color','w','Name',sprintf('REM bout %d – burst locked',bout_idx));
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % (A) REM bout EMG envelope + bursts overview
    nexttile;
    plot(t_seg, emg_env,'k'); hold on;
    yline(thr,'r--','Threshold');
    if ~isempty(b_start)
        plot(t_seg(b_start), emg_env(b_start), 'mo', 'MarkerFaceColor','m');
    end
    xlabel('Time within REM bout (s)');
    ylabel('EMG env');
    title(sprintf('REM bout %d (%.1f s) – EMG envelope & bursts', ...
                  bout_idx, dur_seg));
    grid on;

    % (B) Burst-locked EEG & EMG (mean ± SEM)
    nexttile;
    yyaxis left;
    fill_between(t_win, eeg_mean-eeg_sem, eeg_mean+eeg_sem, [0.8 0.8 1]);
    hold on;
    plot(t_win, eeg_mean, 'b','LineWidth',1.5);
    ylabel('EEG (a.u.)');

    yyaxis right;
    fill_between(t_win, emg_mean-emg_sem, emg_mean+emg_sem, [1 0.8 0.8]);
    hold on;
    plot(t_win, emg_mean,'r','LineWidth',1.5);
    ylabel('Rectified EMG');

    xline(0,'k--','Burst onset');
    xlabel('Time from EMG burst (s)');
    title(sprintf('Burst-locked averages (n = %d bursts)', nUsed));
    grid on;

    % (C) Mean peri-burst spectrogram
    nexttile;
    imagesc(T0, F0, S_mean);
    axis xy;
    xlabel('Time from EMG burst (s)');
    ylabel('Frequency (Hz)');
    title('Mean EEG spectrogram aligned to EMG bursts');
    c = colorbar; 
    c.Label.String = 'Power (dB)';
    xline(0,'k--','Burst onset');
    ylim([0 fmax_spec]);

    title(tl, sprintf('REM bout %d – EMG bursts', bout_idx),'FontWeight','bold');

    % ==========================================================
    % 5) Pre vs post burst PSD comparison (within this bout)
    % ==========================================================
    fprintf('Computing pre/post-burst PSDs...\n');

    pre_win  = [-2 0];  % seconds
    post_win = [0  2];

    % convert to indices in t_win (burst-locked time vector)
    [~, pre_idx1]  = min(abs(t_win - pre_win(1)));
    [~, pre_idx2]  = min(abs(t_win - pre_win(2)));
    [~, post_idx1] = min(abs(t_win - post_win(1)));
    [~, post_idx2] = min(abs(t_win - post_win(2)));

    Ppre  = [];
    Ppost = [];

    % one PSD per burst for pre and post
    for k = 1:nUsed
        eeg_pre  = eeg_win_all(pre_idx1:pre_idx2,  k);
        eeg_post = eeg_win_all(post_idx1:post_idx2, k);

        [Ppre(:,k), Fpsd]  = pwelch(eeg_pre,  [],[],[], Fs, 'onesided');
        [Ppost(:,k), ~]    = pwelch(eeg_post, [],[],[], Fs, 'onesided');
    end

    % average across bursts
    Ppre_mean  = mean(Ppre,  2);
    Ppost_mean = mean(Ppost, 2);

    fmax_plot = 50;
    idx_f = Fpsd <= fmax_plot;

    figure('Color','w','Name',sprintf('Bout %d pre/post PSD',bout_idx));
    plot(Fpsd(idx_f), 10*log10(Ppre_mean(idx_f)),  'k','LineWidth',1.5); hold on;
    plot(Fpsd(idx_f), 10*log10(Ppost_mean(idx_f)), 'r','LineWidth',1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(sprintf('REM bout %d – mean PSD pre vs post EMG burst', bout_idx));
    legend({'Pre (-2–0 s)','Post (0–2 s)'},'Location','best');
    grid on;



        % ==========================================================
    % 6) Band-wise relative power change (burst vs baseline)
    % ==========================================================
    fprintf('Computing band-wise relative power changes...\n');

    % power (not dB) from S_stack: [freq x time x bursts]
    P_stack = abs(S_stack).^2;

    % ---- time windows (in seconds) relative to burst ----
    burst_win    = [-0.5 0.5];   % around burst
    baseline_win = [-2   -0.5];  % pre-burst baseline

    idxT_burst = (T0 >= burst_win(1))    & (T0 <= burst_win(2));
    idxT_base  = (T0 >= baseline_win(1)) & (T0 <= baseline_win(2));

    if ~any(idxT_burst) || ~any(idxT_base)
        warning('Time windows for burst/baseline do not match T0 range.');
    end

    % ---- define frequency bands (Hz) ----
    band_names = { ...
        'Delta 0.5-4', ...
        'Theta 6-10', ...
        'Alpha 10-15', ...
        'Beta 15-30', ...
        'Gamma 30-50'};

    band_edges = [ ...
        0.5   4;  ...
        6    10;  ...
        10   15;  ...
        15   30;  ...
        30   50];

    nBands = size(band_edges,1);
    P_burst  = zeros(nBands,1);
    P_base   = zeros(nBands,1);

    % total power 0.5–50 Hz for normalization  (SUM, not mean)
    idxF_total = (F0 >= 0.5) & (F0 <= 50);

    P_burst_total = sum(P_stack(idxF_total, idxT_burst, :), 'all');
    P_base_total  = sum(P_stack(idxF_total, idxT_base,  :), 'all');

    for b = 1:nBands
        f1 = band_edges(b,1);
        f2 = band_edges(b,2);
        idxF_band = (F0 >= f1) & (F0 < f2);

        % sum over freq, time, bursts
        P_burst(b) = sum(P_stack(idxF_band, idxT_burst, :), 'all');
        P_base(b)  = sum(P_stack(idxF_band, idxT_base,  :), 'all');
    end

    % convert to relative power (% of 0.5–50 Hz)
    rel_burst = 100 * (P_burst / P_burst_total);
    rel_base  = 100 * (P_base  / P_base_total);
    rel_diff  = rel_burst - rel_base;   % percentage-point change

    % ---- print table ----
    fprintf('\nBand-wise relative power (%% of 0.5–50 Hz):\n');
    fprintf('  %-15s   Base(%%)   Burst(%%)   ΔBurst-Base(%%)\n','Band');
    for b = 1:nBands
        fprintf('  %-15s   %6.2f    %6.2f      %+6.2f\n', ...
                band_names{b}, rel_base(b), rel_burst(b), rel_diff(b));
    end
    fprintf('\n');

    % ---- bar plot of change ----
    figure('Color','w','Name',sprintf('Bout %d band-wise change', bout_idx));
    bar(rel_diff);
    set(gca,'XTick',1:nBands,'XTickLabel',band_names,'XTickLabelRotation',45);
    ylabel('\Delta relative power (percentage points)');
    title(sprintf('REM bout %d – burst vs baseline (0.5–50 Hz normalized)', ...
                  bout_idx));
    grid on;



end  % <-- end of main function


function h = fill_between(x, y1, y2, col)
% helper: shaded area between y1 and y2
    h = fill([x(:); flipud(x(:))], [y1(:); flipud(y2(:))], col, ...
             'EdgeColor','none','FaceAlpha',0.4);
end
