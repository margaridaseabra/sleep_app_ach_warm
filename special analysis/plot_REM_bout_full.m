function plot_REM_bout_full(mat_file, REM_CODE, bout_idx)
% plot_REM_bout_full(mat_file, REM_CODE, bout_idx)
%
% For ONE REM bout, plot:
%   1) EEG + EMG traces for the whole bout
%   2) EEG spectrogram
%   3) EEG PSD (dB) for that bout (0–100 Hz)
%
% EXPECTED VARIABLES IN MAT:
%   eeg           : [1 x N] or [N x 1], EEG
%   emg           : [1 x N] or [N x 1], EMG
%   eeg_frequency : scalar, Fs (Hz)
%   sleep_scores  : [1 x T] or [T x 1], one score per second (0/1/2, etc.)
%
% INPUTS
%   mat_file : .mat file path
%   REM_CODE : value in sleep_scores that corresponds to REM (default = 2)
%   bout_idx : index of REM bout to plot (1-based)

    % ---------- inputs ----------
    if nargin < 1 || isempty(mat_file)
        [f,p] = uigetfile('*.mat','Select .mat file');
        if isequal(f,0); return; end
        mat_file = fullfile(p,f);
    end
    if nargin < 2 || isempty(REM_CODE)
        REM_CODE = 2;   % <-- change here if REM is coded differently
    end
    if nargin < 3
        error('You must provide bout_idx (which REM bout to plot).');
    end

    % ---------- load data ----------
    S = load(mat_file);
    EEG          = S.eeg(:)';          % force row
    EMG          = S.emg(:)';
    Fs           = S.eeg_frequency;
    sleep_scores = S.sleep_scores(:)';

    Nsamples = numel(EEG);
    if numel(EMG) ~= Nsamples
        error('EEG and EMG must have same length.');
    end
    Nsec  = numel(sleep_scores);
    t_eeg = (0:Nsamples-1)/Fs;

    % ---------- map 1-s scores to per-sample labels ----------
    est_N = round(Nsec * Fs);
    fprintf('EEG samples = %d, Nsec*Fs ≈ %d (diff = %d)\n', ...
            Nsamples, est_N, Nsamples-est_N);

    score_per_sample = zeros(1,Nsamples);
    for s = 1:Nsec
        i1 = floor((s-1)*Fs) + 1;
        i2 = min(round(s*Fs), Nsamples);
        score_per_sample(i1:i2) = sleep_scores(s);
    end

    isREM = (score_per_sample == REM_CODE);

    % ---------- find REM bouts ----------
    d_rem    = diff([0 isREM 0]);
    rem_start = find(d_rem ==  1);   % sample indices
    rem_end   = find(d_rem == -1) - 1;

    nREM = numel(rem_start);
    if nREM == 0
        error('No REM bouts found (REM_CODE = %d).', REM_CODE);
    end
    if bout_idx < 1 || bout_idx > nREM
        error('bout_idx must be between 1 and %d (asked for %d).', ...
              nREM, bout_idx);
    end

    i1 = rem_start(bout_idx);
    i2 = rem_end(bout_idx);
    seg_len = i2 - i1 + 1;
    dur_seg = seg_len / Fs;

    fprintf('Total REM bouts: %d | Plotting bout %d (%.1f s)\n', ...
            nREM, bout_idx, dur_seg);

    % segment
    eeg_seg = EEG(i1:i2);
    emg_seg = EMG(i1:i2);
    t_rel   = (0:seg_len-1)/Fs;   % time relative to bout start

    % ---------- PSD of this bout ----------
    win   = min(seg_len, round(4*Fs));   % 4-s window or shorter
    nover = floor(0.5*win);
    nfft  = max(2^nextpow2(win), win);

    [Pxx,F] = pwelch(eeg_seg, win, nover, nfft, Fs, 'onesided');
    fmax_plot = 100;
    idx_f = F <= fmax_plot;

    % ---------- Spectrogram ----------
    spec_win   = round(2*Fs);      % 2-s window
    spec_nover = round(1*Fs);      % 50% overlap
    spec_nfft  = max(2^nextpow2(spec_win), spec_win);

    [Sspec,Fspec,Tspec] = spectrogram(eeg_seg, spec_win, spec_nover, ...
                                      spec_nfft, Fs, 'yaxis');
    Sspec_db = 10*log10(abs(Sspec).^2 + eps);
    % shift Tspec to start at 0
    Tspec = Tspec + t_rel(1);
    % limit freqs for display
    fmax_spec = 50;
    idxFs = Fspec <= fmax_spec;

    % ---------- plot ----------
    figure('Color','w','Name',sprintf('REM bout %d', bout_idx));
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % 1) EEG + EMG traces
    nexttile;
    yyaxis left;
    plot(t_rel, eeg_seg);
    ylabel('EEG (a.u.)');
    yyaxis right;
    plot(t_rel, emg_seg);
    ylabel('EMG (a.u.)');
    xlabel('Time within REM bout (s)');
    title(sprintf('REM bout %d – duration %.1f s', bout_idx, dur_seg));
    grid on;

    % 2) Spectrogram
    nexttile;
    imagesc(Tspec, Fspec(idxFs), Sspec_db(idxFs,:));
    axis xy;
    xlabel('Time within REM bout (s)');
    ylabel('Frequency (Hz)');
    title('EEG spectrogram');
    colorbar;
    c = colorbar;
    c.Label.String = 'Power (dB)';
    ylim([0 fmax_spec]);

    % 3) PSD curve (like your group plot but for this bout)
    nexttile;
    plot(F(idx_f), 10*log10(Pxx(idx_f)),'LineWidth',1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('EEG PSD during this REM bout');
    grid on;
    xlim([0 fmax_plot]);

    title(tl, sprintf('REM bout %d', bout_idx), 'FontWeight','bold');

end
