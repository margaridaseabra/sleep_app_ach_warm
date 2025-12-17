function qc = plot_one_episodes_spectrogram(eeg_mat_file, scores_csv, varargin)
% plot_one_episodes_spectrogram
% -------------------------------------------------------------------------
% Interactive QC of REM bouts with:
%   state, EEG, EMG, STFT spectrogram, Morlet spectrogram, PSD + metrics.
%
% Usage:
%   qc = plot_one_episodes_spectrogram( ...
%           '20251001_baseline_mouse1_APP.mat', ...
%           '20251001_baseline_mouse1_APP_scored_scores_1Hz.csv', ...
%           'eegVar','eeg', ...
%           'emgVar','emg', ...
%           'fsVar','eeg_frequency');
%
% Keys inside the viewer:
%   g = good REM (flag = 1)
%   b = bad / mis-scored (flag = 0)
%   m = maybe / unsure   (flag = -1)
%   t = trim start/end of this episode (absolute seconds)
%   p = previous episode
%   q = quit (save what you have)
%
%   → / space / n / enter : go to NEXT episode (flag unchanged)
%   ←                      : go to PREVIOUS episode (flag unchanged)
%
% Output struct qc:
%   qc.rem_st_s   : start times (s) AFTER trimming
%   qc.rem_en_s   : end times   (s) AFTER trimming
%   qc.flags      : 1=good, 0=bad, -1=maybe, NaN=unrated
%
% A MAT file '<eegBase>_REM_QC.mat' is also saved next to eeg_mat_file.
% -------------------------------------------------------------------------

    % ---------- Parse inputs ----------
    p = inputParser;
    addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
    addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
    addParameter(p,'context_s',20,@(x)isscalar(x)&&x>=0);   % pre/post context
    addParameter(p,'eegVar','eeg',@ischar);
    addParameter(p,'emgVar','emg',@ischar);
    addParameter(p,'fsVar','eeg_frequency',@ischar);
    parse(p,varargin{:});
    S = p.Results;
    C = S.codes;

    % ---------- Load scores ----------
    T = readtable(scores_csv);
    assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
        'Scores CSV must contain columns time_s and score');

    score   = double(T.score(:));    % state code per second / epoch
    n_epoch = numel(score);

    % find REM runs (logical mask)
    isREM = (score == C.REM);
    [rem_st_idx, rem_en_idx] = runs_from_codes(isREM);
    if isempty(rem_st_idx)
        warning('No REM bouts found in %s', scores_csv);
        qc = struct('rem_st_s',[],'rem_en_s',[],'flags',[]);
        return;
    end

    % convert epoch indices to seconds
    epoch_sec = S.epoch_sec;
    rem_st_s  = (rem_st_idx - 1) * epoch_sec;      % start at epoch boundary
    rem_en_s  = rem_en_idx * epoch_sec;            % end at end of last epoch
    nREM      = numel(rem_st_s);

    % ---------- Load EEG / EMG ----------
    Smat = load(eeg_mat_file);
    fns  = fieldnames(Smat);

    if ~isfield(Smat, S.eegVar)
        error('EEG variable "%s" not found in %s. Available: %s', ...
              S.eegVar, eeg_mat_file, strjoin(fns, ', '));
    end
    if ~isfield(Smat, S.emgVar)
        error('EMG variable "%s" not found in %s. Available: %s', ...
              S.emgVar, eeg_mat_file, strjoin(fns, ', '));
    end
    if ~isfield(Smat, S.fsVar)
        error('Sampling-rate variable "%s" not found in %s. Available: %s', ...
              S.fsVar, eeg_mat_file, strjoin(fns, ', '));
    end

    eeg = double(Smat.(S.eegVar)(:));
    emg = double(Smat.(S.emgVar)(:));
    fs  = double(Smat.(S.fsVar));

    n_eeg = numel(eeg);
    t_eeg = (0:n_eeg-1)'/fs;

    % ---------- QC loop ----------
    flags = NaN(nREM,1);

    fprintf('\nREM QC for %s\n', eeg_mat_file);
    fprintf('Controls inside figure:\n');
    fprintf('  g = good,  b = bad,  m = maybe,  t = trim,  p = previous,  q = quit\n');
    fprintf('  → / space / n / enter = next episode,  ← = previous episode\n\n');

    i = 1;
    while i <= nREM
        this_st = rem_st_s(i);
        this_en = rem_en_s(i);

        % context window around this episode
        win_pre  = S.context_s;
        win_post = S.context_s;
        win_st   = this_st - win_pre;
        win_en   = this_en + win_post;

        idx_win = t_eeg >= win_st & t_eeg <= win_en;
        t_win   = t_eeg(idx_win);
        eeg_win = eeg(idx_win);
        emg_win = emg(idx_win);

        % expand 1-Hz scoring to continuous samples in window
        sec_idx = floor(t_win / epoch_sec) + 1;
        sec_idx(sec_idx < 1)       = 1;
        sec_idx(sec_idx > n_epoch) = n_epoch;
        state_win = score(sec_idx);

        % show viewer for this REM episode
        fig = show_episode_viewer( ...
            t_win, eeg_win, emg_win, state_win, ...
            this_st, this_en, i, nREM, fs, ...
            'REM');

        % wait for key press (stored in fig.UserData)
        waitfor(fig, 'UserData');
        if ~ishandle(fig)
            fprintf('Figure closed; stopping QC.\n');
            break;
        end
        choice = get(fig,'UserData');
        close(fig);

        switch choice
            case 'g'
                flags(i) = 1;  i = i + 1;
            case 'b'
                flags(i) = 0;  i = i + 1;
            case 'm'
                flags(i) = -1; i = i + 1;
            case 't'
                [rem_st_s(i), rem_en_s(i)] = ask_trim(rem_st_s(i), rem_en_s(i));
                fprintf('Episode %d trimmed to [%.1f, %.1f] s.\n', ...
                        i, rem_st_s(i), rem_en_s(i));
            case {'p','prev'}
                i = max(1, i-1);
            case 'next'
                i = min(nREM, i+1);
            case 'q'
                fprintf('User requested quit; saving partial QC.\n');
                break;
        end
    end

    % ---------- Save QC ----------
    qc = struct('rem_st_s',rem_st_s, ...
                'rem_en_s',rem_en_s, ...
                'flags',   flags);

    [pth,base] = fileparts(eeg_mat_file);
    qc_file = fullfile(pth,[base '_REM_QC.mat']);
    save(qc_file,'qc');
    fprintf('REM QC saved to %s\n', qc_file);
end


% =====================================================================
% Helper: contiguous runs from logical vector
% =====================================================================
function [st, en] = runs_from_codes(mask)
    mask = mask(:) ~= 0;
    if isempty(mask)
        st = []; en = [];
        return;
    end
    d = diff([false; mask; false]);
    st = find(d == 1);
    en = find(d == -1) - 1;
end


% =====================================================================
% Helper: ask user for trimming (absolute seconds)
% =====================================================================
function [new_st, new_en] = ask_trim(old_st, old_en)
    prompt = { ...
        sprintf('New START (s, absolute time). Current = %.2f', old_st), ...
        sprintf('New END (s, absolute time).   Current = %.2f', old_en)};
    dlgtitle = 'Trim REM episode (absolute seconds)';
    dims = [1 60];
    definput = {sprintf('%.2f',old_st), sprintf('%.2f',old_en)};
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer)
        new_st = old_st;
        new_en = old_en;
        return;
    end

    new_st = str2double(answer{1});
    new_en = str2double(answer{2});

    if isnan(new_st) || isnan(new_en) || new_st >= new_en
        warning('Invalid trimming; keeping original bounds.');
        new_st = old_st;
        new_en = old_en;
    end
end


% =====================================================================
% Helper: episode viewer (state, EEG, EMG, STFT, Morlet, PSD)
% =====================================================================
function fig = show_episode_viewer(t, eeg, emg, state_win, ...
                                   st_s, en_s, idx, nEp, fs, ...
                                   label)

    fig = figure('Color','w','Name',[label ' QC viewer'], ...
                 'KeyPressFcn',@keyPressCallback);
    set(fig,'UserData','');

    rel_t  = t - st_s;        % time relative to episode start
    ep_dur = en_s - st_s;

    % 1) State
    ax1 = subplot(6,1,1);
    plot(rel_t, state_win,'k','LineWidth',1); hold on;
    yl = ylim;
    patch([0 ep_dur ep_dur 0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 1], 'EdgeColor','none','FaceAlpha',0.4);
    ylim(yl);
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('state');
    title(sprintf('%s episode %d/%d', label, idx, nEp));

    % 2) EEG
    ax2 = subplot(6,1,2);
    plot(rel_t, eeg); hold on;
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('EEG');
    title('Raw EEG');

    % 3) EMG
    ax3 = subplot(6,1,3);
    plot(rel_t, emg); hold on;
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('EMG');
    title('Raw EMG');

            % 4) EEG spectrogram (STFT, dB) + theta/delta line
    ax4 = subplot(6,1,4);
    try
        % ---- STFT parameters ----
        win_sec  = 4;                       % window length in seconds
        win_len  = round(fs * win_sec);
        if numel(eeg) < win_len
            win_len = numel(eeg);
        end
        noverlap = round(0.9 * win_len);    % 90% overlap
        nfft     = 2^nextpow2(win_len*2);   % lots of freq bins

        eeg_detr = detrend(eeg);
        [S,F,Tspec] = spectrogram(eeg_detr, win_len, noverlap, nfft, fs);

        % ----- power: linear + dB -----
        P_lin = abs(S).^2;                  % for theta/delta
        P     = 10*log10(P_lin + eps);      % for colour map

        % restrict to 0.5–30 Hz
        f_idx = F >= 0.5 & F <= 30;
        F     = F(f_idx);
        P     = P(f_idx,:);
        P_lin = P_lin(f_idx,:);

        % times relative to episode start
        t_spec_rel = rel_t(1) + Tspec;

        % ---- upsample grid for smooth image ----
        t_up = linspace(t_spec_rel(1), t_spec_rel(end), numel(t_spec_rel)*4);
        f_up = linspace(F(1),          F(end),          numel(F)*4);
        [Tgrid,Fgrid] = meshgrid(t_spec_rel, F);
        [Tq,Fq]       = meshgrid(t_up,       f_up);
        Pq            = interp2(Tgrid, Fgrid, P, Tq, Fq, 'linear');

        % ---- Gaussian smoothing ONLY in time (thin freq bands) ----
        sigma_t = 2.0;                       % time smoothing (pixels)
        w_t     = ceil(3*sigma_t);
        ct      = -w_t:w_t;                  % 1-D time axis
        K_t     = exp(-ct.^2/(2*sigma_t^2)); % 1-D Gaussian
        K_t     = K_t / sum(K_t);
        Pq_smooth = conv2(Pq, K_t, 'same');  % smooth along time only

        % ---- colour scaling: make band pop ----
        p_lo = prctile(Pq_smooth(:), 60);
        p_hi = prctile(Pq_smooth(:), 98);

        imagesc(t_up, f_up, Pq_smooth);
        set(gca,'YDir','normal');
        axis tight;
        ylim([0 30]);
        xlim([rel_t(1) rel_t(end)]);
        colormap(ax4, parula);
        caxis([p_lo p_hi]);
        colorbar;
        ylabel('Hz');
        hold on;
        xline(0,'w:','LineWidth',1.2);
        xline(ep_dur,'w:','LineWidth',1.2);
        yline(4,'w:');  % delta upper
        yline(6,'w:');  % theta lower
        yline(9,'w:');  % theta upper

        % ---- theta/delta ratio line (from original, non-upsampled STFT) ----
        % use same F subset as above
        F_sub = F;                          % 0.5–30 Hz

        delta_band = F_sub >= 0.5 & F_sub < 4;
        theta_band = F_sub >= 6   & F_sub < 9;

        Pdelta_t = sum(P_lin(delta_band,:), 1);
        Ptheta_t = sum(P_lin(theta_band,:), 1);
        ratio_td = Ptheta_t ./ max(Pdelta_t, eps);  % theta/delta per time bin

        % smooth ratio over time bins to make a nice line
        ratio_td = movmean(ratio_td, 5);

        % resample ratio to t_up so it matches image x-axis
        ratio_up = interp1(t_spec_rel, ratio_td, t_up, 'linear', 'extrap');

        % normalise ratio to a narrow range around 1 (like GUI)
        r_min = prctile(ratio_up, 5);
        r_max = prctile(ratio_up, 95);
        ratio_norm = (ratio_up - r_min) / (r_max - r_min);
        ratio_norm = 0.98 + 0.08 * ratio_norm;    % map to roughly [0.98, 1.06]

        % plot line on right axis
        yyaxis right;
        plot(t_up, ratio_norm, 'w', 'LineWidth', 1.5);
        ylim([0.96 1.08]);
        ylabel('Theta/Delta');
        yyaxis left;   % go back to left axis for any further plotting

        title('EEG spectrogram (STFT, dB) + \theta/\delta');
        hold off;

    catch ME
        warning(ME.identifier,'STFT spectrogram failed: %s',ME.message);
        text(0.5,0.5,'STFT spectrogram unavailable','Units','normalized', ...
             'HorizontalAlignment','center');
        axis off;
    end

    % 5) Morlet time–frequency (0–30 Hz), z-scored per freq
    ax5 = subplot(6,1,5);
    try
        [wt, f_wt] = cwt(eeg, fs, 'amor', 'FrequencyLimits',[0.5 30]);
        power = abs(wt).^2;

        mu = mean(power, 2);
        sd = std(power, 0, 2) + eps;
        power_z = (power - mu) ./ sd;

        imagesc(rel_t, f_wt, power_z);
        axis xy;
        ylim([0 30]);
        xlim([rel_t(1) rel_t(end)]);
        colormap(ax5, parula);
        caxis([-2 4]);
        colorbar;
        hold on;
        xline(0,'w:','LineWidth',1.2);
        xline(ep_dur,'w:','LineWidth',1.2);
        yline(4,'w:');
        yline(6,'w:');
        yline(9,'w:');
        hold off;
        xlabel('Time rel episode start (s)');
        ylabel('Hz');
        title('Morlet time–frequency (z)');
    catch ME
        warning(ME.identifier, 'cwt failed: %s', ME.message);
        text(0.5,0.5,'cwt unavailable','Units','normalized', ...
             'HorizontalAlignment','center');
        axis off;
    end

    % 6) PSD inside episode core + metrics
    ax6 = subplot(6,1,6);
    core_mask = rel_t >= 0 & rel_t <= ep_dur;
    eeg_core  = eeg(core_mask);
    emg_core  = emg(core_mask);

    if numel(eeg_core) < fs * 2
        text(0.5,0.5,sprintf('%s core too short for PSD',label), ...
            'Units','normalized', 'HorizontalAlignment','center');
        axis off;
    else
        win_len = round(fs * 4);   % 4 s window
        if numel(eeg_core) < win_len
            win_len = numel(eeg_core);
        end
        [pxx, f_psd] = pwelch(detrend(eeg_core), ...
                              hamming(win_len), [], [], fs);

        delta_idx = f_psd >= 0.5 & f_psd < 4;
        theta_idx = f_psd >= 6   & f_psd < 9;
        beta_idx  = f_psd >= 15  & f_psd < 30;

        Pdelta = trapz(f_psd(delta_idx), pxx(delta_idx));
        Ptheta = trapz(f_psd(theta_idx), pxx(theta_idx));
        Pbeta  = trapz(f_psd(beta_idx),  pxx(beta_idx));
        Ptot   = trapz(f_psd,            pxx);

        theta_delta = Ptheta / max(Pdelta, eps);
        theta_rel   = Ptheta / max(Ptot,   eps);
        beta_rel    = Pbeta  / max(Ptot,   eps);
        emg_rms     = rms(emg_core);

        plot(f_psd, 10*log10(pxx));
        xlim([0 30]);
        xlabel('Hz');
        ylabel('Power (dB)');
        title(sprintf('PSD inside %s core (%.1f s)', ...
              label, numel(eeg_core)/fs));
        grid on;

        txt = sprintf('\\theta/\\delta = %.2f   |  \\theta_{rel} = %.2f   |  \\beta_{rel} = %.2f   |  EMG RMS = %.2g', ...
                      theta_delta, theta_rel, beta_rel, emg_rms);
        yl = ylim;
        text(0.5, yl(2)-2, txt, ...
             'Units','data', ...
             'VerticalAlignment','top', ...
             'FontSize',8, ...
             'FontWeight','bold', ...
             'BackgroundColor',[1 1 1 0.7]);
    end

    linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
    xlabel(ax5,'Time rel episode start (s)');
    sgtitle(sprintf('%s QC – g:good  b:bad  m:maybe  t:trim  p:prev  q:quit',label));

    % nested callback to record key
    function keyPressCallback(src,event)
        key = event.Key;
        switch key
            case {'g','b','m','t','p','q'}
                set(src,'UserData',key);
            case {'rightarrow','space','n','return'}
                set(src,'UserData','next');
            case {'leftarrow'}
                set(src,'UserData','prev');
        end
    end
end
