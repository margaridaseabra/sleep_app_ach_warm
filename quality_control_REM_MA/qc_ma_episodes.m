function qc = qc_ma_episodes(eeg_mat_file, scores_csv, varargin)
% qc_ma_episodes  Interactive QC of micro-arousal (MA) bouts.
%
% Same viewer and key controls as qc_rem_episodes, but for MA:
%
%   g = good MA
%   b = bad / mis-scored
%   m = maybe
%   t = trim (absolute seconds)
%   p = previous
%   q = quit
%   → / space / n / enter : next episode
%   ←                      : previous episode
%
% Output struct qc:
%   qc.ma_st_s, qc.ma_en_s, qc.flags
%
% Saved to '<eegBase>_MA_QC.mat'.

    % ---------- Parse inputs ----------
    p = inputParser;
    addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
    addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
    addParameter(p,'context_s',20,@(x)isscalar(x)&&x>=0);
    addParameter(p,'eegVar','eeg',@ischar);
    addParameter(p,'emgVar','emg',@ischar);
    addParameter(p,'fsVar','eeg_frequency',@ischar);
    % for MA we want to see up to high beta
    addParameter(p,'fmax',40,@(x)isscalar(x)&&x>0);  

    parse(p,varargin{:});
    S = p.Results;
    C = S.codes;

    % ---------- Load scores ----------
    T = readtable(scores_csv);
    assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
        'Scores CSV must contain columns time_s and score');

    time_s = double(T.time_s(:));
    score  = double(T.score(:));
    n_epoch = numel(score);

    % find MA runs
    isMA = (score == C.MA);
    [ma_st_idx, ma_en_idx] = runs_from_codes(isMA);
    if isempty(ma_st_idx)
        warning('No MA bouts found in %s', scores_csv);
        qc = struct('ma_st_s',[],'ma_en_s',[],'flags',[]);
        return;
    end

    epoch_sec = S.epoch_sec;
    ma_st_s   = (ma_st_idx - 1) * epoch_sec;
    ma_en_s   = ma_en_idx * epoch_sec;
    nMA       = numel(ma_st_s);

    % ---------- Load EEG / EMG ----------
    Smat = load(eeg_mat_file);

    if ~isfield(Smat, S.eegVar)
        error('EEG variable "%s" not found in %s', S.eegVar, eeg_mat_file);
    end
    if ~isfield(Smat, S.emgVar)
        error('EMG variable "%s" not found in %s', S.emgVar, eeg_mat_file);
    end
    if ~isfield(Smat, S.fsVar)
        error('Sampling-rate variable "%s" not found in %s', S.fsVar, eeg_mat_file);
    end

    eeg = double(Smat.(S.eegVar)(:));
    emg = double(Smat.(S.emgVar)(:));
    fs  = double(Smat.(S.fsVar));

    n_eeg = numel(eeg);
    t_eeg = (0:n_eeg-1)'/fs;

    % ---------- QC loop ----------
    flags = NaN(nMA,1);

    fprintf('\nMA QC for %s\n', eeg_mat_file);
    fprintf('Controls: g/b/m/t/p/q,  → next,  ← prev\n\n');

    i = 1;
    while i <= nMA
        this_st = ma_st_s(i);
        this_en = ma_en_s(i);

        win_pre  = S.context_s;
        win_post = S.context_s;
        win_st   = this_st - win_pre;
        win_en   = this_en + win_post;

        idx_win = t_eeg >= win_st & t_eeg <= win_en;
        t_win   = t_eeg(idx_win);
        eeg_win = eeg(idx_win);
        emg_win = emg(idx_win);

        % map each EEG sample to its 1-Hz epoch for hypnogram
        sec_idx = floor(t_win / epoch_sec) + 1;
        sec_idx(sec_idx < 1)       = 1;
        sec_idx(sec_idx > n_epoch) = n_epoch;
        state_win = score(sec_idx);

        fig = plot_one_episode( ...
            t_win, eeg_win, emg_win, state_win, ...
            this_st, this_en, i, nMA, fs, ...
            'MA', S.fmax);

        waitfor(fig,'UserData');
        if ~ishandle(fig)
            fprintf('Figure closed; stopping MA QC.\n');
            break;
        end
        choice = get(fig,'UserData');
        close(fig);

        switch choice
            case 'g'
                flags(i) = 1;
                i = i + 1;
            case 'b'
                flags(i) = 0;
                i = i + 1;
            case 'm'
                flags(i) = -1;
                i = i + 1;
            case 't'
                [ma_st_s(i), ma_en_s(i)] = ask_trim(ma_st_s(i), ma_en_s(i));
                fprintf('MA %d trimmed to [%.1f, %.1f] s.\n', ...
                        i, ma_st_s(i), ma_en_s(i));
            case {'p','prev'}
                i = max(1, i-1);
            case 'next'
                i = min(nMA, i+1);
            case 'q'
                fprintf('User requested quit; saving partial QC.\n');
                break;
            otherwise
                % ignore
        end
    end

    qc = struct('ma_st_s',ma_st_s, ...
                'ma_en_s',ma_en_s, ...
                'flags',   flags);

    [pth,base] = fileparts(eeg_mat_file);
    qc_file = fullfile(pth,[base '_MA_QC.mat']);
    save(qc_file,'qc');
    fprintf('MA QC saved to %s\n', qc_file);
end


% ---------- shared helpers ----------------------------------------------
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

function [new_st, new_en] = ask_trim(old_st, old_en)
    prompt = { ...
        sprintf('New START (s, absolute time). Current = %.2f', old_st), ...
        sprintf('New END (s, absolute time).   Current = %.2f', old_en)};
    dlgtitle = 'Trim MA episode (absolute seconds)';
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

function fig = plot_one_episode(t, eeg, emg, state_win, ...
                                st_s, en_s, idx, nEp, fs, ...
                                label, fmax)

    if nargin < 11 || isempty(fmax)
        fmax = 40;
    end
    % don't exceed Nyquist
    fmax = min(fmax, fs/2 - 1);

    fig = figure('Color','w','Name',[label ' QC viewer'], ...
                 'KeyPressFcn',@keyPressCallback);
    set(fig,'UserData','');

    rel_t  = t - st_s;
    ep_dur = en_s - st_s;

    % --- 1) State code ---------------------------------------------------
    ax1 = subplot(5,1,1);
    plot(rel_t, state_win,'k','LineWidth',1); hold on;
    yl = ylim;
    patch([0 ep_dur ep_dur 0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 1], 'EdgeColor','none','FaceAlpha',0.4);
    ylim(yl);
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('state');
    title(sprintf('%s episode %d/%d', label, idx, nEp));

    % --- 2) EEG ----------------------------------------------------------
    ax2 = subplot(5,1,2);
    plot(rel_t, eeg); hold on;
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('EEG');
    title('Raw EEG');

    % --- 3) EMG ----------------------------------------------------------
    ax3 = subplot(5,1,3);
    plot(rel_t, emg); hold on;
    xline(0,'k--'); xline(ep_dur,'k--');
    ylabel('EMG');
    title('Raw EMG');

    % --- 4) Morlet time–frequency, tuned for MA -------------------------
    ax4 = subplot(5,1,4);
    try
        % Analyse 0.5–fmax Hz; we will *display* up to ~30 Hz
        f_range = [0.5 fmax];
        [wt, f_wt] = cwt(eeg, fs, 'amor', 'FrequencyLimits',f_range);
        power = abs(wt).^2;

        % z-score along time for each frequency -> bursts pop out
        mu = mean(power, 2);
        sd = std(power, 0, 2) + eps;
        power_z = (power - mu) ./ sd;

        imagesc(rel_t, f_wt, power_z);
        axis xy;
        ylim([0 min(30,fmax)]);             % MA: see delta..beta
        xlim([rel_t(1) rel_t(end)]);
        colormap(turbo);
        caxis([-2 4]);
        colorbar;
        hold on;
        xline(0,'w:','LineWidth',1.2);
        xline(ep_dur,'w:','LineWidth',1.2);

        % band boundaries helpful for MA: delta/theta, theta/sigma, sigma/beta
        if 4 <= fmax,  yline(4,'w:');  end      % delta upper
        if 8 <= fmax,  yline(8,'w:');  end      % theta upper
        if 15 <= fmax, yline(15,'w:'); end      % sigma/beta-ish
        hold off;
        xlabel('Time rel episode start (s)');
        ylabel('Hz');
        title('Morlet power (0.5–30 Hz, z per freq)');
    catch ME
        warning(ME.identifier, 'cwt failed: %s', ME.message);
        text(0.5,0.5,'cwt unavailable','Units','normalized', ...
             'HorizontalAlignment','center');
        axis off;
    end

    % --- 5) PSD inside MA core + metrics --------------------------------
    ax5 = subplot(5,1,5);

    rel_t = t - st_s;                % recompute to be safe
    core_mask = rel_t >= 0 & rel_t <= ep_dur;
    eeg_core  = eeg(core_mask);
    emg_core  = emg(core_mask);

    if numel(eeg_core) < fs * 2
        text(0.5,0.5,sprintf('%s core too short for PSD',label), ...
            'Units','normalized', 'HorizontalAlignment','center');
        axis off;
    else
        win_len = round(fs * 4);
        if numel(eeg_core) < win_len
            win_len = numel(eeg_core);
        end
        [pxx, f_psd] = pwelch(detrend(eeg_core), ...
                              hamming(win_len), [], [], fs);

        % same metrics as REM viewer: delta, theta, beta
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

        emg_rms = rms(emg_core);

        plot(f_psd, 10*log10(pxx));
        xlim([0 fmax]);
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

    linkaxes([ax1 ax2 ax3 ax4],'x');
    xlabel(ax4,'Time rel episode start (s)');
    sgtitle(sprintf('%s QC – g:good  b:bad  m:maybe  t:trim  p:prev  q:quit',label));

    % record key pressed so outer loop can act
    function keyPressCallback(src,event)
        key = event.Key;
        switch key
            case {'g','b','m','t','p','q'}
                set(src,'UserData',key);
            case {'rightarrow','space','n','return'}
                set(src,'UserData','next');
            case {'leftarrow'}
                set(src,'UserData','prev');
            otherwise
        end
    end
end
