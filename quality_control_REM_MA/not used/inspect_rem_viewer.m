function EPISODES = inspect_rem_viewer(eeg_mat_file, scores_csv, varargin)
% inspect_rem_viewer
% -------------------------------------------------------------------------
% Interactive viewer for REM episodes.
%
% Shows, for each REM bout:
%   (1) State code (1 Hz hypnogram fragment)
%   (2) Raw EEG
%   (3) Raw EMG
%   (4) Morlet time–frequency power (EEG, 0–15 Hz, z-scored per frequency)
%   (5) PSD inside the REM core + numeric features:
%           - theta/delta power ratio (6–9 / 0.5–4 Hz)
%           - theta relative power (theta / total 0.5–30 Hz)
%           - beta relative power (15–30 / total)
%           - EMG RMS
%
% You can classify each REM as:
%   g  : good REM
%   b  : bad / mis-scored
%   u  : unsure
%   n, RightArrow : next episode
%   p, LeftArrow  : previous episode
%   q, Esc        : close viewer
%
% Usage:
%   EP = inspect_rem_viewer('20251001_baseline_mouse1_APP.mat', ...
%                           '20251001_baseline_mouse1_APP_scored_scores_1Hz.csv');
%
% Options (name/value):
%   'codes'      : struct with fields WK, NREM, REM, MA (defaults 0,1,2,15)
%   'epoch_sec'  : epoch length for score CSV (default 1)
%   'context_s'  : seconds of context before/after REM bout (default 20)
%   'min_rem_s'  : minimum REM bout duration to include (default 10)
%
% Returns:
%   EPISODES : struct array with fields
%       .start_s, .end_s, .dur_s
%       .idx_start, .idx_end      (indices in 1 Hz code vector)
%       .flag   (NaN = unclassified, 1=good, 0=bad)
%
% Requirements:
%   - Wavelet Toolbox (for cwt)
%   - eeg_mat_file must contain EEG & EMG + sampling rate:
%       EEG candidates: 'eeg','EEG','EEG_rawtrace','Data_EEG'
%       EMG candidates: 'emg','EMG','EMG_rawtrace'
%       fs candidates : 'eeg_frequency','fs_eeg','fs','Fs','sampling_freq'
%
% -------------------------------------------------------------------------

% -------- Parse inputs ---------------------------------------------------
p = inputParser;
addRequired(p,'eeg_mat_file',@ischar);
addRequired(p,'scores_csv',@ischar);

codes_default = struct('WK',0,'NREM',1,'REM',2,'MA',15);
addParameter(p,'codes',codes_default,@isstruct);
addParameter(p,'epoch_sec',1,@(x)isscalar(x)&&x>0);
addParameter(p,'context_s',20,@(x)isscalar(x)&&x>=0);
addParameter(p,'min_rem_s',10,@(x)isscalar(x)&&x>=0);

parse(p, eeg_mat_file, scores_csv, varargin{:});
S = p.Results;
C = S.codes;

% -------- Load scores (1 Hz) --------------------------------------------
T = readtable(scores_csv);
assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
    'Scores CSV must have columns time_s and score.');

time_s = double(T.time_s(:));   % seconds
code   = double(T.score(:));    % state codes

epoch_sec = S.epoch_sec;
n_epoch   = numel(code);
dur_scores = n_epoch * epoch_sec;

fprintf('Loaded scores: %d epochs (%.1f min)\n', n_epoch, dur_scores/60);

% -------- Find REM bouts -------------------------------------------------
isREM = (code == C.REM);

[st_idx, en_idx] = runs_from_codes(isREM);
dur_idx = (en_idx - st_idx + 1);              % in epochs
dur_s   = dur_idx * epoch_sec;

keep = dur_s >= S.min_rem_s;
st_idx = st_idx(keep);
en_idx = en_idx(keep);
dur_s  = dur_s(keep);

n_rem = numel(st_idx);
if n_rem == 0
    error('No REM bouts ≥ %.1f s found.', S.min_rem_s);
end

EP = struct( ...
    'idx_start', num2cell(st_idx), ...
    'idx_end',   num2cell(en_idx), ...
    'start_s',   num2cell((st_idx-1)*epoch_sec), ...
    'end_s',     num2cell(en_idx*epoch_sec), ...
    'dur_s',     num2cell(dur_s), ...
    'flag',      num2cell(NaN(size(st_idx))) );

fprintf('Found %d REM bouts (min dur %.1f s).\n', n_rem, S.min_rem_s);

% -------- Load EEG + EMG + fs -------------------------------------------
EEG_struct = load(eeg_mat_file);

cand_eeg = {'eeg','EEG','EEG_rawtrace','Data_EEG'};
cand_emg = {'emg','EMG','EMG_rawtrace'};
cand_fs  = {'eeg_frequency','fs_eeg','fs','Fs','sampling_freq'};

[eeg, eeg_name] = pick_field(EEG_struct, cand_eeg);
[emg, emg_name] = pick_field(EEG_struct, cand_emg);
[fs, fs_name]   = pick_field(EEG_struct, cand_fs);

eeg = double(eeg(:));
emg = double(emg(:));
fs  = double(fs);

fprintf('EEG: %s (n = %d), EMG: %s (n = %d), fs = %.3f Hz (%s)\n', ...
    eeg_name, numel(eeg), emg_name, numel(emg), fs, fs_name);

n_eeg = numel(eeg);
t_eeg = (0:n_eeg-1)'/fs;      % seconds from start

% 1-Hz time for codes (assume aligned with EEG start)
t_code = ((0:n_epoch-1)' * epoch_sec);

% -------- Create figure & interactive viewer -----------------------------
currentIdx = 1;

fig = figure('Name','REM episode viewer', ...
             'Color','w', ...
             'Units','normalized', ...
             'Position',[0.05 0.05 0.9 0.85]);

set(fig,'KeyPressFcn',@onKeyPress);

plotEpisode(currentIdx);

% wait until figure is closed
waitfor(fig);

if nargout > 0
    EPISODES = EP;
end

% ==================== Nested functions ==================================

    function plotEpisode(k)
        if ~ishandle(fig); return; end

        clf(fig);

        ep = EP(k);
        core_start = ep.start_s;
        core_end   = ep.end_s;
        core_dur   = ep.dur_s;

        % Window with context
        t0 = core_start - S.context_s;
        t1 = core_end   + S.context_s;

        % Clip to actual recording
        t0 = max(t0, 0);
        t1 = min(t1, t_eeg(end));

        win_mask = t_eeg >= t0 & t_eeg <= t1;

        t_win  = t_eeg(win_mask);
        eeg_w  = eeg(win_mask);
        emg_w  = emg(win_mask);

        % Time relative to REM onset (0 = start of REM)
        t_rel = t_win - core_start;

        % ---- State code fragment (1 Hz) ---------------------------------
        code_mask = t_code >= t0 & t_code <= t1;
        t_code_w  = t_code(code_mask);
        code_w    = code(code_mask);
        t_code_rel = t_code_w - core_start;

        subplot(5,1,1);
        hold on;
        % shade REM core
        yl = [min(code_w)-0.5, max(code_w)+0.5];
        patch([0 core_dur core_dur 0], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.9 0.9 1], 'EdgeColor','none','FaceAlpha',0.5);
        stairs(t_code_rel, code_w, 'k','LineWidth',1);
        xlabel('Time rel REM onset (s)');
        ylabel('state code');
        xlim([t_rel(1) t_rel(end)]);
        ylim(yl);
        yline(C.NREM,'Color',[0.6 0.6 0.6],'LineStyle',':');
        yline(C.REM, 'Color',[0 0.5 1],'LineStyle',':');
        box off;

        flag_str = flagToStr(ep.flag);
        title(sprintf('REM episode %d/%d   (flag = %s)', k, n_rem, flag_str));

        % ---- Raw EEG ----------------------------------------------------
        subplot(5,1,2);
        plot(t_rel, eeg_w);
        ylabel('EEG');
        xlim([t_rel(1) t_rel(end)]);
        yline(0,'Color',[0.7 0.7 0.7],'LineStyle',':');
        hold on;
        xline(0,'k:');
        xline(core_dur,'k:');
        title('Raw EEG');

        % ---- Raw EMG ----------------------------------------------------
        subplot(5,1,3);
        plot(t_rel, emg_w);
        ylabel('EMG');
        xlim([t_rel(1) t_rel(end)]);
        yline(0,'Color',[0.7 0.7 0.7],'LineStyle',':');
        hold on;
        xline(0,'k:');
        xline(core_dur,'k:');
        title('Raw EMG');

        % ---- Morlet time–frequency (0–15 Hz), z-scored per freq --------
        subplot(5,1,4);

        try
            % Continuous wavelet transform
            [wt, f_wt] = cwt(eeg_w, fs, 'amor', 'FrequencyLimits',[0.5 15]);

            power = abs(wt).^2;

            % z-score along time for each frequency
            mu = mean(power, 2);
            sd = std(power, 0, 2) + eps;
            power_z = (power - mu) ./ sd;

            imagesc(t_rel, f_wt, power_z);
            axis xy;
            ylim([0 15]);
            xlim([t_rel(1) t_rel(end)]);
            colormap(turbo);
            caxis([-2 4]);
            colorbar;
            hold on;
            xline(0,'w:','LineWidth',1.2);
            xline(core_dur,'w:','LineWidth',1.2);
            % band boundaries
            yline(4,'w:');  % delta upper
            yline(6,'w:');  % theta lower
            yline(9,'w:');  % theta upper
            xlabel('Time rel REM onset (s)');
            ylabel('Hz');
            title('Morlet time–frequency power (z per freq)');

        catch ME
            warning('cwt failed: %s', ME.message);
            text(0.5,0.5,'cwt unavailable','Units','normalized', ...
                'HorizontalAlignment','center');
        end

        % ---- PSD inside REM core + metrics -----------------------------
        subplot(5,1,5);

        core_mask = t_rel >= 0 & t_rel <= core_dur;
        eeg_core = eeg_w(core_mask);
        emg_core = emg_w(core_mask);

        if numel(eeg_core) < fs * 2
            % too short to get a good PSD
            text(0.5,0.5,'REM core too short for PSD','Units','normalized', ...
                'HorizontalAlignment','center');
            axis off;
            return;
        end

        % PSD (Welch) 0.5–30 Hz
        win_len = round(fs * 4);   % 4 s window
        if numel(eeg_core) < win_len
            win_len = numel(eeg_core);
        end
        [pxx, f_psd] = pwelch(detrend(eeg_core), ...
                              hamming(win_len), [], [], fs);

        % Metrics: band powers
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

        % EMG RMS
        emg_rms = rms(emg_core);

        plot(f_psd, 10*log10(pxx));
        xlim([0 30]);
        xlabel('Hz');
        ylabel('Power (dB)');
        title('PSD inside REM core');

        % Display metrics text in top-left
        txt = sprintf('\\theta/\\delta = %.2f   |  \\theta_{rel} = %.2f   |  \\beta_{rel} = %.2f   |  EMG RMS = %.2g', ...
                      theta_delta, theta_rel, beta_rel, emg_rms);

        ylim_current = ylim;
        y_top = ylim_current(2);
        x_left = 0.02 * 30;

        text(1, y_top - 2, txt, ...
            'Units','data', ...
            'VerticalAlignment','top', ...
            'FontSize',8, ...
            'FontWeight','bold', ...
            'BackgroundColor',[1 1 1 0.7]);

        % Link x-axes of the first 4 subplots
        linkaxes(findall(fig,'Type','axes'), 'x');
    end

    % -------- Key press handler -----------------------------------------
    function onKeyPress(~, event)
        switch event.Key
            case {'rightarrow','n'}
                currentIdx = min(currentIdx + 1, n_rem);
                plotEpisode(currentIdx);

            case {'leftarrow','p'}
                currentIdx = max(currentIdx - 1, 1);
                plotEpisode(currentIdx);

            case {'g','G'}
                EP(currentIdx).flag = 1;
                plotEpisode(currentIdx);

            case {'b','B'}
                EP(currentIdx).flag = 0;
                plotEpisode(currentIdx);

            case {'u','U'}
                EP(currentIdx).flag = NaN;
                plotEpisode(currentIdx);

            case {'q','escape'}
                if ishandle(fig)
                    close(fig);
                end
        end
    end

end  % main function


% ==================== Helper functions (outside) ========================

function [data, name] = pick_field(S, candidates)
% Return the first field in S that exists from the candidates list.
for i = 1:numel(candidates)
    if isfield(S, candidates{i})
        data = S.(candidates{i});
        name = candidates{i};
        return;
    end
end
error('None of the candidate fields found: %s', strjoin(candidates, ', '));
end

function [st,en] = runs_from_codes(mask)
% Given a logical vector, return start/end indices of contiguous true runs.
mask = mask(:);
if isempty(mask)
    st = []; en = [];
    return;
end
d = diff([false; mask; false]);
st = find(d == 1);
en = find(d == -1) - 1;
end

function s = flagToStr(flag)
if isnan(flag)
    s = 'NaN';
elseif flag == 1
    s = 'good';
elseif flag == 0
    s = 'bad';
else
    s = num2str(flag);
end
end
