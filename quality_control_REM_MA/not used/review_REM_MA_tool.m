function review_REM_MA_tool(eeg_mat_file, scores_csv, episode_type, varargin)
% review_REM_MA_from_files
% --------------------------------------------------------------
% Interactive viewer to QC REM or MA episodes directly from:
%   - a .mat with EEG/EMG + sampling rate
%   - a CSV with 1-Hz scoring (time_s, score)
%
% INPUTS:
%   eeg_mat_file : string, path to .mat with EEG/EMG
%   scores_csv   : string, path to CSV with columns time_s, score
%   episode_type : 'REM' or 'MA'
%
% OPTIONS (name-value):
%   'codes'      : struct with fields WK,NREM,REM,MA (numeric codes)
%                  default: struct('WK',0,'NREM',1,'REM',2,'MA',15)
%   'pre_sec'    : seconds of context before episode onset (default 20)
%   'post_sec'   : seconds of context after episode end   (default 20)
%   'min_dur_sec': minimum episode duration to include (default 5)
%   'out_file'   : filename to save flags (default '<episode_type>_flags.mat')
%
% KEYBOARD CONTROLS:
%   → or space : next episode
%   ←          : previous episode
%   g          : mark episode as GOOD (flag = 1)
%   b          : mark episode as BAD  (flag = 0)
%   u          : mark episode as UNDECIDED (NaN)
%   s          : save flags to out_file
%
% Saved struct 'ep' has fields:
%   .type, .start_sec, .end_sec, .flag
% --------------------------------------------------------------

%% --- Parse options ---
p = inputParser;
addRequired(p,'eeg_mat_file',@ischar);
addRequired(p,'scores_csv',@ischar);
addRequired(p,'episode_type',@(s)ischar(s) || isstring(s));

codes_default = struct('WK',0,'NREM',1,'REM',2,'MA',15);

addParameter(p,'codes',codes_default,@isstruct);
addParameter(p,'pre_sec',20,@isnumeric);
addParameter(p,'post_sec',20,@isnumeric);
addParameter(p,'min_dur_sec',5,@isnumeric);
addParameter(p,'out_file','',@ischar);

parse(p,eeg_mat_file,scores_csv,episode_type,varargin{:});
OPT = p.Results;

episode_type = upper(string(episode_type));
C = OPT.codes;

if isempty(OPT.out_file)
    OPT.out_file = sprintf('%s_flags.mat', episode_type);
end

%% --- Load EEG/EMG/FS from MAT file ---
EEG_struct = load(eeg_mat_file);

% candidate EEG variable names
cand_eeg = {'eeg','EEG','EEG_rawtrace','EEG_raw'};
eeg = [];
eeg_name = '';
for k = 1:numel(cand_eeg)
    if isfield(EEG_struct, cand_eeg{k})
        eeg = EEG_struct.(cand_eeg{k});
        eeg_name = cand_eeg{k};
        break;
    end
end
if isempty(eeg)
    error('Could not find EEG variable in %s', eeg_mat_file);
end
eeg = double(eeg(:));

% candidate EMG variable names
cand_emg = {'emg','EMG','EMG_rawtrace','EMG_raw'};
emg = [];
emg_name = '';
for k = 1:numel(cand_emg)
    if isfield(EEG_struct, cand_emg{k})
        emg = EEG_struct.(cand_emg{k});
        emg_name = cand_emg{k};
        break;
    end
end
if isempty(emg)
    warning('No EMG variable found in %s. EMG plot will be empty.', eeg_mat_file);
    emg = zeros(size(eeg));
end
emg = double(emg(:));

if numel(emg) ~= numel(eeg)
    warning('EMG length != EEG length; truncating to min length.');
    L = min(numel(eeg), numel(emg));
    eeg = eeg(1:L);
    emg = emg(1:L);
end

% candidate sampling rate variable names
cand_fs = {'eeg_frequency','fs_eeg','fs','Fs','sampling_freq'};
fs = [];
fs_name = '';
for k = 1:numel(cand_fs)
    if isfield(EEG_struct, cand_fs{k})
        fs = EEG_struct.(cand_fs{k});
        fs_name = cand_fs{k};
        break;
    end
end
if isempty(fs)
    error('Could not find sampling rate in %s', eeg_mat_file);
end
fs = double(fs);

fprintf('Loaded EEG (%s) and EMG (%s) from %s, fs = %.2f Hz (%s)\n', ...
    eeg_name, emg_name, eeg_mat_file, fs, fs_name);

total_dur_sec = numel(eeg)/fs;
t_eeg = (0:numel(eeg)-1)'/fs;

%% --- Load scores from CSV ---
T = readtable(scores_csv);
assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
    'CSV must have columns time_s and score.');

t_scores = double(T.time_s(:));
score    = double(T.score(:));

% infer epoch length from time_s (assume mostly regular)
if numel(t_scores) > 1
    dt = median(diff(t_scores));
else
    dt = 1;
end
epoch_len_sec = dt;
fprintf('Scores loaded from %s, epoch_len_sec ≈ %.3f s\n', scores_csv, epoch_len_sec);

%% --- Choose target code ---
switch episode_type
    case "REM"
        target_code = C.REM;
    case "MA"
        target_code = C.MA;
    otherwise
        error('episode_type must be ''REM'' or ''MA''.');
end

%% --- Find contiguous episodes in 1-Hz scoring ---
is_target = (score == target_code);
d = diff([0; is_target; 0]);
start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

episodes = struct('type',{},'start_sec',{},'end_sec',{},'flag',{});

for k = 1:numel(start_idx)
    s_i = start_idx(k);
    e_i = end_idx(k);

    % start/end times; we treat each score as covering [t, t+epoch_len_sec)
    start_sec = t_scores(s_i);
    end_sec   = t_scores(e_i) + epoch_len_sec;
    dur_sec   = end_sec - start_sec;

    if dur_sec >= OPT.min_dur_sec
        episodes(end+1).type      = char(episode_type); %#ok<AGROW>
        episodes(end).start_sec   = start_sec;
        episodes(end).end_sec     = end_sec;
        episodes(end).flag        = NaN;
    end
end

if isempty(episodes)
    error('No %s episodes found with duration >= %.1f s.', episode_type, OPT.min_dur_sec);
end

nEpisodes = numel(episodes);
fprintf('Found %d %s episodes (>= %.1f s).\n', nEpisodes, episode_type, OPT.min_dur_sec);

%% --- Prepare interactive figure ---
curIdx = 1;

fig = figure('Name', sprintf('%s episode review', episode_type), ...
             'Color','w', ...
             'NumberTitle','off', ...
             'KeyPressFcn', @onKeyPress);

plotCurrentEpisode();

    function onKeyPress(~, event)
        switch event.Key
            case {'rightarrow','space'}
                curIdx = min(curIdx + 1, nEpisodes);
                plotCurrentEpisode();
            case 'leftarrow'
                curIdx = max(curIdx - 1, 1);
                plotCurrentEpisode();
            case 'g'
                episodes(curIdx).flag = 1;
                fprintf('Episode %d marked GOOD.\n', curIdx);
                plotCurrentEpisode();
            case 'b'
                episodes(curIdx).flag = 0;
                fprintf('Episode %d marked BAD.\n', curIdx);
                plotCurrentEpisode();
            case 'u'
                episodes(curIdx).flag = NaN;
                fprintf('Episode %d marked UNDECIDED.\n', curIdx);
                plotCurrentEpisode();
            case 's'
                saveFlags();
            otherwise
                % ignore
        end
    end

    function saveFlags()
        ep = episodes; %#ok<NASGU>
        save(OPT.out_file, 'ep');
        fprintf('Saved %d episodes + flags to %s\n', nEpisodes, OPT.out_file);
    end

    function plotCurrentEpisode()
        if ~isvalid(fig); return; end
        clf(fig);

        ep = episodes(curIdx);
        core_start = ep.start_sec;
        core_end   = ep.end_sec;

        % context window
        win_start = max(0, core_start - OPT.pre_sec);
        win_end   = min(total_dur_sec, core_end + OPT.post_sec);

        % sample indices for EEG/EMG
        i0 = max(1, floor(win_start * fs) + 1);
        i1 = min(numel(eeg), ceil(win_end * fs));

        t_win = t_eeg(i0:i1);
        t_rel = t_win - core_start;   % relative to episode onset
        eeg_w = eeg(i0:i1);
        emg_w = emg(i0:i1);

        % core indices (for PSD)
        core_i0 = max(i0, floor(core_start*fs) + 1);
        core_i1 = min(i1, ceil(core_end*fs));
        eeg_core = eeg(core_i0:core_i1);

        % 1) Hypnogram in window
        subplot(5,1,1); hold on;

        % select scores within window (+ a bit of margin)
        in_win = (t_scores >= win_start) & (t_scores <= win_end);
        t_scores_win = t_scores(in_win) - core_start;
        score_win    = score(in_win);

        if ~isempty(t_scores_win)
            stairs(t_scores_win, score_win, 'k','LineWidth',1);
            ylimits = [min(score_win)-0.5, max(score_win)+0.5];
        else
            ylimits = [-0.5 0.5];
        end
        ylim(ylimits);

        % highlight core episode window
        patch([0, core_end-core_start, core_end-core_start, 0], ...
              [ylimits(1), ylimits(1), ylimits(2), ylimits(2)], ...
              [0.9 0.9 1.0], 'FaceAlpha',0.4,'EdgeColor','none');

        xline(0,'b--','LineWidth',1);
        xline(core_end-core_start,'b--','LineWidth',1);

        ylabel('state code');
        title(sprintf('%s episode %d/%d  (flag = %s)', ...
            ep.type, curIdx, nEpisodes, flagStr(ep.flag)));
        grid on;

        % 2) EEG trace
        subplot(5,1,2);
        plot(t_rel, eeg_w);
        xline(0,'b--','LineWidth',1);
        xline(core_end-core_start,'b--','LineWidth',1);
        ylabel('EEG');
        grid on;
        title('Raw EEG');

        % 3) EMG trace
        subplot(5,1,3);
        plot(t_rel, emg_w);
        xline(0,'b--','LineWidth',1);
        xline(core_end-core_start,'b--','LineWidth',1);
        ylabel('EMG');
        grid on;
        title('Raw EMG');

        % 4) Morlet time–frequency map
        subplot(5,1,4);
        try
            [wt, f] = cwt(double(eeg_w), fs, 'amor', 'FrequencyLimits',[0.5 30]);
            powTF = abs(wt).^2;
            imagesc(t_rel, f, powTF);
            ylim([0 15]);
            axis xy;
            colormap(turbo);
            colorbar;
            xline(0,'w--','LineWidth',1);
            xline(core_end-core_start,'w--','LineWidth',1);
            ylabel('Hz');
            title('Morlet time–frequency power (EEG)');
        catch
            text(0.5,0.5,'cwt unavailable (Wavelet Toolbox?)','HorizontalAlignment','center');
            axis off;
        end

        % 5) PSD inside core episode
        subplot(5,1,5);
        if numel(eeg_core) > fs  % at least 1 s
            [Pxx, f_psd] = pwelch(double(eeg_core), hamming(round(fs*2)), [], [], fs);
            plot(f_psd, 10*log10(Pxx));
            xlim([0.5 40]);
            xlabel('Hz');
            ylabel('Power (dB)');
            grid on;
            title('PSD inside episode');
        else
            text(0.5,0.5,'Episode too short for PSD','HorizontalAlignment','center');
            axis off;
        end

        linkaxes(findall(fig,'Type','axes'),'x');
        xlim([t_rel(1), t_rel(end)]);
    end

end

function s = flagStr(f)
    if isnan(f)
        s = 'NaN';
    elseif f == 1
        s = 'GOOD';
    elseif f == 0
        s = 'BAD';
    else
        s = sprintf('%.2f', f);
    end
end
