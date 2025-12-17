function OUT_TRANS = ach_peri_transitions_single(mat_file, scores_csv, codes, varargin)
% ACH_PERI_TRANSITIONS_SINGLE
%   Event-triggered ACh analysis for ONE recording.
%
%   Uses *canonical, stable* transitions:
%     - Wake:  sleep (NREM or REM)  -> Wake
%     - NREM:  Wake                 -> NREM
%     - REM:   NREM                 -> REM
%
%   Each transition must have at least MIN_PRE seconds of the
%   previous state before 0, and MIN_POST seconds of the new
%   state after 0 (defaults: 10 s pre, 20 s post).
%
% REQUIRED
%   mat_file   : path to ACh .mat
%   scores_csv : path to 1-Hz scores CSV
%   codes      : struct with fields .WK .NREM .REM .MA
%
% OPTIONAL name/value
%   't_pre'         : seconds before onset  (default 30)
%   't_post'        : seconds after onset   (default 60)
%   'min_pre'       : min stable pre-state duration (s, default 10)
%   'min_post'      : min stable post-state duration (s, default 20)
%   'out_prefix'    : used in filename      (default '')
%   'out_dir'       : folder to save figs   (default = fileparts(mat_file))
%   'session_label' : string for titles     (default = out_prefix)
%   'zscore_rows'   : z-score each event row in heatmap (default false)
%   'save_fig'      : whether to save PNG   (default true)

p = inputParser;
p.addParameter('t_pre', 30);
p.addParameter('t_post', 60);
p.addParameter('min_pre', 10);
p.addParameter('min_post', 20);
p.addParameter('out_prefix', '');
p.addParameter('out_dir', fileparts(mat_file));
p.addParameter('session_label', '');
p.addParameter('zscore_rows', false);
p.addParameter('save_fig', true);
p.parse(varargin{:});
opt = p.Results;

if isempty(opt.session_label)
    opt.session_label = opt.out_prefix;
end
if ~exist(opt.out_dir, 'dir')
    mkdir(opt.out_dir);
end

% -----------------------------------------------------------
% 1) Load ACh signal + sampling rate
% -----------------------------------------------------------
info = whos('-file', mat_file);
names = {info.name};

ach_candidates = {'ach','ACh','ne','dff','dFF','dff_ach'};
fs_candidates  = {'ach_frequency','ne_frequency','fs_ach','Fs_ach'};

ach_name = '';
for i = 1:numel(ach_candidates)
    if any(strcmp(ach_candidates{i}, names))
        ach_name = ach_candidates{i};
        break;
    end
end
fs_name = '';
for i = 1:numel(fs_candidates)
    if any(strcmp(fs_candidates{i}, names))
        fs_name = fs_candidates{i};
        break;
    end
end

if isempty(ach_name)
    error('Could not find ACh variable in %s', mat_file);
end
if isempty(fs_name)
    error('Could not find ACh sampling frequency in %s', mat_file);
end

S = load(mat_file, ach_name, fs_name);
ach    = S.(ach_name);
fs_ach = S.(fs_name);
ach = ach(:);   % ensure column

% -----------------------------------------------------------
% 2) Load 1-Hz scores
% -----------------------------------------------------------
M = readmatrix(scores_csv);

if size(M,2) == 1
    score     = M(:,1);
    epoch_sec = 1;
else
    time_col  = M(:,1);
    score     = M(:,2);
    dt = diff(time_col);
    epoch_sec = mode(dt(~isnan(dt)));
end
score = score(:);

% -----------------------------------------------------------
% 3) Find canonical, stable transitions
% -----------------------------------------------------------
CODE_WAKE = codes.WK;
CODE_NREM = codes.NREM;
CODE_REM  = codes.REM;

minPreBins  = round(opt.min_pre  / epoch_sec);
minPostBins = round(opt.min_post / epoch_sec);

% --- Wake onset: sleep (NREM or REM) -> Wake
idx_W_fromN = find_stable_transitions(score, CODE_NREM, CODE_WAKE, ...
    minPreBins, minPostBins);
idx_W_fromR = find_stable_transitions(score, CODE_REM,  CODE_WAKE, ...
    minPreBins, minPostBins);
idx_Wake    = sort([idx_W_fromN; idx_W_fromR]);

% --- NREM onset: Wake -> NREM
idx_NREM = find_stable_transitions(score, CODE_WAKE, CODE_NREM, ...
    minPreBins, minPostBins);

% --- REM onset: NREM -> REM
idx_REM  = find_stable_transitions(score, CODE_NREM, CODE_REM, ...
    minPreBins, minPostBins);

onset_idx_struct.Wake = idx_Wake;
onset_idx_struct.NREM = idx_NREM;
onset_idx_struct.REM  = idx_REM;

% Convert to seconds
onset_times_struct = struct();
onset_times_struct.Wake = (idx_Wake - 1) * epoch_sec;
onset_times_struct.NREM = (idx_NREM - 1) * epoch_sec;
onset_times_struct.REM  = (idx_REM  - 1) * epoch_sec;

% -----------------------------------------------------------
% 4) Build peri-event matrices per state
% -----------------------------------------------------------
OUT_TRANS = struct();
OUT_TRANS.mat_file   = mat_file;
OUT_TRANS.scores_csv = scores_csv;
OUT_TRANS.codes      = codes;
OUT_TRANS.t_pre      = opt.t_pre;
OUT_TRANS.t_post     = opt.t_post;
OUT_TRANS.min_pre    = opt.min_pre;
OUT_TRANS.min_post   = opt.min_post;

t_rel_common = [];

stateNames = {'Wake','NREM','REM'};

for sIdx = 1:numel(stateNames)
    stName = stateNames{sIdx};
    onset_times = onset_times_struct.(stName);

    if isempty(onset_times)
        OUT_TRANS.(stName).M           = [];
        OUT_TRANS.(stName).event_times = [];
        OUT_TRANS.(stName).onset_idx   = [];
        continue;
    end

    [Mperi, t_rel, event_times_valid] = make_peri_event_matrix( ...
        ach, fs_ach, onset_times, opt.t_pre, opt.t_post);

    OUT_TRANS.(stName).M           = Mperi;
    OUT_TRANS.(stName).event_times = event_times_valid;
    OUT_TRANS.(stName).onset_idx   = onset_idx_struct.(stName);  % in score bins

    if isempty(t_rel_common)
        t_rel_common = t_rel;
    end
end

OUT_TRANS.t_rel = t_rel_common;

% -----------------------------------------------------------
% 5) Plot per-state heatmaps + mean traces (if requested)
% -----------------------------------------------------------
if (opt.save_fig || nargout == 0) && ~isempty(t_rel_common)
    nStates = numel(stateNames);
    fig = figure('Color','w','Position',[100 100 900 250*nStates]);

    for sIdx = 1:nStates
        stName = stateNames{sIdx};
        Mperi  = OUT_TRANS.(stName).M;

        if isempty(Mperi)
            continue;
        end

        if opt.zscore_rows
            Mplot = (Mperi - mean(Mperi,2,'omitnan')) ./ ...
                     std(Mperi,[],2,'omitnan');
        else
            Mplot = Mperi;
        end

        mean_trace = mean(Mperi, 1, 'omitnan');
        sem_trace  = std(Mperi, [], 1, 'omitnan') ./ sqrt(size(Mperi,1));

        row = sIdx;

        % --- Heatmap ---
        subplot(nStates, 2, 2*row-1);
        imagesc(t_rel_common, 1:size(Mplot,1), Mplot);
        axis xy;
        xlabel('Time from onset (s)');
        ylabel('Event #');
        title(sprintf('%s onset (%d events)', stName, size(Mplot,1)));
        colormap(gca,'hot');
        colorbar;

        % --- Mean ± SEM trace ---
        subplot(nStates, 2, 2*row);
        hold on;
        patch([t_rel_common fliplr(t_rel_common)], ...
              [mean_trace - sem_trace, fliplr(mean_trace + sem_trace)], ...
              0.8*[1 1 1], 'EdgeColor','none','FaceAlpha',0.5);
        plot(t_rel_common, mean_trace, 'k', 'LineWidth',2);
        plot([0 0], ylim, 'r--');
        xlabel('Time from onset (s)');
        ylabel('\DeltaF/F ACh');
        title(sprintf('%s onset mean \x00b1 SEM', stName));
        grid on;
    end

    if ~isempty(opt.session_label)
        sgtitle(strrep(opt.session_label,'_','\_'));
    end

    if opt.save_fig
        if isempty(opt.out_prefix)
            [~,base,~] = fileparts(mat_file);
            opt.out_prefix = base;
        end
        figname = fullfile(opt.out_dir, [opt.out_prefix '_ACh_periOnsets_canonical.png']);
        saveas(fig, figname);
    end
end

end

% ===================================================================
% Helper: find transitions with stable pre/post state
% ===================================================================
function idx_keep = find_stable_transitions(score, prev_code, next_code, ...
                                            minPreBins, minPostBins)
% FIND_STABLE_TRANSITIONS
%   Find indices i where score(i-1)=prev_code, score(i)=next_code,
%   and there are >= minPreBins of prev_code before i and
%   >= minPostBins of next_code starting at i.

n = numel(score);
idx_candidates = find(score(2:end) == next_code & ...
                      score(1:end-1) == prev_code) + 1;

idx_keep = [];

for k = 1:numel(idx_candidates)
    i = idx_candidates(k);

    % pre window [i - minPreBins, i-1]
    i_pre_start = max(1, i - minPreBins);
    pre_segment = score(i_pre_start:i-1);

    % post window [i, i + minPostBins - 1]
    i_post_end  = min(n, i + minPostBins - 1);
    post_segment = score(i:i_post_end);

    if numel(pre_segment)   < minPreBins,  continue; end
    if numel(post_segment)  < minPostBins, continue; end

    if all(pre_segment  == prev_code) && ...
       all(post_segment == next_code)
        idx_keep(end+1,1) = i; %#ok<AGROW>
    end
end
end
