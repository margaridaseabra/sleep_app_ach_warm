function [COND, STATS] = ach_compare_conditions_onsets(OUT_list, cond_labels, t_pre, t_post, dt)
% ACH_COMPARE_CONDITIONS_ONSETS
% -------------------------------------------------------------
% For one mouse, compare ACh onset dynamics across CONDITIONS
% (e.g. baseline / ambtemp / drugs).
%
% INPUTS
%   OUT_list    : array or cell array of OUT structs returned by
%                 plot_ach_eeg_segment (same mouse).
%   cond_labels : cellstr, one condition per OUT
%                 (e.g. {'baseline','baseline','ambtemp','drugs',...})
%   t_pre       : seconds before onset for peri traces (default 20)
%   t_post      : seconds after onset for peri traces (default 40)
%   dt          : step for common time axis (default 0.2 s)
%
% OUTPUTS
%   COND  : struct per condition and state:
%             COND(c).name
%                    .state(s).name
%                               .t_rel
%                               .traces     (events x time)
%                               .mean       (mean trace)
%                               .sem        (SEM trace)
%                               .nEvents
%                               .deltaF_all
%                               .slope_all
%
%   STATS : struct with ANOVA p-values across conditions:
%             STATS.(state).deltaF.p
%             STATS.(state).slope.p
%
% NOTE: This function does **no plotting**. All figures are made by
%       ach_plot_condition_summary_allStates.m

if nargin < 3 || isempty(t_pre),  t_pre  = 20; end
if nargin < 4 || isempty(t_post), t_post = 40; end
if nargin < 5 || isempty(dt),     dt     = 0.2; end

% Allow struct array OR cell array
if ~iscell(OUT_list)
    OUT_list = num2cell(OUT_list);
end
nSeg = numel(OUT_list);

if numel(cond_labels) ~= nSeg
    error('cond_labels must have one entry per OUT struct.');
end

% Unique conditions (order of appearance)
[condNames,~,condIdx_all] = unique(cond_labels,'stable');
nCond = numel(condNames);

% States to analyse
states = {'Wake','NREM','REM'};

% Common peri-onset time axis
t_rel_common = -t_pre:dt:t_post;

% -------------------------------------------------------------
% Initialise COND struct
% -------------------------------------------------------------
COND = struct();
for c = 1:nCond
    COND(c).name = condNames{c};
    for s = 1:numel(states)
        COND(c).state(s).name       = states{s};
        COND(c).state(s).t_rel      = t_rel_common;
        COND(c).state(s).traces     = [];   % events x time
        COND(c).state(s).deltaF_all = [];
        COND(c).state(s).slope_all  = [];
        COND(c).state(s).mean       = [];
        COND(c).state(s).sem        = [];
        COND(c).state(s).nEvents    = 0;
    end
end

% -------------------------------------------------------------
% Helper: stable transitions
% -------------------------------------------------------------
    function idx_keep = find_stable_transitions_local(score, prev_code, next_code, ...
                                                      minPreBins, minPostBins)
        n = numel(score);
        idx_candidates = find(score(2:end) == next_code & ...
                              score(1:end-1) == prev_code) + 1;
        idx_keep = [];
        for kk = 1:numel(idx_candidates)
            i = idx_candidates(kk);
            i_pre_start  = max(1, i - minPreBins);
            pre_segment  = score(i_pre_start:i-1);
            i_post_end   = min(n, i + minPostBins - 1);
            post_segment = score(i:i_post_end);
            if numel(pre_segment)  < minPreBins,  continue; end
            if numel(post_segment) < minPostBins, continue; end
            if all(pre_segment == prev_code) && all(post_segment == next_code)
                idx_keep(end+1,1) = i; %#ok<AGROW>
            end
        end
    end

% -------------------------------------------------------------
% Build peri-onset traces + feature pools by condition
% -------------------------------------------------------------
for seg = 1:nSeg
    OUTe   = OUT_list{seg};
    condId = condIdx_all(seg);

    % ---- get file info ----
    mat_file   = OUTe.mat_file;
    scores_csv = OUTe.scores_csv;
    t_start    = OUTe.t_start;
    t_end      = OUTe.t_end;
    codes      = OUTe.codes;

    % ---- load ACh ----
    info  = whos('-file', mat_file);
    names = {info.name};
    pick  = @(cands) cands{find(ismember(cands,names),1,'first')};

    ach_name    = pick({'ach','ACh','ne','dff','dFF','dff_ach'});
    fs_ach_name = pick({'ach_frequency','ne_frequency','fs_ach','Fs_ach'});
    if isempty(ach_name) || isempty(fs_ach_name)
        warning('No ACh or fs in %s, skipping.', mat_file);
        continue;
    end
    S     = load(mat_file, ach_name, fs_ach_name);
    ach   = S.(ach_name)(:);
    fs_ach= S.(fs_ach_name);

    % ---- load scores ----
    Msc = readmatrix(scores_csv);
    if size(Msc,2) == 1
        score     = Msc(:,1);
        epoch_sec = 1;
        t_scores  = (0:numel(score)-1)' * epoch_sec;
    else
        t_scores  = Msc(:,1);
        score     = Msc(:,2);
        dt_s      = diff(t_scores);
        epoch_sec = mode(dt_s(~isnan(dt_s)));
    end
    score    = score(:);
    t_scores = t_scores(:);

    CODE_WAKE = codes.WK;
    CODE_NREM = codes.NREM;
    CODE_REM  = codes.REM;

    minPreBins  = round(10 / epoch_sec);
    minPostBins = round(20 / epoch_sec);

    % ---- loop states ----
    for s = 1:numel(states)
        stName = states{s};

        switch stName
            case 'Wake'
                idx1 = find_stable_transitions_local(score, CODE_NREM, CODE_WAKE, minPreBins, minPostBins);
                idx2 = find_stable_transitions_local(score, CODE_REM,  CODE_WAKE, minPreBins, minPostBins);
                idx_on = sort([idx1; idx2]);
            case 'NREM'
                idx_on = find_stable_transitions_local(score, CODE_WAKE, CODE_NREM, minPreBins, minPostBins);
            case 'REM'
                idx_on = find_stable_transitions_local(score, CODE_NREM, CODE_REM,  minPreBins, minPostBins);
        end

        if isempty(idx_on), continue; end

        t_onsets = (idx_on-1) * epoch_sec;
        t_onsets = t_onsets(t_onsets >= t_start & t_onsets <= t_end);
        if isempty(t_onsets)
            continue;
        end

        % ---- peri-event matrix for this segment/state ----
        [Mperi, t_rel] = make_peri_event_matrix(ach, fs_ach, t_onsets, t_pre, t_post);
        if isempty(Mperi), continue; end

        % interpolate each event to common axis
        M_common = nan(size(Mperi,1), numel(t_rel_common));
        for ev = 1:size(Mperi,1)
            M_common(ev,:) = interp1(t_rel, Mperi(ev,:), t_rel_common, 'linear', NaN);
        end

        % append traces
        COND(condId).state(s).traces = [COND(condId).state(s).traces; M_common];

        % append features from OUT (already computed in plot_ach_eeg_segment)
        if isfield(OUTe,'features') && isfield(OUTe.features, stName)
            F = OUTe.features.(stName);
            if isfield(F,'deltaF_all') && ~isempty(F.deltaF_all)
                COND(condId).state(s).deltaF_all = [COND(condId).state(s).deltaF_all; F.deltaF_all(:)];
            end
            if isfield(F,'slope_all') && ~isempty(F.slope_all)
                COND(condId).state(s).slope_all  = [COND(condId).state(s).slope_all;  F.slope_all(:)];
            end
        end
    end
end

% -------------------------------------------------------------
% Convert traces to mean/SEM and compute stats
% -------------------------------------------------------------
STATS = struct();

for s = 1:numel(states)
    stName = states{s};

    % mean/SEM per condition
    for c = 1:nCond
        M = COND(c).state(s).traces;
        if ~isempty(M)
            COND(c).state(s).mean    = mean(M,1,'omitnan');
            COND(c).state(s).sem     = std(M,[],1,'omitnan') ./ sqrt(size(M,1));
            COND(c).state(s).nEvents = size(M,1);
        else
            COND(c).state(s).mean    = nan(1,numel(t_rel_common));
            COND(c).state(s).sem     = nan(1,numel(t_rel_common));
            COND(c).state(s).nEvents = 0;
        end
    end

    % pooled features for stats
    all_dF    = [];
    all_slope = [];
    groups    = [];

    for c = 1:nCond
        dF = COND(c).state(s).deltaF_all;
        sl = COND(c).state(s).slope_all;
        all_dF    = [all_dF;    dF(:)];
        all_slope = [all_slope; sl(:)];
        groups    = [groups;    repmat(c,numel(dF),1)];
    end

    st = struct();
    if isempty(all_dF) || numel(unique(groups)) < 2
        st.deltaF.p = NaN;
        st.slope.p  = NaN;
    else
        % one-way ANOVA across whatever conditions actually have events
        st.deltaF.p = anova1(all_dF, groups, 'off');
        st.slope.p  = anova1(all_slope, groups, 'off');
    end

    STATS.(stName) = st;
end
end
