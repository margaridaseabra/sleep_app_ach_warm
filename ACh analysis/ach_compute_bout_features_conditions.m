function [BOUT, STATS] = ach_compute_bout_features_conditions(OUT_list, cond_labels, varargin)
% ACH_COMPUTE_BOUT_FEATURES_CONDITIONS
% -------------------------------------------------------------
% Bout-level ACh + PLV features per CONDITION and STATE.
%
% INPUTS
%   OUT_list    : struct array or cell array of OUT from plot_ach_eeg_segment
%                 (must contain fields: .mat_file, .scores_csv, .t_start,
%                 .t_end, .codes)
%   cond_labels : cellstr, one per OUT (e.g. {'baseline','baseline',
%                                            'ambtemp','ambtemp','drugs',...})
%
% OPTIONAL NAME/VALUE:
%   'plv_eeg_band' : [f1 f2] Hz for EEG phase     (default [5 10])
%   'plv_ach_band' : [f1 f2] Hz for ACh phase     (default [0.05 0.5])
%   'minBoutDur'   : minimum bout duration (s)    (default 2)
%
% OUTPUTS
%   BOUT(c).name            - condition name (e.g. 'baseline')
%        .state(s).name     - 'Wake','NREM','REM'
%                 .dur_all       - bout durations (s)
%                 .meanACh_all   - mean ACh in bout
%                 .slope_all     - linear slope of ACh in bout (ΔF/F/s)
%                 .amp_all       - max-min ACh in bout
%                 .plv_all       - EEG–ACh PLV in the specified bands
%                 .context_prev  - cellstr, bout context wrt previous state
%                 .context_next  - cellstr, bout context wrt next state
%
%   STATS.(state).meanACh.p
%         .slope.p
%         .amp.p
%         .dur.p
%         .plv.p     - one-way ANOVA p across conditions (NaN if not enough data)

% ---------- parse inputs ----------
p = inputParser;
p.addParameter('plv_eeg_band', [5 10]);      % theta-ish
p.addParameter('plv_ach_band', [0.05 0.5]);  % slow ACh
p.addParameter('minBoutDur',   2);           % s
p.parse(varargin{:});
opt = p.Results;

if ~iscell(OUT_list)
    OUT_list = num2cell(OUT_list);
end
nSeg = numel(OUT_list);

if numel(cond_labels) ~= nSeg
    error('cond_labels must have one entry per OUT struct.');
end

[condNames, ~, condIdx_all] = unique(cond_labels,'stable');
nCond = numel(condNames);

states = {'Wake','NREM','REM'};

% ---------- initialise BOUT struct ----------
BOUT = struct([]);
for c = 1:nCond
    BOUT(c).name = condNames{c};
    for s = 1:numel(states)
        BOUT(c).state(s).name         = states{s};
        BOUT(c).state(s).dur_all      = [];
        BOUT(c).state(s).meanACh_all  = [];
        BOUT(c).state(s).slope_all    = [];
        BOUT(c).state(s).amp_all      = [];
        BOUT(c).state(s).plv_all      = [];
        BOUT(c).state(s).context_prev = {};   % cell arrays of strings
        BOUT(c).state(s).context_next = {};
    end
end

% =============================================================
%  Loop over segments / episodes
% =============================================================
for seg = 1:nSeg
    OUTe   = OUT_list{seg};
    condId = condIdx_all(seg);

    if ~isfield(OUTe,'mat_file') || ~isfield(OUTe,'scores_csv')
        warning('OUT(%d) missing mat_file or scores_csv, skipping.', seg);
        continue;
    end

    mat_file   = OUTe.mat_file;
    scores_csv = OUTe.scores_csv;

    if isfield(OUTe,'t_start'), t_start = OUTe.t_start; else, t_start = 0;   end
    if isfield(OUTe,'t_end'),   t_end   = OUTe.t_end;   else, t_end   = Inf; end

    if ~isfield(OUTe,'codes')
        error('OUT(%d) has no .codes field. Re-run plot_ach_eeg_segment so it stores CODES.', seg);
    end
    codes = OUTe.codes;

    CODE_WAKE = codes.WK;
    CODE_NREM = codes.NREM;
    CODE_REM  = codes.REM;

    % ---------- load ACh & EEG ----------
    info  = whos('-file', mat_file);
    names = {info.name};
    pick  = @(cands) cands{find(ismember(cands,names),1,'first')};

    ach_name    = pick_if_present(pick, {'ach','ACh','ne','dff','dFF','dff_ach'});
    fs_ach_name = pick_if_present(pick, {'ach_frequency','ne_frequency','fs_ach','Fs_ach'});
    eeg_name    = pick_if_present(pick, {'eeg','EEG','eeg1','Eeg','eeg_filt'});
    fs_eeg_name = pick_if_present(pick, {'eeg_frequency','fs_eeg','Fs_eeg','EEG_frequency'});

    if isempty(ach_name) || isempty(fs_ach_name)
        warning('No ACh or fs_ach in %s, skipping segment %d.', mat_file, seg);
        continue;
    end

    S  = load(mat_file, ach_name, fs_ach_name);
    ach    = S.(ach_name)(:);
    fs_ach = S.(fs_ach_name);

    have_eeg = ~(isempty(eeg_name) || isempty(fs_eeg_name));
    eeg = [];
    fs_eeg = [];
    if have_eeg
        S2     = load(mat_file, eeg_name, fs_eeg_name);
        eeg    = S2.(eeg_name)(:);
        fs_eeg = S2.(fs_eeg_name);
        eeg    = apply_50Hz_notch_local(eeg, fs_eeg, mat_file);
    end

    % ---------- load scores ----------
    M = readmatrix(scores_csv);
    if size(M,2) == 1
        score     = M(:,1);
        epoch_sec = 1;
        t_scores  = (0:numel(score)-1)' * epoch_sec;
    else
        t_scores  = M(:,1);
        score     = M(:,2);
        dt        = diff(t_scores);
        epoch_sec = mode(dt(~isnan(dt)));
    end
    score    = score(:);
    t_scores = t_scores(:);

    % ---------- build bouts over entire scoring ----------
    svec       = score;
    run_starts = [1; find(diff(svec)~=0)+1];
    run_ends   = [run_starts(2:end)-1; numel(svec)];

    for j = 1:numel(run_starts)
        code = svec(run_starts(j));

        if ~(code == CODE_WAKE || code == CODE_NREM || code == CODE_REM)
            continue;
        end

        prevCode = NaN;
        nextCode = NaN;
        if j > 1,  prevCode = svec(run_starts(j)-1);    end
        if j < numel(run_starts), nextCode = svec(run_starts(j+1)); end

        t0 = t_scores(run_starts(j));
        t1 = t_scores(run_ends(j)) + epoch_sec;

        % intersect with desired window
        if t1 <= t_start || t0 >= t_end
            continue;
        end
        t0c = max(t0, t_start);
        t1c = min(t1, t_end);
        dur = t1c - t0c;
        if dur < opt.minBoutDur
            continue;
        end

        % which state index?
        if code == CODE_WAKE
            stateIdx = 1;
        elseif code == CODE_NREM
            stateIdx = 2;
        else
            stateIdx = 3;
        end

        % ---- ACh segment for this bout ----
        iA0 = max(1, floor(t0c*fs_ach) + 1);
        iA1 = min(numel(ach), floor(t1c*fs_ach));
        if iA1 <= iA0 + 1
            continue;
        end
        segA = ach(iA0:iA1);
        if all(isnan(segA)), continue; end

        dur_s = (iA1 - iA0) / fs_ach;
        t_rel = (0:(numel(segA)-1))' / fs_ach;

        meanA = mean(segA,'omitnan');
        amp   = max(segA) - min(segA);

        X      = [t_rel ones(numel(t_rel),1)];
        b      = X \ double(segA(:));
        slope  = b(1);

        % ---- PLV for this bout (EEG–ACh) ----
        plv_val = NaN;
        if have_eeg
            iE0 = max(1, floor(t0c*fs_eeg) + 1);
            iE1 = min(numel(eeg), floor(t1c*fs_eeg));
            if iE1 > iE0 + 10
                segE = eeg(iE0:iE1);
                plv_val = compute_plv_simple(segE, fs_eeg, segA, fs_ach, ...
                                             opt.plv_eeg_band, opt.plv_ach_band);
            end
        end

        % ---- context labelling ----
        context_prev = 'other';
        context_next = 'other';

        if code == CODE_NREM
            % who comes after NREM?
            if nextCode == CODE_REM
                context_next = 'preREM';
            elseif nextCode == CODE_WAKE
                context_next = 'preWake';
            end
            % who was before NREM?
            if prevCode == CODE_WAKE
                context_prev = 'postWake';
            elseif prevCode == CODE_REM
                context_prev = 'postREM';
            end

        elseif code == CODE_WAKE
            if prevCode == CODE_NREM
                context_prev = 'postNREM';
            elseif prevCode == CODE_REM
                context_prev = 'postREM';
            end
            if nextCode == CODE_NREM
                context_next = 'preNREM';
            elseif nextCode == CODE_REM
                context_next = 'preREM';
            end

        elseif code == CODE_REM
            if prevCode == CODE_NREM
                context_prev = 'postNREM';
            elseif prevCode == CODE_WAKE
                context_prev = 'postWake';
            end
            if nextCode == CODE_WAKE
                context_next = 'preWake';
            elseif nextCode == CODE_NREM
                context_next = 'preNREM';
            end
        end

        % ---- append to BOUT struct ----
        Sstate = BOUT(condId).state(stateIdx);

        Sstate.dur_all(end+1,1)     = dur_s;
        Sstate.meanACh_all(end+1,1) = meanA;
        Sstate.slope_all(end+1,1)   = slope;
        Sstate.amp_all(end+1,1)     = amp;
        Sstate.plv_all(end+1,1)     = plv_val;
        Sstate.context_prev{end+1,1} = context_prev;
        Sstate.context_next{end+1,1} = context_next;

        BOUT(condId).state(stateIdx) = Sstate;
    end
end

% =============================================================
%  Stats across conditions: one-way ANOVA per state & metric
% =============================================================
STATS = struct();
for s = 1:numel(states)
    stName = states{s};

    all_mean = []; g_mean = [];
    all_slope= []; g_slope= [];
    all_amp  = []; g_amp  = [];
    all_dur  = []; g_dur  = [];
    all_plv  = []; g_plv  = [];

    for c = 1:nCond
        S = BOUT(c).state(s);

        all_mean = [all_mean; S.meanACh_all(:)];
        g_mean   = [g_mean;   c*ones(numel(S.meanACh_all),1)];

        all_slope= [all_slope; S.slope_all(:)];
        g_slope  = [g_slope;   c*ones(numel(S.slope_all),1)];

        all_amp  = [all_amp; S.amp_all(:)];
        g_amp    = [g_amp;   c*ones(numel(S.amp_all),1)];

        all_dur  = [all_dur; S.dur_all(:)];
        g_dur    = [g_dur;   c*ones(numel(S.dur_all),1)];

        all_plv  = [all_plv; S.plv_all(:)];
        g_plv    = [g_plv;   c*ones(numel(S.plv_all),1)];
    end

    st = struct();

    % helper for ANOVA
    st.meanACh.p = local_anova(all_mean, g_mean);
    st.slope.p   = local_anova(all_slope, g_slope);
    st.amp.p     = local_anova(all_amp,   g_amp);
    st.dur.p     = local_anova(all_dur,   g_dur);
    st.plv.p     = local_anova(all_plv,   g_plv);

    STATS.(stName) = st;
end
end

% =============================================================
%  Helper functions (local)
% =============================================================

function name = pick_if_present(pick, cands)
    try
        name = pick(cands);
    catch
        name = [];
    end
end

function p = local_anova(vals, groups)
    vals   = vals(:);
    groups = groups(:);
    mask = ~isnan(vals);
    vals   = vals(mask);
    groups = groups(mask);

    if numel(vals) < 2
        p = NaN;
        return;
    end
    if numel(unique(groups)) < 2
        p = NaN;
        return;
    end
    p = anova1(vals, groups, 'off');
end

function eeg_out = apply_50Hz_notch_local(eeg_in, fs_eeg, mat_file)
    eeg_out = eeg_in;
    if fs_eeg <= 120
        return;
    end

    [~, base, ext] = fileparts(mat_file);
    fname = [base ext];

    skip_list = { ...
        '20251001_baseline_mouse1_APP.mat', ...
        '20251002_baseline_mouse2_WT.mat', ...
        '20251003_ambtemp_mouse1_APP.mat', ...
        '20251005_baseline_mouse8_WT.mat', ...
        '20251006_baseline_mouse4_WT.mat' ...
        };

    if ismember(fname, skip_list)
        return;
    end

    wo = 50 / (fs_eeg/2);
    bw = wo / 35;    % ~Q=35
    [bN,aN] = iirnotch(wo, bw);
    eeg_out = filtfilt(bN,aN,double(eeg_in));
end

function plv = compute_plv_simple(eeg_seg, fs_eeg, ach_seg, fs_ach, eeg_band, ach_band)
    plv = NaN;

    eeg_seg = double(eeg_seg(:));
    ach_seg = double(ach_seg(:));

    if numel(eeg_seg) < 50 || numel(ach_seg) < 50
        return;
    end

    % ---- EEG band-pass ----
    WnE = eeg_band / (fs_eeg/2);
    WnE(WnE <= 0) = 0.01;
    WnE(WnE >= 1) = 0.99;
    [bE,aE] = butter(4, WnE, 'bandpass');
    eeg_bp  = filtfilt(bE,aE,eeg_seg);

    % ---- ACh band-pass ----
    WnA = ach_band / (fs_ach/2);
    WnA(WnA <= 0) = 0.01;
    WnA(WnA >= 1) = 0.99;
    [bA,aA] = butter(4, WnA, 'bandpass');
    ach_bp  = filtfilt(bA,aA,ach_seg);

    % ---- Hilbert phase & PLV ----
    phiE = angle(hilbert(eeg_bp));
    phiA = angle(hilbert(ach_bp));

    n = min(numel(phiE), numel(phiA));
    if n < 20
        plv = NaN;
        return;
    end

    dphi = phiE(1:n) - phiA(1:n);
    plv  = abs(mean(exp(1i*dphi)));
end
