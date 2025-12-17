function STATS = compare_ach_segments_onsets_features(OUT_list, labels, t_pre_peri, t_post_peri)
% COMPARE_ACH_SEGMENTS_ONSETS_FEATURES
%   Compare ACh onset dynamics + features across multiple segments.
%
%   OUT_list : array or cell array of OUT structs returned by
%              plot_ach_eeg_segment (we use .mat_file, .scores_csv,
%              .t_start, .t_end and .features).
%   labels   : cellstr, one label per segment (e.g. {'BL-early','BL-late',...})
%   t_pre_peri, t_post_peri : peri-onset window for traces (s, default 20 / 40)
%
% For each state (Wake, NREM, REM):
%   - computes peri-onset mean traces for each segment
%   - overlays them in one figure
%   - makes boxplots of ΔF/F and slope per segment
%   - runs t-test (2 segments) or ANOVA (>2 segments)
%
% STATS.(state).deltaF / .slope contain p-values etc.

if nargin < 3 || isempty(t_pre_peri),  t_pre_peri  = 20; end
if nargin < 4 || isempty(t_post_peri), t_post_peri = 40; end

if ~iscell(OUT_list)
    OUT_list = num2cell(OUT_list);
end
nSeg = numel(OUT_list);

if nargin < 2 || isempty(labels)
    labels = arrayfun(@(k) sprintf('Seg%d',k), 1:nSeg, 'uni',0);
end

states = {'Wake','NREM','REM'};
STATS  = struct();

for s = 1:numel(states)
    stName = states{s};

    % -----------------------------------------------------------------
    % 1) Collect features (deltaF, slope) from OUT.features
    % -----------------------------------------------------------------
    all_dF    = [];
    all_slope = [];
    grp       = [];

    for e = 1:nSeg
        F = OUT_list{e}.features.(stName);
        if isempty(F) || F.n == 0, continue; end
        all_dF    = [all_dF;    F.deltaF_all(:)];
        all_slope = [all_slope; F.slope_all(:)];
        grp       = [grp;       repmat(e, F.n, 1)];
    end

    % -----------------------------------------------------------------
    % 2) Build peri-onset mean traces for each segment
    % -----------------------------------------------------------------
    t_rel_common = [];
    mu_traces    = nan(nSeg, 1);  % will resize
    haveTrace    = false(nSeg,1);

    for e = 1:nSeg
        OUTe = OUT_list{e};

        % reload ACh + scores (we reuse logic from plot_ach_eeg_segment)
        mat_file   = OUTe.mat_file;
        scores_csv = OUTe.scores_csv;
        t_start    = OUTe.t_start;
        t_end      = OUTe.t_end;

        % --- load signals ---
        info = whos('-file', mat_file);
        names = {info.name};
        pick = @(cands) cands{find(ismember(cands, names),1,'first')};
        ach_name    = pick({'ach','ACh','ne','dff','dFF','dff_ach'});
        fs_ach_name = pick({'ach_frequency','ne_frequency','fs_ach','Fs_ach'});
        if isempty(ach_name) || isempty(fs_ach_name), continue; end
        S = load(mat_file, ach_name, fs_ach_name);
        ach    = S.(ach_name)(:);
        fs_ach = S.(fs_ach_name);

        % --- scores ---
        M = readmatrix(scores_csv);
        if size(M,2) == 1
            score     = M(:,1);
            epoch_sec = 1;
            t_scores  = (0:numel(score)-1)' * epoch_sec;
        else
            t_scores  = M(:,1);
            score     = M(:,2);
            dt = diff(t_scores);
            epoch_sec = mode(dt(~isnan(dt)));
        end
        score   = score(:);
        t_scores = t_scores(:);

        % --- canonical, stable transitions ---
        CODES = OUTe.codes;  % save this field in plot_ach_eeg_segment if not already
        CODE_WAKE = CODES.WK;
        CODE_NREM = CODES.NREM;
        CODE_REM  = CODES.REM;

        minPreBins  = round(10 / epoch_sec);
        minPostBins = round(20 / epoch_sec);

        switch stName
            case 'Wake'
                idx1 = find_stable_transitions_local(score, CODE_NREM, CODE_WAKE, minPreBins, minPostBins);
                idx2 = find_stable_transitions_local(score, CODE_REM,  CODE_WAKE, minPreBins, minPostBins);
                idx  = sort([idx1; idx2]);
            case 'NREM'
                idx  = find_stable_transitions_local(score, CODE_WAKE, CODE_NREM, minPreBins, minPostBins);
            case 'REM'
                idx  = find_stable_transitions_local(score, CODE_NREM, CODE_REM,  minPreBins, minPostBins);
        end

        t_onsets = (idx-1)*epoch_sec;
        t_onsets = t_onsets(t_onsets >= t_start & t_onsets <= t_end);

        if isempty(t_onsets)
            continue;
        end

        % --- peri matrix for this segment ---
        [Mperi, t_rel] = make_peri_event_matrix(ach, fs_ach, t_onsets, t_pre_peri, t_post_peri);

        if isempty(t_rel_common)
            % define common time axis (0.1 s step)
            t_rel_common = -t_pre_peri:0.1:t_post_peri;
        end

        mu_evt = mean(Mperi,1,'omitnan');
        mu_interp = interp1(t_rel, mu_evt, t_rel_common, 'linear', NaN);

        if size(mu_traces,2) ~= numel(t_rel_common)
            mu_traces = nan(nSeg, numel(t_rel_common));
        end
        mu_traces(e,:) = mu_interp;
        haveTrace(e) = true;
    end

    % -----------------------------------------------------------------
    % 3) Plot onset traces + feature boxplots
    % -----------------------------------------------------------------
    if any(haveTrace)
        figure('Color','w','Position',[200 200 900 500]);

        % ---- onset traces ----
        subplot(2,2,[1 2]); hold on;
        cols = lines(nSeg);
        leg = {};
        for e = 1:nSeg
            if ~haveTrace(e), continue; end
            plot(t_rel_common, mu_traces(e,:), 'LineWidth',2, 'Color', cols(e,:));
            leg{end+1} = labels{e}; %#ok<AGROW>
        end
        plot([0 0], ylim, 'k--');
        xlabel('Time from onset (s)');
        ylabel('\DeltaF/F ACh');
        title(sprintf('%s onset: mean traces', stName));
        if ~isempty(leg)
            legend(leg, 'Location','bestoutside');
        end
        grid on;

        % ---- ΔF/F boxplot ----
        subplot(2,2,3);
        if ~isempty(all_dF)
            boxplot(all_dF, grp, 'Labels', labels);
            ylabel('\DeltaF/F');
            title(sprintf('%s: \\DeltaF/F at onset', stName));
        else
            text(0.5,0.5,'No events','Units','normalized','HorizontalAlignment','center');
        end

        % ---- slope boxplot ----
        subplot(2,2,4);
        if ~isempty(all_slope)
            boxplot(all_slope, grp, 'Labels', labels);
            ylabel('Slope (\DeltaF/F per s)');
            title(sprintf('%s: slopes at onset', stName));
        else
            text(0.5,0.5,'No events','Units','normalized','HorizontalAlignment','center');
        end

        sgtitle(sprintf('%s comparison across %d segments', stName, nSeg), ...
                'FontWeight','bold');
    end

    % -----------------------------------------------------------------
    % 4) Stats on features
    % -----------------------------------------------------------------
    ST = struct();

    % ΔF/F
    if numel(unique(grp)) > 1
        if nSeg == 2
            g1 = all_dF(grp==1);
            g2 = all_dF(grp==2);
            [~,p,~,stats] = ttest2(g1, g2, 'Vartype','unequal');
            ST.deltaF.test  = 'ttest2';
            ST.deltaF.p     = p;
            ST.deltaF.stats = stats;
        else
            [p,tbl,astats] = anova1(all_dF, grp, 'off');
            ST.deltaF.test  = 'anova1';
            ST.deltaF.p     = p;
            ST.deltaF.tbl   = tbl;
            ST.deltaF.stats = astats;
        end
    else
        ST.deltaF.test = 'none'; ST.deltaF.p = NaN; ST.deltaF.stats = [];
    end

    % slopes
    if numel(unique(grp)) > 1
        if nSeg == 2
            g1 = all_slope(grp==1);
            g2 = all_slope(grp==2);
            [~,p,~,stats] = ttest2(g1, g2, 'Vartype','unequal');
            ST.slope.test  = 'ttest2';
            ST.slope.p     = p;
            ST.slope.stats = stats;
        else
            [p,tbl,astats] = anova1(all_slope, grp, 'off');
            ST.slope.test  = 'anova1';
            ST.slope.p     = p;
            ST.slope.tbl   = tbl;
            ST.slope.stats = astats;
        end
    else
        ST.slope.test = 'none'; ST.slope.p = NaN; ST.slope.stats = [];
    end

    STATS.(stName) = ST;
end
end

% --- local stable-transition helper (stand-alone version) ---------------
function idx_keep = find_stable_transitions_local(score, prev_code, next_code, ...
                                            minPreBins, minPostBins)
n = numel(score);
idx_candidates = find(score(2:end) == next_code & ...
                      score(1:end-1) == prev_code) + 1;
idx_keep = [];
for k = 1:numel(idx_candidates)
    i = idx_candidates(k);
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
