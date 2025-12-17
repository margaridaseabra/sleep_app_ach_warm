function STATS = ach_plv_plot_condition_summary(PLV_list, cond_labels, titleStr)
% ACH_PLV_PLOT_CONDITION_SUMMARY
% -------------------------------------------------------------
% Summarise PLV across conditions for ONE transition type:
%   - mean ± SEM PLV(t) per condition
%   - scalar PLV per episode (mean in [0, +20] s) with ANOVA
%
% INPUTS
%   PLV_list    : struct array OR cell array of structs, one per episode
%                 fields expected:
%                   .t_rel      : 1×T time axis (s, 0 = transition)
%                   .plv_time   : 1×T PLV(t)  [or .plv]
%                   .nEvents    : # of transitions in that episode
%   cond_labels : cellstr, same length as PLV_list
%                 e.g. {'baseline','baseline','ambtemp',...}
%   titleStr    : string for figure title
%
% OUTPUT
%   STATS : struct with fields:
%             .condNames
%             .t_rel
%             .PLV_mean
%             .PLV_sem
%             .nEvents
%             .scalar   (per-episode PLV)
%             .groups   (condition index per episode)
%             .p_anova

    % ---------- empty check ----------
    if isempty(PLV_list)
        warning('PLV_list is empty, nothing to plot.');
        STATS = struct();
        return;
    end

    % ---------- allow cell array OR struct array ----------
    if iscell(PLV_list)
        % concatenate cell contents into 1×N struct array
        PLV_list = [PLV_list{:}];
    end

    if ~isstruct(PLV_list)
        error('PLV_list must be a struct array or a cell array of structs.');
    end

    % ---------- condition labels ----------
    if ~iscell(cond_labels)
        cond_labels = cellstr(cond_labels(:));
    end
    cond_labels = cond_labels(:);
    if numel(cond_labels) ~= numel(PLV_list)
        error('cond_labels must have one entry per PLV struct.');
    end

    nEp = numel(PLV_list);

    % ---------- unique conditions ----------
    [condNames, ~, condIdx] = unique(cond_labels, 'stable');
    nCond = numel(condNames);

    % ---------- common time axis ----------
    t_ref = PLV_list(1).t_rel(:)';   % row
    T = numel(t_ref);

    % Make sure every PLV struct has plv_time on t_ref
    for k = 1:nEp
        P = PLV_list(k);
        t_k = P.t_rel(:)';

        plv_k = get_plv_time(P);    % helper below, returns 1×T row

        if numel(t_k) ~= T || any(abs(t_k - t_ref) > 1e-6)
            % interpolate onto reference axis
            plv_k = interp1(t_k, plv_k, t_ref, 'linear', 'extrap');
            PLV_list(k).t_rel    = t_ref;
            PLV_list(k).plv_time = plv_k;
        else
            % ensure stored as row
            PLV_list(k).plv_time = plv_k;
        end
    end

    % ---------- condition-wise mean ± SEM ----------
    PLV_mean    = nan(nCond, T);
    PLV_sem     = nan(nCond, T);
    nEventsTot  = zeros(1, nCond);
    nEp_perCond = zeros(1, nCond);

    for c = 1:nCond
        idx_c = find(condIdx == c);
        nEp_perCond(c) = numel(idx_c);
        if isempty(idx_c), continue; end

        M    = nan(numel(idx_c), T);
        nEvC = zeros(numel(idx_c),1);

        for ii = 1:numel(idx_c)
            P = PLV_list(idx_c(ii));
            M(ii,:) = P.plv_time(:)';
            if isfield(P,'nEvents')
                nEvC(ii) = P.nEvents;
            else
                nEvC(ii) = NaN;
            end
        end

        PLV_mean(c,:) = mean(M, 1, 'omitnan');

        n_valid = sum(~isnan(M),1);
        sd_c    = std(M, 0, 1, 'omitnan');
        PLV_sem(c,:) = sd_c ./ max(1, sqrt(n_valid));

        nEventsTot(c) = nansum(nEvC);
    end

    % ---------- scalar PLV per episode (mean 0–20 s) ----------
    scalarPLV = nan(nEp,1);
    win0 = 0;
    win1 = 20;

    for k = 1:nEp
        P   = PLV_list(k);
        t   = P.t_rel(:)';
        plv = P.plv_time(:)';

        idxWin = (t >= win0) & (t <= win1);
        if any(idxWin)
            scalarPLV(k) = mean(plv(idxWin), 'omitnan');
        end
    end

    % ---------- ANOVA across conditions ----------
    good = ~isnan(scalarPLV);
    if sum(good) >= 2 && numel(unique(condIdx(good))) >= 2
        p_anova = anova1(scalarPLV(good), condIdx(good), 'off');
    else
        p_anova = NaN;
    end

    % ---------- pack stats ----------
    STATS = struct();
    STATS.condNames = condNames;
    STATS.t_rel     = t_ref;
    STATS.PLV_mean  = PLV_mean;
    STATS.PLV_sem   = PLV_sem;
    STATS.nEvents   = nEventsTot;
    STATS.scalar    = scalarPLV;
    STATS.groups    = condIdx;
    STATS.p_anova   = p_anova;

    % ---------- plotting ----------
    figure('Color','w','Position',[100 100 900 500]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % colours for conditions (baseline, ambtemp, drugs style)
    cols_cond = [ ...
        0.85 0.85 0.90;   % baseline - pale grey/blue
        0.20 0.45 0.85;   % ambtemp - blue
        0.96 0.80 0.60;   % drugs   - beige
        ];

    if nCond > size(cols_cond,1)
        cols_cond = lines(nCond);
    end

    % pretty names
    niceNames = condNames;
    for c = 1:nCond
        name = lower(condNames{c});
        switch name
            case 'baseline'
                niceNames{c} = 'Baseline';
            case 'ambtemp'
                niceNames{c} = 'Amb. temp';
            case 'drugs'
                niceNames{c} = 'Drugs';
            otherwise
                niceNames{c} = regexprep(condNames{c}, '(^.)','${upper($1)}');
        end
    end

    % ----- top: PLV(t) mean ± SEM -----
    nexttile([1 2]); hold on;
    leg = {};

    for c = 1:nCond
        mu = PLV_mean(c,:);
        se = PLV_sem(c,:);
        if all(isnan(mu)), continue; end

        col = cols_cond(c,:);
        t   = t_ref;

        patch([t fliplr(t)], [mu-se fliplr(mu+se)], ...
              col, 'FaceAlpha',0.25, 'EdgeColor','none');
        plot(t, mu, 'Color', col, 'LineWidth', 2);

        leg{end+1} = sprintf('%s (events = %d)', ...
                             niceNames{c}, nEventsTot(c));
    end

    plot([0 0], ylim, 'k--');
    xlabel('Time from transition (s)');
    ylabel('PLV');
    title(titleStr, 'Interpreter','none');
    if ~isempty(leg)
        legend(leg, 'Location','northwest','Box','off');
    end
    grid on;

    yl = ylim;
    maxPLV = max(PLV_mean(:), [], 'omitnan');
    if ~isnan(maxPLV) && maxPLV > 0
        ylim([0 maxPLV*1.2]);
    else
        ylim([0 0.01]);
    end

    % ----- bottom-left: scalar PLV per episode -----
    nexttile; hold on;
    means = nan(1,nCond);
    x = 1:nCond;

    for c = 1:nCond
        vals = scalarPLV(condIdx == c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        means(c) = mean(vals, 'omitnan');
        bar(c, means(c), 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');

        jitter = (rand(size(vals))*0.4 - 0.2);
        plot(c + jitter, vals, 'k.', 'MarkerSize', 8);
    end

    set(gca,'XTick',x,'XTickLabel',niceNames);
    ylabel('Scalar PLV (0–20 s)');
    if ~isnan(p_anova)
        txt = sprintf('ANOVA p = %.3f', p_anova);
    else
        txt = 'ANOVA p = n/a';
    end
    yl = ylim;
    text(0.5, yl(2)*0.95, txt, 'FontWeight','bold');
    box on;

    % ----- bottom-right: leave blank for future metrics -----
    nexttile; axis off;

    sgtitle(titleStr, 'FontWeight','bold');
end

% ===== helper: robustly get PLV(t) field ================================
function plv = get_plv_time(P)
    if isfield(P,'plv_time')
        plv = P.plv_time;
    elseif isfield(P,'plv')
        plv = P.plv;
    else
        error('PLV struct has no field "plv_time" or "plv".');
    end
    plv = plv(:)';   % row vector
end
