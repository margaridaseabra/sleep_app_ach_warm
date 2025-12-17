function ach_plot_condition_summary_allStates(COND, STATS)
% ACH_PLOT_CONDITION_SUMMARY_ALLSTATES
% -------------------------------------------------------------
% For each state (Wake, NREM, REM) make a summary figure:
%   - mean ± SEM onset traces per condition
%   - bar + scatter of peak ΔF/F with ANOVA p-value
%   - bar + scatter of slopes with ANOVA p-value
%
% INPUTS
%   COND, STATS : outputs from ach_compare_conditions_onsets

if isempty(COND)
    warning('COND is empty'); 
    return;
end

% Condition names (prettified)
nCond = numel(COND);
condNames_raw = cell(1,nCond);
for c = 1:nCond
    condNames_raw{c} = COND(c).name;
end

condPretty = cell(1,nCond);
for c = 1:nCond
    raw = lower(strtrim(condNames_raw{c}));
    switch raw
        case 'baseline'
            condPretty{c} = 'Baseline';
        case 'ambtemp'
            condPretty{c} = 'Amb. temp';
        case 'drugs'
            condPretty{c} = 'Drugs';
        otherwise
            r = strrep(condNames_raw{c},'_',' ');
            if ~isempty(r), r(1) = upper(r(1)); end
            condPretty{c} = r;
    end
end

% State list from first condition
stateNames = {COND(1).state.name};
targetStates = {'Wake','NREM','REM'};

cols_cond = lines(nCond);

for sTarget = 1:numel(targetStates)
    stName = targetStates{sTarget};

    % index of this state inside COND(c).state
    idxState = find(strcmpi(stateNames, stName), 1);
    if isempty(idxState)
        warning('State "%s" not found, skipping.', stName);
        continue;
    end

    % common time axis
    t_rel = COND(1).state(idxState).t_rel;
    if isempty(t_rel)
        continue;
    end

    % Gather feature pools across conditions
    all_dF    = [];
    all_slope = [];
    g_dF      = [];
    g_slope   = [];

    mean_dF    = nan(1,nCond);
    mean_slope = nan(1,nCond);

    nEvents_cond = zeros(1,nCond);

    for c = 1:nCond
        S  = COND(c).state(idxState);
        dF = S.deltaF_all(:);
        sl = S.slope_all(:);

        nEvents_cond(c) = numel(dF);

        mean_dF(c)    = mean(dF,'omitnan');
        mean_slope(c) = mean(sl,'omitnan');

        all_dF    = [all_dF;    dF];
        g_dF      = [g_dF;      c*ones(numel(dF),1)];

        all_slope = [all_slope; sl];
        g_slope   = [g_slope;   c*ones(numel(sl),1)];
    end

    if all(isnan(mean_dF)) && all(isnan(mean_slope))
        % no events for this state in any condition
        continue;
    end

    % ---------- Figure for this state ----------
    figure('Color','w','Position',[80 80 1100 600]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('%s ONSETS – ACh features', upper(stName)));

    %% (1) Top-left: mean ± SEM onset traces per condition
    ax1 = nexttile(tl,1); 
    hold(ax1,'on');
    leg = {};

    for c = 1:nCond
        S  = COND(c).state(idxState);
        mu = S.mean;
        se = S.sem;
        nE = S.nEvents;

        if isempty(mu) || all(isnan(mu)) || nE == 0
            continue;
        end

        col = cols_cond(c,:);

        % SEM shading
        patch(ax1, [t_rel fliplr(t_rel)], [mu-se fliplr(mu+se)], ...
              col, 'FaceAlpha',0.15, 'EdgeColor','none');
        plot(ax1, t_rel, mu, 'Color', col, 'LineWidth', 2);

        leg{end+1} = sprintf('%s (events = %d)', condPretty{c}, nE);
    end

    plot(ax1, [0 0], ylim(ax1), 'k--');
    xlabel(ax1, 'Time from transition (s)');
    ylabel(ax1, 'ACh (\DeltaF/F)');
    if ~isempty(leg)
        legend(ax1, leg, 'Location','best','Box','off');
    end
    grid(ax1,'on');
    title(ax1, sprintf('%s onset – mean ACh trace', stName));

    %% (2) Top-right: peak ΔF/F, bar + scatter
    ax2 = nexttile(tl,2); 
    hold(ax2,'on');

    x = 1:nCond;
    bar(ax2, x, mean_dF, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');

    for c = 1:nCond
        S  = COND(c).state(idxState);
        dF = S.deltaF_all(:);
        if isempty(dF), continue; end
        jitter = (rand(size(dF))*0.4 - 0.2);
        plot(ax2, c + jitter, dF, 'k.', 'MarkerSize',8);
    end

    xlim(ax2,[0.5 nCond+0.5]);
    set(ax2,'XTick',x,'XTickLabel',condPretty);
    ylabel(ax2,'Peak ACh (\DeltaF/F)');
    box(ax2,'on');
    title(ax2, sprintf('%s onset – peak response', stName));

    % pick p-value from STATS.(state).deltaF.p if available
    p_txt = 'p = n/a';
    if isfield(STATS, stName) && isfield(STATS.(stName),'deltaF')
        p = STATS.(stName).deltaF.p;
        if ~isnan(p)
            p_txt = sprintf('ANOVA p = %.3g', p);
        end
    end
    yl = ylim(ax2);
    text(ax2, 0.6, yl(2)*0.95, p_txt, 'FontWeight','bold');

    %% (3) Bottom-left: slope, bar + scatter
    ax3 = nexttile(tl,3); 
    hold(ax3,'on');

    bar(ax3, x, mean_slope, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');

    for c = 1:nCond
        S  = COND(c).state(idxState);
        sl = S.slope_all(:);
        if isempty(sl), continue; end
        jitter = (rand(size(sl))*0.4 - 0.2);
        plot(ax3, c + jitter, sl, 'k.', 'MarkerSize',8);
    end

    xlim(ax3,[0.5 nCond+0.5]);
    set(ax3,'XTick',x,'XTickLabel',condPretty);
    ylabel(ax3,'Slope (\DeltaF/F/s)');
    box(ax3,'on');
    title(ax3, sprintf('%s onset – slopes', stName));

    p_txt2 = 'p = n/a';
    if isfield(STATS, stName) && isfield(STATS.(stName),'slope')
        p2 = STATS.(stName).slope.p;
        if ~isnan(p2)
            p_txt2 = sprintf('ANOVA p = %.3g', p2);
        end
    end
    yl = ylim(ax3);
    text(ax3, 0.6, yl(2)*0.95, p_txt2, 'FontWeight','bold');

    %% (4) Bottom-right: leave empty or for future features
    ax4 = nexttile(tl,4);
    axis(ax4,'off');
end
end
