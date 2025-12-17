function ach_plot_bout_features_allStates(BOUT, STATS)
% ACH_PLOT_BOUT_FEATURES_ALLSTATES
% -------------------------------------------------------------
% For each state (Wake/NREM/REM) make a clean summary figure:
%   - Bout duration (s)
%   - Mean ACh in bout
%   - ACh slope in bout (ΔF/F/s)
%   - EEG–ACh PLV in chosen band
%
% INPUTS
%   BOUT, STATS from ach_compute_bout_features_conditions

if isempty(BOUT)
    warning('BOUT is empty, nothing to plot.');
    return;
end

condNames = {BOUT.name};
nCond     = numel(condNames);
cols_cond = lines(nCond);

stateNamesAll = {BOUT(1).state.name};
targetStates  = {'Wake','NREM','REM'};

for sTarget = 1:numel(targetStates)
    stName = targetStates{sTarget};
    idxState = find(strcmpi(stateNamesAll, stName), 1);
    if isempty(idxState)
        continue;
    end

    % ---- gather metrics ----
    dur      = cell(1,nCond);
    meanACh  = cell(1,nCond);
    slope    = cell(1,nCond);
    plv      = cell(1,nCond);

    mean_dur     = nan(1,nCond);
    mean_meanACh = nan(1,nCond);
    mean_slope   = nan(1,nCond);
    mean_plv     = nan(1,nCond);

    nBouts = zeros(1,nCond);

    for c = 1:nCond
        S = BOUT(c).state(idxState);
        dur{c}     = S.dur_all(:);
        meanACh{c} = S.meanACh_all(:);
        slope{c}   = S.slope_all(:);
        plv{c}     = S.plv_all(:);

        nBouts(c)      = numel(dur{c});
        mean_dur(c)    = mean(dur{c},     'omitnan');
        mean_meanACh(c)= mean(meanACh{c}, 'omitnan');
        mean_slope(c)  = mean(slope{c},   'omitnan');
        mean_plv(c)    = mean(plv{c},     'omitnan');
    end

    % skip if absolutely no bouts
    if all(nBouts == 0)
        continue;
    end

    % ---- p-values for this state ----
    p_dur   = NaN; p_meanA = NaN; p_slope = NaN; p_plv = NaN;
    if isfield(STATS, stName)
        if isfield(STATS.(stName),'dur'),     p_dur   = STATS.(stName).dur.p;     end
        if isfield(STATS.(stName),'meanACh'), p_meanA = STATS.(stName).meanACh.p; end
        if isfield(STATS.(stName),'slope'),   p_slope = STATS.(stName).slope.p;   end
        if isfield(STATS.(stName),'plv'),     p_plv   = STATS.(stName).plv.p;     end
    end

    % =========================================================
    %  Figure for this state
    % =========================================================
    fig = figure('Color','w','Position',[80 80 1100 600]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Bout-level features – %s (n bouts: %s)', ...
        upper(stName), sprintf('%d ', nBouts)), ...
        'FontWeight','bold');

    x = 1:nCond;

    % ---- 1) Duration ----
    nexttile(tl); hold on;
    bar(x, mean_dur, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','none');
    for c = 1:nCond
        vals = dur{c};
        if isempty(vals), continue; end
        jitter = (rand(size(vals))*0.5 - 0.25);
        plot(c + jitter, vals, '.', 'Color', cols_cond(c,:), 'MarkerSize',6);
    end
    set(gca,'XTick', x, 'XTickLabel', condNames);
    ylabel('Bout duration (s)');
    box on; grid on;
    txt = sprintf('ANOVA p = %s', pretty_p(p_dur));
    yl  = ylim;
    text(0.5, yl(2)*0.95, txt, 'FontWeight','bold');
    title(sprintf('%s – duration', stName));

    % ---- 2) Mean ACh ----
    nexttile(tl); hold on;
    bar(x, mean_meanACh, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','none');
    for c = 1:nCond
        vals = meanACh{c};
        if isempty(vals), continue; end
        jitter = (rand(size(vals))*0.5 - 0.25);
        plot(c + jitter, vals, '.', 'Color', cols_cond(c,:), 'MarkerSize',6);
    end
    set(gca,'XTick', x, 'XTickLabel', condNames);
    ylabel('Mean ACh (\DeltaF/F)');
    box on; grid on;
    txt = sprintf('ANOVA p = %s', pretty_p(p_meanA));
    yl  = ylim;
    text(0.5, yl(2)*0.95, txt, 'FontWeight','bold');
    title(sprintf('%s – mean ACh', stName));

    % ---- 3) Slope ----
    nexttile(tl); hold on;
    bar(x, mean_slope, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','none');
    for c = 1:nCond
        vals = slope{c};
        if isempty(vals), continue; end
        jitter = (rand(size(vals))*0.5 - 0.25);
        plot(c + jitter, vals, '.', 'Color', cols_cond(c,:), 'MarkerSize',6);
    end
    set(gca,'XTick', x, 'XTickLabel', condNames);
    ylabel('Slope (\DeltaF/F/s)');
    box on; grid on;
    txt = sprintf('ANOVA p = %s', pretty_p(p_slope));
    yl  = ylim;
    text(0.5, yl(2)*0.95, txt, 'FontWeight','bold');
    title(sprintf('%s – ACh slope', stName));

    % ---- 4) PLV ----
    nexttile(tl); hold on;
    bar(x, mean_plv, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','none');
    for c = 1:nCond
        vals = plv{c};
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        jitter = (rand(size(vals))*0.5 - 0.25);
        plot(c + jitter, vals, '.', 'Color', cols_cond(c,:), 'MarkerSize',6);
    end
    set(gca,'XTick', x, 'XTickLabel', condNames);
    ylabel('PLV (EEG–ACh)');
    box on; grid on;
    txt = sprintf('ANOVA p = %s', pretty_p(p_plv));
    yl  = ylim;
    text(0.5, yl(2)*0.95, txt, 'FontWeight','bold');
    title(sprintf('%s – whole-bout PLV', stName));

end
end

function str = pretty_p(p)
    if isnan(p)
        str = 'n/a';
    elseif p < 0.001
        str = '< 0.001';
    else
        str = sprintf('%.3f', p);
    end
end
