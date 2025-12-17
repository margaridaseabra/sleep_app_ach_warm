function plot_APP_WT_metric(M_APP, M_WT, conditions, colors, ylab, ttl)
% M_APP, M_WT: [nAnimals x nCond] matrices (can be 1 x nCond)
% conditions : cellstr, ideally length = nCond
% colors     : struct with .APP (nCond x 3), .WT (nCond x 3)
% ylab       : string for y-label
% ttl        : string for figure title

    % ---- basic checks ----
    if isempty(M_APP) || isempty(M_WT)
        warning('plot_APP_WT_metric: empty APP or WT matrix, skipping plot (%s).', ttl);
        return;
    end

    % Make sure both are 2D and not weirdly squeezed
    M_APP = squeeze(M_APP);
    M_WT  = squeeze(M_WT);

    % If they're vectors, make them 1 x nCond row vectors
    if isvector(M_APP), M_APP = reshape(M_APP, 1, []); end
    if isvector(M_WT),  M_WT  = reshape(M_WT,  1, []); end

    % Determine number of conditions from the data, not from conditions{}
    nCond_app = size(M_APP, 2);
    nCond_wt  = size(M_WT,  2);

    if nCond_app ~= nCond_wt
        error('plot_APP_WT_metric: APP and WT have different number of columns (%d vs %d).', ...
              nCond_app, nCond_wt);
    end

    nCond = nCond_app;
    x     = 1:nCond;

    % Make sure conditions and colours have at least nCond entries
    if numel(conditions) < nCond
        warning('plot_APP_WT_metric: conditions has only %d labels, but data has %d columns. Truncating.', ...
                numel(conditions), nCond);
    end
    cond_plot = conditions(1:nCond);

    if size(colors.APP,1) < nCond || size(colors.WT,1) < nCond
        error('plot_APP_WT_metric: colours do not have enough rows for %d conditions.', nCond);
    end

    % ---- means and SEM ----
    mean_APP = mean(M_APP, 1, 'omitnan');
    mean_WT  = mean(M_WT,  1, 'omitnan');

    sem_APP = std(M_APP, 0, 1, 'omitnan') ./ max(1, sqrt(sum(~isnan(M_APP),1)));
    sem_WT  = std(M_WT,  0, 1, 'omitnan') ./ max(1, sqrt(sum(~isnan(M_WT),1)));

    % ---- create figure ----
    figure('Color','w','Position',[200 200 600 450]); hold on;

    % Grouped bars: nCond rows × 2 columns (APP, WT)
    Y = [mean_APP; mean_WT]';     % [nCond x 2]
    b = bar(x, Y, 'grouped');     % b(1)=APP, b(2)=WT

    % Assign per-condition colours
    b(1).FaceColor = 'flat';
    b(2).FaceColor = 'flat';
    b(1).CData     = colors.APP(1:nCond, :);
    b(2).CData     = colors.WT(1:nCond, :);
    b(1).EdgeColor = 'k';
    b(2).EdgeColor = 'k';

    % ---- errorbars ----
    for c = 1:nCond
        x_APP = b(1).XEndPoints(c);
        x_WT  = b(2).XEndPoints(c);

        errorbar(x_APP, mean_APP(c), sem_APP(c), 'k', ...
                 'LineStyle','none', 'LineWidth',1);
        errorbar(x_WT,  mean_WT(c),  sem_WT(c),  'k', ...
                 'LineStyle','none', 'LineWidth',1);
    end

    % ---- scatter individual animals on top (dots) ----
    jitter = 0.05;
    for c = 1:nCond
        x_APP = b(1).XEndPoints(c);
        x_WT  = b(2).XEndPoints(c);

        for i = 1:size(M_APP,1)
            plot(x_APP + (rand-0.5)*2*jitter, M_APP(i,c), 'o', ...
                 'MarkerFaceColor', colors.APP(c,:), ...
                 'MarkerEdgeColor', 'k');
        end
        for i = 1:size(M_WT,1)
            plot(x_WT + (rand-0.5)*2*jitter, M_WT(i,c), 'o', ...
                 'MarkerFaceColor', colors.WT(c,:), ...
                 'MarkerEdgeColor', 'k');
        end
    end

    set(gca,'XTick',x,'XTickLabel',cond_plot,'FontSize',12);
    ylabel(ylab,'FontSize',13);
    legend({'APP','WT'},'Location','best');
    title(ttl,'FontSize',14);
    grid on; box on;
end
