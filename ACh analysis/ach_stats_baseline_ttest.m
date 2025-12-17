function ach_stats_baseline_ttest(GROUP, out_dir)
% Baseline: compare WT vs APP with unpaired t-tests for a few key metrics.
%
% GROUP : output of run_ach_batch_auto, already subset to baseline.
% out_dir: where to save the figure (can be [] to skip saving).

    if isfield(GROUP, 'metrics_tbl')
        M = GROUP.metrics_tbl;
    else
        M = struct2table(GROUP.metrics);
    end

    M.geno = string(M.geno);

    % Which group is WT, which is APP-like
    if any(M.geno == "WT")
        genoWT  = "WT";
        genoMut = unique(M.geno(M.geno ~= "WT"));
    else
        % Fallback if genotype labels are different
        g_all   = unique(M.geno);
        genoWT  = g_all(1);
        genoMut = g_all(2:end);
    end
    genoMut = genoMut(1);  % assume single mutant group

    % Pick a few important metrics
    metrics      = {'NREM_power', 'NREM_peakHz', 'slope_NREM', 'WakeOn_peak_mean'};
    metricLabels = {'NREM ACh power', 'NREM ACh peak freq (Hz)', ...
                    'NREM ACh slope', 'Wake-onset peak dF/F'};

    f = figure('Name','Baseline: WT vs APP (t-tests)','Color','w');
    nM = numel(metrics);

    for i = 1:nM
        subplot(2,2,i); hold on;

        yWT  = M.(metrics{i})(M.geno == genoWT);
        yAPP = M.(metrics{i})(M.geno == genoMut);

        yWT  = yWT(~isnan(yWT));
        yAPP = yAPP(~isnan(yAPP));

        % Bars
        COL_WT  = [0.6 0.6 0.6];
        COL_APP = [0.392 0.584 0.929];

        bar(1, mean(yWT),  0.6, 'FaceColor', COL_WT);
        bar(2, mean(yAPP), 0.6, 'FaceColor', COL_APP);

        % Scatter data points with jitter
        if ~isempty(yWT)
            scatter(1 + (rand(size(yWT))-0.5)*0.1,  yWT,  25, [0.2 0.2 0.2], 'filled');
        end
        if ~isempty(yAPP)
            scatter(2 + (rand(size(yAPP))-0.5)*0.1, yAPP, 25, [0.2 0.2 0.2], 'filled');
        end

        xlim([0.5 2.5]);
        set(gca,'XTick',[1 2], 'XTickLabel',{char(genoWT), char(genoMut)});
        ylabel(metricLabels{i});
        title(metricLabels{i});
        box off;

        % t-test
        if numel(yWT) > 1 && numel(yAPP) > 1
            [~,p] = ttest2(yWT, yAPP);

            yMax   = max([yWT; yAPP]);
            yMin   = min([yWT; yAPP]);
            yrange = max(yMax - yMin, eps);

            yStar = yMax + 0.15*yrange;
            plot([1 2], [yStar yStar], 'k-', 'LineWidth', 1);
            text(1.5, yStar + 0.03*yrange, p_to_star(p), ...
                 'HorizontalAlignment','center', 'FontSize',10);
            text(1.5, yStar + 0.15*yrange, sprintf('p = %.3g', p), ...
                 'HorizontalAlignment','center', 'FontSize',8);
        end
    end

    sgtitle('Baseline ACh metrics: WT vs APP');

    if nargin > 1 && ~isempty(out_dir)
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        saveas(f, fullfile(out_dir, 'Baseline_WT_vs_APP_ttests.png'));
    end
end
