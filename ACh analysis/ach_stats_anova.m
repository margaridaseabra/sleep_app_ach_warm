function ach_stats_anova(GROUP, out_dir, analysis_name)
% Two-way ANOVA (condition x genotype) for group metrics, with group barplots.
%
% GROUP        : output from run_ach_batch_auto for the relevant folder(s),
%                optionally subset with subset_group_by_cond.
% out_dir      : directory to save figure (can be []).
% analysis_name: string to put in the figure title / filename.

    if isfield(GROUP, 'metrics_tbl')
        M = GROUP.metrics_tbl;
    else
        M = struct2table(GROUP.metrics);
    end

    M.cond = string(M.cond);
    M.geno = string(M.geno);

    conds_all = unique(M.cond, 'stable');
    genos_all = unique(M.geno, 'stable');

    metrics      = {'NREM_power', 'NREM_peakHz', 'slope_NREM', 'WakeOn_peak_mean'};
    metricLabels = {'NREM ACh power', 'NREM ACh peak freq (Hz)', ...
                    'NREM ACh slope', 'Wake-onset peak dF/F'};

    f = figure('Name',['ANOVA: ' analysis_name],'Color','w');
    nM = numel(metrics);

    for mIdx = 1:nM
        subplot(2,2,mIdx); hold on;

        metricName = metrics{mIdx};

        % Draw grouped bars + scatter and keep factor indices for ANOVA
        [y_all, cond_idx, geno_idx] = simple_group_barplot(M, metricName, conds_all, genos_all);

        % Remove NaNs for stats
        ok  = ~isnan(y_all);
        yv  = y_all(ok);
        cf  = cond_idx(ok);
        gf  = geno_idx(ok);

        if numel(yv) > 0
            [p,~,~] = anovan(yv, {cf, gf}, ...
                             'model','interaction', ...
                             'display','off', ...
                             'varnames', {'cond','geno'});
            p_cond = p(1);
            p_geno = p(2);
            p_int  = p(3);
        else
            p_cond = NaN; p_geno = NaN; p_int = NaN;
        end

        title(metricLabels{mIdx});
        ylabel(metricLabels{mIdx});
        box off;

        % Show p-values in the top-left corner of the subplot
        txt = sprintf('p_{cond} = %.3g\np_{geno} = %.3g\np_{int} = %.3g', ...
                      p_cond, p_geno, p_int);
        text(0.02, 0.98, txt, 'Units','normalized', ...
             'VerticalAlignment','top', 'FontSize',8);

        if mIdx == 1
            % Legend only once (WT vs APP color)
            COL_WT  = [0.6 0.6 0.6];
            COL_APP = [0.392 0.584 0.929];
            hold on;
            hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
            hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
            legend([hWT hAPP], {'WT','APP'}, ...
                   'Location','northoutside','Orientation','horizontal','Box','off');
        end
    end

    sgtitle(['ACh metrics – ' analysis_name ' (two-way ANOVA cond × genotype)']);

    if nargin > 1 && ~isempty(out_dir)
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        saveas(f, fullfile(out_dir, ['ANOVA_' analysis_name '.png']));
    end
end
