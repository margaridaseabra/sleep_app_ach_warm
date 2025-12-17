function plot_REM_cluster_types(ALL_SUMMARY, out_dir)
% plot_REM_cluster_types
% -------------------------------------------------------------
% Uses the per-mouse summary proportions of cluster types.

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end

ALL_SUMMARY.geno = string(ALL_SUMMARY.geno);
ALL_SUMMARY.cond = string(ALL_SUMMARY.cond);

conds = unique(ALL_SUMMARY.cond,'stable');
genos = unique(ALL_SUMMARY.geno,'stable');

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

cluster_metrics = {'prop_cluster_shortonly','prop_cluster_shortlong','prop_cluster_longonly'};
cluster_labels  = {'Short-only','Short→Long','Long-only'};

for iM = 1:numel(cluster_metrics)
    metric = cluster_metrics{iM};
    
    figure('Color','w'); hold on;
    for iC = 1:numel(conds)
        for iG = 1:numel(genos)
            sel = (ALL_SUMMARY.cond==conds(iC)) & (ALL_SUMMARY.geno==genos(iG));
            y = ALL_SUMMARY.(metric)(sel);
            if isempty(y), continue; end
            
            xbar = (iC-1)*3 + iG;
            m = mean(y,'omitnan');
            s = std(y,'omitnan')/sqrt(numel(y));
            
            if genos(iG)=="WT"
                col = COL_WT;
            else
                col = COL_APP;
            end
            
            bar(xbar, m, 'FaceColor', col, 'EdgeColor','none'); hold on;
            errorbar(xbar, m, s, 'k','LineStyle','none','LineWidth',1);
            
            xj = xbar + (rand(size(y))-0.5)*0.25;
            plot(xj, y, 'ko','MarkerFaceColor','w');
        end
    end
    
    set(gca,'XTick', (1.5:3:(numel(conds)*3)), ...
            'XTickLabel', cellstr(conds));
    ylabel(metric);
    title(cluster_labels{iM});
    box off;
    legend(genos,'Location','bestoutside');
    
    saveas(gcf, fullfile(out_dir, ['REM_clusters_' metric '.png']));
end
