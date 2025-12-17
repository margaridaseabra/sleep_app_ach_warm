function [y_all, cond_idx_all, geno_idx_all] = simple_group_barplot(M, metricName, conds, genos)
% Draw APP vs WT bars for each condition and return data + factor indices.

    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.392 0.584 0.929];

    % Genotype order: APP-like first, WT second if present
    if any(genos == "WT")
        genoWT  = "WT";
        genoMut = genos(genos ~= "WT");
        if isempty(genoMut)
            geno_order = genos;
        else
            geno_order = [genoMut(1); genoWT];
        end
    else
        geno_order = genos;
    end

    nCond = numel(conds);
    nGeno = numel(geno_order);

    y_all         = [];
    cond_idx_all  = [];
    geno_idx_all  = [];

    for gIdx = 1:nGeno
        g = geno_order(gIdx);
        if strcmpi(g, 'WT')
            col = COL_WT;
        else
            col = COL_APP;
        end

        for cIdx = 1:nCond
            c = conds(cIdx);

            mask = (M.cond == c) & (M.geno == g);
            vals = M.(metricName)(mask);

            mu = mean(vals, 'omitnan');
            se = std(vals, 0, 'omitnan') / max(1, sqrt(sum(~isnan(vals))));

            % Cluster position: APP left, WT right within each condition
            x = cIdx + (gIdx - (nGeno+1)/2)*0.25;

            % Bar + errorbar
            if ~isnan(mu)
                bar(x, mu, 0.22, 'FaceColor', col);
                errorbar(x, mu, se, 'k', 'LineStyle','none', 'CapSize',8);
            end

            % Scatter points with jitter
            if ~isempty(vals)
                scatter(x + (rand(size(vals))-0.5)*0.15, vals, ...
                        20, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha',0.6);
            end

            k = numel(vals);
            y_all        = [y_all; vals(:)];
            cond_idx_all = [cond_idx_all; repmat(cIdx, k, 1)];
            geno_idx_all = [geno_idx_all; repmat(gIdx, k, 1)];
        end
    end

    xlim([0.5 nCond+0.5]);
    set(gca, 'XTick', 1:nCond, 'XTickLabel', conds);
end
