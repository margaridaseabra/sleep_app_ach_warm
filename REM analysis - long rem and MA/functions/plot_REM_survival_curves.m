function plot_REM_survival_curves(ALL_REM, out_dir)
% plot_REM_survival_curves
% -------------------------------------------------------------------------
% For each condition, plot REM bout survival curves:
%   P(remaining in REM) vs time (s), with separate curves per genotype.

if isempty(ALL_REM)
    warning('ALL_REM is empty, skipping survival plots.');
    return;
end

conds = unique(ALL_REM.cond(:));
genos = unique(ALL_REM.geno(:));

for iC = 1:numel(conds)
    cond_i = conds(iC);

    figure('Color','w','Units','normalized','Position',[0.1 0.1 0.5 0.5]);
    hold on;

    for iG = 1:numel(genos)
        geno_i = genos(iG);

        mask = (ALL_REM.cond == cond_i) & (ALL_REM.geno == geno_i);
        d = ALL_REM.dur_s(mask);

        if numel(d) < 2
            continue; % too few bouts to plot
        end

        d = sort(d(:));
        t = unique(d);

        % Survival function S(t) = P(D >= t)
        N = numel(d);
        S = zeros(size(t));
        for k = 1:numel(t)
            S(k) = sum(d >= t(k)) / N;
        end

        plot(t, S, 'LineWidth', 2, 'DisplayName', char(geno_i));
    end

    xlabel('REM bout duration (s)');
    ylabel('P(remaining in REM)');
    title(sprintf('REM bout survival – %s', char(cond_i)), 'Interpreter','none');
    legend('Location','southwest');
    grid on;
    set(gca,'FontSize',12,'LineWidth',1.2);

    fname = fullfile(out_dir, sprintf('REM_survival_%s.png', char(cond_i)));
    exportgraphics(gcf, fname, 'Resolution', 300);
    close(gcf);
end
end
