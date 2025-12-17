function plot_REM_bout_histograms(ALL_REM, out_dir)
% plot_REM_bout_histograms
% -------------------------------------------------------------------------
% For each condition, plot histogram of REM bout durations with
% overlaid distributions for each genotype.

if isempty(ALL_REM)
    warning('ALL_REM is empty, skipping REM histograms.');
    return;
end

conds = unique(ALL_REM.cond(:));
genos = unique(ALL_REM.geno(:));

for iC = 1:numel(conds)
    cond_i = conds(iC);

    figure('Color','w','Units','normalized','Position',[0.1 0.1 0.5 0.5]);
    hold on;

    % Decide on bin edges based on all bouts in this condition
    dur_all = ALL_REM.dur_s(ALL_REM.cond == cond_i);
    if numel(dur_all) < 2
        continue;
    end
    maxDur = prctile(dur_all, 99); % avoid crazy outliers
    edges = 0:10:ceil(maxDur/10)*10;  % 10-second bins

    for iG = 1:numel(genos)
        geno_i = genos(iG);
        mask = (ALL_REM.cond == cond_i) & (ALL_REM.geno == geno_i);
        d = ALL_REM.dur_s(mask);

        if isempty(d)
            continue;
        end

        h = histogram(d, edges, ...
            'Normalization','probability', ...
            'DisplayName', char(geno_i), ...
            'FaceAlpha', 0.4);
        % Let MATLAB pick colors; we just adjust style.
    end

    xlabel('REM bout duration (s)');
    ylabel('Probability');
    title(sprintf('REM bout length distribution – %s', char(cond_i)), ...
          'Interpreter','none');
    legend('Location','northeast');
    grid on;
    set(gca,'FontSize',12,'LineWidth',1.2);

    fname = fullfile(out_dir, sprintf('REM_hist_%s.png', char(cond_i)));
    exportgraphics(gcf, fname, 'Resolution', 300);
    close(gcf);
end
end
