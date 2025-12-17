function plot_psd(f, psd, color)
m = mean(psd, 2);
se = std(psd, 0, 2) / sqrt(size(psd, 2));
fill([f; flipud(f)], [m-se; flipud(m+se)], color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot(f, m, 'Color', color, 'LineWidth', 1.8);
end

function plot_six_group_bar(data, mice, genoWT, genoMut, ...
    COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
    yLabel, show_ids)

hold on;

x_pos = [1, 1.5, 2.5, 3, 4, 4.5];
bar_width = 0.35;

% Plot bars
bar(x_pos(1), mean(data.WT_BASE, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
bar(x_pos(2), mean(data.APP_BASE, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
bar(x_pos(3), mean(data.WT_AMB, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
bar(x_pos(4), mean(data.APP_AMB, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);
bar(x_pos(5), mean(data.WT_DRUGS, 'omitnan'), bar_width, 'FaceColor', COL_WT_DRUGS);
bar(x_pos(6), mean(data.APP_DRUGS, 'omitnan'), bar_width, 'FaceColor', COL_APP_DRUGS);

% Scatter points
jitter = 0.05;
plot_scatter_pts(x_pos(1), data.WT_BASE, show_ids, mice.WT_BASE);
plot_scatter_pts(x_pos(2), data.APP_BASE, show_ids, mice.APP_BASE);
plot_scatter_pts(x_pos(3), data.WT_AMB, show_ids, mice.WT_AMB);
plot_scatter_pts(x_pos(4), data.APP_AMB, show_ids, mice.APP_AMB);
plot_scatter_pts(x_pos(5), data.WT_DRUGS, show_ids, mice.WT_DRUGS);
plot_scatter_pts(x_pos(6), data.APP_DRUGS, show_ids, mice.APP_DRUGS);

xlim([0.5 5]);
set(gca, 'XTick', [1.25, 2.75, 4.25], 'XTickLabel', {'Baseline', 'Ambtemp', 'Drugs'});
ylabel(yLabel);
box off;

% 2-way ANOVA
all_y = [data.WT_BASE; data.APP_BASE; data.WT_AMB; data.APP_AMB; data.WT_DRUGS; data.APP_DRUGS];
cond_f = [ones(size(data.WT_BASE)); ones(size(data.APP_BASE)); ...
          2*ones(size(data.WT_AMB)); 2*ones(size(data.APP_AMB)); ...
          3*ones(size(data.WT_DRUGS)); 3*ones(size(data.APP_DRUGS))];
geno_f = [ones(size(data.WT_BASE)); 2*ones(size(data.APP_BASE)); ...
          ones(size(data.WT_AMB)); 2*ones(size(data.APP_AMB)); ...
          ones(size(data.WT_DRUGS)); 2*ones(size(data.APP_DRUGS))];

ok = ~isnan(all_y);
if sum(ok) > 5
    [p, ~, ~] = anovan(all_y(ok), {cond_f(ok), geno_f(ok)}, ...
        'model', 'interaction', 'display', 'off');
    txt = sprintf('C:%s G:%s I:%s', p_to_star(p(1)), p_to_star(p(2)), p_to_star(p(3)));
    text(0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', 'FontSize', 7, 'FontWeight', 'bold');
end
end