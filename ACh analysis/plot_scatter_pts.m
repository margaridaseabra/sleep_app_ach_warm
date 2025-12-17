function plot_scatter_pts(xpos, y_vals, show_ids, mouse_ids)
if isempty(y_vals); return; end
jitter = 0.05;
xpos_j = xpos + (rand(size(y_vals)) - 0.5) * 2 * jitter;
scatter(xpos_j, y_vals, 20, [0.2 0.2 0.2], 'filled');
if show_ids && ~isempty(mouse_ids)
    for i = 1:numel(y_vals)
        if i <= numel(mouse_ids)
            text(xpos_j(i) + 0.05, y_vals(i), mouse_ids{i}, ...
                'FontSize', 5, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', 'Interpreter', 'none');
        end
    end
end
end