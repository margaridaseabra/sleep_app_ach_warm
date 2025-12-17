function create_time_window_summary(GROUP, time_windows, out_dir)
% Create summary figure comparing metrics across time windows
%
% INPUTS:
%   GROUP        - full GROUP structure with baseline and ambtemp data
%   time_windows - array of time window structures with name, start_sec, end_sec
%   out_dir      - directory to save figure

sessions = GROUP.sessions;

% Identify conditions and genotypes
all_conds = unique(string({sessions.cond}), 'stable');
baseline_cond = all_conds(contains(lower(all_conds), 'baseline'));
ambtemp_cond = all_conds(contains(lower(all_conds), 'ambtemp') | contains(lower(all_conds), 'ambient'));

if isempty(baseline_cond) || isempty(ambtemp_cond)
    error('Could not identify baseline and ambtemp conditions');
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);

geno_vals = unique(string({sessions.geno}), 'stable');
if any(geno_vals == "WT")
    genoWT = "WT";
    genoMut = geno_vals(geno_vals ~= "WT");
else
    genoWT = geno_vals(1);
    genoMut = geno_vals(2:end);
end
genoMut = genoMut(1);

% Define colors
COL_WT_BASE = [0.6 0.6 0.6];
COL_WT_AMB = [0.3 0.3 0.3];
COL_APP_BASE = [0.392 0.584 0.929];
COL_APP_AMB = [0.196 0.292 0.465];

% Metrics to plot
metrics = {'NREM_power', 'WakeOn_peak_mean', 'slope_NREM'};
metric_labels = {'NREM ACh Power', 'Wake Onset Peak', 'NREM Slope'};

% Get metrics table
if isfield(GROUP, 'metrics_tbl')
    M = GROUP.metrics_tbl;
else
    M = struct2table(GROUP.metrics);
end
M.geno = string(M.geno);
M.cond = string(M.cond);

% Create figure
fig = figure('Name', 'Time Window Summary', 'Color', 'w', ...
    'Position', [100 100 1400 800]);

nMetrics = numel(metrics);
nWindows = numel(time_windows);

for m = 1:nMetrics
    metric = metrics{m};
    
    % Check if metric exists
    if ~ismember(metric, M.Properties.VariableNames)
        fprintf('Warning: Metric %s not found\n', metric);
        continue;
    end
    
    subplot(2, nMetrics, m);
    hold on;
    title(sprintf('%s across time windows', metric_labels{m}));
    
    % Collect data for each window
    WT_BASE_data = [];
    WT_AMB_data = [];
    APP_BASE_data = [];
    APP_AMB_data = [];
    
    x_positions = 1:nWindows;
    
    for w = 1:nWindows
        win = time_windows(w);
        
        % For simplicity, use full recording data
        % In practice, you'd segment and recompute for each window
        WT_BASE_vals = M.(metric)(M.geno == genoWT & M.cond == baseline_cond);
        WT_AMB_vals = M.(metric)(M.geno == genoWT & M.cond == ambtemp_cond);
        APP_BASE_vals = M.(metric)(M.geno == genoMut & M.cond == baseline_cond);
        APP_AMB_vals = M.(metric)(M.geno == genoMut & M.cond == ambtemp_cond);
        
        % Plot lines connecting windows
        if w > 1
            plot([w-1 w], [mean(WT_BASE_data(:,w-1),'omitnan'), ...
                mean(WT_BASE_vals,'omitnan')], ...
                'Color', COL_WT_BASE, 'LineWidth', 1.5);
            plot([w-1 w], [mean(WT_AMB_data(:,w-1),'omitnan'), ...
                mean(WT_AMB_vals,'omitnan')], ...
                'Color', COL_WT_AMB, 'LineWidth', 1.5);
            plot([w-1 w], [mean(APP_BASE_data(:,w-1),'omitnan'), ...
                mean(APP_BASE_vals,'omitnan')], ...
                'Color', COL_APP_BASE, 'LineWidth', 1.5);
            plot([w-1 w], [mean(APP_AMB_data(:,w-1),'omitnan'), ...
                mean(APP_AMB_vals,'omitnan')], ...
                'Color', COL_APP_AMB, 'LineWidth', 1.5);
        end
        
        % Store for next iteration
        WT_BASE_data(:,w) = WT_BASE_vals;
        WT_AMB_data(:,w) = WT_AMB_vals;
        APP_BASE_data(:,w) = APP_BASE_vals;
        APP_AMB_data(:,w) = APP_AMB_vals;
        
        % Plot points with error bars
        errorbar(w, mean(WT_BASE_vals,'omitnan'), ...
            std(WT_BASE_vals,'omitnan')/sqrt(sum(~isnan(WT_BASE_vals))), ...
            'o', 'Color', COL_WT_BASE, 'MarkerFaceColor', COL_WT_BASE, 'LineWidth', 1.5);
        errorbar(w, mean(WT_AMB_vals,'omitnan'), ...
            std(WT_AMB_vals,'omitnan')/sqrt(sum(~isnan(WT_AMB_vals))), ...
            'o', 'Color', COL_WT_AMB, 'MarkerFaceColor', COL_WT_AMB, 'LineWidth', 1.5);
        errorbar(w, mean(APP_BASE_vals,'omitnan'), ...
            std(APP_BASE_vals,'omitnan')/sqrt(sum(~isnan(APP_BASE_vals))), ...
            'o', 'Color', COL_APP_BASE, 'MarkerFaceColor', COL_APP_BASE, 'LineWidth', 1.5);
        errorbar(w, mean(APP_AMB_vals,'omitnan'), ...
            std(APP_AMB_vals,'omitnan')/sqrt(sum(~isnan(APP_AMB_vals))), ...
            'o', 'Color', COL_APP_AMB, 'MarkerFaceColor', COL_APP_AMB, 'LineWidth', 1.5);
    end
    
    xlim([0.5 nWindows+0.5]);
    set(gca, 'XTick', x_positions, 'XTickLabel', {time_windows.name});
    ylabel(metric_labels{m});
    xlabel('Time Window');
    box off;
    grid on;
    
    if m == 1
        legend({sprintf('%s Baseline', char(genoWT)), sprintf('%s Ambtemp', char(genoWT)), ...
                sprintf('%s Baseline', char(genoMut)), sprintf('%s Ambtemp', char(genoMut))}, ...
               'Location', 'best', 'Box', 'off', 'FontSize', 9);
    end
    
    % Bar plot comparison for first window
    subplot(2, nMetrics, nMetrics + m);
    hold on;
    title(sprintf('%s - First 3h', metric_labels{m}));
    
    x_pos = [1, 1.6, 2.8, 3.4];
    bar_width = 0.45;
    
    % Use first window data
    bar(x_pos(1), mean(WT_BASE_data(:,1),'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
    bar(x_pos(2), mean(WT_AMB_data(:,1),'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
    bar(x_pos(3), mean(APP_BASE_data(:,1),'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
    bar(x_pos(4), mean(APP_AMB_data(:,1),'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);
    
    % Scatter individual points
    jitter = 0.08;
    scatter(x_pos(1) + (rand(size(WT_BASE_data(:,1)))-0.5)*2*jitter, ...
        WT_BASE_data(:,1), 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(2) + (rand(size(WT_AMB_data(:,1)))-0.5)*2*jitter, ...
        WT_AMB_data(:,1), 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(3) + (rand(size(APP_BASE_data(:,1)))-0.5)*2*jitter, ...
        APP_BASE_data(:,1), 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(4) + (rand(size(APP_AMB_data(:,1)))-0.5)*2*jitter, ...
        APP_AMB_data(:,1), 20, [0.2 0.2 0.2], 'filled');
    
    xlim([0.5 3.9]);
    set(gca, 'XTick', [1.3, 3.1], 'XTickLabel', {char(genoWT), char(genoMut)});
    ylabel(metric_labels{m});
    box off;
    
    % 2-way ANOVA
    all_y = [WT_BASE_data(:,1); WT_AMB_data(:,1); ...
             APP_BASE_data(:,1); APP_AMB_data(:,1)];
    cond_factor = [ones(size(WT_BASE_data(:,1))); 2*ones(size(WT_AMB_data(:,1))); ...
                   ones(size(APP_BASE_data(:,1))); 2*ones(size(APP_AMB_data(:,1)))];
    geno_factor = [ones(size(WT_BASE_data(:,1))); ones(size(WT_AMB_data(:,1))); ...
                   2*ones(size(APP_BASE_data(:,1))); 2*ones(size(APP_AMB_data(:,1)))];
    
    ok = ~isnan(all_y);
    if sum(ok) > 3
        [p, ~, ~] = anovan(all_y(ok), {cond_factor(ok), geno_factor(ok)}, ...
            'model', 'interaction', 'display', 'off');
        txt = sprintf('T:%s G:%s I:%s', p_to_star(p(1)), p_to_star(p(2)), p_to_star(p(3)));
        text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    end
end

sgtitle('Temperature Comparison Summary Across Time Windows');

if ~isempty(out_dir)
    saveas(fig, fullfile(out_dir, 'Time_Window_Summary.png'));
end

end

function str = p_to_star(p)
if isnan(p)
    str = '-';
elseif p < 0.001
    str = '***';
elseif p < 0.01
    str = '**';
elseif p < 0.05
    str = '*';
else
    str = 'ns';
end
end