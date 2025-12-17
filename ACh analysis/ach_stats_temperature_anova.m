function ach_stats_temperature_anova(GROUP, out_dir, analysis_name, time_label)
% Two-way ANOVA comparing baseline vs ambient temperature (ambtemp) conditions
% Compares genotype (WT vs APP) and temperature condition
%
% GROUP        : output from run_ach_batch_auto
% out_dir      : directory to save figure (can be [])
% analysis_name: string for figure title/filename
% time_label   : string describing time period (optional)

if nargin < 4 || isempty(time_label)
    time_label = '';
end

if isfield(GROUP, 'metrics_tbl')
    M = GROUP.metrics_tbl;
else
    M = struct2table(GROUP.metrics);
end

M.cond = string(M.cond);
M.geno = string(M.geno);

% Identify baseline and ambtemp conditions
all_conds = unique(M.cond, 'stable');
baseline_cond = all_conds(contains(lower(all_conds), 'baseline'));
ambtemp_cond = all_conds(contains(lower(all_conds), 'ambtemp') | contains(lower(all_conds), 'ambient'));

if isempty(baseline_cond) || isempty(ambtemp_cond)
    error('Could not identify baseline and ambtemp conditions. Conditions found: %s', ...
        strjoin(all_conds, ', '));
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);

% Filter to only these two conditions
M = M(M.cond == baseline_cond | M.cond == ambtemp_cond, :);
conds_temp = [baseline_cond, ambtemp_cond];

% Identify genotypes
genos_all = unique(M.geno, 'stable');
if any(genos_all == "WT")
    genoWT = "WT";
    genoMut = genos_all(genos_all ~= "WT");
else
    genoWT = genos_all(1);
    genoMut = genos_all(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes');
end
genoMut = genoMut(1);
genos_temp = [genoWT, genoMut];

fprintf('\n=== Temperature Comparison ANOVA ===\n');
fprintf('Baseline condition: %s\n', baseline_cond);
fprintf('Ambient temp condition: %s\n', ambtemp_cond);
fprintf('WT genotype: %s\n', genoWT);
fprintf('Mutant genotype: %s\n', genoMut);

% Check sample sizes
for c = 1:numel(conds_temp)
    for g = 1:numel(genos_temp)
        n = sum(M.cond == conds_temp(c) & M.geno == genos_temp(g));
        fprintf('  %s + %s: n=%d\n', conds_temp(c), genos_temp(g), n);
    end
end
fprintf('====================================\n\n');

% Metrics to analyze
metrics      = {'NREM_power', 'NREM_peakHz', 'slope_NREM', 'WakeOn_peak_mean'};
metricLabels = {'NREM ACh power', 'NREM ACh peak freq (Hz)', ...
                'NREM ACh slope', 'Wake-onset peak dF/F'};

f = figure('Name', ['Temperature ANOVA: ' analysis_name], 'Color', 'w', ...
    'Position', [100 100 1200 800]);
nM = numel(metrics);

for mIdx = 1:nM
    subplot(2, 2, mIdx); hold on;
    metricName = metrics{mIdx};
    
    % Draw grouped bars + scatter and keep factor indices for ANOVA
    [y_all, cond_idx, geno_idx] = temperature_group_barplot(M, metricName, ...
        baseline_cond, ambtemp_cond, genoWT, genoMut);
    
    % Remove NaNs for stats
    ok  = ~isnan(y_all);
    yv  = y_all(ok);
    cf  = cond_idx(ok);
    gf  = geno_idx(ok);
    
    if numel(yv) > 3 && numel(unique(cf)) > 1 && numel(unique(gf)) > 1
        [p, tbl, stats] = anovan(yv, {cf, gf}, ...
            'model', 'interaction', ...
            'display', 'off', ...
            'varnames', {'Temperature', 'Genotype'});
        p_temp = p(1);  % Temperature effect
        p_geno = p(2);  % Genotype effect
        p_int  = p(3);  % Interaction
        
        fprintf('%s:\n', metricLabels{mIdx});
        fprintf('  Temperature: p=%.4f %s\n', p_temp, p_to_star(p_temp));
        fprintf('  Genotype: p=%.4f %s\n', p_geno, p_to_star(p_geno));
        fprintf('  Interaction: p=%.4f %s\n\n', p_int, p_to_star(p_int));
    else
        p_temp = NaN; p_geno = NaN; p_int = NaN;
        fprintf('%s: Insufficient data for ANOVA\n\n', metricLabels{mIdx});
    end
    
    title(metricLabels{mIdx});
    ylabel(metricLabels{mIdx});
    box off;
    
    % Show p-values
    txt = sprintf('p_{temp} = %s\np_{geno} = %s\np_{int} = %s', ...
        format_pval(p_temp), format_pval(p_geno), format_pval(p_int));
    text(0.02, 0.98, txt, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 8);
    
    if mIdx == 1
        % Legend showing all four groups
        COL_WT_BASE = [0.6 0.6 0.6];
        COL_WT_AMB = [0.3 0.3 0.3];
        COL_APP_BASE = [0.392 0.584 0.929];
        COL_APP_AMB = [0.196 0.292 0.465];
        
        h1 = bar(nan, nan, 0.22, 'FaceColor', COL_WT_BASE);
        h2 = bar(nan, nan, 0.22, 'FaceColor', COL_WT_AMB);
        h3 = bar(nan, nan, 0.22, 'FaceColor', COL_APP_BASE);
        h4 = bar(nan, nan, 0.22, 'FaceColor', COL_APP_AMB);
        
        legend([h1 h2 h3 h4], ...
            {[char(genoWT) ' Baseline'], [char(genoWT) ' Ambtemp'], ...
             [char(genoMut) ' Baseline'], [char(genoMut) ' Ambtemp']}, ...
            'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);
    end
end

if ~isempty(time_label)
    sgtitle(['ACh metrics – ' analysis_name ' (two-way ANOVA) – ' time_label]);
else
    sgtitle(['ACh metrics – ' analysis_name ' (two-way ANOVA cond × genotype)']);
end

% Add time period label at bottom
if ~isempty(time_label)
    annotation('textbox', [0.02, 0.02, 0.3, 0.02], 'String', sprintf('Time Period: %s', time_label), ...
        'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

if nargin > 1 && ~isempty(out_dir)
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    saveas(f, fullfile(out_dir, ['ANOVA_Temperature_' analysis_name '.png']));
end
end

% ======================================================================
% Helper function for temperature comparison barplot
% ======================================================================
function [y_all, cond_idx, geno_idx] = temperature_group_barplot(M, metricName, ...
    baseline_cond, ambtemp_cond, genoWT, genoMut)
% Create grouped barplot comparing baseline vs ambtemp for WT and APP

% Define colors
COL_WT_BASE = [0.6 0.6 0.6];      % grey - WT baseline
COL_WT_AMB = [0.3 0.3 0.3];       % darker grey - WT ambient
COL_APP_BASE = [0.392 0.584 0.929]; % cornflower blue - APP baseline
COL_APP_AMB = [0.196 0.292 0.465]; % darker blue - APP ambient

% Check if metric exists
if ~ismember(metricName, M.Properties.VariableNames)
    warning('Metric %s not found in table', metricName);
    y_all = [];
    cond_idx = [];
    geno_idx = [];
    return;
end

% Extract data for each group
WT_base_vals = M.(metricName)(M.cond == baseline_cond & M.geno == genoWT);
WT_amb_vals = M.(metricName)(M.cond == ambtemp_cond & M.geno == genoWT);
APP_base_vals = M.(metricName)(M.cond == baseline_cond & M.geno == genoMut);
APP_amb_vals = M.(metricName)(M.cond == ambtemp_cond & M.geno == genoMut);

% Positions for bars (grouped by genotype)
x_positions = [1, 1.6, 2.8, 3.4]; % WT_base, WT_amb, APP_base, APP_amb
bar_width = 0.45;

% Plot bars
bar(x_positions(1), mean(WT_base_vals, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
bar(x_positions(2), mean(WT_amb_vals, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
bar(x_positions(3), mean(APP_base_vals, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
bar(x_positions(4), mean(APP_amb_vals, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);

% Add error bars (SEM)
jitter = 0.08;

% WT baseline
if ~isempty(WT_base_vals) && any(~isnan(WT_base_vals))
    sem = std(WT_base_vals, 'omitnan') / sqrt(sum(~isnan(WT_base_vals)));
    errorbar(x_positions(1), mean(WT_base_vals, 'omitnan'), sem, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1);
    xpos = x_positions(1) + (rand(size(WT_base_vals)) - 0.5) * 2 * jitter;
    scatter(xpos, WT_base_vals, 20, COL_WT_BASE, 'filled', 'MarkerFaceAlpha', 0.6);
end

% WT ambient
if ~isempty(WT_amb_vals) && any(~isnan(WT_amb_vals))
    sem = std(WT_amb_vals, 'omitnan') / sqrt(sum(~isnan(WT_amb_vals)));
    errorbar(x_positions(2), mean(WT_amb_vals, 'omitnan'), sem, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1);
    xpos = x_positions(2) + (rand(size(WT_amb_vals)) - 0.5) * 2 * jitter;
    scatter(xpos, WT_amb_vals, 20, COL_WT_AMB, 'filled', 'MarkerFaceAlpha', 0.6);
end

% APP baseline
if ~isempty(APP_base_vals) && any(~isnan(APP_base_vals))
    sem = std(APP_base_vals, 'omitnan') / sqrt(sum(~isnan(APP_base_vals)));
    errorbar(x_positions(3), mean(APP_base_vals, 'omitnan'), sem, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1);
    xpos = x_positions(3) + (rand(size(APP_base_vals)) - 0.5) * 2 * jitter;
    scatter(xpos, APP_base_vals, 20, COL_APP_BASE, 'filled', 'MarkerFaceAlpha', 0.6);
end

% APP ambient
if ~isempty(APP_amb_vals) && any(~isnan(APP_amb_vals))
    sem = std(APP_amb_vals, 'omitnan') / sqrt(sum(~isnan(APP_amb_vals)));
    errorbar(x_positions(4), mean(APP_amb_vals, 'omitnan'), sem, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1);
    xpos = x_positions(4) + (rand(size(APP_amb_vals)) - 0.5) * 2 * jitter;
    scatter(xpos, APP_amb_vals, 20, COL_APP_AMB, 'filled', 'MarkerFaceAlpha', 0.6);
end

% Set x-axis
xlim([0.5 3.9]);
set(gca, 'XTick', [1.3, 3.1], 'XTickLabel', {char(genoWT), char(genoMut)});

% Prepare output for ANOVA
% Condition: 1=baseline, 2=ambtemp
% Genotype: 1=WT, 2=APP
y_all = [WT_base_vals; WT_amb_vals; APP_base_vals; APP_amb_vals];
cond_idx = [ones(size(WT_base_vals)); 2*ones(size(WT_amb_vals)); ...
            ones(size(APP_base_vals)); 2*ones(size(APP_amb_vals))];
geno_idx = [ones(size(WT_base_vals)); ones(size(WT_amb_vals)); ...
            2*ones(size(APP_base_vals)); 2*ones(size(APP_amb_vals))];
end

% ======================================================================
% Helper functions
% ======================================================================
function str = format_pval(p)
% Format p-value for display
if isnan(p)
    str = 'N/A';
elseif p < 0.001
    str = '<0.001';
elseif p < 0.01
    str = sprintf('%.3f', p);
else
    str = sprintf('%.2f', p);
end
end

function star = p_to_star(p)
% Convert p-value to significance stars
if isnan(p)
    star = '';
elseif p < 0.001
    star = '***';
elseif p < 0.01
    star = '**';
elseif p < 0.05
    star = '*';
else
    star = '';
end
end