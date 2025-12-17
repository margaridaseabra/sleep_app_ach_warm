function ach_stats_three_conditions(GROUP, out_dir, analysis_name, time_label)
% Two-way ANOVA for three conditions (baseline vs ambtemp vs drug)

if nargin < 4 || isempty(time_label); time_label = ''; end

if isfield(GROUP, 'metrics_tbl'); M = GROUP.metrics_tbl;
else; M = struct2table(GROUP.metrics); end

M.cond = string(M.cond);
M.geno = string(M.geno);

% Identify conditions
all_conds = unique(M.cond, 'stable');
baseline_cond = all_conds(contains(lower(all_conds), 'baseline'));
ambtemp_cond = all_conds(contains(lower(all_conds), 'ambtemp'));
drugs_cond = all_conds(contains(lower(all_conds), 'drugs') | contains(lower(all_conds), 'drug'));

if isempty(baseline_cond) || isempty(ambtemp_cond) || isempty(drugs_cond)
    error('Could not identify all three conditions. Found: %s', strjoin(all_conds, ', '));
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);
drugs_cond = drugs_cond(1);

% Identify genotypes
genos = unique(M.geno, 'stable');
genoWT = genos(genos == "WT"); if isempty(genoWT); genoWT = genos(1); end
genoMut = genos(genos ~= genoWT); genoMut = genoMut(1);

% Filter to these three conditions using ismember
M = M(ismember(M.cond, [baseline_cond, ambtemp_cond, drugs_cond]), :);

% Define colors
COL_WT_BASE = [0.7 0.7 0.7];
COL_APP_BASE = [0.5 0.65 0.85];
COL_WT_AMB = [0.5 0.5 0.5];
COL_APP_AMB = [0.3 0.5 0.75];
COL_WT_DRUGS = [0.3 0.3 0.3];
COL_APP_DRUGS = [0.2 0.35 0.6];

metrics = {'NREM_power', 'NREM_peakHz', 'slope_NREM', 'WakeOn_peak_mean'};
metricLabels = {'NREM ACh power', 'NREM ACh peak freq (Hz)', 'NREM ACh slope', 'Wake-onset peak dF/F'};

f = figure('Name', ['Three Conditions: ' analysis_name], 'Color', 'w', 'Position', [100 100 1400 800]);

for mIdx = 1:numel(metrics)
    subplot(2, 2, mIdx); hold on;
    metricName = metrics{mIdx};
    
    % Extract data for all 6 groups
    y_WT_BASE = M.(metricName)(M.cond == baseline_cond & M.geno == genoWT);
    y_APP_BASE = M.(metricName)(M.cond == baseline_cond & M.geno == genoMut);
    y_WT_AMB = M.(metricName)(M.cond == ambtemp_cond & M.geno == genoWT);
    y_APP_AMB = M.(metricName)(M.cond == ambtemp_cond & M.geno == genoMut);
    y_WT_DRUGS = M.(metricName)(M.cond == drugs_cond & M.geno == genoWT);
    y_APP_DRUGS = M.(metricName)(M.cond == drugs_cond & M.geno == genoMut);
    
    % Bar positions
    x_pos = [1, 1.6, 2.8, 3.4, 4.6, 5.2];
    bar_width = 0.45;
    
    % Plot bars
    bar(x_pos(1), mean(y_WT_BASE, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
    bar(x_pos(2), mean(y_APP_BASE, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
    bar(x_pos(3), mean(y_WT_AMB, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
    bar(x_pos(4), mean(y_APP_AMB, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);
    bar(x_pos(5), mean(y_WT_DRUGS, 'omitnan'), bar_width, 'FaceColor', COL_WT_DRUGS);
    bar(x_pos(6), mean(y_APP_DRUGS, 'omitnan'), bar_width, 'FaceColor', COL_APP_DRUGS);
    
    % Scatter points
    jitter = 0.06;
    scatter(x_pos(1) + (rand(size(y_WT_BASE))-0.5)*2*jitter, y_WT_BASE, 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(2) + (rand(size(y_APP_BASE))-0.5)*2*jitter, y_APP_BASE, 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(3) + (rand(size(y_WT_AMB))-0.5)*2*jitter, y_WT_AMB, 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(4) + (rand(size(y_APP_AMB))-0.5)*2*jitter, y_APP_AMB, 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(5) + (rand(size(y_WT_DRUGS))-0.5)*2*jitter, y_WT_DRUGS, 20, [0.2 0.2 0.2], 'filled');
    scatter(x_pos(6) + (rand(size(y_APP_DRUGS))-0.5)*2*jitter, y_APP_DRUGS, 20, [0.2 0.2 0.2], 'filled');
    
    xlim([0.5 5.7]);
    set(gca, 'XTick', [1.3, 3.1, 4.9], 'XTickLabel', {'Baseline', 'Ambtemp', 'Drugs'});
    ylabel(metricLabels{mIdx});
    title(metricLabels{mIdx});
    box off;
    
    % 2-way ANOVA
    all_y = [y_WT_BASE; y_APP_BASE; y_WT_AMB; y_APP_AMB; y_WT_DRUGS; y_APP_DRUGS];
    cond_f = [ones(size(y_WT_BASE)); ones(size(y_APP_BASE)); ...
              2*ones(size(y_WT_AMB)); 2*ones(size(y_APP_AMB)); ...
              3*ones(size(y_WT_DRUGS)); 3*ones(size(y_APP_DRUGS))];
    geno_f = [ones(size(y_WT_BASE)); 2*ones(size(y_APP_BASE)); ...
              ones(size(y_WT_AMB)); 2*ones(size(y_APP_AMB)); ...
              ones(size(y_WT_DRUGS)); 2*ones(size(y_APP_DRUGS))];
    
    ok = ~isnan(all_y);
    if sum(ok) > 5
        [p, ~, ~] = anovan(all_y(ok), {cond_f(ok), geno_f(ok)}, ...
            'model', 'interaction', 'display', 'off');
        txt = sprintf('C:%s G:%s I:%s', p_to_star(p(1)), p_to_star(p(2)), p_to_star(p(3)));
        text(0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', 'FontSize', 8, 'FontWeight', 'bold');
    end
    
    if mIdx == 1
        h1 = bar(nan, nan, bar_width, 'FaceColor', COL_WT_BASE);
        h2 = bar(nan, nan, bar_width, 'FaceColor', COL_APP_BASE);
        h3 = bar(nan, nan, bar_width, 'FaceColor', COL_WT_AMB);
        h4 = bar(nan, nan, bar_width, 'FaceColor', COL_APP_AMB);
        h5 = bar(nan, nan, bar_width, 'FaceColor', COL_WT_DRUGS);
        h6 = bar(nan, nan, bar_width, 'FaceColor', COL_APP_DRUGS);
        legend([h1 h2 h3 h4 h5 h6], {'WT Base', 'APP Base', 'WT Amb', 'APP Amb', 'WT Drugs', 'APP Drugs'}, ...
            'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8);
    end
end

if ~isempty(time_label)
    sgtitle(['ACh metrics – ' analysis_name ' – ' time_label]);
else
    sgtitle(['ACh metrics – ' analysis_name ' (Three Conditions)']);
end

if ~isempty(out_dir)
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    saveas(f, fullfile(out_dir, ['ANOVA_ThreeConditions_' analysis_name '.png']));
end
end

function str = p_to_star(p)
if p < 0.001; str = '***';
elseif p < 0.01; str = '**';
elseif p < 0.05; str = '*';
else; str = 'ns';
end
end