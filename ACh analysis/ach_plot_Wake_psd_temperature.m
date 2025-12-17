function ach_plot_Wake_psd_temperature(GROUP, out_dir, show_mouse_ids)
% Plot Wake ACh PSD comparing baseline vs ambtemp (WT vs APP) with 2-way ANOVA
%
% INPUTS:
%   GROUP          - structure with sessions and metrics
%   out_dir        - directory to save figures (optional)
%   show_mouse_ids - logical, if true show mouse IDs on scatter plots (default: false)

if nargin < 2
    out_dir = [];
end
if nargin < 3 || isempty(show_mouse_ids)
    show_mouse_ids = false;
end

sessions = GROUP.sessions;

% Identify baseline and ambtemp conditions
all_conds = unique(string({sessions.cond}), 'stable');
baseline_cond = all_conds(contains(lower(all_conds), 'baseline'));
ambtemp_cond = all_conds(contains(lower(all_conds), 'ambtemp'));

if isempty(baseline_cond) || isempty(ambtemp_cond)
    error('Could not identify baseline and ambtemp conditions. Conditions found: %s', ...
        strjoin(all_conds, ', '));
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);

% Figure out genotypes
geno_vals = unique(string({sessions.geno}), 'stable');
if any(geno_vals == "WT")
    genoWT  = "WT";
    genoMut = geno_vals(geno_vals ~= "WT");
else
    genoWT  = geno_vals(1);
    genoMut = geno_vals(2:end);
end
genoMut = genoMut(1);

% Define colors
COL_WT_BASE = [0.6 0.6 0.6];      % grey - WT baseline
COL_WT_AMB = [0.3 0.3 0.3];       % darker grey - WT ambient
COL_APP_BASE = [0.392 0.584 0.929]; % cornflower blue - APP baseline
COL_APP_AMB = [0.196 0.292 0.465]; % darker blue - APP ambient

fprintf('\n=== Temperature PSD Comparison ===\n');
fprintf('Baseline condition: %s\n', baseline_cond);
fprintf('Ambient temp condition: %s\n', ambtemp_cond);
fprintf('WT genotype: %s\n', genoWT);
fprintf('Mutant genotype: %s\n', genoMut);
fprintf('==================================\n\n');

% ---------- collect PSD curves ----------
f_ref = [];
WT_BASE_psd  = [];
WT_AMB_psd   = [];
APP_BASE_psd = [];
APP_AMB_psd  = [];

for k = 1:numel(sessions)
    OUT = GROUP.out{k};
    if ~isfield(OUT, 'psd') || ~isfield(OUT.psd, 'Wake'); continue; end
    
    P = [];
    if isfield(OUT.psd.Wake, 'Pxx')
        P = OUT.psd.Wake.Pxx(:);
    elseif isfield(OUT.psd.Wake, 'psd')
        P = OUT.psd.Wake.psd(:);
    elseif isfield(OUT.psd.Wake, 'power')
        P = OUT.psd.Wake.power(:);
    end
    if isempty(P) || ~isfield(OUT.psd.Wake, 'f'); continue; end
    
    f_this = OUT.psd.Wake.f(:);
    
    % reference frequency grid
    if isempty(f_ref)
        f_ref = f_this;
    elseif numel(f_this) ~= numel(f_ref) || any(abs(f_this - f_ref) > 1e-6)
        P = interp1(f_this, P, f_ref, 'linear', 'extrap');
    end
    
    g = string(sessions(k).geno);
    c = string(sessions(k).cond);
    
    if g == genoWT && c == baseline_cond
        WT_BASE_psd = [WT_BASE_psd P]; %#ok<AGROW>
    elseif g == genoWT && c == ambtemp_cond
        WT_AMB_psd = [WT_AMB_psd P]; %#ok<AGROW>
    elseif g == genoMut && c == baseline_cond
        APP_BASE_psd = [APP_BASE_psd P]; %#ok<AGROW>
    elseif g == genoMut && c == ambtemp_cond
        APP_AMB_psd = [APP_AMB_psd P]; %#ok<AGROW>
    end
end

fprintf('PSD curves collected:\n');
fprintf('  WT Baseline: %d\n', size(WT_BASE_psd, 2));
fprintf('  WT Ambient: %d\n', size(WT_AMB_psd, 2));
fprintf('  APP Baseline: %d\n', size(APP_BASE_psd, 2));
fprintf('  APP Ambient: %d\n', size(APP_AMB_psd, 2));

% ---------- metrics from GROUP.metrics ----------
if isfield(GROUP, 'metrics_tbl')
    M = GROUP.metrics_tbl;
else
    M = struct2table(GROUP.metrics);
end
M.geno = string(M.geno);
M.cond = string(M.cond);

% Extract metrics for each group
y_WT_BASE_power = M.Wake_power(M.geno == genoWT & M.cond == baseline_cond);
y_WT_AMB_power = M.Wake_power(M.geno == genoWT & M.cond == ambtemp_cond);
y_APP_BASE_power = M.Wake_power(M.geno == genoMut & M.cond == baseline_cond);
y_APP_AMB_power = M.Wake_power(M.geno == genoMut & M.cond == ambtemp_cond);

y_WT_BASE_peakHz = M.Wake_peakHz(M.geno == genoWT & M.cond == baseline_cond);
y_WT_AMB_peakHz = M.Wake_peakHz(M.geno == genoWT & M.cond == ambtemp_cond);
y_APP_BASE_peakHz = M.Wake_peakHz(M.geno == genoMut & M.cond == baseline_cond);
y_APP_AMB_peakHz = M.Wake_peakHz(M.geno == genoMut & M.cond == ambtemp_cond);

y_WT_BASE_peakAmp = M.Wake_peakAmp(M.geno == genoWT & M.cond == baseline_cond);
y_WT_AMB_peakAmp = M.Wake_peakAmp(M.geno == genoWT & M.cond == ambtemp_cond);
y_APP_BASE_peakAmp = M.Wake_peakAmp(M.geno == genoMut & M.cond == baseline_cond);
y_APP_AMB_peakAmp = M.Wake_peakAmp(M.geno == genoMut & M.cond == ambtemp_cond);

% ACh cycle frequency (if available)
if ismember('Wake_cycleHz', M.Properties.VariableNames)
    y_WT_BASE_cycleHz = M.Wake_cycleHz(M.geno == genoWT & M.cond == baseline_cond);
    y_WT_AMB_cycleHz = M.Wake_cycleHz(M.geno == genoWT & M.cond == ambtemp_cond);
    y_APP_BASE_cycleHz = M.Wake_cycleHz(M.geno == genoMut & M.cond == baseline_cond);
    y_APP_AMB_cycleHz = M.Wake_cycleHz(M.geno == genoMut & M.cond == ambtemp_cond);
else
    y_WT_BASE_cycleHz = []; y_WT_AMB_cycleHz = [];
    y_APP_BASE_cycleHz = []; y_APP_AMB_cycleHz = [];
end

% Wake slope (if available)
if ismember('slope_Wake', M.Properties.VariableNames)
    y_WT_BASE_slope = M.slope_Wake(M.geno == genoWT & M.cond == baseline_cond);
    y_WT_AMB_slope = M.slope_Wake(M.geno == genoWT & M.cond == ambtemp_cond);
    y_APP_BASE_slope = M.slope_Wake(M.geno == genoMut & M.cond == baseline_cond);
    y_APP_AMB_slope = M.slope_Wake(M.geno == genoMut & M.cond == ambtemp_cond);
else
    y_WT_BASE_slope = []; y_WT_AMB_slope = [];
    y_APP_BASE_slope = []; y_APP_AMB_slope = [];
end

% Mouse IDs for each group
mice_WT_BASE = M.mouse(M.geno == genoWT & M.cond == baseline_cond);
mice_WT_AMB = M.mouse(M.geno == genoWT & M.cond == ambtemp_cond);
mice_APP_BASE = M.mouse(M.geno == genoMut & M.cond == baseline_cond);
mice_APP_AMB = M.mouse(M.geno == genoMut & M.cond == ambtemp_cond);

% ---------- figure ----------
fig = figure('Name', 'Temperature Wake ACh PSD', 'Color', 'w', 'Position', [100 100 1680 300]);

% PSD curves (all 4 groups overlaid)
subplot(1, 6, 1); hold on;
if ~isempty(WT_BASE_psd)
    m = mean(WT_BASE_psd, 2);
    se = std(WT_BASE_psd, 0, 2) / sqrt(size(WT_BASE_psd, 2));
    fill_between(f_ref, m - se, m + se, COL_WT_BASE, 0.2);
    plot(f_ref, m, 'Color', COL_WT_BASE, 'LineWidth', 1.8);
end
if ~isempty(WT_AMB_psd)
    m = mean(WT_AMB_psd, 2);
    se = std(WT_AMB_psd, 0, 2) / sqrt(size(WT_AMB_psd, 2));
    fill_between(f_ref, m - se, m + se, COL_WT_AMB, 0.2);
    plot(f_ref, m, 'Color', COL_WT_AMB, 'LineWidth', 1.8);
end
if ~isempty(APP_BASE_psd)
    m = mean(APP_BASE_psd, 2);
    se = std(APP_BASE_psd, 0, 2) / sqrt(size(APP_BASE_psd, 2));
    fill_between(f_ref, m - se, m + se, COL_APP_BASE, 0.2);
    plot(f_ref, m, 'Color', COL_APP_BASE, 'LineWidth', 1.8);
end
if ~isempty(APP_AMB_psd)
    m = mean(APP_AMB_psd, 2);
    se = std(APP_AMB_psd, 0, 2) / sqrt(size(APP_AMB_psd, 2));
    fill_between(f_ref, m - se, m + se, COL_APP_AMB, 0.2);
    plot(f_ref, m, 'Color', COL_APP_AMB, 'LineWidth', 1.8);
end
xlabel('Frequency (Hz)');
ylabel('Power (A.U.)');
title('Wake ACh PSD');
box off;
xlim([min(f_ref) max(f_ref)]);
legend({sprintf('%s Base', genoWT), sprintf('%s Amb', genoWT), ...
        sprintf('%s Base', genoMut), sprintf('%s Amb', genoMut)}, ...
       'Box', 'off', 'Location', 'northeast', 'FontSize', 7);

% Bar plots with 2-way ANOVA
subplot(1, 6, 2);
plot_four_group_bar(y_WT_BASE_power, y_WT_AMB_power, y_APP_BASE_power, y_APP_AMB_power, ...
    genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
    'Wake power (A.U.)', show_mouse_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB);

subplot(1, 6, 3);
plot_four_group_bar(y_WT_BASE_peakHz, y_WT_AMB_peakHz, y_APP_BASE_peakHz, y_APP_AMB_peakHz, ...
    genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
    'Wake peak freq (Hz)', show_mouse_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB);

subplot(1, 6, 4);
plot_four_group_bar(y_WT_BASE_peakAmp, y_WT_AMB_peakAmp, y_APP_BASE_peakAmp, y_APP_AMB_peakAmp, ...
    genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
    'Wake peak amp (A.U.)', show_mouse_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB);

subplot(1, 6, 5);
if ~isempty([y_WT_BASE_cycleHz; y_WT_AMB_cycleHz; y_APP_BASE_cycleHz; y_APP_AMB_cycleHz])
    plot_four_group_bar(y_WT_BASE_cycleHz, y_WT_AMB_cycleHz, y_APP_BASE_cycleHz, y_APP_AMB_cycleHz, ...
        genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
        'ACh cycles (Hz)', show_mouse_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB);
else
    text(0.5, 0.5, 'No cycle data', 'HorizontalAlignment', 'center');
    axis off;
end

subplot(1, 6, 6);
if ~isempty([y_WT_BASE_slope; y_WT_AMB_slope; y_APP_BASE_slope; y_APP_AMB_slope])
    plot_four_group_bar(y_WT_BASE_slope, y_WT_AMB_slope, y_APP_BASE_slope, y_APP_AMB_slope, ...
        genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
        'Wake slope (ΔF/F/s)', show_mouse_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB);
else
    text(0.5, 0.5, 'No slope data', 'HorizontalAlignment', 'center');
    axis off;
end

sgtitle('Temperature Comparison – Wake ACh PSD & metrics (Baseline vs Ambtemp)');

if ~isempty(out_dir)
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    saveas(fig, fullfile(out_dir, 'Temperature_Wake_ACh_PSD.png'));
end
end

% ----------------------------------------------------------------------
function plot_four_group_bar(y_WT_BASE, y_WT_AMB, y_APP_BASE, y_APP_AMB, ...
    genoWT, genoMut, COL_WT_BASE, COL_WT_AMB, COL_APP_BASE, COL_APP_AMB, ...
    yLabel, show_ids, mice_WT_BASE, mice_WT_AMB, mice_APP_BASE, mice_APP_AMB)

hold on;

% Bar positions
x_positions = [1, 1.6, 2.8, 3.4];
bar_width = 0.45;

% Plot bars
bar(x_positions(1), mean(y_WT_BASE, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
bar(x_positions(2), mean(y_WT_AMB, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
bar(x_positions(3), mean(y_APP_BASE, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
bar(x_positions(4), mean(y_APP_AMB, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);

% Scatter points with jitter
jitter = 0.08;
plot_scatter_group(x_positions(1), y_WT_BASE, show_ids, mice_WT_BASE, [0.2 0.2 0.2]);
plot_scatter_group(x_positions(2), y_WT_AMB, show_ids, mice_WT_AMB, [0.2 0.2 0.2]);
plot_scatter_group(x_positions(3), y_APP_BASE, show_ids, mice_APP_BASE, [0.2 0.2 0.2]);
plot_scatter_group(x_positions(4), y_APP_AMB, show_ids, mice_APP_AMB, [0.2 0.2 0.2]);

xlim([0.5 3.9]);
set(gca, 'XTick', [1.3, 3.1], 'XTickLabel', {char(genoWT), char(genoMut)});
ylabel(yLabel);
box off;

% Perform 2-way ANOVA
all_y = [y_WT_BASE; y_WT_AMB; y_APP_BASE; y_APP_AMB];
cond_factor = [ones(size(y_WT_BASE)); 2*ones(size(y_WT_AMB)); ...
               ones(size(y_APP_BASE)); 2*ones(size(y_APP_AMB))];
geno_factor = [ones(size(y_WT_BASE)); ones(size(y_WT_AMB)); ...
               2*ones(size(y_APP_BASE)); 2*ones(size(y_APP_AMB))];

ok = ~isnan(all_y);
if sum(ok) > 3 && numel(unique(cond_factor(ok))) > 1 && numel(unique(geno_factor(ok))) > 1
    [p, ~, ~] = anovan(all_y(ok), {cond_factor(ok), geno_factor(ok)}, ...
        'model', 'interaction', 'display', 'off', 'varnames', {'Temp', 'Geno'});
    
    p_temp = p(1);
    p_geno = p(2);
    p_int = p(3);
    
    % Display p-values on plot
    yMax = max(all_y, [], 'omitnan');
    yMin = min(all_y, [], 'omitnan');
    yrange = max(yMax - yMin, eps);
    
    txt = sprintf('T:%s G:%s I:%s', p_to_star(p_temp), p_to_star(p_geno), p_to_star(p_int));
    text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7);
end
end

% ----------------------------------------------------------------------
function plot_scatter_group(xpos, y_vals, show_ids, mouse_ids, color)
if isempty(y_vals); return; end

jitter = 0.08;
xpos_jittered = xpos + (rand(size(y_vals)) - 0.5) * 2 * jitter;
scatter(xpos_jittered, y_vals, 25, color, 'filled');

if show_ids && ~isempty(mouse_ids)
    for i = 1:numel(y_vals)
        if i <= numel(mouse_ids)
            text(xpos_jittered(i), y_vals(i), sprintf(' %s', mouse_ids{i}), ...
                'FontSize', 7, 'HorizontalAlignment', 'left');
        end
    end
end
end

% ----------------------------------------------------------------------
function fill_between(x, y1, y2, col, alphaVal)
x = x(:);
y1 = y1(:);
y2 = y2(:);
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
    'EdgeColor', 'none', 'FaceAlpha', alphaVal);
end

% ----------------------------------------------------------------------
function str = p_to_star(p)
% Convert p-value to star notation for significance
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