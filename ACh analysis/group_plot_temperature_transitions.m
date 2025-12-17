function group_plot_temperature_transitions(GROUP, out_dir, show_mouse_ids, time_label)
% Compare baseline vs ambient temperature transitions with 2-way ANOVA
% 
% INPUTS:
%   GROUP - structure with sessions and transitions
%   out_dir - directory to save figures (optional)
%   show_mouse_ids - logical, if true show mouse IDs on scatter plots (default: false)
%   time_label - string describing time period (e.g., '0-3h', 'Full Recording') (optional)

if nargin < 2
    out_dir = [];
end
if nargin < 3 || isempty(show_mouse_ids)
    show_mouse_ids = false;
end
if nargin < 4 || isempty(time_label)
    time_label = 'Full Recording';
    warning('time_label not provided, using default: Full Recording');
end

trans_list  = {'Wake_onset','NREM_onset','REM_onset'};
title_list  = {'WAKE ONSETS','NREM ONSETS','REM ONSETS'};

% === DIAGNOSTIC: Check what transitions are available ===
fprintf('\n=== Transition Data Availability ===\n');
sessions = GROUP.sessions;
for k = 1:numel(sessions)
    fprintf('Session %d (%s, %s, %s): ', k, sessions(k).mouse, sessions(k).geno, sessions(k).cond);
    OUT = GROUP.out{k};
    if isfield(OUT, 'transitions') && ~isempty(OUT.transitions)
        trans_names = {OUT.transitions.name};
        fprintf('%s\n', strjoin(trans_names, ', '));
    else
        fprintf('NO TRANSITIONS\n');
    end
end
fprintf('====================================\n\n');

% Plot each transition type
for i = 1:numel(trans_list)
    % IMPORTANT: Pass ALL parameters including time_label
    fig = plot_one_temperature_transition(GROUP, trans_list{i}, title_list{i}, show_mouse_ids, time_label);
    
    if ~isempty(out_dir)
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        fname = fullfile(out_dir, sprintf('Transitions_%s_%s.png', trans_list{i}, strrep(time_label, ' ', '_')));
        saveas(fig, fname);
    end
end
end

% ======================================================================
function fig = plot_one_temperature_transition(GROUP, trans_name, big_title, show_mouse_ids, time_label)

sessions = GROUP.sessions;

% Identify baseline and ambient temperature conditions
conds = unique(string({sessions.cond}),'stable');
baseline_cond = conds(contains(lower(conds), 'baseline'));
ambient_cond = conds(contains(lower(conds), 'ambtemp') | contains(lower(conds), 'ambient'));

if isempty(baseline_cond) || isempty(ambient_cond)
    error('Could not identify baseline and ambient temperature conditions. Conditions found: %s', strjoin(conds, ', '));
end

baseline_cond = baseline_cond(1);
ambient_cond = ambient_cond(1);

% figure out genotype labels (WT first if present)
geno_vals = unique(string({sessions.geno}),'stable');
if any(geno_vals=="WT")
    genoWT  = "WT";
    genoMut = geno_vals(geno_vals~="WT");
else
    genoWT  = geno_vals(1);
    genoMut = geno_vals(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes (WT + mutant).');
end
genoMut = genoMut(1);

% Define colors
COL_WT_BASE = [0.6 0.6 0.6];      % grey - WT baseline
COL_WT_AMB = [0.3 0.3 0.3];       % darker grey - WT ambient
COL_APP_BASE = [0.392 0.584 0.929]; % cornflower blue - APP baseline
COL_APP_AMB = [0.196 0.292 0.465]; % darker blue - APP ambient

% === DIAGNOSTIC ===
fprintf('\n=== %s - Temperature Comparison ===\n', trans_name);
fprintf('Baseline condition: %s\n', baseline_cond);
fprintf('Ambient condition: %s\n', ambient_cond);
fprintf('WT genotype: %s\n', genoWT);
fprintf('Mutant genotype: %s\n', genoMut);

% Check distribution
for c = 1:2
    cond_arr = [baseline_cond, ambient_cond];
    cond = cond_arr(c);
    n_wt = sum(string({sessions.cond}) == cond & string({sessions.geno}) == genoWT);
    n_app = sum(string({sessions.cond}) == cond & string({sessions.geno}) == genoMut);
    fprintf(' %s: WT n=%d, %s n=%d\n', cond, n_wt, genoMut, n_app);
end
fprintf('=====================================\n\n');

fig = figure('Name',sprintf('%s – %s (%s)',big_title,trans_name,time_label),'Color','w', ...
    'Position',[100 100 1400 400]);

% Collect data for all groups
[WT_BASE_traces, WT_BASE_peaks, WT_BASE_mice] = collect_data(sessions, GROUP, baseline_cond, genoWT, trans_name);
[WT_AMB_traces, WT_AMB_peaks, WT_AMB_mice] = collect_data(sessions, GROUP, ambient_cond, genoWT, trans_name);
[APP_BASE_traces, APP_BASE_peaks, APP_BASE_mice] = collect_data(sessions, GROUP, baseline_cond, genoMut, trans_name);
[APP_AMB_traces, APP_AMB_peaks, APP_AMB_mice] = collect_data(sessions, GROUP, ambient_cond, genoMut, trans_name);

% Get common time axis
t_ref = get_time_axis(sessions, GROUP, trans_name);

fprintf('%s: WT_BASE=%d, WT_AMB=%d, APP_BASE=%d, APP_AMB=%d traces\n', ...
    trans_name, size(WT_BASE_traces,2), size(WT_AMB_traces,2), ...
    size(APP_BASE_traces,2), size(APP_AMB_traces,2));

% ========== LEFT: Mean traces ± SEM ==================
subplot(1,2,1); hold on;

% Plot all four groups
if ~isempty(WT_BASE_traces)
    plot_trace_with_sem(t_ref, WT_BASE_traces, COL_WT_BASE);
end
if ~isempty(WT_AMB_traces)
    plot_trace_with_sem(t_ref, WT_AMB_traces, COL_WT_AMB);
end
if ~isempty(APP_BASE_traces)
    plot_trace_with_sem(t_ref, APP_BASE_traces, COL_APP_BASE);
end
if ~isempty(APP_AMB_traces)
    plot_trace_with_sem(t_ref, APP_AMB_traces, COL_APP_AMB);
end

plot([0 0], ylim, 'k--', 'LineWidth', 1);
xlabel('Time from transition (s)');
ylabel('ACh (ΔF/F)');
title(sprintf('%s - Temperature Comparison', big_title));
box off;

legend({[char(genoWT) ' Baseline'], [char(genoWT) ' Ambtemp'], ...
        [char(genoMut) ' Baseline'], [char(genoMut) ' Ambtemp']}, ...
       'Box', 'off', 'Location', 'best', 'FontSize', 9);

% ========== RIGHT: Peak comparison with 2-way ANOVA ==================
subplot(1,2,2); hold on;

% Prepare data for 2-way ANOVA and plotting
all_peaks = [WT_BASE_peaks; WT_AMB_peaks; APP_BASE_peaks; APP_AMB_peaks];
genotype_factor = [repmat({'WT'}, numel(WT_BASE_peaks) + numel(WT_AMB_peaks), 1);
                   repmat({char(genoMut)}, numel(APP_BASE_peaks) + numel(APP_AMB_peaks), 1)];
condition_factor = [repmat({'Baseline'}, numel(WT_BASE_peaks), 1);
                    repmat({'Ambtemp'}, numel(WT_AMB_peaks), 1);
                    repmat({'Baseline'}, numel(APP_BASE_peaks), 1);
                    repmat({'Ambtemp'}, numel(APP_AMB_peaks), 1)];

% Bar plots with grouped layout
x_positions = [1, 1.8, 3, 3.8]; % WT_BASE, WT_AMB, APP_BASE, APP_AMB
bar_width = 0.6;

bar(x_positions(1), mean(WT_BASE_peaks, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
bar(x_positions(2), mean(WT_AMB_peaks, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
bar(x_positions(3), mean(APP_BASE_peaks, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
bar(x_positions(4), mean(APP_AMB_peaks, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);

% Scatter individual points with jitter
jitter = 0.08;
plot_scatter_with_labels(x_positions(1), WT_BASE_peaks, COL_WT_BASE, jitter, show_mouse_ids, WT_BASE_mice);
plot_scatter_with_labels(x_positions(2), WT_AMB_peaks, COL_WT_AMB, jitter, show_mouse_ids, WT_AMB_mice);
plot_scatter_with_labels(x_positions(3), APP_BASE_peaks, COL_APP_BASE, jitter, show_mouse_ids, APP_BASE_mice);
plot_scatter_with_labels(x_positions(4), APP_AMB_peaks, COL_APP_AMB, jitter, show_mouse_ids, APP_AMB_mice);

xlim([0.5 4.3]);
set(gca,'XTick', [1.4, 3.4], 'XTickLabel', {char(genoWT), char(genoMut)});
ylabel('Peak ACh (ΔF/F)');
title('Peak Response - 2-way ANOVA');
box off;

% Perform 2-way ANOVA
if numel(all_peaks) > 3 && numel(unique(genotype_factor)) > 1 && numel(unique(condition_factor)) > 1
    [p_values, tbl, stats] = anovan(all_peaks, {genotype_factor, condition_factor}, ...
        'model', 'interaction', 'varnames', {'Genotype', 'Condition'}, 'display', 'off');
    
    % Display results on plot
    ymax = max(all_peaks(:));
    ypos = ymax * 1.15;
    
    % Format p-values
    text_str = {sprintf('Temperature: p=%s', format_pval(p_values(1))), ...
                sprintf('Genotype: p=%s', format_pval(p_values(2))), ...
                sprintf('Interaction: p=%s', format_pval(p_values(3)))};
    
    text(0.7, ypos, text_str, 'FontSize', 8, 'VerticalAlignment', 'top');
    
    fprintf('\n2-way ANOVA results for %s:\n', trans_name);
    fprintf('  Temperature effect: p = %.4f %s\n', p_values(1), p_to_star(p_values(1)));
    fprintf('  Genotype effect: p = %.4f %s\n', p_values(2), p_to_star(p_values(2)));
    fprintf('  Interaction: p = %.4f %s\n', p_values(3), p_to_star(p_values(3)));
end

% Add time period label
annotation('textbox', [0.02, 0.98, 0.3, 0.02], 'String', sprintf('Time Period: %s', time_label), ...
    'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

end

% ======================================================================
% Helper functions
% ======================================================================

function [traces, peaks, mouse_ids] = collect_data(sessions, GROUP, cond, geno, trans_name)
% Collect traces and peaks for a specific condition and genotype
traces = [];
peaks = [];
mouse_ids = {};
t_ref = [];

for k = 1:numel(sessions)
    if string(sessions(k).cond) ~= cond || string(sessions(k).geno) ~= geno
        continue;
    end
    
    OUT = GROUP.out{k};
    idxT = [];
    if isfield(OUT,'transitions') && ~isempty(OUT.transitions)
        idxT = find(strcmp({OUT.transitions.name}, trans_name));
    end
    
    if isempty(idxT)
        continue;
    end
    
    Ttr = OUT.transitions(idxT);
    
    if ~isfield(Ttr, 'mean') || isempty(Ttr.mean) || all(isnan(Ttr.mean))
        continue;
    end
    
    % Get time axis
    t_this = get_time_vector(Ttr);
    
    % Set or match reference time axis
    if isempty(t_ref)
        t_ref = t_this;
    elseif numel(t_this) ~= numel(t_ref) || any(abs(t_this - t_ref) > 1e-6)
        Ttr.mean = interp1(t_this, Ttr.mean, t_ref, 'linear', NaN);
    end
    
    traces = [traces Ttr.mean(:)]; %#ok<AGROW>
    
    if isfield(Ttr, 'peak_mean') && ~isnan(Ttr.peak_mean)
        peaks = [peaks; Ttr.peak_mean]; %#ok<AGROW>
        mouse_ids{end+1} = sessions(k).mouse; %#ok<AGROW>
    elseif isfield(Ttr, 'peak') && ~isnan(Ttr.peak)
        % Try alternative field name
        peaks = [peaks; Ttr.peak]; %#ok<AGROW>
        mouse_ids{end+1} = sessions(k).mouse; %#ok<AGROW>
    else
        % Compute peak from mean trace if not stored
        peak_val = max(Ttr.mean);
        if ~isnan(peak_val)
            peaks = [peaks; peak_val]; %#ok<AGROW>
            mouse_ids{end+1} = sessions(k).mouse; %#ok<AGROW>
        end
    end
end
end

function t_ref = get_time_axis(sessions, GROUP, trans_name)
% Get a reference time axis from the first valid transition
t_ref = [];
for k = 1:numel(sessions)
    OUT = GROUP.out{k};
    if ~isfield(OUT,'transitions') || isempty(OUT.transitions)
        continue;
    end
    idxT = find(strcmp({OUT.transitions.name}, trans_name));
    if ~isempty(idxT)
        Ttr = OUT.transitions(idxT);
        if isfield(Ttr, 'mean') && ~isempty(Ttr.mean)
            t_ref = get_time_vector(Ttr);
            return;
        end
    end
end
% Default if nothing found
t_ref = linspace(-100, 100, 100)';
end

function t_vec = get_time_vector(Ttr)
% Extract time vector from transition structure
if isfield(Ttr, 't_axis')
    t_vec = Ttr.t_axis(:);
elseif isfield(Ttr, 'time')
    t_vec = Ttr.time(:);
elseif isfield(Ttr, 't')
    t_vec = Ttr.t(:);
else
    if isfield(Ttr, 'win_sec')
        n_pts = numel(Ttr.mean);
        t_vec = linspace(Ttr.win_sec(1), Ttr.win_sec(2), n_pts)';
    else
        n_pts = numel(Ttr.mean);
        t_vec = linspace(-100, 100, n_pts)';
    end
end
end

function plot_trace_with_sem(t_ref, traces, color)
% Plot mean trace with SEM shading
if isempty(traces) || all(isnan(traces(:)))
    return;
end

m = mean(traces, 2, 'omitnan');
se = std(traces, 0, 2, 'omitnan') / sqrt(size(traces,2));

% Shaded SEM
fill_between(t_ref, m - se, m + se, color, 0.2);

% Mean line
plot(t_ref, m, 'Color', color, 'LineWidth', 1.5);
end

function fill_between(x, y1, y2, col, alphaVal)
% Helper for shaded SEM
x = x(:); y1 = y1(:); y2 = y2(:);
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
    'EdgeColor','none','FaceAlpha',alphaVal);
end

function plot_scatter_with_labels(xpos, peaks, color, jitter, show_labels, mouse_ids)
% Scatter plot with optional mouse ID labels
if isempty(peaks)
    return;
end

xpos_jittered = xpos + (rand(size(peaks))-0.5)*2*jitter;
scatter(xpos_jittered, peaks, 25, color, 'filled');

if show_labels && ~isempty(mouse_ids)
    for i = 1:numel(peaks)
        if i <= numel(mouse_ids)
            text(xpos_jittered(i), peaks(i), sprintf(' %s', mouse_ids{i}), ...
                'FontSize', 7, 'HorizontalAlignment', 'left');
        end
    end
end
end

function str = format_pval(p)
% Format p-value for display
if p < 0.001
    str = '<0.001';
elseif p < 0.01
    str = sprintf('=%.3f', p);
else
    str = sprintf('=%.2f', p);
end
end

function star = p_to_star(p)
% Convert p-value to significance stars
if p < 0.001
    star = '***';
elseif p < 0.01
    star = '**';
elseif p < 0.05
    star = '*';
else
    star = '';
end
end