function group_plot_three_conditions_transitions(GROUP, out_dir, show_mouse_ids, time_label)
% Compare baseline vs ambtemp vs drug transitions with mixed-effects ANOVA
% Colors match: WT baseline (light grey), APP baseline (light blue),
%               WT ambtemp (medium grey), APP ambtemp (medium blue),
%               WT drug (dark grey), APP drug (dark blue)

if nargin < 2; out_dir = []; end
if nargin < 3 || isempty(show_mouse_ids); show_mouse_ids = false; end
if nargin < 4 || isempty(time_label); time_label = 'Full Recording'; end

trans_list = {'Wake_onset','NREM_onset','REM_onset'};
title_list = {'WAKE ONSETS','NREM ONSETS','REM ONSETS'};

for i = 1:numel(trans_list)
    fig = plot_one_three_cond_transition(GROUP, trans_list{i}, title_list{i}, show_mouse_ids, time_label);
    
    if ~isempty(out_dir)
        if ~exist(out_dir,'dir'); mkdir(out_dir); end
        fname = fullfile(out_dir, sprintf('Transitions_%s_%s.png', trans_list{i}, strrep(time_label, ' ', '_')));
        saveas(fig, fname);
    end
end
end

function fig = plot_one_three_cond_transition(GROUP, trans_name, big_title, show_mouse_ids, time_label)

sessions = GROUP.sessions;
conds = unique(string({sessions.cond}),'stable');

% Identify conditions
baseline_cond = conds(contains(lower(conds), 'baseline'));
ambtemp_cond = conds(contains(lower(conds), 'ambtemp') | contains(lower(conds), 'ambient'));
drugs_cond = conds(contains(lower(conds), 'drugs') | contains(lower(conds), 'drug'));

if isempty(baseline_cond) || isempty(ambtemp_cond) || isempty(drugs_cond)
    error('Could not identify all three conditions. Found: %s', strjoin(conds, ', '));
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);
drugs_cond = drugs_cond(1);

% Identify genotypes
geno_vals = unique(string({sessions.geno}),'stable');
if any(geno_vals=="WT")
    genoWT = "WT";
    genoMut = geno_vals(geno_vals~="WT");
else
    genoWT = geno_vals(1);
    genoMut = geno_vals(2:end);
end
genoMut = genoMut(1);

% Define colors matching your image
COL_WT_BASE = [0.7 0.7 0.7];        % Light grey - WT baseline
COL_APP_BASE = [0.5 0.65 0.85];     % Light blue - APP baseline
COL_WT_AMB = [0.5 0.5 0.5];         % Medium grey - WT ambtemp
COL_APP_AMB = [0.3 0.5 0.75];       % Medium blue - APP ambtemp
COL_WT_DRUGS = [0.3 0.3 0.3];       % Dark grey - WT drugs
COL_APP_DRUGS = [0.2 0.35 0.6];     % Dark blue - APP drugs

fig = figure('Name',sprintf('%s – %s (%s)',big_title,trans_name,time_label),'Color','w', ...
    'Position',[100 100 1600 400]);

% Collect data for all 6 groups
[WT_BASE_tr, WT_BASE_pk, WT_BASE_mice] = collect_data(sessions, GROUP, baseline_cond, genoWT, trans_name);
[APP_BASE_tr, APP_BASE_pk, APP_BASE_mice] = collect_data(sessions, GROUP, baseline_cond, genoMut, trans_name);
[WT_AMB_tr, WT_AMB_pk, WT_AMB_mice] = collect_data(sessions, GROUP, ambtemp_cond, genoWT, trans_name);
[APP_AMB_tr, APP_AMB_pk, APP_AMB_mice] = collect_data(sessions, GROUP, ambtemp_cond, genoMut, trans_name);
[WT_DRUGS_tr, WT_DRUGS_pk, WT_DRUGS_mice] = collect_data(sessions, GROUP, drugs_cond, genoWT, trans_name);
[APP_DRUGS_tr, APP_DRUGS_pk, APP_DRUGS_mice] = collect_data(sessions, GROUP, drugs_cond, genoMut, trans_name);

t_ref = get_time_axis(sessions, GROUP, trans_name);

% ========== LEFT: Mean traces ± SEM ==================
subplot(1,2,1); hold on;

if ~isempty(WT_BASE_tr); plot_trace_with_sem(t_ref, WT_BASE_tr, COL_WT_BASE); end
if ~isempty(APP_BASE_tr); plot_trace_with_sem(t_ref, APP_BASE_tr, COL_APP_BASE); end
if ~isempty(WT_AMB_tr); plot_trace_with_sem(t_ref, WT_AMB_tr, COL_WT_AMB); end
if ~isempty(APP_AMB_tr); plot_trace_with_sem(t_ref, APP_AMB_tr, COL_APP_AMB); end
if ~isempty(WT_DRUGS_tr); plot_trace_with_sem(t_ref, WT_DRUGS_tr, COL_WT_DRUGS); end
if ~isempty(APP_DRUGS_tr); plot_trace_with_sem(t_ref, APP_DRUGS_tr, COL_APP_DRUGS); end

plot([0 0], ylim, 'k--', 'LineWidth', 1);
xlabel('Time from transition (s)');
ylabel('ACh (ΔF/F)');
title(sprintf('%s - Three Condition Comparison', big_title));
box off;

legend({[char(genoWT) ' Baseline'], [char(genoMut) ' Baseline'], ...
        [char(genoWT) ' Ambtemp'], [char(genoMut) ' Ambtemp'], ...
        [char(genoWT) ' Drugs'], [char(genoMut) ' Drugs']}, ...
       'Box', 'off', 'Location', 'best', 'FontSize', 8);

% ========== RIGHT: Peak bars ==================
subplot(1,2,2); hold on;

x_positions = [1, 1.6, 2.8, 3.4, 4.6, 5.2];
bar_width = 0.45;

bar(x_positions(1), mean(WT_BASE_pk, 'omitnan'), bar_width, 'FaceColor', COL_WT_BASE);
bar(x_positions(2), mean(APP_BASE_pk, 'omitnan'), bar_width, 'FaceColor', COL_APP_BASE);
bar(x_positions(3), mean(WT_AMB_pk, 'omitnan'), bar_width, 'FaceColor', COL_WT_AMB);
bar(x_positions(4), mean(APP_AMB_pk, 'omitnan'), bar_width, 'FaceColor', COL_APP_AMB);
bar(x_positions(5), mean(WT_DRUGS_pk, 'omitnan'), bar_width, 'FaceColor', COL_WT_DRUGS);
bar(x_positions(6), mean(APP_DRUGS_pk, 'omitnan'), bar_width, 'FaceColor', COL_APP_DRUGS);

jitter = 0.06;
plot_scatter(x_positions(1), WT_BASE_pk, COL_WT_BASE, jitter, show_mouse_ids, WT_BASE_mice);
plot_scatter(x_positions(2), APP_BASE_pk, COL_APP_BASE, jitter, show_mouse_ids, APP_BASE_mice);
plot_scatter(x_positions(3), WT_AMB_pk, COL_WT_AMB, jitter, show_mouse_ids, WT_AMB_mice);
plot_scatter(x_positions(4), APP_AMB_pk, COL_APP_AMB, jitter, show_mouse_ids, APP_AMB_mice);
plot_scatter(x_positions(5), WT_DRUGS_pk, COL_WT_DRUGS, jitter, show_mouse_ids, WT_DRUGS_mice);
plot_scatter(x_positions(6), APP_DRUGS_pk, COL_APP_DRUGS, jitter, show_mouse_ids, APP_DRUGS_mice);

xlim([0.5 5.7]);
set(gca,'XTick', [1.3, 3.1, 4.9], 'XTickLabel', {'Baseline', 'Ambtemp', 'Drugs'});
ylabel('Peak ACh (ΔF/F)');
title('Peak Response - 2-way ANOVA');
box off;

% Two-way ANOVA
all_peaks = [WT_BASE_pk; APP_BASE_pk; WT_AMB_pk; APP_AMB_pk; WT_DRUGS_pk; APP_DRUGS_pk];
cond_factor = [repmat(1, numel(WT_BASE_pk), 1); repmat(1, numel(APP_BASE_pk), 1);
               repmat(2, numel(WT_AMB_pk), 1); repmat(2, numel(APP_AMB_pk), 1);
               repmat(3, numel(WT_DRUGS_pk), 1); repmat(3, numel(APP_DRUGS_pk), 1)];
geno_factor = [ones(numel(WT_BASE_pk), 1); 2*ones(numel(APP_BASE_pk), 1);
               ones(numel(WT_AMB_pk), 1); 2*ones(numel(APP_AMB_pk), 1);
               ones(numel(WT_DRUGS_pk), 1); 2*ones(numel(APP_DRUGS_pk), 1)];

ok = ~isnan(all_peaks);
if sum(ok) > 5
    [p, ~, ~] = anovan(all_peaks(ok), {cond_factor(ok), geno_factor(ok)}, ...
        'model', 'interaction', 'display', 'off', 'varnames', {'Condition', 'Genotype'});
    
    txt = sprintf('Cond: p=%s\nGeno: p=%s\nInt: p=%s', ...
        format_pval(p(1)), format_pval(p(2)), format_pval(p(3)));
    text(0.98, 0.98, txt, 'Units', 'normalized', 'FontSize', 8, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
end

annotation('textbox', [0.02, 0.98, 0.3, 0.02], 'String', sprintf('Time Period: %s', time_label), ...
    'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Helper functions
function [traces, peaks, mouse_ids] = collect_data(sessions, GROUP, cond, geno, trans_name)
% Collect traces and peaks for a specific condition and genotype
% Updated to compute peaks the same way as group_plot_temperature_transitions

traces = [];
peaks = [];
mouse_ids = {};
t_ref = [];

for k = 1:numel(sessions)

    % --- Match condition & genotype ---
    if string(sessions(k).cond) ~= cond || string(sessions(k).geno) ~= geno
        continue;
    end

    OUT = GROUP.out{k};
    if ~isfield(OUT,'transitions'); continue; end

    idxT = find(strcmp({OUT.transitions.name}, trans_name));
    if isempty(idxT); continue; end

    Ttr = OUT.transitions(idxT);

    % --- Must have a trace ---
    if ~isfield(Ttr, 'mean') || isempty(Ttr.mean) || all(isnan(Ttr.mean))
        continue;
    end

    % --- Time axis ---
    t_this = get_time_vector(Ttr);

    % --- Set / align to reference axis ---
    if isempty(t_ref)
        t_ref = t_this(:);
    elseif numel(t_this) ~= numel(t_ref) || any(abs(t_this(:) - t_ref) > 1e-6)
        try
            Ttr.mean = interp1(t_this(:), Ttr.mean(:), t_ref, 'linear', NaN);
        catch
            fprintf('Warning: Could not interpolate trace for session %d\n', k);
            continue;
        end
    end

    traces = [traces Ttr.mean(:)]; %#ok<AGROW>

    % --------------------------
    %     PEAK COMPUTATION
    % --------------------------
    if isfield(Ttr, 'peak_mean') && ~isnan(Ttr.peak_mean)
        peak_val = Ttr.peak_mean;

    elseif isfield(Ttr, 'peak') && ~isnan(Ttr.peak)
        peak_val = Ttr.peak;

    else
        % Compute peak from mean trace — SAME as in temperature version
        peak_val = max(Ttr.mean(:), [], 'omitnan');
    end

    if ~isnan(peak_val)
        peaks = [peaks; peak_val]; %#ok<AGROW>
        mouse_ids{end+1} = sessions(k).mouse; %#ok<AGROW>
    end

end
end

function t_ref = get_time_axis(sessions, GROUP, trans_name)
t_ref = linspace(-100, 100, 100)';
for k = 1:numel(sessions)
    OUT = GROUP.out{k};
    if ~isfield(OUT,'transitions'); continue; end
    idxT = find(strcmp({OUT.transitions.name}, trans_name));
    if ~isempty(idxT) && isfield(OUT.transitions(idxT), 'mean')
        t_ref = get_time_vector(OUT.transitions(idxT)); return;
    end
end
end

function t_vec = get_time_vector(Ttr)
if isfield(Ttr, 't_axis'); t_vec = Ttr.t_axis(:);
elseif isfield(Ttr, 'time'); t_vec = Ttr.time(:);
elseif isfield(Ttr, 't'); t_vec = Ttr.t(:);
elseif isfield(Ttr, 'win_sec')
    t_vec = linspace(Ttr.win_sec(1), Ttr.win_sec(2), numel(Ttr.mean))';
else
    t_vec = linspace(-100, 100, numel(Ttr.mean))';
end
end

function plot_trace_with_sem(t_ref, traces, color)
if isempty(traces); return; end

% Ensure traces match t_ref length
if size(traces, 1) ~= numel(t_ref)
    fprintf('Warning: Trace length (%d) does not match time axis (%d), skipping plot\n', ...
        size(traces, 1), numel(t_ref));
    return;
end

m = mean(traces, 2, 'omitnan');
se = std(traces, 0, 2, 'omitnan') / sqrt(size(traces,2));

% Ensure all are column vectors of same length
m = m(:); 
se = se(:); 
t_ref = t_ref(:);

% Check dimensions before fill
if numel(m) ~= numel(t_ref) || numel(se) ~= numel(t_ref)
    fprintf('Warning: Dimension mismatch - t_ref: %d, m: %d, se: %d\n', ...
        numel(t_ref), numel(m), numel(se));
    return;
end

fill([t_ref; flipud(t_ref)], [m-se; flipud(m+se)], color, 'EdgeColor','none','FaceAlpha',0.2);
plot(t_ref, m, 'Color', color, 'LineWidth', 1.5);
end

function plot_scatter(xpos, peaks, color, jitter, show_ids, mice)
if isempty(peaks); return; end
xpos_j = xpos + (rand(size(peaks))-0.5)*2*jitter;
scatter(xpos_j, peaks, 25, color, 'filled');
if show_ids && ~isempty(mice)
    for i = 1:numel(peaks)
        if i <= numel(mice)
            text(xpos_j(i), peaks(i), sprintf(' %s', mice{i}), 'FontSize', 6);
        end
    end
end
end

function str = format_pval(p)
if p < 0.001; str = '<0.001';
elseif p < 0.01; str = sprintf('%.3f', p);
else; str = sprintf('%.2f', p);
end
end