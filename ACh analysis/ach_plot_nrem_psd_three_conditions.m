function ach_plot_nrem_psd_three_conditions(GROUP, out_dir, show_mouse_ids, time_label)
% Plot NREM ACh PSD comparing baseline vs ambtemp vs DRUGS (WT vs APP) with 2-way ANOVA

if nargin < 2; out_dir = []; end
if nargin < 3 || isempty(show_mouse_ids); show_mouse_ids = false; end
if nargin < 4 || isempty(time_label); time_label = 'Full Recording'; end

sessions = GROUP.sessions;

% Identify conditions
all_conds = unique(string({sessions.cond}), 'stable');
baseline_cond = all_conds(contains(lower(all_conds), 'baseline'));
ambtemp_cond = all_conds(contains(lower(all_conds), 'ambtemp'));
drugs_cond = all_conds(contains(lower(all_conds), 'drugs') | contains(lower(all_conds), 'drug'));

if isempty(baseline_cond) || isempty(ambtemp_cond) || isempty(drugs_cond)
    error('Could not identify all three conditions');
end

baseline_cond = baseline_cond(1);
ambtemp_cond = ambtemp_cond(1);
drugs_cond = drugs_cond(1);

% Identify genotypes
geno_vals = unique(string({sessions.geno}), 'stable');
if any(geno_vals == "WT"); genoWT = "WT"; genoMut = geno_vals(geno_vals ~= "WT");
else; genoWT = geno_vals(1); genoMut = geno_vals(2:end); end
genoMut = genoMut(1);

% Define colors
COL_WT_BASE = [0.7 0.7 0.7];
COL_APP_BASE = [0.5 0.65 0.85];
COL_WT_AMB = [0.5 0.5 0.5];
COL_APP_AMB = [0.3 0.5 0.75];
COL_WT_DRUGS = [0.3 0.3 0.3];
COL_APP_DRUGS = [0.2 0.35 0.6];

% Collect PSD curves
f_ref = [];
[WT_BASE_psd, WT_AMB_psd, WT_DRUGS_psd] = deal([]);
[APP_BASE_psd, APP_AMB_psd, APP_DRUGS_psd] = deal([]);

for k = 1:numel(sessions)
    OUT = GROUP.out{k};
    if ~isfield(OUT, 'psd') || ~isfield(OUT.psd, 'NREM'); continue; end
    
    P = [];
    if isfield(OUT.psd.NREM, 'Pxx'); P = OUT.psd.NREM.Pxx(:);
    elseif isfield(OUT.psd.NREM, 'psd'); P = OUT.psd.NREM.psd(:);
    elseif isfield(OUT.psd.NREM, 'power'); P = OUT.psd.NREM.power(:);
    end
    if isempty(P) || ~isfield(OUT.psd.NREM, 'f'); continue; end
    
    f_this = OUT.psd.NREM.f(:);
    if isempty(f_ref); f_ref = f_this;
    elseif numel(f_this) ~= numel(f_ref)
        P = interp1(f_this, P, f_ref, 'linear', 'extrap');
    end
    
    g = string(sessions(k).geno);
    c = string(sessions(k).cond);
    
    if g == genoWT && c == baseline_cond; WT_BASE_psd = [WT_BASE_psd P];
    elseif g == genoWT && c == ambtemp_cond; WT_AMB_psd = [WT_AMB_psd P];
    elseif g == genoWT && c == drugs_cond; WT_DRUGS_psd = [WT_DRUGS_psd P];
    elseif g == genoMut && c == baseline_cond; APP_BASE_psd = [APP_BASE_psd P];
    elseif g == genoMut && c == ambtemp_cond; APP_AMB_psd = [APP_AMB_psd P];
    elseif g == genoMut && c == drugs_cond; APP_DRUGS_psd = [APP_DRUGS_psd P];
    end
end

% Get metrics
if isfield(GROUP, 'metrics_tbl'); M = GROUP.metrics_tbl;
else; M = struct2table(GROUP.metrics); end
M.geno = string(M.geno);
M.cond = string(M.cond);

% Extract metrics for each group
metrics = {'NREM_power', 'NREM_peakHz', 'NREM_peakAmp', 'NREM_cycleHz', 'slope_NREM'};
data = struct();
mice = struct();

for i = 1:numel(metrics)
    metric = metrics{i};
    if ~ismember(metric, M.Properties.VariableNames); continue; end
    
    data.(metric).WT_BASE = M.(metric)(M.geno == genoWT & M.cond == baseline_cond);
    data.(metric).APP_BASE = M.(metric)(M.geno == genoMut & M.cond == baseline_cond);
    data.(metric).WT_AMB = M.(metric)(M.geno == genoWT & M.cond == ambtemp_cond);
    data.(metric).APP_AMB = M.(metric)(M.geno == genoMut & M.cond == ambtemp_cond);
    data.(metric).WT_DRUGS = M.(metric)(M.geno == genoWT & M.cond == drugs_cond);
    data.(metric).APP_DRUGS = M.(metric)(M.geno == genoMut & M.cond == drugs_cond);
    
    mice.(metric).WT_BASE = M.mouse(M.geno == genoWT & M.cond == baseline_cond);
    mice.(metric).APP_BASE = M.mouse(M.geno == genoMut & M.cond == baseline_cond);
    mice.(metric).WT_AMB = M.mouse(M.geno == genoWT & M.cond == ambtemp_cond);
    mice.(metric).APP_AMB = M.mouse(M.geno == genoMut & M.cond == ambtemp_cond);
    mice.(metric).WT_DRUGS = M.mouse(M.geno == genoWT & M.cond == drugs_cond);
    mice.(metric).APP_DRUGS = M.mouse(M.geno == genoMut & M.cond == drugs_cond);
end

% Create figure
fig = figure('Name', 'NREM PSD - Three Conditions', 'Color', 'w', 'Position', [100 100 1680 300]);

% PSD curves
subplot(1, 6, 1); hold on;
if ~isempty(WT_BASE_psd); plot_psd(f_ref, WT_BASE_psd, COL_WT_BASE); end
if ~isempty(APP_BASE_psd); plot_psd(f_ref, APP_BASE_psd, COL_APP_BASE); end
if ~isempty(WT_AMB_psd); plot_psd(f_ref, WT_AMB_psd, COL_WT_AMB); end
if ~isempty(APP_AMB_psd); plot_psd(f_ref, APP_AMB_psd, COL_APP_AMB); end
if ~isempty(WT_DRUGS_psd); plot_psd(f_ref, WT_DRUGS_psd, COL_WT_DRUGS); end
if ~isempty(APP_DRUGS_psd); plot_psd(f_ref, APP_DRUGS_psd, COL_APP_DRUGS); end
xlabel('Frequency (Hz)'); ylabel('Power (A.U.)'); title('NREM ACh PSD'); box off;
xlim([min(f_ref) max(f_ref)]);
legend({'WT Base', 'APP Base', 'WT Amb', 'APP Amb', 'WT Drugs', 'APP Drugs'}, ...
    'Box', 'off', 'Location', 'northeast', 'FontSize', 7);

% Bar plots
subplot(1, 6, 2);
if isfield(data, 'NREM_power')
    plot_six_group_bar(data.NREM_power, mice.NREM_power, genoWT, genoMut, ...
        COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
        'NREM power (A.U.)', show_mouse_ids);
end

subplot(1, 6, 3);
if isfield(data, 'NREM_peakHz')
    plot_six_group_bar(data.NREM_peakHz, mice.NREM_peakHz, genoWT, genoMut, ...
        COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
        'NREM peak freq (Hz)', show_mouse_ids);
end

subplot(1, 6, 4);
if isfield(data, 'NREM_peakAmp')
    plot_six_group_bar(data.NREM_peakAmp, mice.NREM_peakAmp, genoWT, genoMut, ...
        COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
        'NREM peak amp (A.U.)', show_mouse_ids);
end

subplot(1, 6, 5);
if isfield(data, 'NREM_cycleHz')
    plot_six_group_bar(data.NREM_cycleHz, mice.NREM_cycleHz, genoWT, genoMut, ...
        COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
        'ACh cycles (Hz)', show_mouse_ids);
else
    text(0.5, 0.5, 'No cycle data', 'HorizontalAlignment', 'center'); axis off;
end

subplot(1, 6, 6);
if isfield(data, 'slope_NREM')
    plot_six_group_bar(data.slope_NREM, mice.slope_NREM, genoWT, genoMut, ...
        COL_WT_BASE, COL_APP_BASE, COL_WT_AMB, COL_APP_AMB, COL_WT_DRUGS, COL_APP_DRUGS, ...
        'NREM slope (ΔF/F/s)', show_mouse_ids);
else
    text(0.5, 0.5, 'No slope data', 'HorizontalAlignment', 'center'); axis off;
end

sgtitle(['Three Conditions – NREM ACh PSD & metrics – ' time_label]);

annotation('textbox', [0.02, 0.98, 0.3, 0.02], 'String', sprintf('Time Period: %s', time_label), ...
    'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

if ~isempty(out_dir)
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    saveas(fig, fullfile(out_dir, sprintf('NREM_PSD_ThreeCond_%s.png', strrep(time_label, ' ', '_'))));
end
end

% =========================================================================

% =========================================================================
% Helper functions
% =========================================================================




