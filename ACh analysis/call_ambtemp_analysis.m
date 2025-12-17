%% 3h cropped analysis: 3h baseline vs 3h ambient temperature
% Produces:
%   - transitions (Wake/NREM/REM) traces + peak bars + t-tests
%   - ANOVA bar plots for group metrics (incl. NREM PSD metrics, slopes…)
%   - NREM PSD traces + "peak oscillation" bars for each condition separately

% -------- EDIT THIS PATH --------
sigDir3h = '/Users/margaridaseabra/24.11_cropped_ambtemp';
% --------------------------------
scoreDir3h = sigDir3h;   % mat + csv in same folder
outDir3h   = fullfile(sigDir3h, 'ACh_group_3h');

if ~exist(outDir3h,'dir'); mkdir(outDir3h); end

%% 1) Run batch on 3h cropped files
GROUP_3h_all = run_ach_batch_auto(sigDir3h, scoreDir3h);

%% Quick diagnostic: Check transition counts
fprintf('\n=== Checking transitions per session ===\n');
for k = 1:numel(GROUP_3h_all.sessions)
    sess = GROUP_3h_all.sessions(k);
    OUT = GROUP_3h_all.out{k};
    fprintf('%s (%s, %s): ', sess.mouse, sess.geno, sess.cond);
    if isfield(OUT, 'transitions')
        for t = 1:numel(OUT.transitions)
            fprintf('%s=%d events, ', OUT.transitions(t).name, OUT.transitions(t).n_events);
        end
        fprintf('\n');
    else
        fprintf('NO transitions\n');
    end
end
fprintf('========================================\n\n');

%% 3) Transitions: Wake / NREM / REM onset
%   - rows = conditions (baseline, ambtemp)
%   - left = mean ± SEM traces
%   - right = peak bars + t-tests
group_plot_temperature_transitions(GROUP_3h_all, outDir3h);

%% 4) Group metrics across conditions (two-way ANOVA: cond x genotype)
% By default ach_stats_anova uses:
%   NREM_power, NREM_peakHz, slope_NREM, WakeOn_peak_mean
% but you can add more metrics inside that function.
ach_stats_temperature_anova(GROUP_3h_all, outDir3h, '3h_baseline_vs_ambtemp');

%% 5) Wake / NREM / REM onset traces + peak bar plots + stats
group_plot_temperature_transitions(GROUP_3h_all, outDir3h,true);

% 6) NREM ACh PSD traces + "peak oscillation" bars (power, freq, amp)
ach_plot_nrem_psd_temperature(GROUP_3h_all, outDir3h,true);

% 7) Wake ACh PSD traces + "peak oscillation" bars (power, freq, amp)
ach_plot_Wake_psd_temperature(GROUP_3h_all, outDir3h,true);





%%Other

%% 5) NREM ACh PSD + "peak oscillation" metrics, per condition
% This gives you NE-style panels: PSD trace + bars for power/freq/amp.
conds3h = {'baseline','ambtemp'};
for i = 1:numel(conds3h)
    c = conds3h{i};
    Gc = subset_group_by_cond(GROUP_3h, c);
    out_c = fullfile(outDir3h, ['PSD_' c]);
    ach_plot_nrem_psd_baseline(Gc, out_c);
end

fprintf('\n[3h analysis] Done. Figures saved to %s\n', outDir3h);
