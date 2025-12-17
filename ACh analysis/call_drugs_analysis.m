%% 6h cropped analysis: baseline vs ambient temperature vs drugs
% Produces:
%   - transitions (Wake/NREM/REM) traces + peak bars + t-tests
%   - ANOVA bar plots for group metrics (incl. NREM PSD metrics, slopes…)
%   - NREM PSD traces + "peak oscillation" bars for each condition separately

% -------- EDIT THIS PATH --------
sigDir6h = '/Users/margaridaseabra/24.11_cropped_3h_washout';
% --------------------------------
scoreDir6h = sigDir6h;   % mat + csv in same folder
outDir6h   = fullfile(sigDir6h, 'ACh_group_6h');

if ~exist(outDir6h,'dir'); mkdir(outDir6h); end

%% 1) Run batch on 6h cropped files
GROUP_6h_all = run_ach_batch_auto(sigDir6h, scoreDir6h);

%% 2) Keep only the three conditions of interest
GROUP_6h = subset_group_by_cond(GROUP_6h_all, ...
                                {'baseline','ambtemp','drugs'});

save(fullfile(outDir6h, 'ACh_GROUP_6h_baseline_ambtemp_drugs.mat'), ...
     'GROUP_6h','-v7.3');

%% 3) Transitions: Wake / NREM / REM onset
%   - rows = baseline, ambtemp, drugs
%   - left = mean ± SEM traces
%   - right = peak bars + t-tests
group_plot_all_transitions(GROUP_6h, outDir6h);

%% 4) Group metrics (two-way ANOVA cond x genotype)
ach_stats_anova(GROUP_6h, outDir6h, '6h_baseline_ambtemp_drugs');

%% 5) NREM ACh PSD + "peak oscillation" metrics, per condition
conds6h = {'baseline','ambtemp','drugs'};
for i = 1:numel(conds6h)
    c = conds6h{i};
    Gc = subset_group_by_cond(GROUP_6h, c);
    out_c = fullfile(outDir6h, ['PSD_' c]);
    ach_plot_nrem_psd_baseline(Gc, out_c);
end

fprintf('\n[6h analysis] Done. Figures saved to %s\n', outDir6h);
