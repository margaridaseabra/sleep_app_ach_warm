%% Three Condition Comparison Analysis: Baseline vs Ambtemp vs Drug (WT vs APP)
% This script compares baseline, ambient temperature, and drug conditions
% across different time windows: 0-3h, 3-6h, and washout

% --- EDIT THESE PATHS ---
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';
outDir   = fullfile(sigDir, 'ACh_three_condition_comparison');
% -------------------------

% Ensure output directory exists
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output directory: %s\n', outDir);
end

%% 1) Run the batch analysis (all conditions present in those folders)
fprintf('\n========================================\n');
fprintf('Step 1: Running batch analysis...\n');
fprintf('========================================\n');
GROUP_all = run_ach_batch_auto(sigDir, scoreDir);

%% 2) Keep only 'baseline', 'ambtemp', and 'drugs' sessions
fprintf('\n========================================\n');
fprintf('Step 2: Filtering for baseline, ambtemp, and drugs conditions...\n');
fprintf('========================================\n');
GROUP_three = subset_group_by_cond(GROUP_all, {'baseline', 'ambtemp', 'drugs'});

fprintf('Sessions after filtering: %d\n', numel(GROUP_three.sessions));
fprintf('  Baseline sessions: %d\n', sum(string({GROUP_three.sessions.cond}) == "baseline"));
fprintf('  Ambtemp sessions: %d\n', sum(string({GROUP_three.sessions.cond}) == "ambtemp"));
fprintf('  Drugs sessions: %d\n', sum(string({GROUP_three.sessions.cond}) == "drugs"));

% 3) Save for safety
fprintf('Saving GROUP_three to: %s\n', fullfile(outDir, 'ACh_GROUP_three_conditions.mat'));
save(fullfile(outDir, 'ACh_GROUP_three_conditions.mat'), ...
     'GROUP_three', '-v7.3');
fprintf('Saved: ACh_GROUP_three_conditions.mat\n');

%% 4) Overall three-condition comparison (full recording)
fprintf('\n========================================\n');
fprintf('Step 3: Running full recording analysis...\n');
fprintf('========================================\n');
%%
% 4a) Summary statistics with 2-way ANOVA
fprintf('  - Statistics...\n');
ach_stats_three_conditions(GROUP_three, outDir, 'Full_Recording', 'Full Recording');
%%
% 4b) State transition plots (Wake/NREM/REM onsets)
fprintf('  - Transitions...\n');
group_plot_three_conditions_transitions(GROUP_three, outDir, false, 'Full Recording');
%%
% 4c) NREM ACh PSD comparison
fprintf('  - NREM PSD...\n');
ach_plot_nrem_psd_three_conditions(GROUP_three, outDir, false, 'Full Recording');
%%
% 4d) Wake ACh PSD comparison
fprintf('  - Wake PSD...\n');
ach_plot_wake_psd_three_conditions(GROUP_three, outDir, false, 'Full Recording');

fprintf('Full recording analysis complete!\n');

%% 5) Time-windowed analysis: 0-3h, 3-6h, washout
fprintf('\n========================================\n');
fprintf('Step 4: Time-windowed analysis...\n');
fprintf('========================================\n');

time_windows = struct([]);
time_windows(1).name = '0-3h';
time_windows(1).start_sec = 0;
time_windows(1).end_sec = 3 * 3600;

time_windows(2).name = '3-6h';
time_windows(2).start_sec = 3 * 3600;
time_windows(2).end_sec = 6 * 3600;

time_windows(3).name = 'washout';
time_windows(3).start_sec = 6 * 3600;
time_windows(3).end_sec = Inf;

fprintf('\nStarting time-windowed analysis...\n');
fprintf('Number of time windows: %d\n', numel(time_windows));

for w = 1:numel(time_windows)
    win = time_windows(w);
    fprintf('\n========================================\n');
    if isinf(win.end_sec)
        fprintf('Processing time window %d/%d: %s (%d-end sec)\n', ...
            w, numel(time_windows), win.name, win.start_sec);
    else
        fprintf('Processing time window %d/%d: %s (%d-%d sec)\n', ...
            w, numel(time_windows), win.name, win.start_sec, win.end_sec);
    end
    fprintf('========================================\n');
    
    win_dir = fullfile(outDir, win.name);
    if ~exist(win_dir, 'dir')
        mkdir(win_dir);
        fprintf('Created directory: %s\n', win_dir);
    end
    
    fprintf('Segmenting data by time window...\n');
    try
        GROUP_win = segment_group_by_time(GROUP_three, win.start_sec, win.end_sec);
        fprintf('Segmentation complete. Sessions: %d\n', numel(GROUP_win.sessions));
    catch ME
        fprintf('ERROR during segmentation: %s\n', ME.message);
        continue;
    end
    
    save(fullfile(win_dir, sprintf('ACh_GROUP_%s.mat', win.name)), ...
         'GROUP_win', '-v7.3');
    fprintf('Saved: %s\n', fullfile(win_dir, sprintf('ACh_GROUP_%s.mat', win.name)));
    
   %% % Run all analyses
    fprintf('\n--- Running statistics ---\n');
    try
        ach_stats_three_conditions(GROUP_win, win_dir, win.name, win.name);
        fprintf('Statistics complete.\n');
    catch ME
        fprintf('ERROR in statistics: %s\n', ME.message);
    end
    %%
    fprintf('\n--- Plotting transitions ---\n');
    try
        group_plot_three_conditions_transitions(GROUP_win, win_dir, false, win.name);
        fprintf('Transitions complete.\n');
    catch ME
        fprintf('ERROR in transitions: %s\n', ME.message);
    end
    %%
    fprintf('\n--- Plotting NREM PSD ---\n');
    try
        ach_plot_nrem_psd_three_conditions(GROUP_win, win_dir, false, win.name);
        fprintf('NREM PSD complete.\n');
    catch ME
        fprintf('ERROR in NREM PSD: %s\n', ME.message);
    end
    %%
    fprintf('\n--- Plotting Wake PSD ---\n');
    try
        ach_plot_wake_psd_three_conditions(GROUP_win, win_dir, false, win.name);
        fprintf('Wake PSD complete.\n');
    catch ME
        fprintf('ERROR in Wake PSD: %s\n', ME.message);
    end
    
    fprintf('Time window %s analysis complete!\n', win.name);
end

fprintf('\n========================================\n');
fprintf('Three condition comparison analysis complete!\n');
fprintf('Results saved to: %s\n', outDir);
fprintf('========================================\n');
fprintf('\nOutput structure:\n');
fprintf('  - ACh_GROUP_three_conditions.mat\n');
fprintf('  - Full recording plots (ANOVA, transitions, PSDs)\n');
fprintf('  - 0-3h/ 3-6h/ washout/ folders\n');
fprintf('========================================\n');