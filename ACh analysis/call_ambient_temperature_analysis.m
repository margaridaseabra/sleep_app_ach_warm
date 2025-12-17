%% Temperature Comparison Analysis: Baseline vs Ambtemp (WT vs APP)
% This script compares baseline and ambient temperature conditions
% across different time windows: 0-3h, 3-6h, and washout

% --- EDIT THESE PATHS ---
sigDir   = '/Users/margaridaseabra/24.11 signalnotscored';
scoreDir = '/Users/margaridaseabra/24.11scores';
outDir   = fullfile(sigDir, 'ACh_temperature_comparison');
% -------------------------

if ~exist(outDir, 'dir'); mkdir(outDir); end

%% 1) Run the batch analysis (all conditions present in those folders)
% This will process all sessions including baseline and ambtemp
fprintf('\n========================================\n');
fprintf('Step 1: Running batch analysis...\n');
fprintf('========================================\n');
GROUP_all = run_ach_batch_auto(sigDir, scoreDir);

%% 2) Keep only 'baseline' and 'ambtemp' sessions
fprintf('\n========================================\n');
fprintf('Step 2: Filtering for baseline and ambtemp conditions...\n');
fprintf('========================================\n');
GROUP_temp = subset_group_by_cond(GROUP_all, {'baseline', 'ambtemp'});

fprintf('Sessions after filtering: %d\n', numel(GROUP_temp.sessions));
fprintf('  Baseline sessions: %d\n', sum(string({GROUP_temp.sessions.cond}) == "baseline"));
fprintf('  Ambtemp sessions: %d\n', sum(string({GROUP_temp.sessions.cond}) == "ambtemp"));

% 3) Save for safety
save(fullfile(outDir, 'ACh_GROUP_temperature_all.mat'), ...
     'GROUP_temp', '-v7.3');
fprintf('Saved: ACh_GROUP_temperature_all.mat\n');

%% 4) Overall temperature comparison (full recording)
% This uses all data from each session

fprintf('\n========================================\n');
fprintf('Step 3: Running full recording analysis...\n');
fprintf('========================================\n');

% 4a) Summary statistics with 2-way ANOVA
fprintf('  - Statistics...\n');
ach_stats_temperature_anova(GROUP_temp, outDir, 'Full_Recording', 'Full Recording');

% 4b) State transition plots (Wake/NREM/REM onsets)
fprintf('  - Transitions...\n');
group_plot_temperature_transitions(GROUP_temp, outDir, true, 'Full Recording');

% 4c) NREM ACh PSD comparison
fprintf('  - NREM PSD...\n');
ach_plot_nrem_psd_temperature(GROUP_temp, outDir, true, 'Full Recording');

% 4d) Wake ACh PSD comparison
fprintf('  - Wake PSD...\n');
ach_plot_wake_psd_temperature2(GROUP_temp, outDir, true, 'Full Recording');

fprintf('Full recording analysis complete!\n');

%% 5) Time-windowed analysis: 0-3h, 3-6h, washout
% Now we'll segment each recording into time windows and analyze separately

fprintf('\n========================================\n');
fprintf('Step 4: Time-windowed analysis...\n');
fprintf('========================================\n');

time_windows = struct([]);
time_windows(1).name = '0-3h';
time_windows(1).start_sec = 0;
time_windows(1).end_sec = 3 * 3600;  % 3 hours in seconds

time_windows(2).name = '3-6h';
time_windows(2).start_sec = 3 * 3600;
time_windows(2).end_sec = 6 * 3600;

time_windows(3).name = 'washout';
time_windows(3).start_sec = 6 * 3600;
time_windows(3).end_sec = Inf;  % rest of recording

% Create time-windowed groups
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
    
    % Create subdirectory for this window
    win_dir = fullfile(outDir, win.name);
    if ~exist(win_dir, 'dir')
        mkdir(win_dir);
        fprintf('Created directory: %s\n', win_dir);
    else
        fprintf('Using existing directory: %s\n', win_dir);
    end
    
    % Segment GROUP by time window
    fprintf('Segmenting data by time window...\n');
    try
        GROUP_win = segment_group_by_time(GROUP_temp, win.start_sec, win.end_sec);
        fprintf('Segmentation complete. Sessions: %d\n', numel(GROUP_win.sessions));
    catch ME
        fprintf('ERROR during segmentation: %s\n', ME.message);
        fprintf('Skipping this time window.\n');
        continue;
    end
    
    % Save windowed group
    save(fullfile(win_dir, sprintf('ACh_GROUP_%s.mat', win.name)), ...
         'GROUP_win', '-v7.3');
    fprintf('Saved: %s\n', fullfile(win_dir, sprintf('ACh_GROUP_%s.mat', win.name)));
    
    % Run all analyses for this time window
    fprintf('\n--- Running statistics ---\n');
    try
        ach_stats_temperature_anova(GROUP_win, win_dir, win.name, win.name);
        fprintf('Statistics complete.\n');
    catch ME
        fprintf('ERROR in statistics: %s\n', ME.message);
        fprintf('Stack: %s\n', ME.stack(1).name);
    end
    
    fprintf('\n--- Plotting transitions ---\n');
    try
        group_plot_temperature_transitions(GROUP_win, win_dir, true, win.name);
        fprintf('Transitions complete.\n');
    catch ME
        fprintf('ERROR in transitions: %s\n', ME.message);
        fprintf('Stack: %s\n', ME.stack(1).name);
    end
    
    fprintf('\n--- Plotting NREM PSD ---\n');
    try
        ach_plot_nrem_psd_temperature(GROUP_win, win_dir, true, win.name);
        fprintf('NREM PSD complete.\n');
    catch ME
        fprintf('ERROR in NREM PSD: %s\n', ME.message);
        fprintf('Stack: %s\n', ME.stack(1).name);
    end
    
    fprintf('\n--- Plotting Wake PSD ---\n');
    try
        ach_plot_wake_psd_temperature2(GROUP_win, win_dir, true, win.name);
        fprintf('Wake PSD complete.\n');
    catch ME
        fprintf('ERROR in Wake PSD: %s\n', ME.message);
        fprintf('Stack: %s\n', ME.stack(1).name);
    end
    
    fprintf('Time window %s analysis complete!\n', win.name);
end

%% 6) Create summary comparison figure across time windows
fprintf('\n========================================\n');
fprintf('Step 5: Creating summary figure...\n');
fprintf('========================================\n');

try
    create_time_window_summary(GROUP_temp, time_windows, outDir);
    fprintf('Summary figure created successfully.\n');
catch ME
    fprintf('ERROR creating summary: %s\n', ME.message);
    fprintf('Stack: %s\n', ME.stack(1).name);
end

fprintf('\n========================================\n');
fprintf('Temperature comparison analysis complete!\n');
fprintf('Results saved to: %s\n', outDir);
fprintf('========================================\n');
fprintf('\nOutput structure:\n');
fprintf('  - ACh_GROUP_temperature_all.mat (full dataset)\n');
fprintf('  - Full recording plots (ANOVA, transitions, PSDs)\n');
fprintf('  - 0-3h/ (first 3 hours analysis)\n');
fprintf('  - 3-6h/ (second 3 hours analysis)\n');
fprintf('  - washout/ (remaining time analysis)\n');
fprintf('  - Time_Window_Summary.png (cross-window comparison)\n');
fprintf('========================================\n');