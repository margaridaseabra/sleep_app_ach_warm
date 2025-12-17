% run_REM_MA_pipeline.m
% -------------------------------------------------------
% High-level REM + MA analysis + plots across all mice.
% Now: EEG .mat in data_dir, scores in separate scores_dir as .csv

clear; clc;

data_dir   = fullfile('/Users/margaridaseabra/24.11notchedsignal');          % where your .mat EEG lives
scores_dir = fullfile('/Users/margaridaseabra/24.11scores');    % where your scoring CSVs live
out_dir    = fullfile(pwd,'REM_MA_out');    % output folder
if ~isfolder(out_dir); mkdir(out_dir); end

addpath(fullfile(pwd,'functions'));         % helper functions folder

files = dir(fullfile(data_dir,'*.mat'));

ALL_REM      = table();
ALL_CLUSTERS = table();
ALL_SUMMARY  = table();

for iFile = 1:numel(files)
    eeg_fname = files(iFile).name;
    fprintf('Processing %s...\n', eeg_fname);
    
    % ---------- 1) Load EEG (notched) from .mat -------------------------
    S = load(fullfile(data_dir, eeg_fname));
    
    % >>> ADAPT THESE NAMES TO YOUR .MAT VARIABLES <<<
    eeg_notched = S.eeg;
    fs = S.eeg_frequency;
    
   
    
    % ---------- 2) Load scores from CSV (NOT in data_dir) --------------
        % --------- 2) Derive the correct scores CSV name from EEG filename ------
    % EEG filename pattern (your .mat):
    %   YYYYMMDD-cond-mouseID_GENO.mat
    %
    % Example:
    %   20251024-drugs-mouse10_APP.mat
    %
    % Scores CSV pattern:
    %   YYYYMMDD_cond_mouseID_GENO_scored_scores_1Hz.csv
    %
    % Example:
    %   20251024_drugs_mouse10_APP_scored_scores_1Hz.csv

    [~, eegbase, ~] = fileparts(eeg_fname);

    % Use regexp to parse date, condition, mouse and genotype
    tokens = regexp(eegbase, ...
        '^(?<date>\d{8})_(?<cond>[^-]+)_(?<mouse>[^-_]+)_(?<geno>[^_.]+)', ...
        'names');

    if isempty(tokens)
        error('Unexpected EEG filename format: %s', eegbase);
    end

    date_str  = string(tokens.date);   % "20251024"
    cond_str  = string(tokens.cond);   % "drugs"
    mouse_str = string(tokens.mouse);  % "mouse10"
    geno_str  = string(tokens.geno);   % "APP"

    % Now build the scores CSV base:
    %   YYYYMMDD_cond_mouse_GENO_scored_scores_1Hz.csv
    scores_base = sprintf('%s_%s_%s_%s_scored_scores_1Hz', ...
                          date_str, cond_str, mouse_str, geno_str);

    scores_file = fullfile(scores_dir, scores_base + ".csv");

    if ~isfile(scores_file)
        error('Scores CSV not found for %s at %s', eeg_fname, scores_file);
    end


    if ~isfile(scores_file)
        error('Scores CSV not found for %s at %s', eeg_fname, scores_file);
    end

    % This still parses states, times, MA, and (redundantly) metadata if you want
    [states, epochs_t, MA_tbl, mouse_id, geno, cond] = read_scores_csv(scores_file);

    
    % ---------- 3) Epoch spectra (optional but handy) -------------------
    EP = compute_epoch_spectra(eeg_notched, fs, states, epochs_t, ...
                               'mouse', mouse_id, ...
                               'geno',  geno, ...
                               'cond',  cond);
    save(fullfile(out_dir, [eegbase '_EP.mat']), 'EP');
    
    % ---------- 4) Bout table from scores ------------------------------
    BOUTS = build_bout_table(states, epochs_t, mouse_id, geno, cond);
    save(fullfile(out_dir, [eegbase '_BOUTS.mat']), 'BOUTS');
    
    % ---------- 5) Annotate REM with MAs -------------------------------
    REM = annotate_REM_with_MA(BOUTS, MA_tbl, ...
                               'pre_window_s', 60, ...
                               'long_quantile', 0.75);
    save(fullfile(out_dir, [eegbase '_REM.mat']), 'REM');
    
    % ---------- 6) Cluster REM bouts -----------------------------------
    CLUST = cluster_REM_bouts(REM, 'max_gap_s', 60);
    save(fullfile(out_dir, [eegbase '_CLUST.mat']), 'CLUST');
    
    % ---------- 7) Per-recording summary -------------------------------
    % (same block as before – no change needed)
    SUMMARY = table;
    SUMMARY.mouse = string(mouse_id);
    SUMMARY.geno  = string(geno);
    SUMMARY.cond  = string(cond);
    
    SUMMARY.n_REM_total = height(REM);
    SUMMARY.n_REM_long  = sum(REM.is_long);
    SUMMARY.propREM_long = SUMMARY.n_REM_long / max(SUMMARY.n_REM_total,1);
    
    SUMMARY.total_REM_dur_s  = sum(REM.dur_s);
    SUMMARY.total_long_dur_s = sum(REM.dur_s(REM.is_long));
    SUMMARY.propREM_time_long = SUMMARY.total_long_dur_s / ...
                                max(SUMMARY.total_REM_dur_s,1);
    
    SUMMARY.mean_n_MA_pre_normal = mean(REM.n_MA_pre(~REM.is_long), 'omitnan');
    SUMMARY.mean_n_MA_pre_long   = mean(REM.n_MA_pre( REM.is_long), 'omitnan');

    % ---- REM fragmentation metrics ----
    % Average REM bout length (s)
    SUMMARY.mean_REM_bout_len_s = SUMMARY.total_REM_dur_s / ...
                                  max(SUMMARY.n_REM_total, 1);

    % Fragmentation index: number of REM bouts / total REM time (per minute of REM)
    SUMMARY.REM_frag_bouts_per_min_REM = ...
        (SUMMARY.n_REM_total / max(SUMMARY.total_REM_dur_s, 1)) * 60;

    
    if ~isempty(CLUST)
        nC = height(CLUST);
        SUMMARY.prop_cluster_shortonly  = sum(CLUST.cluster_type=="ShortOnly")  / nC;
        SUMMARY.prop_cluster_shortlong  = sum(CLUST.cluster_type=="ShortLong")  / nC;
        SUMMARY.prop_cluster_longonly   = sum(CLUST.cluster_type=="LongOnly")   / nC;
        SUMMARY.mean_n_short_before_long = mean( ...
            CLUST.n_short_before_long(CLUST.cluster_type=="ShortLong"), ...
            'omitnan');
    else
        SUMMARY.prop_cluster_shortonly  = NaN;
        SUMMARY.prop_cluster_shortlong  = NaN;
        SUMMARY.prop_cluster_longonly   = NaN;
        SUMMARY.mean_n_short_before_long = NaN;
    end
    
    ALL_REM      = [ALL_REM; REM];
    ALL_CLUSTERS = [ALL_CLUSTERS; CLUST];
    ALL_SUMMARY  = [ALL_SUMMARY; SUMMARY];
end

% ---------- Add REM size categories (short / medium / long) ----------
[ALL_REM, rem_size_edges_s] = add_REM_size_category(ALL_REM);
save(fullfile(out_dir, 'REM_size_edges.mat'), 'rem_size_edges_s');

% ---------- Save group-level tables ----------
save(fullfile(out_dir, 'ALL_REM.mat'),      'ALL_REM');
save(fullfile(out_dir, 'ALL_CLUSTERS.mat'), 'ALL_CLUSTERS');
save(fullfile(out_dir, 'ALL_SUMMARY.mat'),  'ALL_SUMMARY');

% ---------- Create group-level plots ----------
plot_REM_long_stats(ALL_SUMMARY, out_dir);        % existing
plot_MA_REM_relationships(ALL_REM, out_dir);      % existing
plot_REM_cluster_types(ALL_SUMMARY, out_dir);     % existing

% New: survival curves + histograms of REM durations
plot_REM_survival_curves(ALL_REM, out_dir);
plot_REM_bout_histograms(ALL_REM, out_dir);

% New: basic APP vs WT stats for REM metrics (writes a CSV)
run_REM_group_stats(ALL_SUMMARY, out_dir);

fprintf('\n✓ Analysis complete. Plots + stats saved in: %s\n', out_dir);
