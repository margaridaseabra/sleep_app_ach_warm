function [ALL_REM, ALL_SUMMARY] = run_REM_MA_pipeline_drugs_6h( ...
                                        data_dir, scores_dir, out_dir)
% run_REM_MA_pipeline_drugs_6h
% -------------------------------------------------------------------------
% Like run_REM_MA_pipeline_ambtemp_3h, but:
%   - only BASELINE + drugs
%   - uses CROPPED scores files: *_scored_scores_1Hz_crop.csv
%   - REM metrics are thus restricted to the 6 h window
%
% INPUTS
%   data_dir   : folder with EEG .mat:
%                  YYYYMMDD_condition_mouseX_GENO.mat
%   scores_dir : folder with CROPPED CSVs:
%                  YYYYMMDD_condition_mouseX_GENO_scored_scores_1Hz_crop.csv
%   out_dir    : where to save per-file intermediate results
%
% OUTPUTS
%   ALL_REM      : all REM bouts across mice/conds (6 h windows only)
%   ALL_SUMMARY  : per-recording REM summary (6 h windows only)
% -------------------------------------------------------------------------

if nargin < 3 || isempty(out_dir)
    % different folder name so ambtemp 3h and drugs 6h never collide
    out_dir = fullfile(scores_dir,'REM_drugs_6h_out');
end
if ~isfolder(out_dir), mkdir(out_dir); end

files = dir(fullfile(data_dir,'*.mat'));

ALL_REM      = table();
ALL_SUMMARY  = table();

for iFile = 1:numel(files)
    eeg_fname = files(iFile).name;
    [~, eegbase, ~] = fileparts(eeg_fname);

    % Expect: YYYYMMDD_condition_mouseX_GENO.mat
    tok = regexp(eegbase, ...
        '^(?<date>\d{8})_(?<cond>[^_]+)_(?<mouse>mouse\d+)_(?<geno>APP|WT)$', ...
        'names');

    if isempty(tok)
        fprintf('Skipping %s (unexpected EEG filename pattern).\n', eeg_fname);
        continue;
    end

    cond = lower(string(tok.cond));

    % --------- HERE we keep BASELINE + drugs instead of ambtemp ----------
    if ~(cond=="baseline" || cond=="drugs")
        % only baseline + drugs in this pipeline
        continue;
    end

    mouse_id = string(tok.mouse);
    geno     = string(tok.geno);

    % Build expected CROPPED scores filename
    % YYYYMMDD_condition_mouseX_GENO_scored_scores_1Hz_crop.csv
    scores_base = sprintf('%s_%s_%s_%s_scored_scores_1Hz_crop', ...
                          tok.date, tok.cond, tok.mouse, tok.geno);
    scores_file = fullfile(scores_dir, [scores_base '.csv']);

    if ~isfile(scores_file)
        fprintf('No cropped scores for %s (%s), skipping.\n', eeg_fname, scores_file);
        continue;
    end

    fprintf('Processing %s with %s (6 h window)...\n', eeg_fname, scores_file);

    % ---------- 1) Load EEG from .mat ----------
    S = load(fullfile(data_dir, eeg_fname));
    eeg_notched = S.eeg;              % adapt if needed
    fs          = S.eeg_frequency;    % adapt if needed

    % ---------- 2) Load states / MA from CROPPED scores ----------
    % If your helper is generic, you can keep using read_scores_csv_ambtemp.
    % If you prefer, make a read_scores_csv_drugs with the same interface.
    [states, epochs_t, MA_tbl, mouse_label, geno_label, cond_label] = ...
        read_scores_csv_ambtemp(scores_file);

    mouse_label = string(mouse_label);
    geno_label  = string(geno_label);
    cond_label  = string(cond_label);

    % ---------- 3) Epoch spectra ----------
    EP = compute_epoch_spectra(eeg_notched, fs, states, epochs_t, ...
                               'mouse', mouse_label, ...
                               'geno',  geno_label, ...
                               'cond',  cond_label);
    save(fullfile(out_dir, [eegbase '_EP_6h.mat']), 'EP');

    % ---------- 4) Bout table from scores ----------
    BOUTS = build_bout_table(states, epochs_t, mouse_label, geno_label, cond_label);
    save(fullfile(out_dir, [eegbase '_BOUTS_6h.mat']), 'BOUTS');

    % ---------- 5) Annotate REM with MAs ----------
    REM = annotate_REM_with_MA(BOUTS, MA_tbl, ...
                               'pre_window_s', 60, ...
                               'long_quantile', 0.75);
    save(fullfile(out_dir, [eegbase '_REM_6h.mat']), 'REM');

    % ---------- 6) (Optional) Cluster REM bouts ----------
    CLUST = cluster_REM_bouts(REM, 'max_gap_s', 60);
    save(fullfile(out_dir, [eegbase '_CLUST_6h.mat']), 'CLUST');

    % ---------- 7) Per-recording summary (6 h window) ----------
    SUMMARY = table;
    SUMMARY.mouse = mouse_label;
    SUMMARY.geno  = geno_label;
    SUMMARY.cond  = cond_label;

    SUMMARY.n_REM_total = height(REM);
    SUMMARY.n_REM_long  = sum(REM.is_long);
    SUMMARY.propREM_long = SUMMARY.n_REM_long / max(SUMMARY.n_REM_total,1);

    SUMMARY.total_REM_dur_s  = sum(REM.dur_s);
    SUMMARY.total_long_dur_s = sum(REM.dur_s(REM.is_long));
    SUMMARY.propREM_time_long = SUMMARY.total_long_dur_s / ...
                                max(SUMMARY.total_REM_dur_s,1);

    SUMMARY.mean_n_MA_pre_normal = mean(REM.n_MA_pre(~REM.is_long), 'omitnan');
    SUMMARY.mean_n_MA_pre_long   = mean(REM.n_MA_pre( REM.is_long), 'omitnan');

    SUMMARY.mean_REM_bout_len_s = SUMMARY.total_REM_dur_s / ...
                                  max(SUMMARY.n_REM_total, 1);

    SUMMARY.REM_frag_bouts_per_min_REM = ...
        (SUMMARY.n_REM_total / max(SUMMARY.total_REM_dur_s, 1)) * 60;

    if ~isempty(CLUST)
        nC = height(CLUST);
        SUMMARY.prop_cluster_shortonly   = sum(CLUST.cluster_type=="ShortOnly")  / nC;
        SUMMARY.prop_cluster_shortlong   = sum(CLUST.cluster_type=="ShortLong")  / nC;
        SUMMARY.prop_cluster_longonly    = sum(CLUST.cluster_type=="LongOnly")   / nC;
        SUMMARY.mean_n_short_before_long = mean( ...
            CLUST.n_short_before_long(CLUST.cluster_type=="ShortLong"), ...
            'omitnan');
    else
        SUMMARY.prop_cluster_shortonly   = NaN;
        SUMMARY.prop_cluster_shortlong   = NaN;
        SUMMARY.prop_cluster_longonly    = NaN;
        SUMMARY.mean_n_short_before_long = NaN;
    end

    ALL_REM     = [ALL_REM; REM];
    ALL_SUMMARY = [ALL_SUMMARY; SUMMARY];
end

% ---------- Add REM size category (short/med/long) ----------
[ALL_REM, rem_size_edges_s] = add_REM_size_category(ALL_REM);
save(fullfile(out_dir, 'REM_size_edges_6h.mat'), 'rem_size_edges_s');

fprintf('\n✓ REM+MA pipeline (6 h baseline + drugs) finished.\n');
fprintf('   n REM bouts: %d,  n recordings: %d\n', height(ALL_REM), height(ALL_SUMMARY));
end
