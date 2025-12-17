function GROUP = run_ach_batch_auto(sigDir, scoreDir)
% run_ach_batch_auto
% -------------------------------------------------------------------------
% Automatically:
%   - finds all ACh .mat signal files in sigDir
%   - finds all 1-Hz score CSVs in scoreDir
%   - matches them by mouse ID and condition
%   - runs ach_analysis() for each pair
%   - collects per-session metrics into GROUP and ACh_group_metrics.csv
%
% Expected file name patterns:
%   Signals (.mat):
%       YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
%       e.g. 20251026-drugs-mouse12_notched50Hz_bw3.mat
%
%   Scores (.csv):
%       YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
%       e.g. 20251012-baseline-mouse12-WT_scored_scores_1Hz.csv
%
%   where:
%       cond  = baseline | ambtemp | drugs (etc., no dashes inside)
%       GENO  = WT | APPPS1 (or whatever label you use)
%
% Usage:
%   GROUP = run_ach_batch_auto('/path/to/mat_folder', ...
%                              '/path/to/csv_folder');
%
% Requires: ach_analysis.m on path.

    if nargin < 1 || isempty(sigDir)
        sigDir = uigetdir(pwd,'Select folder with ACh .mat files');
        if sigDir == 0, error('No signal folder selected.'); end
    end
    if nargin < 2 || isempty(scoreDir)
        scoreDir = uigetdir(pwd,'Select folder with 1-Hz score CSVs');
        if scoreDir == 0, error('No score folder selected.'); end
    end

    fprintf('\nSignals dir: %s\nScores dir : %s\n', sigDir, scoreDir);

    % ---------------------------------------------------------------------
    % 1) Find and parse all score files
    % ---------------------------------------------------------------------
  
    scoreFiles = dir(fullfile(scoreDir, '*_scores_1Hz*.csv'));

    scoreMeta  = struct([]);     % let first entry define the fields

    for k = 1:numel(scoreFiles)
        fn = scoreFiles(k).name;
        info = parse_score_filename(fn);
        if isempty(info)
            fprintf('  [WARN] Could not parse score filename: %s\n', fn);
            continue;
        end
        info.file = fullfile(scoreDir, fn);

        if isempty(scoreMeta)
            scoreMeta = info;
        else
            scoreMeta(end+1) = info; 
        end
    end

    if isempty(scoreMeta)
        error('No valid score CSVs found in %s', scoreDir);
    end

    % ---------------------------------------------------------------------
    % 2) Find and parse all signal (.mat) files, match to scores
    % ---------------------------------------------------------------------
    sigFiles = dir(fullfile(sigDir, '*.mat'));
    SESS = struct([]);           % again, let first entry define fields

    for k = 1:numel(sigFiles)
        fn = sigFiles(k).name;
        sInfo = parse_signal_filename(fn);
        if isempty(sInfo)
            fprintf('  [WARN] Could not parse signal filename: %s\n', fn);
            continue;
        end

        % Build key for matching: lower(cond) & mouse ID
        condKey  = lower(sInfo.cond);
        mouseKey = sInfo.mouse;           % number as string

        % find first score file with same cond + mouse
        idx = find(strcmpi({scoreMeta.cond}, condKey) & ...
                   strcmp({scoreMeta.mouse}, mouseKey), 1, 'first');

        if isempty(idx)
            fprintf('  [WARN] No score CSV for %s (mouse %s, cond %s)\n', ...
                    fn, mouseKey, condKey);
            continue;
        end

        sm = scoreMeta(idx);

        sess.mouse = sprintf('mouse%s', mouseKey); % e.g. 'mouse12'
        sess.geno  = sm.geno;                      % WT / APPPS1
        sess.cond  = condKey;                      % baseline / ambtemp / drugs
        sess.mat   = fullfile(sigDir, fn);
        sess.csv   = sm.file;

        if isempty(SESS)
            SESS = sess;
        else
            SESS(end+1) = sess; %#ok<AGROW>
        end
    end

    if isempty(SESS)
        error('No matched signal+score pairs found.');
    end

    fprintf('\nMatched %d sessions (signal + scores).\n', numel(SESS));

    % ---------------------------------------------------------------------
    % 3) Run ach_analysis for each session and collect metrics
    % ---------------------------------------------------------------------
    CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);   % adjust if needed

    nSess = numel(SESS);
    GROUP.sessions = SESS;
    GROUP.out      = cell(nSess,1);
    GROUP.metrics  = struct([]);

    for k = 1:nSess
        s = SESS(k);
        fprintf('\n=== Processing %s | %s | %s ===\n', ...
                s.mouse, s.geno, s.cond);

        out_prefix = sprintf('%s_%s_%s', s.mouse, s.geno, s.cond);

        OUT = ach_analysis_new( ...
            s.mat, s.csv, ...
            'codes', CODES, ...
            'mouse_id', s.mouse, ...
            'session', s.cond, ...
            'out_prefix', out_prefix, ...
            'verbose', true);   % set false if you prefer silence
             close all;
        GROUP.out{k} = OUT;

        fprintf('\n=== DEBUG: Session %d ===\n', k);
        if isfield(OUT,'psd') && isfield(OUT.psd,'Wake')
            fprintf('Wake PSD fields:\n');
            disp(fieldnames(OUT.psd.Wake));
        else
            fprintf('No Wake PSD found\n');
        end

        % ---------- Light per-session summary for group stats ----------
        GROUP.metrics(k).mouse = s.mouse;
        GROUP.metrics(k).geno  = s.geno;
        GROUP.metrics(k).cond  = s.cond;

        % NREM PSD metrics (if present)
        if isfield(OUT,'psd') && isfield(OUT.psd,'NREM') && ~isempty(OUT.psd.NREM.f)
            GROUP.metrics(k).NREM_power   = OUT.psd.NREM.band_power;
            GROUP.metrics(k).NREM_peakHz  = OUT.psd.NREM.peak_freq;
            GROUP.metrics(k).NREM_peakAmp = OUT.psd.NREM.peak_amp;
        else
            GROUP.metrics(k).NREM_power   = NaN;
            GROUP.metrics(k).NREM_peakHz  = NaN;
            GROUP.metrics(k).NREM_peakAmp = NaN;
        end

        % Wake PSD metrics (if present)
        if isfield(OUT,'psd') && isfield(OUT.psd,'Wake') && ~isempty(OUT.psd.Wake.f)
            GROUP.metrics(k).Wake_power   = OUT.psd.Wake.band_power;
            GROUP.metrics(k).Wake_peakHz  = OUT.psd.Wake.peak_freq;
            GROUP.metrics(k).Wake_peakAmp = OUT.psd.Wake.peak_amp;
            
            % Cycle frequency
            if isfield(OUT.psd.Wake, 'cycle_freq')
                GROUP.metrics(k).Wake_cycleHz = OUT.psd.Wake.cycle_freq;
            else
                GROUP.metrics(k).Wake_cycleHz = NaN;
            end
        else
            GROUP.metrics(k).Wake_power   = NaN;
            GROUP.metrics(k).Wake_peakHz  = NaN;
            GROUP.metrics(k).Wake_peakAmp = NaN;
            GROUP.metrics(k).Wake_cycleHz = NaN;
        end
        
        % NREM cycle frequency
        if isfield(OUT,'psd') && isfield(OUT.psd,'NREM') && ~isempty(OUT.psd.NREM.f)
            if isfield(OUT.psd.NREM, 'cycle_freq')
                GROUP.metrics(k).NREM_cycleHz = OUT.psd.NREM.cycle_freq;
            else
                GROUP.metrics(k).NREM_cycleHz = NaN;
            end
        else
            GROUP.metrics(k).NREM_cycleHz = NaN;
        end

        % Compute Wake cycle frequency from raw signal
        GROUP.metrics(k).Wake_cycleHz = compute_cycle_freq(OUT, 'Wake');
        
        % Similarly for NREM if you want
        GROUP.metrics(k).NREM_cycleHz = compute_cycle_freq(OUT, 'NREM');

        % State-wise ACh slopes (if you implemented OUT.state_slopes)
        if isfield(OUT,'state_slopes')
            GROUP.metrics(k).slope_Wake = get_state_slope(OUT,'Wake');
            GROUP.metrics(k).slope_NREM = get_state_slope(OUT,'NREM');
            GROUP.metrics(k).slope_REM  = get_state_slope(OUT,'REM');
        else
            GROUP.metrics(k).slope_Wake = NaN;
            GROUP.metrics(k).slope_NREM = NaN;
            GROUP.metrics(k).slope_REM  = NaN;
        end

        % % ---------- Transition metrics (Wake, NREM, REM onsets) ----------
        % GROUP.metrics(k) = add_transition_metric(GROUP.metrics(k), OUT, ...
        %                                          'Wake_onset','WakeOn');
        % GROUP.metrics(k) = add_transition_metric(GROUP.metrics(k), OUT, ...
        %                                          'NREM_onset','NREMOn');
        % GROUP.metrics(k) = add_transition_metric(GROUP.metrics(k), OUT, ...
        %                                          'REM_onset','REMOn');
                % ---------- Transition metrics (Wake, NREM, REM onsets) ----------
        % ---- WAKE onset ----
        idxW = [];
        if isfield(OUT,'transitions')
            idxW = find(strcmp({OUT.transitions.name}, 'Wake_onset'));
        end
        if ~isempty(idxW)
            Tw = OUT.transitions(idxW);
            if isfield(Tw,'peaks') && ~isempty(Tw.peaks)
                GROUP.metrics(k).WakeOn_peak_mean = mean(Tw.peaks, 'omitnan');
            else
                GROUP.metrics(k).WakeOn_peak_mean = NaN;
            end
            if isfield(Tw,'slopes') && ~isempty(Tw.slopes)
                GROUP.metrics(k).WakeOn_slope_mean = mean(Tw.slopes, 'omitnan');
            else
                GROUP.metrics(k).WakeOn_slope_mean = NaN;
            end
        else
            GROUP.metrics(k).WakeOn_peak_mean  = NaN;
            GROUP.metrics(k).WakeOn_slope_mean = NaN;
        end

        % ---- NREM onset ----
        idxN = [];
        if isfield(OUT,'transitions')
            idxN = find(strcmp({OUT.transitions.name}, 'NREM_onset'));
        end
        if ~isempty(idxN)
            Tn = OUT.transitions(idxN);
            if isfield(Tn,'peaks') && ~isempty(Tn.peaks)
                GROUP.metrics(k).NREMOn_peak_mean = mean(Tn.peaks, 'omitnan');
            else
                GROUP.metrics(k).NREMOn_peak_mean = NaN;
            end
            if isfield(Tn,'slopes') && ~isempty(Tn.slopes)
                GROUP.metrics(k).NREMOn_slope_mean = mean(Tn.slopes, 'omitnan');
            else
                GROUP.metrics(k).NREMOn_slope_mean = NaN;
            end
        else
            GROUP.metrics(k).NREMOn_peak_mean  = NaN;
            GROUP.metrics(k).NREMOn_slope_mean = NaN;
        end

        % ---- REM onset ----
        idxR = [];
        if isfield(OUT,'transitions')
            idxR = find(strcmp({OUT.transitions.name}, 'REM_onset'));
        end
        if ~isempty(idxR)
            Tr = OUT.transitions(idxR);
            if isfield(Tr,'peaks') && ~isempty(Tr.peaks)
                GROUP.metrics(k).REMOn_peak_mean = mean(Tr.peaks, 'omitnan');
            else
                GROUP.metrics(k).REMOn_peak_mean = NaN;
            end
            if isfield(Tr,'slopes') && ~isempty(Tr.slopes)
                GROUP.metrics(k).REMOn_slope_mean = mean(Tr.slopes, 'omitnan');
            else
                GROUP.metrics(k).REMOn_slope_mean = NaN;
            end
        else
            GROUP.metrics(k).REMOn_peak_mean  = NaN;
            GROUP.metrics(k).REMOn_slope_mean = NaN;
        end

    end
    % ---------------------------------------------------------------------
    % 4) Save group metrics table to CSV + .mat
    % ---------------------------------------------------------------------
    GROUP.metrics_tbl = struct2table(GROUP.metrics);
    csv_out = fullfile(sigDir, 'ACh_group_metrics.csv');
    mat_out = fullfile(sigDir, 'ACh_GROUP_all.mat');

    writetable(GROUP.metrics_tbl, csv_out);
    save(mat_out, 'GROUP');

    fprintf('\nSaved:\n  %s\n  %s\n', csv_out, mat_out);
    fprintf('Done.\n\n');
end

% =======================================================================
% --------------------- Helper functions --------------------------------
% =======================================================================

function info = parse_signal_filename(fname)
% Parse signal MAT names in either of these forms:
%   1) YYYYMMDD-cond-mouseNN_something.mat
%   2) YYYYMMDD_cond_mouseNN_something.mat
%
% Examples:
%   20251001-baseline-mouse1_notched50Hz_bw3.mat
%   20251001_baseline_mouse1_APP.mat
%   20251001_baseline_mouse1_APP_REM_QC.mat
%   20251022_drugs_mouse2_WT.mat

    % Allow both '-' and '_' as separators
    pat = '^(?<date>\d{8})[-_](?<cond>[^-_]+)[-_]mouse(?<mouse>\d+).*\.mat$';

    m = regexp(fname, pat, 'names');
    if isempty(m)
        info = [];
    else
        info = struct();
        info.date  = m.date;
        info.cond  = lower(m.cond);   % use lowercase for matching
        info.mouse = m.mouse;         % just the number as string
    end
end

function cycleHz = compute_cycle_freq(OUT, state)
% Compute ACh cycle frequency from peak-to-peak intervals in raw signal
% 
% INPUTS:
%   OUT   - output structure from ach_analysis
%   state - 'Wake', 'NREM', or 'REM'
%
% OUTPUT:
%   cycleHz - average cycle frequency in Hz (NaN if cannot compute)

    cycleHz = NaN;
    
    % Check if we have the state-segmented signal
    if ~isfield(OUT, 'sig_state') || ~isfield(OUT.sig_state, state)
        return;
    end
    
    sig = OUT.sig_state.(state);
    if isempty(sig) || numel(sig) < 10
        return;
    end
    
    % Get sampling frequency
    if isfield(OUT, 'fs')
        fs = OUT.fs;
    elseif isfield(OUT, 'Fs')
        fs = OUT.Fs;
    else
        warning('No sampling frequency found, assuming 1000 Hz');
        fs = 1000;
    end
    
    % Detrend and normalize the signal
    sig = detrend(sig(:));
    sig = sig / std(sig);
    
    % Find peaks with reasonable constraints
    % MinPeakHeight: at least 0.5 std above mean
    % MinPeakDistance: at least 0.2s between peaks (max 5 Hz)
    % MinPeakProminence: peak must stand out by at least 0.3 std
    [pks, locs] = findpeaks(sig, ...
        'MinPeakHeight', 0.5, ...
        'MinPeakDistance', round(fs * 0.2), ...
        'MinPeakProminence', 0.3);
    
    % Need at least 3 peaks to compute reliable frequency
    if numel(locs) < 3
        return;
    end
    
    % Compute inter-peak intervals in seconds
    interpeak_intervals = diff(locs) / fs;
    
    % Remove outliers (intervals more than 3 std from median)
    med_interval = median(interpeak_intervals);
    std_interval = std(interpeak_intervals);
    valid_idx = abs(interpeak_intervals - med_interval) < 3 * std_interval;
    
    if sum(valid_idx) < 2
        return;
    end
    
    % Average period and convert to frequency
    avg_period = mean(interpeak_intervals(valid_idx));
    cycleHz = 1 / avg_period;
    
    % Sanity check: typical ACh cycles are 0.01-5 Hz
    if cycleHz < 0.01 || cycleHz > 5
        cycleHz = NaN;
    end
end

