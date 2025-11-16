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
    scoreFiles = dir(fullfile(scoreDir, '*_scores_1Hz.csv'));
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
            scoreMeta(end+1) = info; %#ok<AGROW>
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

        OUT = ach_analysis( ...
            s.mat, s.csv, ...
            'codes', CODES, ...
            'mouse_id', s.mouse, ...
            'session', s.cond, ...
            'out_prefix', out_prefix, ...
            'verbose', true);   % set false if you prefer silence

        GROUP.out{k} = OUT;

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

        % Wake-onset transition peaks / slopes (if present)
        idxW = [];
        if isfield(OUT,'transitions')
            idxW = find(strcmp({OUT.transitions.name},'Wake_onset'));
        end
        if ~isempty(idxW)
            Tw = OUT.transitions(idxW);
            GROUP.metrics(k).WakeOn_peak_mean  = mean(Tw.peaks,  'omitnan');
            if isfield(Tw,'slopes') && ~isempty(Tw.slopes)
                GROUP.metrics(k).WakeOn_slope_mean = mean(Tw.slopes,'omitnan');
            else
                GROUP.metrics(k).WakeOn_slope_mean = NaN;
            end
        else
            GROUP.metrics(k).WakeOn_peak_mean  = NaN;
            GROUP.metrics(k).WakeOn_slope_mean = NaN;
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
% Parse: YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
pat = '^(?<date>\d{8})-(?<cond>[^-]+)-mouse(?<mouse>\d+).*\.mat$';
m = regexp(fname, pat, 'names');
if isempty(m)
    info = [];
else
    info = struct();
    info.date  = m.date;
    info.cond  = lower(m.cond);  % use lowercase for matching
    info.mouse = m.mouse;        % just the number as string
end
end

function info = parse_score_filename(fname)
% Parse: YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
pat = ['^(?<date>\d{8})-(?<cond>[^-]+)-mouse(?<mouse>\d+)-' ...
       '(?<geno>[^-_]+)_scored_scores_1Hz\.csv$'];
m = regexp(fname, pat, 'names');
if isempty(m)
    info = [];
else
    info = struct();
    info.date  = m.date;
    info.cond  = lower(m.cond);      % lowercase for matching
    info.mouse = m.mouse;            % number as string
    info.geno  = m.geno;             % keep genotype label as-is
end
end

function val = get_state_slope(OUT, stateName)
% Gracefully get mean slope for a state from OUT.state_slopes
if isfield(OUT,'state_slopes') && isfield(OUT.state_slopes,stateName)
    val = OUT.state_slopes.(stateName).mean;
else
    val = NaN;
end
end
