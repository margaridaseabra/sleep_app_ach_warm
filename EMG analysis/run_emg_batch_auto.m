function EMG_GROUP = run_emg_batch_auto(sigDir, scoreDir)
% RUN_EMG_BATCH_AUTO
% -------------------------------------------------------------------------
% Batch EMG burst analysis for all mice / conditions.
%
% Uses same filename conventions as run_ach_batch_auto:
%   Signals (.mat): YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
%   Scores  (.csv): YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
%
% OUTPUT:
%   EMG_GROUP.sessions    : metadata per session (mouse, geno, cond, files)
%   EMG_GROUP.out         : cell array of per-session OUT structs
%   EMG_GROUP.metrics_tbl : table with per-session state metrics suitable
%                           for group plots (bursts/min & peak amp per state)

if nargin < 1 || isempty(sigDir)
    sigDir = uigetdir(pwd,'Select folder with .mat signal files');
    if sigDir == 0, error('No signal folder selected.'); end
end
if nargin < 2 || isempty(scoreDir)
    scoreDir = uigetdir(pwd,'Select folder with 1-Hz score CSVs');
    if scoreDir == 0, error('No score folder selected.'); end
end

fprintf('\nSignals dir: %s\nScores dir : %s\n', sigDir, scoreDir);

% -------------------------------------------------------------------------
% 1) Parse score filenames
% -------------------------------------------------------------------------
scoreFiles = dir(fullfile(scoreDir, '*_scores_1Hz.csv'));
scoreMeta  = struct([]);

for k = 1:numel(scoreFiles)
    fn = scoreFiles(k).name;
    info = parse_score_filename_emg(fn);
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

% -------------------------------------------------------------------------
% 2) Parse signal filenames & match to scores
% -------------------------------------------------------------------------
sigFiles = dir(fullfile(sigDir, '*.mat'));
SESS = struct([]);

for k = 1:numel(sigFiles)
    fn = sigFiles(k).name;
    sInfo = parse_signal_filename_emg(fn);
    if isempty(sInfo)
        fprintf('  [WARN] Could not parse signal filename: %s\n', fn);
        continue;
    end

    condKey  = lower(sInfo.cond);
    mouseKey = sInfo.mouse;

    idx = find(strcmpi({scoreMeta.cond}, condKey) & ...
               strcmp({scoreMeta.mouse}, mouseKey), 1, 'first');

    if isempty(idx)
        fprintf('  [WARN] No score CSV for %s (mouse %s, cond %s)\n', ...
                fn, mouseKey, condKey);
        continue;
    end

    sm = scoreMeta(idx);

    sess.mouse = sprintf('mouse%s', mouseKey);
    sess.geno  = sm.geno;
    sess.cond  = condKey;
    sess.mat   = fullfile(sigDir, fn);
    sess.csv   = sm.file;

    if isempty(SESS)
        SESS = sess;
    else
        SESS(end+1) = sess; %#ok<AGROW>
    end
end

if isempty(SESS)
    error('No matched EMG signal+score pairs found.');
end

fprintf('\nMatched %d sessions for EMG burst analysis.\n', numel(SESS));

% -------------------------------------------------------------------------
% 3) Run emg_burst_analysis on each session and collect metrics
% -------------------------------------------------------------------------
CODES = struct('WK',0,'NREM',1,'REM',2,'MA',15);

nSess = numel(SESS);
EMG_GROUP.sessions = SESS;
EMG_GROUP.out      = cell(nSess,1);
EMG_GROUP.metrics  = struct([]);

for k = 1:nSess
    s = SESS(k);
    fprintf('\n=== EMG session %d/%d: %s | %s | %s ===\n', ...
            k, nSess, s.mouse, s.geno, s.cond);

    out_prefix = sprintf('%s_%s_%s', s.mouse, s.geno, s.cond);

    % condition-specific notch: ambtemp has strong 50 Hz + harmonics
    if strcmpi(s.cond, 'ambtemp')
        notch_freqs = [50 100 150];  % from PSD inspection
    else
        notch_freqs = [];
    end

    OUT = emg_burst_analysis( ...
        s.mat, s.csv, ...
        'codes', CODES, ...
        'mouse_id', s.mouse, ...
        'session', s.cond, ...
        'out_prefix', out_prefix, ...
        'notch_freqs', notch_freqs, ...
        'notch_Q', 30, ...
        'verbose', true);


    EMG_GROUP.out{k} = OUT;

    % per-session state metrics
    ss = OUT.state_summary;
    % helper to fetch by state label
    fm = @(lab, field) ss.(field)(strcmp(ss.state, lab));

    EMG_GROUP.metrics(k).mouse = s.mouse;
    EMG_GROUP.metrics(k).geno  = s.geno;
    EMG_GROUP.metrics(k).cond  = s.cond;

    for st = ["Wake","NREM","REM"]
        st = char(st);
        suf = ['_' st];
        EMG_GROUP.metrics(k).(['bursts_per_min' suf]) = fm(st,'bursts_per_min_mean');
        EMG_GROUP.metrics(k).(['peak_amp' suf])       = fm(st,'peak_amp_mean');
    end
end

EMG_GROUP.metrics_tbl = struct2table(EMG_GROUP.metrics);

csv_out = fullfile(sigDir, 'EMG_group_metrics.csv');
mat_out = fullfile(sigDir, 'EMG_GROUP_all.mat');
writetable(EMG_GROUP.metrics_tbl, csv_out);
save(mat_out, 'EMG_GROUP');

fprintf('\nSaved:\n  %s\n  %s\n', csv_out, mat_out);
fprintf('Done EMG batch.\n\n');

end

% ---------------- filename parsers (same logic as ACh) ------------------
function info = parse_signal_filename_emg(fname)
% Parse: YYYYMMDD-cond-mouseNN_notched50Hz_bw3.mat
pat = '^(?<date>\d{8})-(?<cond>[^-]+)-mouse(?<mouse>\d+).*\.mat$';
m = regexp(fname, pat, 'names');
if isempty(m)
    info = [];
else
    info = struct();
    info.date  = m.date;
    info.cond  = lower(m.cond);
    info.mouse = m.mouse;
end
end

function info = parse_score_filename_emg(fname)
% Parse: YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
pat = ['^(?<date>\d{8})-(?<cond>[^-]+)-mouse(?<mouse>\d+)-' ...
       '(?<geno>[^-_]+)_scored_scores_1Hz\.csv$'];
m = regexp(fname, pat, 'names');
if isempty(m)
    info = [];
else
    info = struct();
    info.date  = m.date;
    info.cond  = lower(m.cond);
    info.mouse = m.mouse;
    info.geno  = m.geno;
end
end
