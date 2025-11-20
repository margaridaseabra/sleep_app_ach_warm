function SUMMARY = run_batch_export_scores_MAincluded(input_dir, varargin)
% Run export_sleep_scores_with_MA over all .mat files in a folder tree.
%
% SUMMARY = run_batch_export_scores('~/data/scored', ...
%     'out_base', '~/data/exports', ...
%     'pattern', '*.mat', ...               % or '*_scored.mat'
%     'keep_subfolders', true, ...
%     'epochSec', [], ...                   % set 1 or 5 to force
%     'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), ...
%     'ma_thresh_sec', 15, ...
%     'include_files_without_scores', false, ... % if true: warn & skip cleanly
%     'skip_if_exists', true, ...           % don't redo if CSV already exists
%     'dry_run', false, ...
%     'workers', 0, ...                     % 0/1 = serial; >1 = parfor
%     'write_log', true);
%
% Returns a table-like struct SUMMARY and also writes:
%   <out_base>/batch_export_summary.csv
%   <out_base>/batch_export_summary.mat

p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*.mat',@ischar);
addParameter(p,'out_base','',@ischar);
addParameter(p,'keep_subfolders',true,@islogical);
addParameter(p,'epochSec',[],@(x) isempty(x)||isscalar(x));
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);
addParameter(p,'include_files_without_scores',false,@islogical);
addParameter(p,'skip_if_exists',true,@islogical);
addParameter(p,'dry_run',false,@islogical);
addParameter(p,'workers',0,@(x)isscalar(x) && x>=0);
addParameter(p,'write_log',true,@islogical);
parse(p, input_dir, varargin{:});
S = p.Results;

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);
if isempty(S.out_base), S.out_base = fullfile(S.input_dir, 'exports'); end
if ~isfolder(S.out_base), mkdir(S.out_base); end

% Collect files (recursive)
files = dir(fullfile(S.input_dir, '**', S.pattern));
files = files(~[files.isdir]);

if isempty(files)
    warning('No files matched pattern "%s" under %s', S.pattern, S.input_dir);
    SUMMARY = struct('n_files',0,'rows',table());
    return;
end

% Prepare jobs
N = numel(files);
rows(N,1) = struct('src','', 'out_dir','', 'out_prefix','', 'success',false, ...
    'epochSec',NaN,'n_MA_bouts',NaN,'MA_sec',NaN,'WK_sec_before',NaN, ...
    'WK_sec_after',NaN,'pct_wake_reclassified',NaN, ...
    'scores_1Hz_csv','', 'bouts_csv','', 'timetable_mat','', 'log','', ...
    'message','');

for i = 1:N
    f = files(i);
    src_path = fullfile(f.folder, f.name);
    rel      = erase(src_path, [S.input_dir filesep]);
    reldir   = fileparts(rel);
    out_dir  = S.out_base;
    if S.keep_subfolders && ~isempty(reldir)
        out_dir = fullfile(S.out_base, reldir);
        if ~S.dry_run && ~isfolder(out_dir), mkdir(out_dir); end
    end
    [~, base, ~] = fileparts(src_path);
    out_prefix   = base;

    rows(i).src        = src_path;
    rows(i).out_dir    = out_dir;
    rows(i).out_prefix = out_prefix;
end

% Optionally open a parallel pool
usePar = S.workers > 1;
if usePar
    try
        pool = gcp('nocreate');
        if isempty(pool) || pool.NumWorkers ~= S.workers
            parpool('local', S.workers);
        end
    catch
        warning('Could not start parallel pool. Running serially.');
        usePar = false;
    end
end

% Process files
fprintf('Batch exporting %d file(s) from %s -> %s\n', N, S.input_dir, S.out_base);

runner = @(i) process_one(rows(i), S);

if usePar
    parfor i = 1:N
        rows(i) = runner(i); %#ok<PFBNS>
    end
else
    for i = 1:N
        rows(i) = runner(i);
    end
end

% Build summary table
T = struct2table(rows);
summary_csv = fullfile(S.out_base, 'batch_export_summary.csv');
summary_mat = fullfile(S.out_base, 'batch_export_summary.mat');
writetable(T, summary_csv);
save(summary_mat, 'T', 'S');

% Return
SUMMARY = struct('n_files', N, 'rows', T, 'summary_csv', summary_csv, 'summary_mat', summary_mat);
fprintf('✅ Batch done. Summary: %s\n', summary_csv);
end

% ------------------- helpers -------------------

function row = process_one(row, S)
% skip if outputs exist
if S.skip_if_exists && ~S.dry_run
    f1 = fullfile(row.out_dir, sprintf('%s_scores_1Hz.csv', row.out_prefix));
    f2 = fullfile(row.out_dir, sprintf('%s_scores_bouts.csv', row.out_prefix));
    f3 = fullfile(row.out_dir, sprintf('%s_scores_tt.mat', row.out_prefix));
    if isfile(f1) && isfile(f2) && isfile(f3)
        row.success = true;
        row.scores_1Hz_csv = f1; row.bouts_csv = f2; row.timetable_mat = f3;
        row.message = 'Skipped (exists)';
        return;
    end
end

% load just to check presence of sleep_scores (cheap)
try
    info = whos('-file', row.src);
    has_scores = any(strcmp({info.name}, 'sleep_scores'));
catch
    row.message = 'Could not read file header';
    return;
end

if ~has_scores
    if S.include_files_without_scores
        row.message = 'No sleep_scores (skipped)';
        return;
    else
        row.message = 'No sleep_scores (skipped)';
        return;
    end
end

% create out_dir
if ~S.dry_run && ~isfolder(row.out_dir)
    try
        mkdir(row.out_dir);
    catch ME
        row.message = sprintf('Cannot create out_dir: %s', ME.message);
        return;
    end
end

% call exporter
try
    if S.dry_run
        row.message = 'Dry-run OK';
        row.success = true;
        return;
    end

    R = export_sleep_scores_with_MA(row.src, row.out_prefix, ...
            'out_dir', row.out_dir, ...
            'epochSec', S.epochSec, ...
            'codes', S.codes, ...
            'ma_thresh_sec', S.ma_thresh_sec, ...
            'write_log', S.write_log);

    row.success = logical(R.success);
    row.epochSec = R.epochSec;
    row.n_MA_bouts = R.n_MA_bouts;
    row.MA_sec = R.MA_sec;
    row.WK_sec_before = R.WK_sec_before;
    row.WK_sec_after = R.WK_sec_after;
    row.pct_wake_reclassified = R.pct_wake_reclassified;
    if isfield(R,'files')
        flds = {'scores_1Hz_csv','bouts_csv','timetable_mat','log'};
        for k = 1:numel(flds)
            fld = flds{k};
            if isfield(R.files,fld), row.(fld) = R.files.(fld); end
        end
    end
    row.message = ternary(row.success,'OK','Exporter returned success=false');

    % -------- NEW STEP: force WK bouts <= ma_thresh_sec to MA in 1-Hz CSV --------
    % This assumes *_scores_1Hz.csv is a 1 Hz time series.
    % Every contiguous run of WK (code S.codes.WK) with duration <= ma_thresh_sec
    % seconds is re-coded to MA (code S.codes.MA).
    try
        if row.success && ~isempty(row.scores_1Hz_csv) && isfile(row.scores_1Hz_csv)
            recode_short_WK_to_MA(row.scores_1Hz_csv, ...
                                  S.codes.WK, S.codes.MA, S.ma_thresh_sec);
        end
    catch ME2
        % Do not fail the job; just append a warning to the message.
        row.message = sprintf('%s | MA recode warning: %s', ...
                              row.message, ME2.message);
    end

catch ME
    row.success = false;
    row.message = ME.message;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

% ------------------- NEW SUBFUNCTION -------------------
function recode_short_WK_to_MA(csv_file, codeWK, codeMA, maxDurSec)
% recode_short_WK_to_MA
%   Operates on a 1-Hz CSV exported by export_sleep_scores_with_MA.
%   Any contiguous run of WK samples (value == codeWK) with duration
%   <= maxDurSec seconds is re-labeled to MA (value == codeMA).
%
%   Assumes the CSV has one row per second.

    T = readtable(csv_file);

    % Try to find the score column in a robust way
    cand = {'state','score','Stage','stage','sleep_state','code'};
    scoreVarName = '';
    for i = 1:numel(cand)
        if ismember(cand{i}, T.Properties.VariableNames)
            scoreVarName = cand{i};
            break;
        end
    end
    if isempty(scoreVarName)
        error('Could not find a score/state column in %s', csv_file);
    end

    s = T.(scoreVarName);
    if ~isnumeric(s)
        s = double(s);
    end

    isWK = (s == codeWK);
    if ~any(isWK)
        % nothing to do
        return;
    end

    % Find contiguous WK runs
    d = diff([0; isWK; 0]);  % +1 at starts, -1 at ends
    startIdx = find(d == 1);
    endIdx   = find(d == -1) - 1;
    len      = endIdx - startIdx + 1;   % length in samples (1 sample = 1 s)

    % WK bouts with duration <= maxDurSec seconds
    shortMask = (len <= maxDurSec);
    if ~any(shortMask)
        % no short WK bouts
        return;
    end

    shortRuns = find(shortMask);
    for k = shortRuns(:)'
        s(startIdx(k):endIdx(k)) = codeMA;
    end

    T.(scoreVarName) = s;
    writetable(T, csv_file);
end
