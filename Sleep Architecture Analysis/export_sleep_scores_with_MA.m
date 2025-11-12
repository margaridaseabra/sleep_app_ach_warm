function R = export_sleep_scores_with_MA(mat_file, out_prefix, varargin)
% export_sleep_scores_with_MA('20251002-baseline-mouse2_scored.mat','mouse2_base', ...
%     'out_dir','C:\data\exports', 'ma_thresh_sec',15)
%
% Outputs in out_dir (or pwd if omitted):
%   <prefix>_scores_1Hz.csv
%   <prefix>_scores_bouts.csv
%   <prefix>_scores_tt.mat
%   <prefix>_export.log.txt   (optional summary log)
%
% Returns struct R with fields:
%   success, out_dir, files (struct), n_MA_bouts, MA_sec, WK_sec_before, WK_sec_after,
%   pct_wake_reclassified, epochSec, t0

p = inputParser;
addRequired(p,'mat_file',@ischar);
addRequired(p,'out_prefix',@ischar);

addParameter(p,'out_dir','',@ischar);                 % output folder (default: pwd)
addParameter(p,'scoreVar','sleep_scores',@ischar);
addParameter(p,'epochSec',[],@(x) isempty(x)||isscalar(x));
addParameter(p,'startTime',0,@isscalar);              % seconds offset if needed
addParameter(p,'fsVar','eeg_frequency',@ischar);      % used to infer epochSec if missing
addParameter(p,'eegVar','eeg',@ischar);
addParameter(p,'write_log',true,@islogical);

% Code scheme (edit to your scoring)
C_default = struct('WK',0,'NREM',1,'REM',0,'MA',15);
addParameter(p,'codes',C_default,@isstruct);

% MA rule
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);  % "less than" threshold
parse(p,mat_file,out_prefix,varargin{:});
S = p.Results; C = S.codes;

% ---- Prepare output dir
if isempty(S.out_dir), S.out_dir = pwd; end
if ~exist(S.out_dir,'dir'), mkdir(S.out_dir); end

% ---- Init return struct
R = struct('success',false,'out_dir',S.out_dir,'files',struct(), ...
           'n_MA_bouts',0,'MA_sec',0,'WK_sec_before',0,'WK_sec_after',0, ...
           'pct_wake_reclassified',0,'epochSec',NaN,'t0',S.startTime);

try
    D = load(S.mat_file);
    assert(isfield(D,S.scoreVar),'Missing %s',S.scoreVar);
    scores_epoch = D.(S.scoreVar)(:).';    % row vector of epoch labels
    Nepoch = numel(scores_epoch);

    % -------- Detect epoch length (sec)
    epochSec = S.epochSec;
    if isempty(epochSec)
        cand = [];
        for k = ["epoch_len_sec","epoch_length","epoch_sec","scoring_epoch_sec"]
            if isfield(D,k), cand = D.(k); break; end
        end
        if ~isempty(cand)
            epochSec = double(cand);
        else
            % Infer from EEG if present
            assert(isfield(D,S.fsVar)&&isfield(D,S.eegVar), ...
                'Provide epochSec or include EEG + fs to infer it.');
            fs = double(D.(S.fsVar)); Neeg = numel(D.(S.eegVar));
            T  = Neeg / fs; epochSec = round(T / Nepoch);
            if abs(epochSec-1)<0.2, epochSec=1; end
            if abs(epochSec-5)<0.5, epochSec=5; end
        end
    end
    R.epochSec = epochSec;

    % -------- Expand to 1 Hz (NO interpolation)
    t0 = S.startTime; R.t0 = t0;
    t_sec      = (0:Nepoch*epochSec-1) + t0;                 % absolute seconds
    score_1Hz  = repelem(uint16(scores_epoch), epochSec);    % 1 label per second
    score_1Hz  = score_1Hz(:); t_sec = t_sec(:);

    % Baseline wake seconds (before reclass)
    WK_sec_before = nnz(score_1Hz == C.WK);

    % -------- Reclassify micro-arousals:
    % Any contiguous WK run with duration < ma_thresh_sec becomes MA
    [wk_st, wk_en] = run_starts_ends(score_1Hz == C.WK);
    wk_dur = wk_en - wk_st + 1;                 % seconds
    ma_runs = find(wk_dur < S.ma_thresh_sec);
    for i = reshape(ma_runs,1,[])
        score_1Hz(wk_st(i):wk_en(i)) = uint16(C.MA);
    end

    % After reclass
    WK_sec_after = nnz(score_1Hz == C.WK);
    [ma_st, ma_en] = run_starts_ends(score_1Hz == C.MA);
    MA_sec  = sum(ma_en - ma_st + 1);
    n_MA    = numel(ma_st);

    pct_wake_reclassified = 100 * (WK_sec_before - WK_sec_after) / max(1, WK_sec_before);

    % -------- Save files
    f1 = fullfile(S.out_dir, sprintf('%s_scores_1Hz.csv', S.out_prefix));
    f2 = fullfile(S.out_dir, sprintf('%s_scores_bouts.csv', S.out_prefix));
    f3 = fullfile(S.out_dir, sprintf('%s_scores_tt.mat', S.out_prefix));
    flog = fullfile(S.out_dir, sprintf('%s_export.log.txt', S.out_prefix));

    % 1Hz CSV
    T1 = table(t_sec, score_1Hz, 'VariableNames',{'time_s','score'});
    writetable(T1, f1);

    % Bouts CSV (after MA reclass)
    idx = [true; diff(score_1Hz)~=0];
    b_st = find(idx);
    b_en = [b_st(2:end)-1; numel(score_1Hz)];
    Tb = table( ...
        score_1Hz(b_st), t_sec(b_st), t_sec(b_en), b_en-b_st+1, ...
        'VariableNames',{'score','start_s','end_s','dur_s'});
    writetable(Tb, f2);

    % Timetable MAT
    TT = timetable(seconds(t_sec), score_1Hz, 'VariableNames',{'score'});
    save(f3, 'TT','epochSec','t0','C');

    % Log (optional)
    if S.write_log
        fid = fopen(flog,'w');
        if fid>0
            fprintf(fid, 'Export summary — %s\n', datestr(now));
            fprintf(fid, 'Source file: %s\n', S.mat_file);
            fprintf(fid, 'Output dir : %s\n', S.out_dir);
            fprintf(fid, 'Epoch size : %d s\n', epochSec);
            fprintf(fid, 'MA rule    : WK runs < %d s -> MA\n', S.ma_thresh_sec);
            fprintf(fid, 'MA bouts   : %d\n', n_MA);
            fprintf(fid, 'MA seconds : %d\n', MA_sec);
            fprintf(fid, 'WK seconds (before): %d\n', WK_sec_before);
            fprintf(fid, 'WK seconds (after) : %d\n', WK_sec_after);
            fprintf(fid, '%% Wake reclassified: %.2f%%\n', pct_wake_reclassified);
            fprintf(fid, 'Files:\n  %s\n  %s\n  %s\n', f1, f2, f3);
            fclose(fid);
        end
    end

    % ---- Populate return + console summary
    R.success = true;
    R.files = struct('scores_1Hz_csv',f1,'bouts_csv',f2,'timetable_mat',f3,'log',flog);
    R.n_MA_bouts = n_MA;
    R.MA_sec = MA_sec;
    R.WK_sec_before = WK_sec_before;
    R.WK_sec_after  = WK_sec_after;
    R.pct_wake_reclassified = pct_wake_reclassified;

    fprintf(['✅ Exported sleep scoring with MA\n' ...
             '📂 Output: %s\n' ...
             '🧩 Epoch size: %d s | MA rule: WK runs < %d s\n' ...
             '🔎 MA bouts: %d | MA seconds: %d | Wake reclassified: %.2f%%\n' ...
             '🗃 Files:\n  %s\n  %s\n  %s\n'], ...
             S.out_dir, epochSec, S.ma_thresh_sec, n_MA, MA_sec, pct_wake_reclassified, f1, f2, f3);

catch ME
    warning('Export failed: %s', ME.message);
    R.success = false;
end
end

% ===== Helper: find contiguous run starts/ends from a logical vector =====
function [starts, ends_] = run_starts_ends(mask)
if isempty(mask), starts = []; ends_ = []; return; end
mask = mask(:) ~= 0;
d = diff([false; mask; false]);
starts = find(d == 1);
ends_  = find(d == -1) - 1;
end
