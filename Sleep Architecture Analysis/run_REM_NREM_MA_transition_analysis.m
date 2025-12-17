function OUT = run_REM_NREM_MA_transition_analysis(input_dir, varargin)
% run_REM_NREM_MA_transition_analysis
% -------------------------------------------------------------------------
% For ALL scored recordings in input_dir:
%   - run sleep_architecture_from_scores to get Runs table
%   - extract all REM bouts and classify them as Short/Medium/Long
%     based on global REM-duration quantiles.
%   - for each REM bout compute:
%         * prev_state, prev_state_dur_s
%         * prev_nrem_dur_s  (if previous state is NREM)
%         * ma_prev_window_sec  (total MA in last ma_window_sec before REM)
%         * next_state, next_state_dur_s
%         * wake_after_rem_dur_s (continuous wake after REM)
%
%   - for baseline only, for each NREM bout compute whether:
%         * NREM -> REM -> Wake  (triad)
%         * NREM -> Wake        (no REM in between)
%     and if triad: which REM type (Short/Medium/Long).
%
% INPUT:
%   input_dir : folder with *_scores_1Hz.csv
%
% OPTIONS (name-value):
%   'pattern'         : file pattern (default '*_scores_1Hz.csv')
%   'codes'           : struct('WK',0,'NREM',1,'REM',2,'MA',15)
%   'quantiles'       : e.g. [33 66] for Short/Med/Long
%   'ma_window_sec'   : window before REM onset to integrate MA (default 120)
%   'wake_include_MA' : if true, treat MA as wake when computing
%                       wake_after_rem_dur_s (default true)
%   'out_dir'         : where to save CSVs (default input_dir/rem_transitions)
%
% OUTPUT:
%   OUT.REM_TRANS   : per-REM table with transition features + rem_type
%   OUT.NREM_TRANS  : per-NREM table with sequence labels
%   OUT.REM_THRESH  : thresholds used for Short/Med/Long
%   OUT.params      : struct of parameters
%   OUT.files       : struct with CSV paths
% -------------------------------------------------------------------------

p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*_scores_1Hz.csv',@ischar);
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'quantiles',[33 66],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'ma_window_sec',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'wake_include_MA',true,@islogical);
addParameter(p,'out_dir','',@ischar);
parse(p, input_dir, varargin{:});
S = p.Results;
C = S.codes;

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);

if isempty(S.out_dir)
    S.out_dir = fullfile(S.input_dir, 'rem_transitions');
end
if ~isfolder(S.out_dir)
    mkdir(S.out_dir);
end

if exist('sleep_architecture_from_scores','file') ~= 2
    error(['sleep_architecture_from_scores.m is not on the MATLAB path.\n' ...
           'Add it, e.g.: addpath(genpath(''%s''))'], pwd);
end

% ---------- discover files ----------
F = dir(fullfile(S.input_dir, '**', S.pattern));
F = F(~[F.isdir]);
assert(~isempty(F), 'No files matched "%s" under %s', S.pattern, S.input_dir);

% ---------- collectors ----------
REM_TRANS = table();   % one row per REM bout
NREM_TRANS = table();  % one row per NREM bout (for sequence analysis)
all_rem_durs = [];     % to build global quantiles

for i = 1:numel(F)
    csv_path = fullfile(F(i).folder, F(i).name);
    info = parse_info_from_filename_relaxed(F(i).name);  % same relaxed parser as before

    if ~info.ok
        warning('Filename not recognized pattern, skipping: %s', F(i).name);
        continue;
    end

    try
        OUT_i = sleep_architecture_from_scores(csv_path, ...
            'codes', C, ...
            'includeMA', true, ...
            'ma_thresh_sec', 15, ...
            'reclassify_short_wake_to_MA', true, ...
            'out_prefix','', 'out_dir','', ...
            'write_log', false, 'verbose', false);
    catch ME
        warning('Failed file %s: %s', F(i).name, ME.message);
        continue;
    end

    Runs = OUT_i.runs;
    if isempty(Runs), continue; end

    % ensure dur_s
    if ~ismember('dur_s', Runs.Properties.VariableNames)
        Runs.dur_s = Runs.end_s - Runs.start_s + 1;
    end

    nB = height(Runs);

    % ---------------------------------------------------------------------
    % 1) PER-REM FEATURES
    % ---------------------------------------------------------------------
    for j = 1:nB
        if Runs.state_code(j) ~= C.REM
            continue;
        end

        this_start = Runs.start_s(j);
        this_end   = Runs.end_s(j);
        this_dur   = Runs.dur_s(j);

        % previous run
        if j > 1
            prev_code = Runs.state_code(j-1);
            prev_dur  = Runs.dur_s(j-1);
        else
            prev_code = NaN;
            prev_dur  = NaN;
        end

        prev_state = code_to_state(prev_code, C);
        prev_nrem_dur = NaN;
        if ~isnan(prev_code) && prev_code == C.NREM
            prev_nrem_dur = prev_dur;
        end

        % MA in window before REM onset
        t0 = this_start - S.ma_window_sec;
        t1 = this_start - 1;
        ma_prev_window = 0;
        if t0 < Runs.start_s(1)
            t0 = Runs.start_s(1);
        end
        if t1 >= t0
            isMA = Runs.state_code == C.MA;
            for k = find(isMA).'
                rs = Runs.start_s(k);
                re = Runs.end_s(k);
                overlap = max(0, min(re, t1) - max(rs, t0) + 1);
                ma_prev_window = ma_prev_window + overlap;
            end
        end

        % next run
        if j < nB
            next_code = Runs.state_code(j+1);
            next_dur  = Runs.dur_s(j+1);
        else
            next_code = NaN;
            next_dur  = NaN;
        end
        next_state = code_to_state(next_code, C);

        % wake after REM: continuous wake starting right after REM
        wake_after = 0;
        if j < nB
            if S.wake_include_MA
                wake_codes = [C.WK, C.MA];
            else
                wake_codes = C.WK;
            end

            k = j+1;
            while k <= nB && ismember(Runs.state_code(k), wake_codes)
                wake_after = wake_after + Runs.dur_s(k);
                k = k+1;
            end
        end

        % Append row (rem_type assigned later, after quantiles)
        REM_TRANS = [REM_TRANS; ...
            table( ...
                string(F(i).name), ...
                string(info.date), ...
                string(info.condition), ...
                string(info.mouse), ...
                string(info.genotype), ...
                this_start, this_end, this_dur, ...
                prev_state, prev_dur, prev_nrem_dur, ma_prev_window, ...
                next_state, next_dur, wake_after, ...
                'VariableNames', {'file','date','condition','mouse','genotype', ...
                                  'rem_start_s','rem_end_s','rem_dur_s', ...
                                  'prev_state','prev_state_dur_s','prev_nrem_dur_s', ...
                                  'ma_prev_window_sec', ...
                                  'next_state','next_state_dur_s', ...
                                  'wake_after_rem_dur_s'})]; %#ok<AGROW>

        all_rem_durs = [all_rem_durs; this_dur]; %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % 2) PER-NREM SEQUENCES (baseline only)
    % ---------------------------------------------------------------------
    is_baseline = (lower(strtrim(info.condition)) == "baseline");
    if ~is_baseline
        continue;
    end

    for j = 1:nB
        if Runs.state_code(j) ~= C.NREM
            continue;
        end

        seq_type = "Other";
        seq_rem_dur = NaN;   % if triad includes REM
        % NREM -> REM -> Wake vs NREM -> Wake

        if j < nB
            code2 = Runs.state_code(j+1);

            if code2 == C.REM
                % NREM -> REM -> ?
                if j+2 <= nB
                    code3 = Runs.state_code(j+2);
                else
                    code3 = NaN;
                end

                % define wake (could include MA)
                if S.wake_include_MA
                    wake_codes = [C.WK, C.MA];
                else
                    wake_codes = C.WK;
                end

                if ~isnan(code3) && ismember(code3, wake_codes)
                    seq_type   = "NREM_REM_WAKE";
                    seq_rem_dur = Runs.dur_s(j+1);
                else
                    seq_type   = "NREM_REM_other";
                    seq_rem_dur = Runs.dur_s(j+1);
                end

            else
                % NREM -> Wake (no REM in between)
                if S.wake_include_MA
                    wake_codes = [C.WK, C.MA];
                else
                    wake_codes = C.WK;
                end

                if ismember(code2, wake_codes)
                    seq_type = "NREM_WAKE";
                else
                    seq_type = "NREM_other";
                end
            end
        end

        NREM_TRANS = [NREM_TRANS; ...
            table( ...
                string(F(i).name), ...
                string(info.date), ...
                string(info.condition), ...
                string(info.mouse), ...
                string(info.genotype), ...
                Runs.start_s(j), Runs.end_s(j), Runs.dur_s(j), ...
                seq_type, seq_rem_dur, ...
                'VariableNames', {'file','date','condition','mouse','genotype', ...
                                  'nrem_start_s','nrem_end_s','nrem_dur_s', ...
                                  'seq_type','seq_rem_dur_s'})]; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
% 3) Global REM quantiles -> Short / Medium / Long
% -------------------------------------------------------------------------
all_rem_durs = double(all_rem_durs);
all_rem_durs = all_rem_durs(~isnan(all_rem_durs) & all_rem_durs > 0);

q = prctile(all_rem_durs, S.quantiles);
q1 = q(1); q2 = q(2);

REM_THRESH = struct('quantiles',S.quantiles,'q1_sec',q1,'q2_sec',q2);

% Assign types to REM_TRANS
rem_type = strings(height(REM_TRANS),1);
rem_type(REM_TRANS.rem_dur_s <  q1) = "Short";
rem_type(REM_TRANS.rem_dur_s >= q1 & REM_TRANS.rem_dur_s < q2) = "Medium";
rem_type(REM_TRANS.rem_dur_s >= q2) = "Long";
REM_TRANS.rem_type = rem_type;

% Assign types to NREM_TRANS seq_rem_dur_s where applicable
nrem_rem_type = strings(height(NREM_TRANS),1);
isTriad = ~isnan(NREM_TRANS.seq_rem_dur_s) & NREM_TRANS.seq_rem_dur_s > 0;
d = NREM_TRANS.seq_rem_dur_s;

nrem_rem_type(isTriad & d <  q1) = "Short";
nrem_rem_type(isTriad & d >= q1 & d < q2) = "Medium";
nrem_rem_type(isTriad & d >= q2)          = "Long";
NREM_TRANS.seq_rem_type = nrem_rem_type;

% -------------------------------------------------------------------------
% 4) Save CSVs
% -------------------------------------------------------------------------
rem_trans_csv  = fullfile(S.out_dir, 'REM_transitions_all.csv');
nrem_trans_csv = fullfile(S.out_dir, 'NREM_sequences_baseline.csv');

if ~isempty(REM_TRANS),  writetable(REM_TRANS,  rem_trans_csv); end
if ~isempty(NREM_TRANS), writetable(NREM_TRANS, nrem_trans_csv); end

% -------------------------------------------------------------------------
% 5) Pack output
% -------------------------------------------------------------------------
OUT = struct();
OUT.REM_TRANS   = REM_TRANS;
OUT.NREM_TRANS  = NREM_TRANS;
OUT.REM_THRESH  = REM_THRESH;
OUT.params      = S;
OUT.files = struct('rem_trans_csv',rem_trans_csv, ...
                   'nrem_trans_csv',nrem_trans_csv);

fprintf('✅ REM–NREM–MA transition analysis finished. Outputs in: %s\n', S.out_dir);
end

% =====================================================================
% Helper: map state_code -> string label using C
% =====================================================================
function st = code_to_state(code, C)
    st = "NONE";
    if isnan(code), return; end
    if isfield(C,'WK')   && code == C.WK,   st = "WK";   return; end
    if isfield(C,'NREM') && code == C.NREM, st = "NREM"; return; end
    if isfield(C,'REM')  && code == C.REM,  st = "REM";  return; end
    if isfield(C,'MA')   && code == C.MA,   st = "MA";   return; end
end

% =====================================================================
% Helper: relaxed filename parser (same idea as in REM_length script)
% =====================================================================
function info = parse_info_from_filename_relaxed(fname)
    info = struct('date','','condition','','mouse','','genotype','','ok',false);
    [~, name, ~] = fileparts(fname);
    parts = regexp(name,'[-_]','split');
    if numel(parts) < 4, return; end
    info.date      = parts{1};
    info.condition = lower(parts{2});
    info.mouse     = parts{3};
    info.genotype  = upper(parts{4});
    info.ok        = true;
end
