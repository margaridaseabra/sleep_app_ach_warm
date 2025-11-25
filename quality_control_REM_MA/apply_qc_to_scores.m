function apply_qc_to_scores(scores_csv, qc_rem_mat, qc_ma_mat, varargin)
% apply_qc_to_scores
% -------------------------------------------------------------------------
% Create a NEW scoring CSV (in a new folder) that includes:
%   - original score
%   - score_qc           : corrected score after REM & MA QC
%   - REM_good_epoch     : logical (1/0) mask per epoch
%   - MA_good_epoch      : logical (1/0) mask per epoch
%
% IMPORTANT:
%   - Original scores_csv is NEVER modified.
%   - A new CSV is written to: <scores_folder>/QC_scores/<basename>_QC.csv
%
% QC logic:
%   - Uses original runs of REM and MA from *original* `score`.
%   - For REM:
%       * flag == 1 (good): use trimmed start/end from qc_rem,
%         recode REM only inside those new boundaries.
%       * flag == 0 (bad): remove REM label for that original bout
%         (set to WK).
%       * flag NaN/-1: keep original REM scoring for that bout.
%   - Same logic for MA bouts with qc_ma.
%
% Usage:
%   apply_qc_to_scores('20251006_baseline_mouse4_WT_scored_scores_1Hz.csv', ...
%                      '20251006-baseline-mouse4_WT_REM_QC.mat', ...
%                      '20251006-baseline-mouse4_WT_MA_QC.mat');
%
% Options:
%   'codes' : struct with fields WK, NREM, REM, MA (default: 0,1,2,15)

    p = inputParser;
    addRequired(p,'scores_csv',@ischar);
    addRequired(p,'qc_rem_mat',@ischar);
    addRequired(p,'qc_ma_mat',@ischar);
    addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
    parse(p,scores_csv,qc_rem_mat,qc_ma_mat,varargin{:});
    S = p.Results;
    C = S.codes;

    %% 1) Load original scoring
    T = readtable(scores_csv);
    assert(all(ismember({'time_s','score'}, T.Properties.VariableNames)), ...
        'scores_csv must contain columns time_s and score');

    time_s = double(T.time_s(:));
    score  = double(T.score(:));
    n = numel(score);

    % New columns
    score_qc      = score;             % start as copy of original
    REM_good_epoch = false(n,1);
    MA_good_epoch  = false(n,1);

    %% 2) Recompute original REM & MA bouts from raw score
    [st_idx, en_idx] = runs_from_codes(score);
    run_code = score(st_idx);

    % REM bouts
    rem_idx   = find(run_code == C.REM);
    rem_st_idx = st_idx(rem_idx);
    rem_en_idx = en_idx(rem_idx);
    rem_st_s0  = time_s(rem_st_idx);
    rem_en_s0  = time_s(rem_en_idx);

    % MA bouts
    ma_idx   = find(run_code == C.MA);
    ma_st_idx = st_idx(ma_idx);
    ma_en_idx = en_idx(ma_idx);
    ma_st_s0  = time_s(ma_st_idx);
    ma_en_s0  = time_s(ma_en_idx);

    %% 3) Load QC files
    qc_rem = [];
    qc_ma  = [];
    if exist(qc_rem_mat,'file')
        tmp = load(qc_rem_mat);
        qc_rem = tmp.qc;
    else
        warning('No REM QC file: %s', qc_rem_mat);
    end
    if exist(qc_ma_mat,'file')
        tmp = load(qc_ma_mat);
        qc_ma = tmp.qc;
    else
        warning('No MA QC file: %s', qc_ma_mat);
    end

    %% 4) Apply REM QC
    if ~isempty(qc_rem)
        if numel(qc_rem.rem_st_s) ~= numel(rem_st_s0)
            warning('REM QC bouts (%d) do not match original REM bouts (%d). Skipping REM QC.', ...
                numel(qc_rem.rem_st_s), numel(rem_st_s0));
        else
            for k = 1:numel(rem_st_s0)
                flag  = qc_rem.flags(k);
                st0_i = rem_st_idx(k);
                en0_i = rem_en_idx(k);

                % clear original REM label for this bout if we'll change it
                if flag == 1 || flag == 0
                    score_qc(st0_i:en0_i) = C.WK;   % remove REM
                end

                if flag == 1
                    % Good REM: apply trimmed boundaries from qc_rem
                    st_new_s = qc_rem.rem_st_s(k);
                    en_new_s = qc_rem.rem_en_s(k);

                    idx_st_new = find(time_s >= st_new_s, 1, 'first');
                    idx_en_new = find(time_s <= en_new_s, 1, 'last');

                    if isempty(idx_st_new) || isempty(idx_en_new) || idx_st_new > idx_en_new
                        warning('Invalid trimmed REM bout %d; skipping.', k);
                        continue;
                    end

                    score_qc(idx_st_new:idx_en_new) = C.REM;
                    REM_good_epoch(idx_st_new:idx_en_new) = true;

                elseif flag == 0
                    % Bad REM: we already cleared the REM label (-> WK)
                    % nothing more to do
                else
                    % NaN or -1: keep original REM scoring
                    REM_good_epoch(st0_i:en0_i) = false;  % not "good"
                end
            end
        end
    end

    %% 5) Apply MA QC
    if ~isempty(qc_ma)
        if numel(qc_ma.ma_st_s) ~= numel(ma_st_s0)
            warning('MA QC bouts (%d) do not match original MA bouts (%d). Skipping MA QC.', ...
                numel(qc_ma.ma_st_s), numel(ma_st_s0));
        else
            for k = 1:numel(ma_st_s0)
                flag  = qc_ma.flags(k);
                st0_i = ma_st_idx(k);
                en0_i = ma_en_idx(k);

                if flag == 1 || flag == 0
                    score_qc(st0_i:en0_i) = C.WK;   % clear original MA label
                end

                if flag == 1
                    st_new_s = qc_ma.ma_st_s(k);
                    en_new_s = qc_ma.ma_en_s(k);

                    idx_st_new = find(time_s >= st_new_s, 1, 'first');
                    idx_en_new = find(time_s <= en_new_s, 1, 'last');

                    if isempty(idx_st_new) || isempty(idx_en_new) || idx_st_new > idx_en_new
                        warning('Invalid trimmed MA bout %d; skipping.', k);
                        continue;
                    end

                    score_qc(idx_st_new:idx_en_new) = C.MA;
                    MA_good_epoch(idx_st_new:idx_en_new) = true;

                elseif flag == 0
                    % bad MA: just removed it
                else
                    % NaN/-1: keep original
                    MA_good_epoch(st0_i:en0_i) = false;
                end
            end
        end
    end

    %% 6) Write NEW CSV in QC_scores folder
    T.score_qc         = score_qc;
    T.REM_good_epoch   = REM_good_epoch;
    T.MA_good_epoch    = MA_good_epoch;

    [pth,base,ext] = fileparts(scores_csv);
    out_dir = fullfile(pth,'QC_scores');
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    out_name = fullfile(out_dir, [base '_QC' ext]);
    writetable(T,out_name);

    fprintf('QC scoring written to:\n   %s\n', out_name);
end

%% ---- helper: contiguous runs from codes ----
function [st,en] = runs_from_codes(code)
    code = code(:);
    if isempty(code)
        st = []; en = [];
        return;
    end
    d = diff(code);
    st = [1; find(d~=0)+1];
    en = [st(2:end)-1; numel(code)];
end
