function [states, epochs_t, MA_tbl, mouse_id, geno, cond] = read_scores_csv_drugs(scores_file)
% read_scores_csv
% -------------------------------------------------------------
% Compatible with CSVs like:
%   time_s, score
% where "score" is a numeric sleep stage code:
%   0 = Wake, 1 = NREM, 2 = REM, 3 = MA  (edit if needed)
%
% Filename format (flexible):
%   YYYYMMDD_cond_mouseID_GENO_scored_scores_1Hz.csv
%   YYYYMMDD-cond-mouseID-GENO_scored_scores_1Hz.csv
%
% Examples:
%   20251001_baseline_mouse1_APP_scored_scores_1Hz.csv
%   20251001-baseline-mouse1-APP_scored_scores_1Hz.csv

    % ---------- A) Parse metadata from filename -------------------------
    [~, fname, ~] = fileparts(scores_file);

    % Allow both '-' and '_' between date/cond/mouse/geno
    tokens = regexp(fname, ...
        '^(?<date>\d{8})[-_](?<cond>[^-_]+)[-_](?<mouse>[^-_]+)[-_](?<geno>[^-_.]+)_scored_scores_1Hz_crop', ...
        'names');

    if isempty(tokens)
        error(['Unexpected scores filename format: %s\n' ...
               'Expected something like: 20251001_baseline_mouse1_APP_scored_scores_1Hz_crop'], ...
               fname);
    end

    cond     = string(tokens.cond);
    mouse_id = string(tokens.mouse);
    geno     = string(tokens.geno);

    % ---------- B) Read the table --------------------------------------
    T = readtable(scores_file);

    % ---------- C) Map columns -----------------------------------------
    if all(ismember({'time_s','score'}, T.Properties.VariableNames))
        epochs_t   = T.time_s;
        score_code = T.score;
    else
        error('Expected columns "time_s" and "score" in %s', scores_file);
    end

    % ---------- D) Convert numeric score → state labels ----------------
    % Adjust these mappings if your scoring codes differ!
    states = strings(height(T),1);
    states(score_code==0) = "WK";      % wake
    states(score_code==1) = "NREM";    % NREM
    states(score_code==2) = "REM";     % REM
    states(score_code==3) = "MA";      % microarousal

    % ---------- E) Build MA table --------------------------------------
    isMA  = (states=="MA");
    MA_tbl = table();
    if any(isMA)
        idx = find(isMA);
        is_new          = [true; diff(idx)>1];
        bout_start_idx  = idx(is_new);
        bout_end_idx    = [idx(find(is_new)-1); idx(end)]; 

        nB      = numel(bout_start_idx);
        start_s = zeros(nB,1);
        end_s   = zeros(nB,1);
        for iB = 1:nB
            start_s(iB) = epochs_t(bout_start_idx(iB));
            end_s(iB)   = epochs_t(bout_end_idx(iB)) + 1; % 1 s step
        end
        MA_tbl = table(start_s, end_s);
        MA_tbl.Properties.VariableNames = {'start_s','end_s'};
    end
end
