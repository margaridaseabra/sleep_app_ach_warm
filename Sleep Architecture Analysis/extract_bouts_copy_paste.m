function extract_bouts_copy_paste(data_dir, state, condition, period_xlsx)
% EXTRACT_BOUTS_COPY_PASTE - Get lists ready to copy-paste into Prism
%
% Usage:
%   extract_bouts_copy_paste('/Users/.../scores', 'REM', 'baseline');
%   extract_bouts_copy_paste('/Users/.../scores', 'REM', 'ambtemp', 'ambtemp_window.xlsx');
%
% For ambtemp, if period_xlsx is provided, only bouts whose *entire*
% duration lies between [Time started (s), Time ended (s)] for that
% mouse/genotype (from the XLSX) are kept.

if nargin < 4
    period_xlsx = [];
end

codes = struct('WK',0,'NREM',1,'REM',2,'MA',15);
state_code = codes.(state);

F = dir(fullfile(data_dir, '**', '*_scores_1Hz.csv'));
F = F(~[F.isdir]);

APP_bouts = [];
WT_bouts  = [];

% ---------- Load ambtemp window info (if applicable) ----------
use_window = strcmpi(condition,'ambtemp') && ~isempty(period_xlsx);
if use_window
    W = readtable(period_xlsx);

    % Mouse column numeric
    if ~isnumeric(W.Mouse)
        W.Mouse = double(W.Mouse);
    end

    % Genotype column (if present)
    if ismember('Genotype', W.Properties.VariableNames)
        W.Genotype = upper(string(W.Genotype));
    else
        W.Genotype = repmat("", height(W), 1);
    end

    vnames = W.Properties.VariableNames;

    % Start column: "Time started (s)" or similar
    start_idx = find(strcmpi(vnames, 'Time started (s)'), 1);
    if isempty(start_idx)
        mask_start = contains(lower(vnames), 'time') & ...
                     contains(lower(vnames), 'start') & ...
                     contains(lower(vnames), 's');
        start_idx = find(mask_start, 1);
    end
    if isempty(start_idx)
        error('Period file does not contain a "Time started (s)"-like column.');
    end
    start_colname = vnames{start_idx};

    % Optional end column: e.g. "Time ended (s)" (if not present, window goes to end)
    mask_end = contains(lower(vnames), 'time') & ...
               contains(lower(vnames), 'end')  & ...
               contains(lower(vnames), 's');
    if any(mask_end)
        end_colname = vnames{find(mask_end, 1)};
    else
        end_colname = '';
    end
end

epoch_sec = 1;  % scores_1Hz → 1 s per index

for i = 1:numel(F)
    [~, base, ~] = fileparts(F(i).name);

    expr = '^(?<date>\d{8})_(?<cond>[^_]+)_(?<mouse>mouse\d+)_(?<geno>APP|WT)';
    m = regexp(base, expr, 'names');
    if isempty(m) || ~strcmp(m.cond, condition)
        continue;
    end

    % ---------- Window for this mouse/genotype (ambtemp only) ----------
    win_start_s = 0;
    win_end_s   = Inf;
    if use_window
        % mouse10 -> 10
        tok = regexp(m.mouse, 'mouse(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            mouseID = str2double(tok{1});
        else
            mouseID = NaN;
        end

        geno_here = upper(string(m.geno));
        idxW = (W.Mouse == mouseID);
        if any(idxW) && ismember('Genotype', W.Properties.VariableNames)
            idxW = idxW & (W.Genotype == geno_here);
        end

        if any(idxW)
            rowW = find(idxW, 1, 'first');
            win_start_s = W{rowW, start_colname};
            if ~isempty(end_colname)
                win_end_s = W{rowW, end_colname};
            else
                win_end_s = Inf;  % if no end column, go until end of recording
            end

            if isfinite(win_end_s)
                fprintf('Using window %.1f–%.1f s for %s_%s\n', ...
                        win_start_s, win_end_s, m.mouse, m.geno);
            else
                fprintf('Using window from %.1f s to end for %s_%s\n', ...
                        win_start_s, m.mouse, m.geno);
            end
        else
            fprintf('No window info for %s_%s in %s → using full recording.\n', ...
                    m.mouse, m.geno, period_xlsx);
        end
    end

    % ---------- Load scores ----------
    T = readtable(fullfile(F(i).folder, F(i).name));

    % Find score column
    if ismember('state', T.Properties.VariableNames)
        scores = T.state;
    elseif ismember('score', T.Properties.VariableNames)
        scores = T.score;
    else
        scores = T{:,1};
    end

    % ---------- Detect bouts (but keep only those fully inside window) ----------
    in_bout = false;
    start_idx = NaN;

    for j = 1:numel(scores)
        if scores(j) == state_code && ~in_bout
            in_bout = true;
            start_idx = j;

        elseif scores(j) ~= state_code && in_bout
            duration = j - start_idx;  % in seconds (1 Hz)
            % Check if this bout is fully inside [win_start_s, win_end_s]
            keep_bout = true;
            if use_window
                bout_start_s = (start_idx - 1) * epoch_sec;
                bout_end_s   = bout_start_s + duration * epoch_sec; % exclusive
                if ~(bout_start_s >= win_start_s && bout_end_s <= win_end_s)
                    keep_bout = false;
                end
            end

            if keep_bout
                if strcmp(m.geno, 'APP')
                    APP_bouts = [APP_bouts; duration]; %#ok<AGROW>
                else
                    WT_bouts = [WT_bouts; duration];  %#ok<AGROW>
                end
            end

            in_bout = false;
        end
    end

    % Bout still ongoing at the end of the file
    if in_bout
        duration = numel(scores) - start_idx + 1;
        keep_bout = true;
        if use_window
            bout_start_s = (start_idx - 1) * epoch_sec;
            bout_end_s   = bout_start_s + duration * epoch_sec;
            if ~(bout_start_s >= win_start_s && bout_end_s <= win_end_s)
                keep_bout = false;
            end
        end

        if keep_bout
            if strcmp(m.geno, 'APP')
                APP_bouts = [APP_bouts; duration]; %#ok<AGROW>
            else
                WT_bouts = [WT_bouts; duration];  %#ok<AGROW>
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('BOUT DURATIONS FOR PRISM\n');
fprintf('========================================\n\n');
fprintf('APP: %d bouts (median = %.1f sec)\n', length(APP_bouts), median(APP_bouts));
fprintf('WT:  %d bouts (median = %.1f sec)\n\n', length(WT_bouts), median(WT_bouts));

% Save to text files
writematrix(APP_bouts, sprintf('%s_%s_APP.txt', state, condition));
writematrix(WT_bouts,  sprintf('%s_%s_WT.txt',  state, condition));
writematrix(ones(size(APP_bouts)), sprintf('%s_%s_APP_Y.txt', state, condition));
writematrix(ones(size(WT_bouts)),  sprintf('%s_%s_WT_Y.txt',  state, condition));

fprintf('FILES SAVED:\n');
fprintf('  %s_%s_APP.txt     <- APP bout durations (X values)\n', state, condition);
fprintf('  %s_%s_APP_Y.txt   <- APP event codes (Y values, all 1s)\n', state, condition);
fprintf('  %s_%s_WT.txt      <- WT bout durations (X values)\n', state, condition);
fprintf('  %s_%s_WT_Y.txt    <- WT event codes (Y values, all 1s)\n\n', state, condition);

fprintf('IN PRISM:\n');
fprintf('1. NEW TABLE & GRAPH → XY → Survival\n');
fprintf('2. In the X column: Copy-paste from %s_%s_APP.txt, then %s_%s_WT.txt below\n', ...
        state, condition, state, condition);
fprintf('3. In Group A (Y): Copy-paste from %s_%s_APP_Y.txt (all 1s)\n', state, condition);
fprintf('4. In Group B (Y): Copy-paste from %s_%s_WT_Y.txt (all 1s)\n', state, condition);
fprintf('5. Analyze → Kaplan-Meier\n');
fprintf('========================================\n\n');
end
