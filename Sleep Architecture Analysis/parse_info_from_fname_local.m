% =====================================================================
% Local helper: parse your filename pattern
% Expected pattern:
%   YYYYMMDD_condition_mouseX_APP_scored_scores_1Hz.csv
%   YYYYMMDD_condition_mouseX_WT_scored_scores_1Hz.csv
% =====================================================================
function info = parse_info_from_fname_local(fname)
    % Returns struct with: date, condition, mouse, genotype, ok

    [~, base, ~] = fileparts(fname);

    info = struct('date',"", 'condition',"", 'mouse',"", ...
                  'genotype',"", 'ok', false);

    % Pattern:
    %   <date>_<cond>_<mouse>_<APP/WT>[_scored]_scores_1Hz
    %
    % Examples:
    %   20251001_baseline_mouse1_APP_scored_scores_1Hz
    %   20251002_baseline_mouse2_WT_scored_scores_1Hz
    %
    expr = ['^(?<date>\d{8})_' ...         % 8-digit date
            '(?<condition>[^_]+)_' ...     % condition up to next underscore
            '(?<mouse>mouse\d+)_' ...      % mouse + digits
            '(?<genotype>APP|WT)' ...      % APP or WT
            '(?<scored>_scored)?' ...      % optional "_scored"
            '_scores_1Hz$'];               % final suffix

    m = regexp(base, expr, 'names');

    if isempty(m)
        % pattern not recognized
        return;
    end

    info.date      = string(m.date);
    info.condition = string(m.condition);
    info.mouse     = string(m.mouse);
    info.genotype  = string(m.genotype);
    info.ok        = true;
end
