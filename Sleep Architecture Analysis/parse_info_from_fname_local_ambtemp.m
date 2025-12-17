function info = parse_info_from_fname_local_ambtemp(fname)
    % Returns struct with: date, condition, mouse, genotype, ok

    [~, base, ~] = fileparts(fname);

    info = struct('date',"", 'condition',"", 'mouse',"", ...
                  'genotype',"", 'ok', false);

    % Expected pattern (no .csv here, we use "base"):
    %   YYYYMMDD_condition_mouseX_APP_scored_scores_1Hz_crop
    %   YYYYMMDD_condition_mouseX_WT_scored_scores_1Hz_crop
    %
    expr = ['^(?<date>\d{8})_' ...         % 8-digit date
            '(?<condition>[^_]+)_' ...     % condition up to next underscore
            '(?<mouse>mouse\d+)_' ...      % "mouse" + digits
            '(?<genotype>APP|WT)' ...      % genotype
            '(?<scored>_scored)?' ...      % optional "_scored"
            '_scores_1Hz_crop$'];          % final suffix, end of string

    m = regexp(base, expr, 'names', 'once');

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
