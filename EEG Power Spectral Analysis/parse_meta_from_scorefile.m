function meta = parse_meta_from_scorefile(scoreFile)
    % Parse date, condition, mouseID, genotype from a score filename
    % Example: 20251025-ambtemp-mouse12-WT_scored_scores_1Hz.csv

    [~, name, ~] = fileparts(scoreFile);

    % Remove suffix starting at '_scored'
    base = regexprep(name, '_scored.*$', '');
    % base = '20251025-ambtemp-mouse12-WT'

    parts = split(base, '-');

    if numel(parts) < 4
        error('Unexpected score filename format: %s', name);
    end

    meta.date      = parts{1};  % '20251025'
    meta.condition = parts{2};  % 'ambtemp', 'drugs', 'baseline', etc.
    meta.mouseID   = parts{3};  % 'mouse12'
    meta.genotype  = parts{4};  % 'WT', 'APP', etc.
end
