function meta = parse_meta_from_scorefile_ambtemp(scoreFile)
% Extract date, condition, mouse ID and genotype from score CSV filename.
%
% Expected filename (no path, no extension):
%   YYYYMMDD_cond_mouseN_GENO_scored_scores_1Hz
%
% Example:
%   20251023_drugs_mouse3_APP_scored_scores_1Hz

    [~, name, ~] = fileparts(scoreFile);  % e.g. 20251023_drugs_mouse3_APP_scored_scores_1Hz

    % Normalise separators just in case some have '-'
    name = strrep(name, '-', '_');

    % Remove the suffix
    name = regexprep(name, '_scored_scores_1Hz_crop', '');

    % Now expect: YYYYMMDD_cond_mouseN_GENO
    tokens = regexp(name, ...
        '^(?<date>\d{8})_(?<cond>[^_]+)_mouse(?<mouse>\d+)_(?<geno>WT|APP)$', ...
        'names');

    if isempty(tokens)
        error('Unexpected score filename format: %s', name);
    end

    meta = struct();
    meta.date      = tokens.date;
    meta.condition = tokens.cond;
    meta.mouseID   = sprintf('mouse%s', tokens.mouse);
    meta.genotype  = tokens.geno;
end
