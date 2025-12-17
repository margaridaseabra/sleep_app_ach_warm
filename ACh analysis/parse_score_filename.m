function info = parse_score_filename(fname)
% Parse score CSV names in either of these forms:
%   1) YYYYMMDD-cond-mouseNN-GENO_scored_scores_1Hz.csv
%   2) YYYYMMDD_cond_mouseNN_GENO_scored_scores_1Hz.csv
%   3) Same as above but with "_crop" before ".csv"
%
% Examples:
%   20251012-baseline-mouse12-WT_scored_scores_1Hz.csv
%   20251002_baseline_mouse2_WT_scored_scores_1Hz_crop.csv

    % Allow both '-' and '_' as separators between fields
    % and make the trailing "_crop" optional.
    pat = ['^(?<date>\d{8})[-_](?<cond>[^-_]+)[-_]mouse(?<mouse>\d+)[-_]' ...
           '(?<geno>[^-_]+)_scored_scores_1Hz(_crop)?\.csv$'];

    m = regexp(fname, pat, 'names');
    if isempty(m)
        info = [];
    else
        info = struct();
        info.date  = m.date;
        info.cond  = lower(m.cond);   % lowercase for matching
        info.mouse = m.mouse;         % number as string
        info.geno  = m.geno;          % keep genotype label as-is
    end
end
