function [rows_clean, OUTLIERS] = remove_outlier_mice_perhour(rows_perhr, states_to_check, varargin)
% remove_outlier_mice_perhour
% -------------------------------------------------------------------------
% For BASELINE data in rows_perhr, and for each requested STATE + GENOTYPE:
%   - Compute mean bouts_per_h per mouse across all hours
%   - Use isoutlier to flag "outlier mice" within that state+genotype
%   - Remove those mouse×state rows from rows_perhr (all hours, baseline)
%
% This is meant as a *preprocessing* step before plotting / ANOVA.
%
% INPUTS
%   rows_perhr      : table from run_group_sleep_architecture (group_per_hour)
%   states_to_check : string/cell array of states, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'Method'        : method for isoutlier ('median','gesd','mean'); default 'median'
%   'Threshold'     : ThresholdFactor for isoutlier (e.g. 3); default 3
%   'Verbose'       : true/false, print info; default true
%
% OUTPUTS
%   rows_clean  : rows_perhr with outlier mouse×state rows removed (baseline only)
%   OUTLIERS    : table listing which mice were flagged, with their mean bouts/h
% -------------------------------------------------------------------------

if nargin < 2 || isempty(states_to_check)
    states_to_check = ["WK","MA","NREM","REM"];
end
states_to_check = string(states_to_check(:)).';

p = inputParser;
addParameter(p,'Method','median',@(x)ischar(x) || isstring(x));
addParameter(p,'Threshold',3,@(x)isscalar(x) && x>0);
addParameter(p,'Verbose',true,@(x)islogical(x) && isscalar(x));
parse(p, varargin{:});

method    = p.Results.Method;
thresh    = p.Results.Threshold;
verbose   = p.Results.Verbose;

rows_clean = rows_perhr;
OUTLIERS   = table([],[],[],[],[], ...
    'VariableNames', {'state','genotype','mouse','mean_bouts_per_h','is_outlier'});

% work only on baseline rows for outlier detection
T = rows_perhr;
cond_str    = lower(strtrim(T.condition));
is_baseline = cond_str == "baseline";
T = T(is_baseline, :);

if isempty(T)
    warning('No baseline rows found in rows_perhr. No outlier detection done.');
    return;
end

for st = states_to_check
    for gen = ["WT","APP"]
        sub = T(T.state == st & T.genotype == gen, :);
        if isempty(sub)
            continue;
        end

        % mean bouts_per_h per mouse
        Gm = groupsummary(sub, 'mouse', 'mean', 'bouts_per_h');
        Gm.Properties.VariableNames{end} = 'mean_bouts_per_h';

        vals = Gm.mean_bouts_per_h;
        if numel(vals) < 3
            % too few mice to sensibly define outliers
            continue;
        end

        is_out = isoutlier(vals, method, 'ThresholdFactor', thresh);

        if any(is_out)
            these = Gm(is_out, :);

            % build OUTLIERS rows
            new_rows = table( ...
                repmat(st, height(these), 1), ...
                repmat(gen, height(these), 1), ...
                these.mouse, ...
                these.mean_bouts_per_h, ...
                is_out(is_out), ...
                'VariableNames', OUTLIERS.Properties.VariableNames);

            OUTLIERS = [OUTLIERS; new_rows]; %#ok<AGROW>

            if verbose
                fprintf('\n[Outlier detection] State=%s, Genotype=%s\n', st, gen);
                disp(these(:, {'mouse','mean_bouts_per_h'}));
            end

             % remove these mouse×state baseline rows from rows_clean
            for i = 1:height(these)
                mID = these.mouse(i);

                % recompute baseline mask on *rows_clean* (current size)
                cond_str_clean    = lower(strtrim(rows_clean.condition));
                is_baseline_clean = cond_str_clean == "baseline";

                mask_rm = (rows_clean.state    == st)   & ...
                          (rows_clean.genotype == gen)  & ...
                          (rows_clean.mouse    == mID)  & ...
                          is_baseline_clean;  % baseline only

                rows_clean(mask_rm, :) = [];
            end

        end
    end
end

if verbose && isempty(OUTLIERS)
    fprintf('[Outlier detection] No outlier mice detected with method=%s, threshold=%.2f.\n', ...
            method, thresh);
end
end
