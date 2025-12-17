function rows_overall_ambtemp = run_group_sleep_architecture_ambtemp_overall(input_dir, varargin)
% run_group_sleep_architecture_ambtemp_overall
% -------------------------------------------------------------------------
% Minimal version for AMBTEMP project:
%   - Takes cropped *_scored_scores_1Hz_crop.csv files
%   - Runs sleep_architecture_from_scores on each
%   - Builds ONLY an OVERALL table:
%
%       state, n_bouts, total_dur_s, mean_bout_dur_s,
%       file, date, condition, mouse, genotype
%
%   - Returns this as rows_overall_ambtemp
%   - Saves:
%       rows_overall_ambtemp.mat
%       rows_overall_ambtemp.csv
%
% No per-hour tables, no ANOVA, no plots.
%
% Assumes filenames like:
%   20251005_baseline_mouse8_WT_scored_scores_1Hz_crop.csv
%   20251014_ambtemp_mouse11_WT_scored_scores_1Hz_crop.csv
%
% i.e. pattern:
%   YYYYMMDD_condition_mouseX_APP_scored_scores_1Hz_crop.csv
%   YYYYMMDD_condition_mouseX_WT_scored_scores_1Hz_crop.csv
%
% OPTIONS (name–value):
%   'pattern' : default '*_scored_scores_1Hz_crop.csv'
%   'codes'   : struct('WK',0,'NREM',1,'REM',2,'MA',15)
%   'includeMA', 'ma_thresh_sec', 'reclassify_short_wake_to_MA'
%   'out_dir' : folder to save group files (default: input_dir/group_overall_ambtemp)
% -------------------------------------------------------------------------

% ----------------- args -----------------
p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*_scored_scores_1Hz_crop.csv',@ischar);
addParameter(p,'codes',struct('WK',0,'NREM',1,'REM',2,'MA',15),@isstruct);
addParameter(p,'includeMA',true,@islogical);
addParameter(p,'ma_thresh_sec',15,@(x)isscalar(x)&&x>0);
addParameter(p,'reclassify_short_wake_to_MA',true,@islogical);
addParameter(p,'out_dir','',@ischar);
parse(p, input_dir, varargin{:});
S = p.Results;

% ---- guard: function must be on path
if exist('sleep_architecture_from_scores','file') ~= 2
    error(['sleep_architecture_from_scores.m is not on the MATLAB path.\n' ...
           'Add it, e.g.: addpath(genpath(''%s''))'], pwd);
end

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);
if isempty(S.out_dir)
    S.out_dir = fullfile(S.input_dir, 'group_overall_ambtemp');
end
if ~isfolder(S.out_dir), mkdir(S.out_dir); end

% ----------------- discover files -----------------
F = dir(fullfile(S.input_dir, '**', S.pattern));
F = F(~[F.isdir]);
assert(~isempty(F), 'No files matched "%s" under %s', S.pattern, S.input_dir);

fprintf('Found %d cropped scoring files.\n', numel(F));

% ----------------- collector schema -----------------
schema_overall_names = {'state','n_bouts','total_dur_s','mean_bout_dur_s', ...
                        'file','date','condition','mouse','genotype'};
schema_overall_types = {'string','double','double','double', ...
                        'string','string','string','string','string'};
rows_overall = table('Size',[0 numel(schema_overall_names)], ...
    'VariableTypes', schema_overall_types, ...
    'VariableNames', schema_overall_names);

% ----------------- loop files -----------------
for i = 1:numel(F)
    csv_path = fullfile(F(i).folder, F(i).name);
    fprintf('Processing %s...\n', F(i).name);

    % Parse date / condition / mouse / genotype from filename
    info = parse_info_from_fname_local(F(i).name);

    if ~info.ok
        warning('  -> Filename not recognized pattern, skipping: %s', F(i).name);
        continue;
    end

    % use filename without suffix as prefix for per-file outputs (if you want logs)
    base_prefix = erase(F(i).name, '_scores_1Hz_crop.csv');

    try
        OUT_i = sleep_architecture_from_scores(csv_path, ...
            'codes', S.codes, ...
            'includeMA', S.includeMA, ...
            'ma_thresh_sec', S.ma_thresh_sec, ...
            'reclassify_short_wake_to_MA', S.reclassify_short_wake_to_MA, ...
            'out_prefix', base_prefix, ...
            'out_dir', S.out_dir, ...
            'write_log', true, ...
            'verbose', false);
    catch ME
        warning('  -> Failed file %s: %s', F(i).name, ME.message);
        continue;
    end

    O = OUT_i.overall;
    if istable(O) && ~isempty(O)
        O2 = table( ...
            string(O.state), ...
            double(O.n_bouts), ...
            double(O.total_dur_s), ...
            double(O.mean_bout_dur_s), ...
            repmat(string(F(i).name),     height(O),1), ...
            repmat(string(info.date),     height(O),1), ...
            repmat(string(info.condition),height(O),1), ...
            repmat(string(info.mouse),    height(O),1), ...
            repmat(string(info.genotype), height(O),1), ...
            'VariableNames', schema_overall_names);
        rows_overall = [rows_overall; O2]; %#ok<AGROW>
    end
end

% ----------------- save & return -----------------
rows_overall_ambtemp = rows_overall;

mat_file  = fullfile(S.out_dir, 'rows_overall_ambtemp.mat');
csv_file  = fullfile(S.out_dir, 'rows_overall_ambtemp.csv');

save(mat_file, 'rows_overall_ambtemp');
if ~isempty(rows_overall_ambtemp)
    writetable(rows_overall_ambtemp, csv_file);
else
    warning('rows_overall_ambtemp is empty; not writing CSV.');
end

fprintf('\n✅ rows_overall_ambtemp created with %d rows.\n', height(rows_overall_ambtemp));
fprintf('   MAT: %s\n', mat_file);
fprintf('   CSV: %s\n', csv_file);
end


% =====================================================================
% Local helper: parse your cropped filename pattern
% Expected pattern (no .csv in BASE because of fileparts):
%   YYYYMMDD_condition_mouseX_APP_scored_scores_1Hz_crop
%   YYYYMMDD_condition_mouseX_WT_scored_scores_1Hz_crop
% =====================================================================
function info = parse_info_from_fname_local(fname)
    [~, base, ~] = fileparts(fname);

    info = struct('date',"", 'condition',"", 'mouse',"", ...
                  'genotype',"", 'ok', false);

    expr = ['^(?<date>\d{8})_' ...         % 8-digit date
            '(?<condition>[^_]+)_' ...     % condition
            '(?<mouse>mouse\d+)_' ...      % mouse + digits
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
