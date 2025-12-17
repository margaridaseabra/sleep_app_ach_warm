function OUT = extract_bouts_for_prism_v3(input_dir, varargin)
% EXTRACT_BOUTS_FOR_PRISM_V3 - FINAL CORRECTED VERSION
% Outputs format for Prism SURVIVAL tables (each column = one group)
%
% Usage:
%   OUT = extract_bouts_for_prism_v3('/path/to/data', ...
%       'state', 'REM', ...
%       'condition', 'baseline', ...
%       'output', 'REM_baseline.csv');

p = inputParser;
addRequired(p, 'input_dir', @ischar);
addParameter(p, 'state', 'REM', @ischar);
addParameter(p, 'condition', '', @ischar);
addParameter(p, 'codes', struct('WK',0,'NREM',1,'REM',2,'MA',15), @isstruct);
addParameter(p, 'pattern', '*_scores_1Hz.csv', @ischar);
addParameter(p, 'output', 'bouts_for_prism.csv', @ischar);
addParameter(p, 'time_unit', 'sec', @(x)ismember(x,{'sec','min'}));
parse(p, input_dir, varargin{:});
S = p.Results;

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);

% Find files
F = dir(fullfile(S.input_dir, '**', S.pattern));
F = F(~[F.isdir]);
assert(~isempty(F), 'No files matched "%s"', S.pattern);

fprintf('Found %d files\n', numel(F));

% Get state code
state_code = S.codes.(S.state);

% Collect bouts by genotype
all_durations_APP = [];
all_durations_WT = [];

for i = 1:numel(F)
    csv_path = fullfile(F(i).folder, F(i).name);
    info = parse_filename_local(F(i).name);
    
    if ~info.ok || (~isempty(S.condition) && ~strcmp(info.condition, S.condition))
        continue;
    end
    
    try
        T = readtable(csv_path);
    catch
        warning('Failed to read %s', F(i).name);
        continue;
    end
    
    score_col = find_score_column(T);
    if isempty(score_col), continue; end
    
    scores = T.(score_col);
    bouts = detect_bouts(scores, state_code);
    
    if isempty(bouts), continue; end
    
    durations = bouts(:,3); % duration in seconds
    if strcmp(S.time_unit, 'min')
        durations = durations / 60;
    end
    
    fprintf('  %s: %d %s bouts\n', F(i).name, length(durations), S.state);
    
    % Separate by genotype
    if strcmp(info.genotype, 'APP')
        all_durations_APP = [all_durations_APP; durations]; %#ok<AGROW>
    else
        all_durations_WT = [all_durations_WT; durations]; %#ok<AGROW>
    end
end

assert(~isempty(all_durations_APP) || ~isempty(all_durations_WT), 'No bouts extracted');

fprintf('\nTotal bouts:\n');
fprintf('  APP: %d\n', length(all_durations_APP));
fprintf('  WT: %d\n', length(all_durations_WT));

% Create output in PRISM SURVIVAL format
% Each column = one group, each row = one bout duration
max_n = max(length(all_durations_APP), length(all_durations_WT));

% Pad with NaN (Prism will ignore empty cells)
APP_padded = [all_durations_APP; NaN(max_n - length(all_durations_APP), 1)];
WT_padded = [all_durations_WT; NaN(max_n - length(all_durations_WT), 1)];

% Create table with APP and WT columns
OUT = table(APP_padded, WT_padded, 'VariableNames', {'APP', 'WT'});

% Write to CSV
writetable(OUT, S.output);

fprintf('\n✓ Output saved: %s\n', S.output);
fprintf('\n=== PRISM IMPORT INSTRUCTIONS ===\n');
fprintf('1. In Prism: NEW TABLE & GRAPH → XY → Survival\n');
fprintf('2. Choose: "Enter or import data into a new table"\n');
fprintf('3. Choose: "Enter elapsed time as number of days (or months...)"\n');
fprintf('4. Click CREATE\n');
fprintf('5. In the table: Data menu → Import from Excel/CSV\n');
fprintf('6. Select: %s\n', S.output);
fprintf('7. Prism will import APP and WT as separate groups\n');
fprintf('8. Analyze → Survival → Kaplan Meier\n');
fprintf('==================================\n\n');

% Print summary stats
fprintf('Summary statistics (duration in %s):\n', S.time_unit);
fprintf('  APP:\n');
fprintf('    n = %d\n', length(all_durations_APP));
fprintf('    mean = %.2f %s\n', mean(all_durations_APP), S.time_unit);
fprintf('    median = %.2f %s\n', median(all_durations_APP), S.time_unit);
fprintf('    min = %.2f %s\n', min(all_durations_APP), S.time_unit);
fprintf('    max = %.2f %s\n', max(all_durations_APP), S.time_unit);
fprintf('  WT:\n');
fprintf('    n = %d\n', length(all_durations_WT));
fprintf('    mean = %.2f %s\n', mean(all_durations_WT), S.time_unit);
fprintf('    median = %.2f %s\n', median(all_durations_WT), S.time_unit);
fprintf('    min = %.2f %s\n', min(all_durations_WT), S.time_unit);
fprintf('    max = %.2f %s\n', max(all_durations_WT), S.time_unit);

end

% =====================================================================
% Helper functions
% =====================================================================

function col = find_score_column(T)
    candidates = {'state', 'State', 'score', 'Score', 'scoring', 'Scoring'};
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            col = candidates{i};
            return;
        end
    end
    for i = 1:width(T)
        if isnumeric(T{:,i})
            col = T.Properties.VariableNames{i};
            return;
        end
    end
    col = '';
end

function bouts = detect_bouts(scores, state_code)
    bouts = [];
    in_bout = false;
    start_idx = NaN;
    
    for i = 1:numel(scores)
        if scores(i) == state_code
            if ~in_bout
                in_bout = true;
                start_idx = i;
            end
        else
            if in_bout
                end_idx = i - 1;
                duration = end_idx - start_idx + 1;
                bouts = [bouts; start_idx, end_idx, duration]; %#ok<AGROW>
                in_bout = false;
            end
        end
    end
    
    if in_bout
        end_idx = numel(scores);
        duration = end_idx - start_idx + 1;
        bouts = [bouts; start_idx, end_idx, duration];
    end
end

function info = parse_filename_local(fname)
    [~, base, ~] = fileparts(fname);
    info = struct('date', '', 'condition', '', 'mouse', '', 'genotype', '', 'ok', false);
    
    expr = ['^(?<date>\d{8})_' ...
            '(?<condition>[^_]+)_' ...
            '(?<mouse>mouse\d+)_' ...
            '(?<genotype>APP|WT)' ...
            '(?<scored>_scored)?' ...
            '_scores_1Hz$'];
    
    m = regexp(base, expr, 'names');
    if isempty(m), return; end
    
    info.date = m.date;
    info.condition = m.condition;
    info.mouse = m.mouse;
    info.genotype = m.genotype;
    info.ok = true;
end