function csv_path = export_state_survival_raw_forPrism(ALL_BOUTS, out_dir, varargin)
% export_state_survival_raw_forPrism
% -------------------------------------------------------------------------
% Build a very simple Prism-ready CSV for survival analysis:
%
%   Columns: WT, APP
%   Rows   : one bout per row (time to event, in seconds)
%
% No censoring column, no extra info.
% Perfect for Prism's "enter time to event for each subject" survival table.
%
% INPUT
%   ALL_BOUTS : table with at least columns:
%               state, geno, cond, dur_s
%               (e.g. from build_state_bout_table_from_scores)
%   out_dir   : folder to save CSV
%
% NAME–VALUE OPTIONS
%   'baseline_cond' : which condition to keep (default: "baseline")
%   'state_label'   : which state to keep (e.g. "REM" or "NREM")
%   'filename_tag'  : extra tag for filename (default: state_label)
%
% OUTPUT
%   csv_path : full path to saved CSV file
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

p = inputParser;
addParameter(p,'baseline_cond',"baseline",@(s)ischar(s) || isstring(s));
addParameter(p,'state_label',"REM",@(s)ischar(s) || isstring(s));
addParameter(p,'filename_tag',"",@(s)ischar(s) || isstring(s));
parse(p, varargin{:});

baseline_cond = string(p.Results.baseline_cond);
state_label   = string(p.Results.state_label);
filename_tag  = string(p.Results.filename_tag);

if filename_tag == ""
    filename_tag = state_label;
end

% ---- Basic checks ----
needed = {'state','geno','cond','dur_s'};
if ~all(ismember(needed, ALL_BOUTS.Properties.VariableNames))
    error('ALL_BOUTS must contain columns: state, geno, cond, dur_s');
end

T = ALL_BOUTS;

% Normalize to string
T.state = string(T.state);
T.geno  = string(T.geno);
T.cond  = string(T.cond);

% ---- Filter by condition and state ----
mask = (T.cond == baseline_cond) & (T.state == state_label);
Tsub = T(mask, :);

if isempty(Tsub)
    warning('No bouts found for cond="%s", state="%s".', ...
            baseline_cond, state_label);
    csv_path = "";
    return;
end

% ---- Extract WT vs APP durations ----
dur_WT  = Tsub.dur_s(Tsub.geno == "WT");
dur_APP = Tsub.dur_s(Tsub.geno == "APP");

dur_WT  = dur_WT(~isnan(dur_WT) & dur_WT > 0);
dur_APP = dur_APP(~isnan(dur_APP) & dur_APP > 0);

nWT  = numel(dur_WT);
nAPP = numel(dur_APP);

if nWT == 0 && nAPP == 0
    warning('No valid durations >0 for cond="%s", state="%s".', ...
            baseline_cond, state_label);
    csv_path = "";
    return;
end

% ---- Build Prism-style table: columns WT, APP ----
nRows = max(nWT, nAPP);
WT_col  = NaN(nRows,1);
APP_col = NaN(nRows,1);

if nWT > 0
    WT_col(1:nWT) = dur_WT;
end
if nAPP > 0
    APP_col(1:nAPP) = dur_APP;
end

Tprism = table(WT_col, APP_col, 'VariableNames', {'WT','APP'});

% ---- Save CSV ----
fname    = sprintf('survival_raw_%s_%s_forPrism.csv', ...
                   char(baseline_cond), char(filename_tag));
csv_path = fullfile(out_dir, fname);

writetable(Tprism, csv_path);
fprintf('✅ Prism raw survival CSV saved to: %s\n', csv_path);
fprintf('   WT bouts: %d, APP bouts: %d\n', nWT, nAPP);
end
