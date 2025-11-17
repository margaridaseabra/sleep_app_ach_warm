function ANOVA = run_sleep_arch_anova_two_way(OVERALL2, out_dir)
% RUN_SLEEP_ARCH_ANOVA_TWO_WAY
% -------------------------------------------------------------------------
% Two-way ANOVA on % time spent in each state.
%
% Factors:
%   - Genotype: WT vs APP
%   - Condition: baseline / ambtemp / drugs
%
% Done separately for each state (WK, NREM, REM, MA if present).
% Uses % of recording per (mouse × condition × genotype group) as response.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

ANOVA = struct();

T = OVERALL2;

if ~ismember('total_dur_s', T.Properties.VariableNames)
    warning('run_sleep_arch_anova_two_way: total_dur_s not found in OVERALL2.');
    return;
end

% Normalize to string
if ismember('mouse', T.Properties.VariableNames)
    subj = string(T.mouse);
elseif ismember('file', T.Properties.VariableNames)
    subj = string(T.file);
else
    subj = string((1:height(T))');
end

T.state     = string(T.state);
T.condition = string(T.condition);
T.genotype  = string(T.genotype);

% Collapse genotype: everything not WT is APP
geno = T.genotype;
geno(geno ~= "WT") = "APP";
T.geno_group = geno;

% --- Compute % of recording per row (same as in plotting helper) -------
[grp, ~, ~, ~] = findgroups(subj, T.condition, T.geno_group);
tot_s = splitapply(@(x) nansum(x), T.total_dur_s, grp);
T.rec_total_s = tot_s(grp);

T.rec_total_s(T.rec_total_s<=0 | isnan(T.rec_total_s)) = NaN;
T.pct_rec = 100 * T.total_dur_s ./ T.rec_total_s;

state_pref    = ["WK","NREM","REM","MA"];
presentStates = state_pref(ismember(state_pref, unique(T.state)));

cond_order = ["baseline","ambtemp","drugs"];

logfile = fullfile(out_dir, 'sleep_arch_two_way_anova.log.txt');
fid = fopen(logfile, 'w');
if fid>0
    fprintf(fid, 'Two-way ANOVA on %% time in each state\n');
    fprintf(fid, 'Factors: Genotype (WT vs APP) × Condition (baseline/ambtemp/drugs)\n\n');
end

for s = 1:numel(presentStates)
    st = presentStates(s);
    Ts = T(T.state == st & ~isnan(T.pct_rec), :);
    if isempty(Ts), continue; end

    % Restrict to baseline/ambtemp/drugs
    presentConds = unique(Ts.condition, 'stable');
    conds = cond_order(ismember(cond_order, presentConds));
    if numel(conds) < 2
        continue;
    end
    Ts = Ts(ismember(Ts.condition, conds), :);

    y      = double(Ts.pct_rec);
    geno_f = cellstr(Ts.geno_group);
    cond_f = cellstr(Ts.condition);

    [p, tbl, stats] = anovan(y, {geno_f, cond_f}, ...
                             'model',    'interaction', ...
                             'varnames', {'Genotype','Condition'}, ...
                             'display',  'off');

    st_key = lower(char(st));  % 'wk','nrem','rem',...
    ANOVA.(st_key).p        = p;
    ANOVA.(st_key).tbl      = tbl;
    ANOVA.(st_key).stats    = stats;
    ANOVA.(st_key).conds    = conds;
    ANOVA.(st_key).statestr = string(st);

    pg = p(1);  % Genotype
    pc = p(2);  % Condition
    pi = p(3);  % Interaction

    if fid>0
        fprintf(fid, 'State: %s\n', st_key);
        fprintf(fid, '  p(Genotype)    = %.4g\n', pg);
        fprintf(fid, '  p(Condition)   = %.4g\n', pc);
        fprintf(fid, '  p(Interaction) = %.4g\n\n', pi);
    end
end

if fid>0
    fclose(fid);
    fprintf('📄 Two-way ANOVA summary written to: %s\n', logfile);
end
end
