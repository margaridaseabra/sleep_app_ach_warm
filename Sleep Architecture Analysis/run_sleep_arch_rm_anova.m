function STATS = run_sleep_arch_rm_anova(OVERALL2, out_dir)
% RUN_SLEEP_ARCH_RM_ANOVA
% -------------------------------------------------------------------------
% Repeated-measures 2-way ANOVA on % time in each state.
%
% Design:
%   Subject = mouse
%   Within-subject factor: Condition (baseline / ambtemp / drugs)
%   Between-subject factor: Genotype (WT vs APP)
%
% Each observation is ONE mouse × condition (pooled across all files).
%
% For each state (WK, NREM, REM, MA if present):
%   1) Pool total seconds per mouse × condition × state.
%   2) Compute % of total recording time for that mouse × condition.
%   3) Build wide table: baseline, ambtemp, drugs as repeated measures.
%   4) Fit repeated-measures model: fitrm + ranova.
%   5) Post-hoc Sidak (Dunn–Sidak) multiple comparisons:
%        - Condition main effect
%        - Condition within genotype
%        - Genotype within condition
%
% RETURNS
%   STATS.<state> struct with fields:
%       .rm, .ranova, .anova_between, .post_cond,
%       .post_cond_by_geno, .post_geno_by_cond, .Twide, .conds
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

STATS = struct();

T = OVERALL2;

if ~ismember('total_dur_s', T.Properties.VariableNames)
    warning('run_sleep_arch_rm_anova: total_dur_s not found in OVERALL2.');
    return;
end

% --------- Normalize text variables -------------------------------------
mouse     = string(T.mouse);
cond      = string(T.condition);
state     = string(T.state);
genotype  = string(T.genotype);

% Collapse genotype: everything not WT is APP
geno_group = genotype;
geno_group(geno_group ~= "WT") = "APP";

% --------- Aggregate per mouse × condition × state × genotype ----------
[Gc, uMouse, uCond, uState, uGeno] = findgroups(mouse, cond, state, geno_group);
state_tot_s = splitapply(@(x) nansum(x), T.total_dur_s, Gc);

Agg = table(uMouse, uCond, uState, uGeno, state_tot_s, ...
    'VariableNames', {'mouse','condition','state','geno_group','state_tot_s'});

% Total recording duration per mouse × condition (sum over states)
[Gr, uMouseR, uCondR] = findgroups(Agg.mouse, Agg.condition);
rec_total_s = splitapply(@(x) nansum(x), Agg.state_tot_s, Gr);

Rec = table(uMouseR, uCondR, rec_total_s, ...
    'VariableNames', {'mouse','condition','rec_total_s'});

% Join back so each state row knows total recording
Agg = outerjoin(Agg, Rec, 'Keys',{'mouse','condition'}, 'MergeKeys',true);

Agg.pct_rec = 100 * Agg.state_tot_s ./ Agg.rec_total_s;
Agg.rec_total_s(Agg.rec_total_s<=0 | isnan(Agg.rec_total_s)) = NaN;

% Clean up types
Agg.state      = string(Agg.state);
Agg.condition  = string(Agg.condition);
Agg.geno_group = string(Agg.geno_group);

% --------- State and condition order ------------------------------------
state_pref    = ["WK","NREM","REM","MA"];
presentStates = state_pref(ismember(state_pref, unique(Agg.state)));

cond_order    = ["baseline","ambtemp","drugs"];

logfile = fullfile(out_dir, 'sleep_arch_rm_two_way_anova.log.txt');
fid = fopen(logfile, 'w');
if fid>0
    fprintf(fid, 'Repeated-measures 2-way ANOVA on %% time in each state\n');
    fprintf(fid, 'Subject = mouse | Within: Condition | Between: Genotype\n');
    fprintf(fid, 'Conditions: baseline, ambtemp, drugs (where available)\n\n');
end

for s = 1:numel(presentStates)
    st = presentStates(s);
    mask_st = Agg.state == st;
    Ts = Agg(mask_st & ~isnan(Agg.pct_rec), :);
    if isempty(Ts), continue; end

    % Restrict to our main 3 conditions, in order
    presentConds = unique(Ts.condition,'stable');
    conds = cond_order(ismember(cond_order, presentConds));
    if numel(conds) < 2
        % Need at least 2 repeated levels
        continue;
    end
    Ts = Ts(ismember(Ts.condition, conds), :);

    % ---------- Wide table: one row per mouse ---------------------------
    Tstate = Ts(:, {'mouse','condition','geno_group','pct_rec'});
    Tstate.geno_group = categorical(Tstate.geno_group);
    Tstate.condition  = categorical(Tstate.condition);

    % unstack: columns = conditions (baseline/ambtemp/drugs present in this state)
    Twide = unstack(Tstate, 'pct_rec', 'condition');
    % Twide now has columns: mouse, geno_group, <cond1>, <cond2>, ...

    % Keep only rows with genotype info
    Twide = Twide(~isundefined(Twide.geno_group), :);

    % Determine actual repeated-measure columns from BOTH cond_order and Twide
    allCols = Twide.Properties.VariableNames;
    % intersection in the desired order (baseline, ambtemp, drugs)
    rm_cond_cols = cellstr(intersect(cond_order, string(allCols), 'stable'));

    % We need at least 2 repeated columns to fit RM-ANOVA
    if numel(rm_cond_cols) < 2
        continue;
    end

    % Build formula like 'baseline,ambtemp,drugs ~ geno_group'
    respList = strjoin(rm_cond_cols, ',');
    formula = sprintf('%s ~ geno_group', respList);

    % Within-subject design table: one row per repeated column
    WithinDesign = table(categorical(rm_cond_cols.'), ...
                         'VariableNames', {'Condition'});

    % Fit repeated-measures model (now #rows(WithinDesign) matches #responses)
    rm = fitrm(Twide, formula, 'WithinDesign', WithinDesign);

    % Within-subject ANOVA (Condition + interaction with Genotype)
    ranovatbl = ranova(rm, 'WithinModel', 'Condition');

    % Between-subject effect (Genotype)
    anovatbl_between = anova(rm);

    % Post-hoc (Dunn–Sidak == Sidak)
    post_cond = multcompare(rm, 'Condition', ...
                            'ComparisonType', 'dunn-sidak');

    post_cond_by_geno = multcompare(rm, 'Condition', ...
                                    'By','geno_group', ...
                                    'ComparisonType','dunn-sidak');

    post_geno_by_cond = multcompare(rm, 'geno_group', ...
                                    'By','Condition', ...
                                    'ComparisonType','dunn-sidak');

    st_key = lower(char(st));

    STATS.(st_key).rm                = rm;
    STATS.(st_key).ranova            = ranovatbl;
    STATS.(st_key).anova_between     = anovatbl_between;
    STATS.(st_key).post_cond         = post_cond;
    STATS.(st_key).post_cond_by_geno = post_cond_by_geno;
    STATS.(st_key).post_geno_by_cond = post_geno_by_cond;
    STATS.(st_key).Twide             = Twide;
    STATS.(st_key).conds             = rm_cond_cols;

    if fid>0
        fprintf(fid, 'State: %s\n', st_key);
        fprintf(fid, '  Repeated measures cols: %s\n', strjoin(rm_cond_cols, ', '));
        fprintf(fid, '  (see ranova / multcompare outputs in STATS.%s)\n\n', st_key);
    end
end

if fid>0
    fclose(fid);
    fprintf('📄 Repeated-measures ANOVA summary written to: %s\n', logfile);
end
end
