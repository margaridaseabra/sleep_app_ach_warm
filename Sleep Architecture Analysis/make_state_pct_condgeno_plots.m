function make_state_pct_condgeno_plots(OVERALL2, out_dir, COL_WT, COL_APP, OPT)
% MAKE_STATE_PCT_CONDGENO_PLOTS
% -------------------------------------------------------------------------
% For each state (WK / NREM / REM / MA if present), make ONE figure with:
%   x-axis: conditions (baseline, ambtemp, drugs)
%   bars:  WT vs APP (2 bars per condition, grouped)
%   y-axis: mean % of recording time in that state
%
% Also, if sleep_arch_rm_two_way_anova.mat exists in out_dir, it reads
% repeated-measures 2-way ANOVA stats (Genotype × Condition) and writes
% a compact summary [Genotype: */** | Condition: n.s. | G×C: *] in the title.
% -------------------------------------------------------------------------

T = OVERALL2;

if ~ismember('total_dur_s', T.Properties.VariableNames)
    warning('make_state_pct_condgeno_plots: total_dur_s not found in OVERALL2.');
    return;
end

% --------- Try to load RM-ANOVA results (optional) ----------------------
STATS = [];
stats_file = fullfile(out_dir, 'sleep_arch_rm_two_way_anova.mat');
if exist(stats_file, 'file')
    S = load(stats_file, 'STATS');
    if isfield(S, 'STATS')
        STATS = S.STATS;
    end
end

% --- Normalize variables & types ----------------------------------------
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

T.state     = categorical(T.state);
T.condition = categorical(T.condition);
T.genotype  = categorical(T.genotype);

% Collapse genotype: everything not WT is APP
geno = string(T.genotype);
geno(geno ~= "WT") = "APP";
T.geno_group = categorical(geno);

% --- Total recording duration per subject × condition × geno_group ------
[grp, ~, ~, ~] = findgroups(subj, string(T.condition), string(T.geno_group));
tot_s = splitapply(@(x) nansum(x), T.total_dur_s, grp);
T.rec_total_s = tot_s(grp);

T.rec_total_s(T.rec_total_s<=0 | isnan(T.rec_total_s)) = NaN;
T.pct_rec = 100 * T.total_dur_s ./ T.rec_total_s;

% --- States & conditions we care about ----------------------------------
state_pref    = {'WK','NREM','REM','MA'};
presentStates = state_pref(ismember(state_pref, cellstr(categories(T.state))));
if isempty(presentStates)
    warning('make_state_pct_condgeno_plots: no WK/NREM/REM states found.');
    return;
end

cond_order = {'baseline','ambtemp','drugs'};
geno_order = {'WT','APP'};

for s = 1:numel(presentStates)
    st = presentStates{s};                  % e.g. 'WK'
    mask_state = (T.state == st);
    Ts = T(mask_state & ~isnan(T.pct_rec), :);
    if isempty(Ts), continue; end

    % Conditions actually present for this state
    presentConds = cellstr(categories(Ts.condition));
    conds = cond_order(ismember(cond_order, presentConds));
    if isempty(conds), continue; end

    nC = numel(conds);
    nG = numel(geno_order);

    M  = NaN(nC, nG);
    SE = NaN(nC, nG);

    for ci = 1:nC
        for gi = 1:nG
            mask = (Ts.condition == conds{ci}) & (Ts.geno_group == geno_order{gi});
            vals = Ts.pct_rec(mask);
            if isempty(vals), continue; end

            m   = mean(vals, 'omitnan');
            sd  = std(vals,  'omitnan');
            n   = sum(~isnan(vals));
            sem = sd ./ max(sqrt(n), 1);   % SEM

            M(ci,gi)  = m;
            SE(ci,gi) = sem;
        end
    end

    if all(isnan(M(:)))
        continue;
    end

    % ----- Make the figure ----------------------------------------------
    fig = figure('Name', sprintf('Pct time in %s by condition×genotype', st), ...
                 'Color','w');
    hold on

    cats = categorical(conds, conds);
    b = bar(cats, M, 'grouped');

    % Colours per genotype
    if nG >= 1, b(1).FaceColor = COL_WT;  end   % WT
    if nG >= 2, b(2).FaceColor = COL_APP; end   % APP

    % Error bars (SEM)
    for gi = 1:nG
        if gi > numel(b), continue; end
        x_e = b(gi).XEndPoints;
        errorbar(x_e, M(:,gi), SE(:,gi), 'k.', 'LineWidth',1);
    end

    % ----- Per-subject dots ---------------------------------------------
    has_out_flag = ismember('total_min_is_outlier', Ts.Properties.VariableNames);
    if has_out_flag
        is_out_all = Ts.total_min_is_outlier;
    else
        is_out_all = false(height(Ts),1);
    end

    for ii = 1:height(Ts)
        cstr = char(string(Ts.condition(ii)));
        gstr = char(string(Ts.geno_group(ii)));

        ci = find(strcmp(conds, cstr), 1);
        gi = find(strcmp(geno_order, gstr), 1);
        if isempty(ci) || isempty(gi) || gi > numel(b)
            continue;
        end

        x0 = b(gi).XEndPoints(ci);
        x  = x0 + (rand-0.5)*2*OPT.jitter;
        y  = Ts.pct_rec(ii);

        if is_out_all(ii)
            plot(x, y, 'r.', 'MarkerSize',10);

            if OPT.label_outliers && ismember('mouse', Ts.Properties.VariableNames)
                raw = char(string(Ts.mouse(ii)));
                num = regexp(raw, '\d+$', 'match', 'once'); % trailing digits
                if ~isempty(num)
                    lbl = ['m' num];
                else
                    lbl = raw;
                end
                text(x, y, [' ' lbl], ...
                    'Color',[0.6 0 0], 'FontSize',7, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom');
            end
        else
            plot(x, y, 'k.', 'MarkerSize',8);
        end
    end

    % ----- Labels & title + RM-ANOVA stars in title ---------------------
    if strcmp(st,'WK')
        st_label = 'Wake';
    elseif strcmp(st,'NREM')
        st_label = 'NREM';
    elseif strcmp(st,'REM')
        st_label = 'REM';
    else
        st_label = st;
    end

    ylim([0 100]);
    ylabel(sprintf('Time in %s (%% of recording)', st_label));
    xlabel('Condition');

    % Base title
    ttl_main = sprintf('%% time in %s – WT (grey) vs APP (cornflower blue)', st_label);
    ttl_anno = '';

    % If RM-ANOVA results exist, pull p-values and convert to stars
    if ~isempty(STATS)
        st_key = lower(st);   % 'wk','nrem','rem',...
        if isfield(STATS, st_key)
            [pG, pC, pI] = get_pvals_from_rmstats(STATS.(st_key));

            if ~isempty(pG)
                starsG = p_to_stars(pG);
            else
                starsG = 'n.s.';
            end
            if ~isempty(pC)
                starsC = p_to_stars(pC);
            else
                starsC = 'n.s.';
            end
            if ~isempty(pI)
                starsI = p_to_stars(pI);
            else
                starsI = 'n.s.';
            end

            ttl_anno = sprintf(' [Genotype: %s | Condition: %s | G×C: %s]', ...
                               starsG, starsC, starsI);
        end
    end

    ttl_full = [ttl_main ttl_anno];   % char row
    title(ttl_full, 'FontWeight','bold');

    legend(geno_order, 'Location','northoutside', ...
           'Orientation','horizontal','Box','off');

    xtickangle(0);
    box off;
    grid on;

    % ----- Save figure ---------------------------------------------------
    fname = sprintf('fig_pcttime_%s_condgeno.png', lower(st));
    saveas(fig, fullfile(out_dir, fname));
end
end

% --------- helper: p → stars --------------------------------------------
function s = p_to_stars(p)
    if ~isscalar(p) || ~isnumeric(p) || isnan(p)
        s = 'n.s.';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end

% --------- helper: extract p-values from STATS.(state) ------------------
function [pG, pC, pI] = get_pvals_from_rmstats(Sst)
% pG : Genotype main effect (between-subjects, from anova(rm))
% pC : Condition main effect (within, from ranova)
% pI : Interaction Condition×Genotype (within, from ranova)

pG = []; pC = []; pI = [];

% --- Genotype main effect (between) -------------------------------------
if isfield(Sst, 'anova_between')
    btbl = Sst.anova_between;
    % Try Term column, then row names
    if ismember('Term', btbl.Properties.VariableNames)
        term = btbl.Term;
        pcol = btbl.pValue;
        for i = 1:numel(term)
            if strcmp(char(term(i)), 'geno_group')
                pG = pcol(i);
                break;
            end
        end
    else
        rn = btbl.Properties.RowNames;
        if ~isempty(rn) && ismember('geno_group', rn)
            idx = find(strcmp(rn, 'geno_group'),1);
            if ismember('pValue', btbl.Properties.VariableNames)
                pG = btbl.pValue(idx);
            end
        end
    end
end

% --- Condition & Interaction (within) -----------------------------------
if isfield(Sst, 'ranova')
    rtbl = Sst.ranova;
    pcol_name = '';
    if ismember('pValue', rtbl.Properties.VariableNames)
        pcol_name = 'pValue';
    elseif ismember('pValueGG', rtbl.Properties.VariableNames)
        pcol_name = 'pValueGG';
    end

    if ~isempty(pcol_name)
        pcol = rtbl.(pcol_name);

        % row names like 'Condition' and 'Condition:geno_group'
        rn = rtbl.Properties.RowNames;

        if ~isempty(rn)
            % Condition main effect
            idxC = find(strcmp(rn, 'Condition'), 1);
            if ~isempty(idxC)
                pC = pcol(idxC);
            end

            % Interaction row often 'Condition:geno_group'
            idxI = find(~cellfun(@isempty, strfind(rn, 'Condition:geno_group')), 1);
            if ~isempty(idxI)
                pI = pcol(idxI);
            end
        end
    end
end
end
