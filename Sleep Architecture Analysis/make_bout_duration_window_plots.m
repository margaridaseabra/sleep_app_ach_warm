function make_bout_duration_window_plots(PERHOUR2, out_dir, COL_WT, COL_APP, OPT)
% MAKE_BOUT_DURATION_WINDOW_PLOTS
% -------------------------------------------------------------------------
% For each state (WK / NREM / REM if present), make ONE figure:
%
%   x-axis: conditions (baseline, ambtemp, drugs)
%   per condition: 6 bars, in this order (tightly grouped):
%       WT 0–3 h, APP 0–3 h, WT 3–6 h, APP 3–6 h, WT >6 h, APP >6 h
%
% Bars = mean ± SEM of mean bout duration (s) per mouse.
%
% Visuals:
%   - Bars within each condition are contiguous; gaps only between conditions.
%   - On each bar, a small label: "0–3 h", "3–6 h", or ">6 h".
%   - For ambtemp and drugs, the 3–6 h bars are shaded, with "manipulation".
%
% Stats:
%   - For each state + condition, 2-way RM ANOVA (Window (within) × Genotype (between)).
%   - Above each condition group:  G:*, Win:n.s., Int:**  (using p_to_stars).
% -------------------------------------------------------------------------

PH = PERHOUR2;

% Optional: silence "ModifiedAndSavedVarnames" warning from fitrm/ranova/anova
ws = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
cleanupObj = onCleanup(@() warning(ws));

need = {'hour_idx','dur_s','bouts_per_h','state','condition','mouse','genotype'};
if ~all(ismember(need, PH.Properties.VariableNames))
    warning('make_bout_duration_window_plots: PERHOUR2 missing required columns.');
    return;
end

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---------- define 3 time windows from hour_idx (0-based) ---------------
h = double(PH.hour_idx);
win_idx = zeros(size(h));           % 1,2,3
win_idx(h <= 2)          = 1;       % 0–3 h
win_idx(h >= 3 & h <= 5) = 2;       % 3–6 h
win_idx(h >= 6)          = 3;       % >6 h
PH.win = win_idx;

PH = PH(PH.win>=1 & PH.win<=3, :);
if isempty(PH)
    warning('make_bout_duration_window_plots: no data in 0–3/3–6/>6 h windows.');
    return;
end

% ---------- normalize text columns --------------------------------------
PH.state     = string(PH.state);
PH.condition = string(PH.condition);
PH.mouse     = string(PH.mouse);
PH.genotype  = string(PH.genotype);

% collapse genotype: everything not WT is APP
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% ---------- per-mouse mean bout duration per state×cond×geno×window -----
[gid, st, cond, geno_g, mouse, win] = findgroups( ...
    PH.state, PH.condition, PH.geno_group, PH.mouse, PH.win);

tot_dur   = splitapply(@(x) sum(x,'omitnan'), double(PH.dur_s),      gid);
tot_bouts = splitapply(@(x) sum(x,'omitnan'), double(PH.bouts_per_h),gid);

mean_bout = tot_dur ./ max(tot_bouts,1);
mean_bout(tot_bouts == 0) = NaN;

Tmouse = table(st, cond, geno_g, mouse, win, mean_bout, ...
    'VariableNames',{'state','condition','geno','mouse','win','mean_bout_dur_win_s'});

% ---------- states & conditions of interest -----------------------------
state_pref = ["WK","NREM","REM"];
cond_pref  = ["baseline","ambtemp","drugs"];

states = state_pref(ismember(state_pref, unique(Tmouse.state)));
conds  = cond_pref(ismember(cond_pref, unique(Tmouse.condition)));

if isempty(states)
    warning('make_bout_duration_window_plots: no WK/NREM/REM in Tmouse.');
    return;
end

nBarsPerCond = 6;                        % WT/APP × 3 windows
win_labels   = ["0–3 h","3–6 h",">6 h"]; % for bar text

for si = 1:numel(states)
    st_name = states(si);
    Ts = Tmouse(Tmouse.state == st_name, :);
    if isempty(Ts), continue; end

    % conditions actually present for this state
    conds_st = conds(ismember(conds, unique(Ts.condition)));
    if isempty(conds_st), continue; end

    nC = numel(conds_st);
    N  = nC*nBarsPerCond;

    means   = NaN(N,1);
    sems    = NaN(N,1);
    xCond   = strings(N,1);
    xGeno   = strings(N,1);
    xWin    = NaN(N,1);

    % ---- manual x positions: tight within condition, big gaps between ----
    group_size = nBarsPerCond;   % 6 bars per condition
    group_gap  = 4;              % increase for more space between conditions
    x          = NaN(N,1);

    k = 0;
    for ci = 1:nC
        cnd = conds_st(ci);

        % leftmost x for this condition block
        base = (ci-1)*(group_size + group_gap);

        % j = 1..6 inside the condition
        for w = 1:3                % window 1..3
            for gname = ["WT","APP"]   % WT then APP
                k = k + 1;

                j = k - (ci-1)*group_size;   % 1..6 within this condition
                x(k) = base + j;             % contiguous bars inside block

                mask = Ts.condition == cnd & Ts.win == w & Ts.geno == gname;
                vals = Ts.mean_bout_dur_win_s(mask);

                if ~isempty(vals)
                    m   = mean(vals, 'omitnan');
                    sd  = std(vals,  'omitnan');
                    n   = sum(~isnan(vals));
                    sem = sd ./ max(sqrt(n),1);   % *** SEM, not SD ***

                    means(k) = m;
                    sems(k)  = sem;
                end

                xCond(k) = cnd;
                xGeno(k) = gname;
                xWin(k)  = w;
            end
        end
    end

    if all(isnan(means))
        warning('No bout-duration data for state %s', st_name);
        continue;
    end

    % ---------- plotting -------------------------------------------------
    f = figure('Color','w','Name', sprintf('Bout duration windows — %s', st_name), ...
               'Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    hold on

    % bars (within-condition bars tightly packed because x increments by 1)
    for i = 1:N
        if isnan(means(i)), continue; end
        if strcmpi(xGeno(i),'WT')
            fc = COL_WT;
        else
            fc = COL_APP;
        end
        bar(x(i), means(i), 0.8, 'FaceColor', fc, 'EdgeColor','none');
    end

    % error bars
    errorbar(x, means, sems, 'k.', 'LineStyle','none','LineWidth',1);

    % y limits
    yl   = ylim;
    yTop = yl(2);
    ylim([0 yTop*1.2]);   % extra room for labels & stats
    yl   = ylim;
    yTop = yl(2);

    % per-bar time-window label ("0–3 h", "3–6 h", ">6 h") near bottom
    yLabel = yTop * 0.03;   % 3% up from zero
    for i = 1:N
        if isnan(means(i)), continue; end
        if isnan(xWin(i)),  continue; end
        w = xWin(i);
        if w>=1 && w<=3
            txtw = win_labels(w);
            text(x(i), yLabel, txtw, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',8);
        end
    end

    % shade 3–6 h window for ambtemp & drugs (+ label "manipulation")
    for ci = 1:nC
        cnd = conds_st(ci);
        if ~(cnd=="ambtemp" || cnd=="drugs"), continue; end

        % indices of this condition in the global vectors
        idx_cond = find(xCond == cnd);
        if numel(idx_cond) ~= nBarsPerCond, continue; end

        % indices for window 2 (3–6 h): within group these are bars 3 & 4
        idx_w2 = idx_cond(3:4);
        x1 = min(x(idx_w2)) - 0.4;
        x2 = max(x(idx_w2)) + 0.4;

        patch([x1 x2 x2 x1], [0 0 yTop*0.9 yTop*0.9], ...
              [0.7 0.7 0.7], 'FaceAlpha',0.6, 'EdgeColor','none');

        % label "manipulation" roughly in middle of patch
        xm = (x1 + x2)/2;
        ym = yTop*0.9;
        text(xm, ym*0.97, 'manipulation', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontSize',8, 'FontAngle','italic');
    end

    % condition labels under group centers
    for ci = 1:nC
        cnd = conds_st(ci);
        idx_cond = find(xCond == cnd);
        xc = mean(x(idx_cond));
        text(xc, -0.03*yTop, upper(char(cnd)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold','FontSize',10, 'Units','data');
    end

    % ticks only at condition centers, not at each bar
    cond_centers = zeros(nC,1);
    for ci = 1:nC
        cnd = conds_st(ci);
        idx_cond = find(xCond == cnd);
        cond_centers(ci) = mean(x(idx_cond));
    end
    set(gca,'XTick',cond_centers,'XTickLabel',upper(string(conds_st)));
    xtickangle(0);

    ylabel('Mean bout duration (s)');

    % pretty state label
    if st_name=="WK"
        st_label = 'Wake';
    elseif st_name=="NREM"
        st_label = 'NREM';
    elseif st_name=="REM"
        st_label = 'REM';
    elseif st_label== "MA"
        st_label = 'MA';
    else
        st_label = char(st_name);
    end

    title(sprintf('Mean bout duration in 3 time windows — %s', st_label));
    box off; grid on

    % ---------- per-condition RM ANOVA (Window × Genotype) overlay -------
    for ci = 1:nC
        cnd = conds_st(ci);

        % subtable for this state+condition
        mask_sc = Ts.condition == cnd & ~isnan(Ts.mean_bout_dur_win_s);
        Tc = Ts(mask_sc, :);
        if height(Tc) < 3, continue; end

        % wide table: rows = mouse, cols = windows
        try
            Tw = unstack(Tc, 'mean_bout_dur_win_s', 'win');  % columns like '1','2','3'
        catch
            continue;
        end

        % get measurement vars (window columns)
        measVars = Tw.Properties.VariableNames;
        measVars = measVars(~ismember(measVars, {'state','condition','mouse','geno'}));
        if numel(measVars) < 2, continue; end

        % remove rows where all window measurements are NaN
        Tw = rmmissing(Tw, 'DataVariables', measVars, 'MinNumMissing', numel(measVars));
        if height(Tw) < 3, continue; end

        % fit repeated-measures model: Window (within) × geno (between)
        formula = sprintf('%s-%s ~ geno', measVars{1}, measVars{end});
        WithinDesign = table((1:numel(measVars))','VariableNames',{'Window'});

        try
            rm = fitrm(Tw, formula, 'WithinDesign', WithinDesign);
        catch
            continue;
        end

        % within-subject: Window main effect & Window×geno interaction
        rtbl = ranova(rm, 'WithinModel','Window');
        rn = rtbl.Properties.RowNames;
        pW = NaN; pInt = NaN;
        rowW   = strcmp(rn,'Window');
        rowInt = strcmp(rn,'Window:geno');
        if any(rowW),   pW   = rtbl.pValue(rowW);   end
        if any(rowInt), pInt = rtbl.pValue(rowInt); end

        % between-subject: geno main effect
        bt = anova(rm);
        pG = NaN;
        if ismember('Term', bt.Properties.VariableNames)
            term = bt.Term;
            if iscellstr(term) || isstring(term)
                rowG = strcmp(string(term),'geno');
            else
                rowG = false(size(term));
            end
            if any(rowG) && ismember('pValue', bt.Properties.VariableNames)
                pG = bt.pValue(rowG);
            end
        else
            rn_bt = bt.Properties.RowNames;
            rowG = strcmp(rn_bt,'geno');
            if any(rowG) && ismember('pValue', bt.Properties.VariableNames)
                pG = bt.pValue(rowG);
            end
        end

        if all(isnan([pG pW pInt])), continue; end

        % where to write stats: center above this condition group
        idx_cond = find(xCond == cnd);
        xc_group = mean(x(idx_cond));
        yl = ylim; yTXT = yl(2)*1.05;

        txt_stats = sprintf('G:%s  Win:%s  Int:%s', ...
            p_to_stars(pG), p_to_stars(pW), p_to_stars(pInt));

        text(xc_group, yTXT, txt_stats, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',8);
    end

    % ---------- save figure ---------------------------------------------
    fname = sprintf('boutdur_3windows_%s.png', lower(char(st_name)));
    saveas(f, fullfile(out_dir, fname));
end
end
