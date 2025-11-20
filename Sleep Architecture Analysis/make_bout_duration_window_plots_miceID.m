function make_bout_duration_window_plots_miceID(PERHOUR, out_dir, COL_WT, COL_APP, OPT)
% MAKE_BOUT_DURATION_WINDOW_PLOTS
% -------------------------------------------------------------------------
% For each state (WK, NREM, REM), plot:
%   "Mean bout duration in 3 time windows" (0–3 h, 3–6 h, >6 h)
% for each condition (baseline, ambtemp, drugs).
%
% Bars: WT (grey) vs APP (cornflower blue) per window
% Dots: each mouse × window, WITH MOUSE ID LABEL on every dot.
%
% Uses PERHOUR table: hour_idx, dur_s, bouts_per_h, state, condition,
% mouse, genotype.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

if nargin < 3 || isempty(COL_WT)
    COL_WT = [0.6 0.6 0.6];
end
if nargin < 4 || isempty(COL_APP)
    COL_APP = [0.392 0.584 0.929];
end
if nargin < 5 || isempty(OPT)
    OPT.label_outliers = true;
    OPT.jitter = 0.12;
end

PH = PERHOUR;

needed = {'hour_idx','dur_s','bouts_per_h','state','condition','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    warning('make_bout_duration_window_plots: PERHOUR missing required columns, skipping.');
    return;
end

% ---- define windows: 0–3h, 3–6h, >6h (same as run_bout_window_stats) ---
h = double(PH.hour_idx);
win_idx = zeros(size(h));
win_idx(h <= 2)          = 1;   % 0–3 h
win_idx(h >= 3 & h <= 5) = 2;   % 3–6 h
win_idx(h >= 6)          = 3;   % >6 h
PH.win = win_idx;

PH = PH(PH.win>=1 & PH.win<=3, :);
if isempty(PH)
    warning('make_bout_duration_window_plots: no data in 0–3/3–6/>6 h windows.');
    return;
end

% Normalise text columns
PH.state     = string(PH.state);
PH.condition = string(PH.condition);
PH.mouse     = string(PH.mouse);
PH.genotype  = string(PH.genotype);

% Collapse all non-WT genotypes to APP
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno = geno;

% ---- per-mouse mean bout duration per state×cond×geno×window ----------
[gid, st, cond, geno_g, mouse, win] = findgroups( ...
    PH.state, PH.condition, PH.geno, PH.mouse, PH.win);

tot_dur   = splitapply(@(x) sum(x,'omitnan'), double(PH.dur_s),      gid);
tot_bouts = splitapply(@(x) sum(x,'omitnan'), double(PH.bouts_per_h),gid);

mean_bout = tot_dur ./ max(tot_bouts,1);
mean_bout(tot_bouts==0) = NaN;

T = table(st, cond, geno_g, mouse, win, mean_bout, ...
    'VariableNames',{'state','condition','geno','mouse','win','mean_bout_dur_win_s'});

% ---- which states / conditions exist? ----------------------------------
state_pref = ["WK","NREM","REM","MA"];
cond_pref  = ["baseline","ambtemp","drugs"];

states = state_pref(ismember(state_pref, unique(T.state)));
conds  = cond_pref(ismember(cond_pref, unique(T.condition)));

if isempty(states) || isempty(conds)
    warning('make_bout_duration_window_plots: no matching states/conditions to plot.');
    return;
end

win_labels = {'0–3 h','3–6 h (manip)','>6 h'};

% ================= main loop: one figure per state ======================
for si = 1:numel(states)
    st_name = states(si);

    Tst = T(T.state == st_name & ~isnan(T.mean_bout_dur_win_s), :);
    if isempty(Tst), continue; end

    conds_st = conds(ismember(conds, unique(Tst.condition)));
    if isempty(conds_st), continue; end

    % One figure per state
    fig = figure('Color','w', ...
        'Name', sprintf('Mean bout duration in 3 time windows — %s', st_name));

    nC = numel(conds_st);
    tiledlayout(1, nC, 'TileSpacing','compact', 'Padding','compact');

    for ci = 1:nC
        cnd = conds_st(ci);
        Tc = Tst(Tst.condition == cnd, :);
        if isempty(Tc), continue; end

        nexttile; hold on;

        % ---- bar means per window × genotype --------------------------
        muWT  = nan(1,3); seWT  = nan(1,3);
        muAPP = nan(1,3); seAPP = nan(1,3);

        for w = 1:3
            rowsW_WT  = (Tc.win==w & Tc.geno=="WT");
            rowsW_APP = (Tc.win==w & Tc.geno=="APP");

            valsWT  = Tc.mean_bout_dur_win_s(rowsW_WT);
            valsAPP = Tc.mean_bout_dur_win_s(rowsW_APP);

            if ~isempty(valsWT) && any(~isnan(valsWT))
                nWT = sum(~isnan(valsWT));
                muWT(w) = mean(valsWT,'omitnan');
                seWT(w) = std(valsWT,'omitnan') / max(sqrt(nWT),1);
            end
            if ~isempty(valsAPP) && any(~isnan(valsAPP))
                nAPP = sum(~isnan(valsAPP));
                muAPP(w) = mean(valsAPP,'omitnan');
                seAPP(w) = std(valsAPP,'omitnan') / max(sqrt(nAPP),1);
            end
        end

        % x positions: WT slightly left, APP slightly right of each integer window
        offs = 0.18;
        for w = 1:3
            xWT  = w - offs;
            xAPP = w + offs;

            if ~isnan(muWT(w))
                bar(xWT, muWT(w), 0.35, 'FaceColor', COL_WT, 'EdgeColor','none');
                if ~isnan(seWT(w))
                    errorbar(xWT, muWT(w), seWT(w), 'k', ...
                        'LineStyle','none', 'CapSize',6, 'LineWidth',1);
                end
            end

            if ~isnan(muAPP(w))
                bar(xAPP, muAPP(w), 0.35, 'FaceColor', COL_APP, 'EdgeColor','none');
                if ~isnan(seAPP(w))
                    errorbar(xAPP, muAPP(w), seAPP(w), 'k', ...
                        'LineStyle','none', 'CapSize',6, 'LineWidth',1);
                end
            end
        end

        % ---- per-mouse dots WITH LABELS on every point ----------------
        umice = unique(Tc.mouse);
        for mk = 1:numel(umice)
            mk_id = umice(mk);
            rowsM = Tc.mouse == mk_id;

            for w = 1:3
                rowsMW = rowsM & Tc.win == w;
                if ~any(rowsMW), continue; end

                vals = Tc.mean_bout_dur_win_s(rowsMW);
                if isempty(vals) || all(isnan(vals)), continue; end

                y = mean(vals,'omitnan');
                gk = unique(Tc.geno(rowsMW));
                if isempty(gk), continue; end

                if gk == "WT"
                    baseX = w - offs;
                else
                    baseX = w + offs;
                end

                x = baseX + (rand-0.5)*2*OPT.jitter;

                % same color for all dots (black) – you can tweak if you want
                plot(x, y, 'o', 'MarkerSize',5, ...
                    'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor','none');

                text(x, y, " "+mouse_short(mk_id), ...
                    'Color',[0 0 0], 'FontSize',7, ...
                    'HorizontalAlignment','left', 'VerticalAlignment','middle');
            end
        end

        xlim([0.5 3.5]);
        xticks(1:3);
        xticklabels(win_labels);
        xlabel('Time window');
        ylabel('Mean bout duration (s)');
        title(char(cnd), 'Interpreter','none');
        box off; grid on;
    end

    sgtitle(sprintf('Mean bout duration in 3 time windows — %s', st_name), ...
        'FontWeight','bold');

    out_file = fullfile(out_dir, ...
        sprintf('fig_mean_bout_3win_%s.png', lower(char(st_name))));
    saveas(fig, out_file);
end

end
