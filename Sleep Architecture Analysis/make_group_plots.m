function make_group_plots(OVERALL, PERHOUR, out_dir, varargin)
% make_group_plots with outlier overlays (per-mouse) and SEM bars.
% Optional name-value:
%   'label_outliers' (false) : add mouse labels next to red dots
%   'jitter'         (0.12)  : horizontal jitter for per-mouse dots

ip = inputParser;
addParameter(ip,'label_outliers',true,@islogical);
addParameter(ip,'jitter',0.12,@(x)isscalar(x)&&x>=0&&x<=0.5);
parse(ip,varargin{:});
OPT = ip.Results;

% Palette per state
COL = struct('WK',[0.55 0.55 0.55], 'NREM',[0.30 0.50 0.85], 'REM',[0.80 0.35 0.35], 'MA',[0.60 0.60 0.20]);

% ---------- Compute outlier flags ----------
OVERALL2 = OVERALL;
OVERALL2.total_min = OVERALL2.total_dur_s/60;
OVERALL2 = add_outlier_cols(OVERALL2, {'state','genotype','condition'}, ...
            {'total_min','n_bouts','mean_bout_dur_s'});

PERHOUR2 = add_outlier_cols(PERHOUR, {'state','genotype','condition','hour_idx'}, ...
            {'bouts_per_h','dur_s','mean_bout_dur_s'});

% Which states exist?
states = cellstr(unique(OVERALL2.state,'stable'))';

% ---------- Fig 1: Total duration (min) ----------
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL2(strcmp(OVERALL2.state, st), :);
    if isempty(sub), continue; end
    G = agg_mean_sem(sub, {'genotype','condition'}, 'total_min');
    if isempty(G), continue; end

    f = figure('Color','w','Name', ['Total ' st ' (min)']); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    xcats = categorical(cats, cats);
    b = bar(xcats, G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(xcats, G.mean, G.sem, 'k.', 'LineWidth',1);

    % per-mouse overlay
    x_mouse = categorical(sub.genotype + " | " + sub.condition, cats);
    xnum = double(x_mouse);
    jit = (rand(height(sub),1)-0.5)*OPT.jitter;
    plot(xnum + jit, sub.total_min, 'k.', 'MarkerSize',10);
    is_out = sub.total_min_is_outlier;
    plot(xnum(is_out)+jit(is_out), sub.total_min(is_out), 'r.', 'MarkerSize',12);
    if OPT.label_outliers
        for ii = find(is_out).'
            text(xnum(ii)+jit(ii), sub.total_min(ii), " "+mouse_short(sub.mouse(ii)), ...
             'Color',[0.6 0 0], 'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end
    end

    ylabel('Total duration (min)'); title(sprintf('Total %s (min) — genotype × condition', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_total_%s.png', lower(st))));
end

% ---------- Fig 2: Number of bouts ----------
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL2(strcmp(OVERALL2.state, st), :);
    if isempty(sub), continue; end
    G = agg_mean_sem(sub, {'genotype','condition'}, 'n_bouts');
    if isempty(G), continue; end

    f = figure('Color','w','Name', ['Bouts ' st]); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    xcats = categorical(cats, cats);
    b = bar(xcats, G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(xcats, G.mean, G.sem, 'k.', 'LineWidth',1);

    x_mouse = categorical(sub.genotype + " | " + sub.condition, cats);
    xnum = double(x_mouse);
    jit = (rand(height(sub),1)-0.5)*OPT.jitter;
    plot(xnum + jit, sub.n_bouts, 'k.', 'MarkerSize',10);
    is_out = sub.n_bouts_is_outlier;
    plot(xnum(is_out)+jit(is_out), sub.n_bouts(is_out), 'r.', 'MarkerSize',12);
    if OPT.label_outliers
        for ii = find(is_out).'
            text(xnum(ii)+jit(ii), sub.n_bouts(ii), " "+mouse_short(sub.mouse(ii)), ...
                 'Color',[0.6 0 0], 'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','middle');

        end
    end

    ylabel('Number of bouts'); title(sprintf('Bouts — %s', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_bouts_%s.png', lower(st))));
end

% ---------- Fig 3: Mean bout duration (s) ----------
for s = 1:numel(states)
    st = states{s};
    sub = OVERALL2(strcmp(OVERALL2.state, st), :);
    if isempty(sub), continue; end
    G = agg_mean_sem(sub, {'genotype','condition'}, 'mean_bout_dur_s');
    if isempty(G), continue; end

    f = figure('Color','w','Name', ['Mean bout dur ' st]); hold on
    cats = strcat(G.genotype, " | ", G.condition);
    xcats = categorical(cats, cats);
    b = bar(xcats, G.mean); set(b,'FaceColor', pick_col(COL, st));
    errorbar(xcats, G.mean, G.sem, 'k.', 'LineWidth',1);

    x_mouse = categorical(sub.genotype + " | " + sub.condition, cats);
    xnum = double(x_mouse);
    jit = (rand(height(sub),1)-0.5)*OPT.jitter;
    plot(xnum + jit, sub.mean_bout_dur_s, 'k.', 'MarkerSize',10);
    is_out = sub.mean_bout_dur_s_is_outlier;
    plot(xnum(is_out)+jit(is_out), sub.mean_bout_dur_s(is_out), 'r.', 'MarkerSize',12);
    if OPT.label_outliers
        for ii = find(is_out).'
            text(xnum(ii)+jit(ii), sub.mean_bout_dur_s(ii), " "+mouse_short(sub.mouse(ii)), ...
                 'Color',[0.6 0 0], 'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','middle');

        end
    end

    ylabel('Mean bout duration (s)'); title(sprintf('Mean bout duration — %s', st));
    xtickangle(30); box off
    saveas(f, fullfile(out_dir, sprintf('fig_mean_bout_%s.png', lower(st))));
end

% ---------- Fig 4: Bouts per hour (mean ± SEM) + outlier points ----------
PH = PERHOUR2; PH = PH(~isnan(PH.bouts_per_h),:);
if ~isempty(PH)
    Tagg = agg_mean_sem(PH, {'state','genotype','condition','hour_idx'}, 'bouts_per_h');
    Tagg.Properties.VariableNames{'hour_idx'} = 'hour';

    for s = 1:numel(states)
        st = states{s};
        sub = Tagg(strcmp(Tagg.state, st), :);
        if isempty(sub), continue; end

        f = figure('Color','w','Name', ['Bouts per hour — ' st]); tiledlayout('flow');
        Ugeno = unique(sub.genotype,'stable'); Ucond = unique(sub.condition,'stable');

        for g = 1:numel(Ugeno)
            for c = 1:numel(Ucond)
                sc = sub(strcmp(sub.genotype,Ugeno{g}) & strcmp(sub.condition,Ucond{c}), :);
                if isempty(sc), continue; end
                nexttile; hold on
                errorbar(double(sc.hour), sc.mean, sc.sem, 'o-','LineWidth',1.25, ...
                         'Color', pick_col(COL, st));

                % per-file overlay (faint)
                pm = PH(strcmp(PH.state, st) & strcmp(PH.genotype,Ugeno{g}) & strcmp(PH.condition,Ucond{c}), :);
                files_u = unique(pm.file,'stable');
                for ff = 1:numel(files_u)
                    p1 = pm(strcmp(pm.file, files_u(ff)), :);
                    plot(double(p1.hour_idx), p1.bouts_per_h, '-', 'Color',[0 0 0 0.15], 'LineWidth',0.75);
                    % highlight outlier points
                    outmask = p1.bouts_per_h_is_outlier;
                    if any(outmask)
                        plot(double(p1.hour_idx(outmask)), p1.bouts_per_h(outmask), 'ro', 'MarkerSize',5, 'LineWidth',1);
                        if OPT.label_outliers
                            ms = mouse_short(p1.mouse(outmask));
                            for ii = 1:numel(ms)
                                text(double(p1.hour_idx(outmask(ii)))+0.05, p1.bouts_per_h(outmask(ii)), " "+ms(ii), ...
                                    'Color',[0.6 0 0], 'FontSize',7, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
                            end
                        end
                    end
                end
                title(sprintf('%s | %s', Ugeno{g}, Ucond{c}));
                xlabel('Hour'); ylabel('Bouts per hour'); grid on; box off
            end
        end
        sgtitle(sprintf('Bouts per hour — %s', st));
        saveas(f, fullfile(out_dir, sprintf('fig_bph_byhour_%s.png', lower(st))));
    end
end
end

% ================= helpers =================

function T2 = add_outlier_cols(T, group_vars, metric_list)
% Adds *_z, *_rz, *_is_outlier per metric within groups defined by group_vars.
T2 = T;

% Build normalized grouping vectors once
gv = cell(1, numel(group_vars));
for k = 1:numel(group_vars)
    v = T2.(group_vars{k});
    if ischar(v) || isstring(v) || iscellstr(v)
        gv{k} = string(v);
    else
        gv{k} = v;
    end
end
G = findgroups(gv{:});   % <- only capture the group index

for m = 1:numel(metric_list)
    metric = metric_list{m};
    if ~ismember(metric, T2.Properties.VariableNames), continue; end

    zname   = [metric '_z'];
    rzname  = [metric '_rz'];
    flagcol = [metric '_is_outlier'];

    T2.(zname)   = NaN(height(T2),1);
    T2.(rzname)  = NaN(height(T2),1);
    T2.(flagcol) = false(height(T2),1);

    x = double(T2.(metric));
    ug = unique(G);
    for gi = reshape(ug(~isnan(ug)),1,[])
        idx = (G == gi);
        xi  = x(idx);

        mu  = mean(xi,'omitnan');
        sd  = std(xi,'omitnan');
        med = median(xi,'omitnan');
        madv = mad_n(xi);

        z  = (xi - mu) ./ max(sd,  eps);
        rz = (xi - med)./ max(madv, eps);

        out = abs(rz) > 3 | abs(z) > 3.5;

        T2.(zname)(idx)   = z;
        T2.(rzname)(idx)  = rz;
        T2.(flagcol)(idx) = out;
    end
end
end


function G = agg_mean_sem(T, groupVars, valueVar)
vals = double(T.(valueVar));
gv = cell(1, numel(groupVars));
for k = 1:numel(groupVars)
    v = T.(groupVars{k});
    if ischar(v) || isstring(v) || iscellstr(v)
        gv{k} = string(v);
    else
        gv{k} = v;
    end
end
keysOut = cell(1, numel(groupVars));
[Gid, keysOut{:}] = findgroups(gv{:});
m  = splitapply(@(x) mean(x,'omitnan'), vals, Gid);
sd = splitapply(@(x) std(x,'omitnan'),  vals, Gid);
n  = splitapply(@(x) sum(~isnan(x)),    vals, Gid);
sem = sd ./ max(sqrt(n), 1);
G = table();
for k = 1:numel(groupVars)
    G.(groupVars{k}) = keysOut{k};
end
G.mean = m; G.sem = sem;
end

function s = mad_n(x)
m = median(x,'omitnan');
s = 1.4826 * median(abs(x - m),'omitnan');
end

function c = pick_col(COL, st)
switch upper(st)
    case 'WK',   c = COL.WK;
    case 'NREM', c = COL.NREM;
    case 'REM',  c = COL.REM;
    case 'MA',   c = COL.MA;
    otherwise,   c = [0.5 0.5 0.5];
end
end
function lbl = mouse_short(s)
% Make compact labels like "m1" from "mouse1"/"Mouse-03"/"03".
s = string(s);
lbl = strings(size(s));
for k = 1:numel(s)
    t = s(k);
    num = regexp(t, '\d+$', 'match', 'once');  % direct trailing digits
    if ~isempty(num)
        lbl(k) = "m" + string(num);
    else
        % fallback: trimmed original, max 8 chars
        t = regexprep(t, '^\s+|\s+$', '');
        if strlength(t) > 8, t = extractBefore(t, 9); end
        lbl(k) = t;
    end
end
end
