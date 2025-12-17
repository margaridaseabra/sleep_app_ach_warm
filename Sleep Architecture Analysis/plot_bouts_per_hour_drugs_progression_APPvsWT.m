function OUT = plot_bouts_per_hour_drugs_progression_APPvsWT( ...
                rows_perhr, meta_tbl, states_to_plot, out_dir, varargin)
% plot_bouts_per_hour_drugs_progression_APPvsWT
% -------------------------------------------------------------------------
% Progression of bouts/hour around drugs start, WT vs APP.
%
% INPUTS
%   rows_perhr : group_per_hour table from run_group_sleep_architecture
%                Needs columns:
%                  hour_idx, hour_start_s, bouts_per_h, state,
%                  condition, mouse, genotype
%
%   meta_tbl   : Excel table with columns (VariableNamingRule 'preserve'):
%                  'Mouse','Genotype','Time started (s)', ...
%
%   states_to_plot : e.g. ["WK","MA","NREM","REM"]
%
%   out_dir    : folder for figures
%
% NAME–VALUE
%   'relHourRange'        : [min max] hours relative to drugs start (default [-3 5])
%   'drugWindow_h'         : duration of the manipulation (for shade), default 3
%   'useFDRforStars'      : use BH-FDR per state (default true)
%   'minNperGroupForStats': min WT and APP per bin for stats (default 3)
% -------------------------------------------------------------------------

OUT = struct('success',false,'states',[], ...
             'fig_files',struct(), ...
             'stats',struct(), ...
             'anova',struct());

% ---------- args ----------
if nargin < 4 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

if nargin < 3 || isempty(states_to_plot)
    states_to_plot = ["WK","MA","NREM","REM"];
end
states_to_plot = string(states_to_plot(:)).';

p = inputParser;
addParameter(p,'relHourRange',[-3 5],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'drugWindow_h',6,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'useFDRforStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p,varargin{:});

relRange   = p.Results.relHourRange;
drugWindow  = p.Results.drugWindow_h;
useFDR     = p.Results.useFDRforStars;
minN       = p.Results.minNperGroupForStats;

% ---------- 1) Normalize rows_perhr ----------
P = rows_perhr;

needVars = {'hour_idx','hour_start_s','bouts_per_h','state', ...
            'condition','mouse','genotype'};
missing = setdiff(needVars, P.Properties.VariableNames);
if ~isempty(missing)
    error('rows_perhr is missing columns: %s', strjoin(missing,', '));
end

P.hour_idx     = double(P.hour_idx);
P.hour_start_s = double(P.hour_start_s);
P.bouts_per_h  = double(P.bouts_per_h);
P.state        = string(P.state);
P.condition    = lower(strtrim(string(P.condition)));
P.mouse        = string(P.mouse);
P.genotype     = string(P.genotype);

% keep only drugs recordings
P = P(P.condition=="drugs",:);
if isempty(P)
    warning('No drugs rows in rows_perhr.');
    return;
end

% ---------- 1b) numeric mouse_id ----------
mouse_id = nan(height(P),1);
for i = 1:height(P)
    mstr = char(P.mouse(i));          % e.g. 'mouse8'
    d    = regexp(mstr,'\d+','match','once');
    mouse_id(i) = str2double(d);
end
P.mouse_id = mouse_id;

% ---------- 2) meta table ----------
if ~ismember("Mouse", meta_tbl.Properties.VariableNames)
    error('meta_tbl must have a "Mouse" column (use VariableNamingRule ''preserve'').');
end
if ~ismember("Time started (s)", meta_tbl.Properties.VariableNames)
    error('meta_tbl must have a "Time started (s)" column.');
end

meta_small = table;
rawMouse = meta_tbl.("Mouse");
if isnumeric(rawMouse)
    meta_small.mouse_id = double(rawMouse);
else
    rawStr = string(rawMouse);
    mid = nan(numel(rawStr),1);
    for i = 1:numel(rawStr)
        d = regexp(char(rawStr(i)),'\d+','match','once');
        mid(i) = str2double(d);
    end
    meta_small.mouse_id = mid;
end
meta_small.t_start_s = double(meta_tbl.("Time started (s)"));

[~, ia] = unique(meta_small.mouse_id);
meta_small = meta_small(ia,:);

% ---------- 3) Join + relative time ----------
J = innerjoin(P, meta_small, 'Keys','mouse_id');

if isempty(J)
    warning('No overlap between rows_perhr and meta_tbl (by numeric mouse_id).');
    return;
end

J.rel_hour = (J.hour_start_s - J.t_start_s) / 3600;

J = J(J.rel_hour >= relRange(1) & J.rel_hour <= relRange(2), :);
if isempty(J)
    warning('After relHourRange filter, no data left.');
    return;
end

J.rel_hour_bin = round(J.rel_hour);   % -3,-2,...,5
all_bins = unique(J.rel_hour_bin);
all_bins = all_bins(:).';
nBins    = numel(all_bins);

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

% ---------- 5) loop states ----------
for s = 1:numel(states_to_plot)
    st = states_to_plot(s);

    Jst = J(J.state==st,:);
    if isempty(Jst)
        warning('No drugs rows for state %s in progression; skipping.', st);
        continue;
    end

    hasWT  = any(Jst.genotype=="WT");
    hasAPP = any(Jst.genotype=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state %s; skipping.', st);
        continue;
    end

    % ---------- 5a) per-mouse x bin matrix ----------
    mice = unique(Jst.mouse);
    nM   = numel(mice);

    dataMat = nan(nM, nBins);
    genoVec = strings(nM,1);

    for iM = 1:nM
        mID = mice(iM);
        maskM = (Jst.mouse==mID);
        g_this = unique(Jst.genotype(maskM));
        if numel(g_this) ~= 1
            g_this = g_this(1);
        end
        genoVec(iM) = g_this;

        for b = 1:nBins
            rb = all_bins(b);
            mask = maskM & (Jst.rel_hour_bin == rb);
            vals = Jst.bouts_per_h(mask);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                dataMat(iM,b) = mean(vals);
            end
        end
    end

    % ---------- 5b) Mixed-effects ANOVA: Time × Genotype ----------
    L = table;
    L.mouse    = categorical(Jst.mouse);
    L.genotype = categorical(Jst.genotype);
    L.Time     = Jst.rel_hour_bin;      % discrete relative hour
    L.bph      = Jst.bouts_per_h;

    L = L(~isnan(L.bph), :);

    pGen = NaN; pTime = NaN; pInt = NaN;
    lme  = [];
    aL   = table();

    try
        if numel(categories(L.genotype)) >= 2 && numel(unique(L.mouse)) >= 2
            lme = fitlme(L, 'bph ~ Time * genotype + (1|mouse)', ...
                           'DummyVarCoding','effects');
            aL  = anova(lme);

            termsA = string(aL.Term);

            rowG = strcmp(termsA,"genotype");
            rowT = strcmp(termsA,"Time");
            rowI = strcmp(termsA,"Time:genotype") | strcmp(termsA,"genotype:Time");

            if any(rowG), pGen  = aL.pValue(rowG); end
            if any(rowT), pTime = aL.pValue(rowT); end
            if any(rowI), pInt  = aL.pValue(rowI); end
        end
    catch ME
        warning('Mixed-effects ANOVA failed for state %s: %s', st, ME.message);
    end

    nameSt = matlab.lang.makeValidName(st);
    OUT.anova.(nameSt).lme   = lme;
    OUT.anova.(nameSt).anova = aL;
    OUT.anova.(nameSt).pGen  = pGen;
    OUT.anova.(nameSt).pTime = pTime;
    OUT.anova.(nameSt).pInt  = pInt;

    fprintf('\n[RM (mixed) ANOVA] %s bouts/hour (rel hours %d..%d)\n', ...
        st, relRange(1), relRange(2));
    fprintf('  Genotype p = %.4g | Time p = %.4g | Time×Genotype p = %.4g\n', ...
        pGen, pTime, pInt);

    % ---------- 5c) per-bin WT vs APP stats ----------
    meanWT  = nan(1,nBins); semWT  = nan(1,nBins);
    meanAPP = nan(1,nBins); semAPP = nan(1,nBins);

    p_t     = nan(1,nBins);
    p_rs    = nan(1,nBins);
    d_crit  = nan(1,nBins);
    nWT_vec = nan(1,nBins);
    nAPP_vec= nan(1,nBins);
    max_y   = nan(1,nBins);

    for b = 1:nBins
        rb = all_bins(b);
        maskBin = (Jst.rel_hour_bin == rb);

        valsWT  = Jst.bouts_per_h(maskBin & Jst.genotype=="WT");
        valsAPP = Jst.bouts_per_h(maskBin & Jst.genotype=="APP");
        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        if ~isempty(valsWT)
            meanWT(b) = mean(valsWT);
            semWT(b)  = std(valsWT)/sqrt(numel(valsWT));
        end
        if ~isempty(valsAPP)
            meanAPP(b) = mean(valsAPP);
            semAPP(b)  = std(valsAPP)/sqrt(numel(valsAPP));
        end

        if ~isempty(valsWT) && ~isempty(valsAPP)
            nWT_vec(b)  = numel(valsWT);
            nAPP_vec(b) = numel(valsAPP);

            [~, p_t(b)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            p_rs(b)     = ranksum(valsWT, valsAPP);

            m1 = mean(valsWT); m2 = mean(valsAPP);
            s1 = std(valsWT);  s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
            d_crit(b) = (m2 - m1)/sp;
        end

        all_vals = [valsWT; valsAPP];
        if ~isempty(all_vals)
            max_y(b) = max(all_vals);
        end
    end

    % ---------- 5d) BH–FDR across bins ----------
    p_t_fdr = nan(size(p_t));
    valid = ~isnan(p_t);
    pv = p_t(valid);
    if ~isempty(pv)
        [spv, idx] = sort(pv(:));
        m = numel(spv);
        adj = spv .* (m./(1:m)');
        for k = m-1:-1:1
            adj(k) = min(adj(k), adj(k+1));
        end
        adj(adj>1) = 1;
        adj_full = nan(size(pv));
        adj_full(idx) = adj;
        p_t_fdr(valid) = adj_full;
    end

    % ---------- 5e) Plot (MEAN LINES + SHADED SEM, NO DOTS) ----------
    fig = figure('Color','w','Position',[150 150 900 450]); hold on;

    % --- shaded SEM bands first (so they go behind the lines) ---
    % WT band
    if hasWT
        validWT = ~isnan(meanWT) & ~isnan(semWT);
        if any(validWT)
            xWT = all_bins(validWT);
            mw  = meanWT(validWT);
            sw  = semWT(validWT);
            yu = mw + sw;
            yl = mw - sw;
            xp = [xWT(:); flipud(xWT(:))];
            yp = [yl(:);  flipud(yu(:))];
            patch(xp, yp, COL_WT, 'FaceAlpha',0.15, 'EdgeColor','none');
        end
    end

    % APP band
    if hasAPP
        validAPP = ~isnan(meanAPP) & ~isnan(semAPP);
        if any(validAPP)
            xAPP = all_bins(validAPP);
            ma   = meanAPP(validAPP);
            sa   = semAPP(validAPP);
            yu = ma + sa;
            yl = ma - sa;
            xp = [xAPP(:); flipud(xAPP(:))];
            yp = [yl(:);   flipud(yu(:))];
            patch(xp, yp, COL_APP, 'FaceAlpha',0.15, 'EdgeColor','none');
        end
    end

    % --- mean ± SEM lines (errorbar) ---
    if hasWT
        errorbar(all_bins, meanWT, semWT, '-o', ...
            'Color', COL_WT, 'MarkerFaceColor', COL_WT, ...
            'LineWidth',1.2,'MarkerSize',5);
    end
    if hasAPP
        errorbar(all_bins, meanAPP, semAPP, '-o', ...
            'Color', COL_APP, 'MarkerFaceColor', COL_APP, ...
            'LineWidth',1.2,'MarkerSize',5);
    end

    % y-limits
    y_global = max([ ...
        max(max_y,[],'omitnan'), ...
        max(meanWT+semWT,[],'omitnan'), ...
        max(meanAPP+semAPP,[],'omitnan') ...
    ],[],'omitnan');
    if isempty(y_global) || isnan(y_global)
        y_global = 1;
    end
    ylim([0, y_global*1.35]);

    % vertical shading for drugs window [0, drugWindow]
    yl = ylim;
    patch([0 drugWindow drugWindow 0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.2);
    uistack(findobj(gca,'Type','patch','-not','FaceAlpha',0.2),'bottom');

    % per-bin stars
    if any(~isnan(max_y))
        for b = 1:nBins
            if useFDR
                p_here = p_t_fdr(b);
            else
                p_here = p_t(b);
            end
            if isnan(p_here) || p_here >= 0.05
                continue;
            end
            if nWT_vec(b) < minN || nAPP_vec(b) < minN
                continue;
            end

            stars = local_p2stars(p_here);
            y_star = max_y(b);
            if isnan(y_star)
                y_star = y_global*0.8;
            end
            y_star = y_star + 0.05*y_global;

            text(all_bins(b), y_star, stars, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',11,'FontWeight','bold');
        end
    end

    % ANOVA stars for G, T, G×T
    starG = local_p2stars(pGen);
    starT = local_p2stars(pTime);
    starI = local_p2stars(pInt);

    txt_anova = sprintf('G %s (p=%.3f)  |  T %s (p=%.3f)  |  G×T %s (p=%.3f)', ...
                        starG, pGen, starT, pTime, starI, pInt);
    text(min(all_bins), y_global*1.30, txt_anova, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'FontSize',9);

    xlabel('Hours relative to drugs start');
    ylabel(sprintf('Bouts/hour (%s)', st));
    title(sprintf('drugs progression: %s bouts/hour (WT vs APP)', st));

    if hasWT && hasAPP
        legend({'WT mean \pm SEM','APP mean \pm SEM'}, ...
               'Location','northoutside','Orientation','horizontal');
    elseif hasWT
        legend({'WT mean \pm SEM'}, 'Location','northoutside');
    elseif hasAPP
        legend({'APP mean \pm SEM'}, 'Location','northoutside');
    end

    set(gca,'Box','off','FontSize',12);

    fname = sprintf('drugs_progression_bouts_per_hour_%s_APPvsWT.png', lower(st));
    out_file = fullfile(out_dir, fname);
    saveas(fig, out_file);
    % NOTE: we **do not** close(fig); so the figure stays open

    OUT.fig_files.(nameSt) = out_file;

    % stats table per bin
    stats_idx = ~isnan(p_t);
    if any(stats_idx)
        statTbl = table( ...
            all_bins(stats_idx)', ...
            nWT_vec(stats_idx)', ...
            nAPP_vec(stats_idx)', ...
            meanWT(stats_idx)', ...
            meanAPP(stats_idx)', ...
            p_t(stats_idx)', ...
            p_t_fdr(stats_idx)', ...
            p_rs(stats_idx)', ...
            d_crit(stats_idx)', ...
            'VariableNames', {'Rel_hour','nWT','nAPP', ...
                              'MeanWT','MeanAPP', ...
                              'p_ttest','p_ttest_FDR','p_ranksum','Cohen_d'});
        OUT.stats.(nameSt) = statTbl;
        fprintf('Per-bin WT vs APP stats for %s:\n', st);
        disp(statTbl);
    else
        OUT.stats.(nameSt) = table();
        fprintf('[Stats] No bins with both WT and APP for %s.\n', st);
    end

    OUT.states = [OUT.states, st];
end

OUT.success = true;
fprintf('✅ drugs progression plots + stats saved in %s\n', out_dir);

end

% --------- local helper for p -> stars ----------
function stars = local_p2stars(p)
    if isnan(p)
        stars = 'n.s.';
    elseif p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end
