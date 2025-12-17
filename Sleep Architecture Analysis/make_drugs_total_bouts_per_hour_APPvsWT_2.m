function OUT = make_drugs_total_bouts_per_hour_APPvsWT_2(rows_perhr, out_dir, states_to_plot, varargin)
% make_drugs_total_bouts_per_hour_APPvsWT_2
% -------------------------------------------------------------------------
% For drugs recordings only (condition=="drugs"):
%   For each requested STATE (WK, MA, NREM, REM), make:
%
%   1) Time-course plot:
%       X: hour_idx (0, 1, 2, ...)
%       Y: bouts_per_h in that state
%       Curves: WT vs APP (mean ± SEM)
%
%   2) Between-mice variability plot:
%       Per mouse: mean over hours -> compare WT vs APP (bars + dots)
%
%   3) Tables per state:
%       a) OUT.hourByMouse.(STATE)  [and CSV]
%            "Hour;APP_mouse1;WT_mouse2;..."
%            rows = hours, columns = mice
%
%       b) OUT.perMouseHour.(STATE)      (mouse x hour, internal)
%          OUT.perMouseHour_long.(STATE) (one row per mouse–hour, internal)
%
%   The "Hour;APP_mouse1;WT_mouse2;..." CSV is:
%       drugs_bouts_per_hour_<state>_APPvsWT_hour_by_mouse.csv
%
% INPUT
%   rows_perhr    : table from run_group_sleep_architecture (group_per_hour)
%                   must have columns: condition, state, genotype, mouse,
%                                      hour_idx, bouts_per_h
%   out_dir       : folder to save figures/CSVs (default: pwd)
%   states_to_plot: string/cell array, default ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'showIDs'              : true/false, show mouse IDs next to dots (default: false)
%   'showStars'            : true/false, draw per-hour stars (default: true)
%   'useFDRforStars'       : true/false, use FDR p for stars (default: true)
%   'minNperGroupForStats' : minimum n per genotype to draw a star (default: 3)
%   'pointStyle'           : 'dots' | 'bands' | 'both' | 'none' (default: 'dots')
%   'bandAlpha'            : transparency of SEM band (default: 0.2)
% -------------------------------------------------------------------------

% ----------- Parse inputs -----------
if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

if nargin < 3 || isempty(states_to_plot)
    states_to_plot = ["WK","MA","NREM","REM"];
end
states_to_plot = string(states_to_plot(:)).';   % row string array

p = inputParser;
addParameter(p, 'showIDs', false, @(x)islogical(x) && isscalar(x));
addParameter(p, 'showStars', true, @(x)islogical(x) && isscalar(x));
addParameter(p, 'useFDRforStars', true, @(x)islogical(x) && isscalar(x));
addParameter(p, 'minNperGroupForStats', 3, @(x)isscalar(x) && x>=1);
addParameter(p, 'pointStyle', 'dots', @(s)ischar(s) || isstring(s));
addParameter(p, 'bandAlpha', 0.2, @(x)isscalar(x) && x>=0 && x<=1);
parse(p, varargin{:});

showIDs              = p.Results.showIDs;
showStars            = p.Results.showStars;
useFDRforStars       = p.Results.useFDRforStars;
minNperGroupForStats = p.Results.minNperGroupForStats;
pointStyle           = lower(string(p.Results.pointStyle));
bandAlpha            = p.Results.bandAlpha;

P = rows_perhr;

% ----------- 1) drugs only -----------
cond_str    = lower(strtrim(P.condition));
is_drugs  = cond_str == "drugs";   % change this if your label differs
P = P(is_drugs, :);

if isempty(P)
    warning('No drugs rows found in rows_perhr (condition=="drugs"). Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no drugs data');
    return;
end

OUT = struct();
OUT.success        = true;
OUT.states_plotted = [];
OUT.files          = struct();
OUT.var_files      = struct();
OUT.stats          = struct();
OUT.varstats       = struct();
OUT.anova          = struct();
OUT.perMouseHour       = struct();  % internal: mouse x hour
OUT.perMouseHour_long  = struct();  % internal: long format
OUT.hourByMouse        = struct();  % NEW: Hour x mouse (format you want)

% Colors
COL_WT      = [0.6 0.6 0.6];       % grey
COL_APP     = [0.39 0.58 0.93];    % cornflower-ish blue
COL_WT_DOT  = [0.3 0.3 0.3];       % darker grey
COL_APP_DOT = [0.1 0.2 0.6];       % darker blue

for s = 1:numel(states_to_plot)
    st = states_to_plot(s);

    % ----------- 2) Filter to this state -----------
    Pst = P(P.state == st, :);
    if isempty(Pst)
        warning('No drugs rows for state "%s". Skipping.', st);
        continue;
    end

    G = Pst;  % shorthand

    hasWT  = any(G.genotype == "WT");
    hasAPP = any(G.genotype == "APP");

    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state "%s". Skipping.', st);
        continue;
    end

    stField = matlab.lang.makeValidName(st);   % struct field name, e.g. 'REM'

    % ----------- 2b) Prepare data per mouse × hour -----------
    all_hours = unique(G.hour_idx);
    all_hours = sort(all_hours);
    nH = numel(all_hours);

    mice = unique(G.mouse);
    nM   = numel(mice);

    dataMat   = nan(nM, nH);   % per-mouse (rows) × hour (cols)
    genoVec   = strings(nM,1); % genotype per mouse

    for i = 1:nM
        mID = mice(i);
        maskMouse = (G.mouse == mID);

        % genotype for this mouse
        g_this = unique(G.genotype(maskMouse));
        if numel(g_this) ~= 1
            warning('Mouse %s has multiple genotypes? Using first.', string(mID));
            genoVec(i) = g_this(1);
        else
            genoVec(i) = g_this;
        end

        for h = 1:nH
            hr = all_hours(h);
            mask = maskMouse & (G.hour_idx == hr);
            vals = G.bouts_per_h(mask);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                dataMat(i,h) = mean(vals);  % in case of duplicates
            end
        end
    end

    % ----------- 2c) WIDE TABLE (per mouse × per hour) - internal --------
    hourNames = arrayfun(@(hr) sprintf('H%d', hr), all_hours, 'UniformOutput', false);
    Twide = table(mice, genoVec, 'VariableNames', {'mouse','genotype'});
    for h = 1:nH
        Twide.(hourNames{h}) = dataMat(:,h);
    end

    OUT.perMouseHour.(stField) = Twide;

    % ----------- 2d) LONG TABLE (mouse, genotype, hour_idx, value) -------
    long_mouse    = [];
    long_genotype = [];
    long_hour     = [];
    long_value    = [];

    for i = 1:nM
        for h = 1:nH
            val = dataMat(i,h);
            if ~isnan(val)
                long_mouse    = [long_mouse;    mice(i)];
                long_genotype = [long_genotype; genoVec(i)];
                long_hour     = [long_hour;     all_hours(h)];
                long_value    = [long_value;    val];
            end
        end
    end

    Tlong = table(long_mouse, long_genotype, long_hour, long_value, ...
        'VariableNames', {'mouse','genotype','hour_idx','bouts_per_h'});

    OUT.perMouseHour_long.(stField) = Tlong;

    % ----------- 2e) HOUR × MOUSE TABLE (the format you want) ------------
    %   Columns: Hour, GENO_mouseID, GENO_mouseID2, ...
    %
    %   We build varnames like "APP_mouse1", "WT_mouse11", etc.
    %
    hourTbl = table(all_hours, 'VariableNames', {'Hour'});

    for i = 1:nM
        % create a nice variable name, e.g. 'APP_mouse1'
        mouseLabel = string(mice(i));
        genoLabel  = string(genoVec(i));
        varNameStr = strcat(genoLabel, "_", mouseLabel);
        varName    = matlab.lang.makeValidName(varNameStr);  % safe for MATLAB

        hourTbl.(varName) = dataMat(i, :)';  % one column per mouse
    end

    OUT.hourByMouse.(stField) = hourTbl;

    % Save this as a CSV with ';' as delimiter, matching your example
    csv_name_hour = sprintf('drugs_bouts_per_hour_%s_APPvsWT_hour_by_mouse.csv', lower(st));
    csv_path_hour = fullfile(out_dir, csv_name_hour);
    try
        writetable(hourTbl, csv_path_hour, 'Delimiter',';');
        fprintf('📄 Saved HOUR×MOUSE table for %s to: %s\n', st, csv_path_hour);
    catch ME
        warning('Could not save hour×mouse CSV for %s: %s', st, ME.message);
    end

    % ----------- 2f) RM ANOVA (Time × Genotype) – robust to 1 genotype ----
    % Drop mice that are NaN for all hours
    validRows   = any(~isnan(dataMat), 2);
    Twide_anova = Twide(validRows, :);
    geno_anova  = genoVec(validRows);

    uniqG = unique(geno_anova);
    if numel(uniqG) < 2
        % Only one genotype present -> cannot fit model "~ genotype"
        warning('State %s (drugs): only one genotype present for RM ANOVA. Skipping fitrm.', st);
        rm         = [];
        withinTbl  = table();
        betweenTbl = table();
    else
        WithinDesign = table(all_hours(:), 'VariableNames', {'Time'});
        measureStr   = sprintf('%s-%s', hourNames{1}, hourNames{end});

        try
            rm = fitrm(Twide_anova, sprintf('%s ~ genotype', measureStr), ...
                       'WithinDesign', WithinDesign);
            withinTbl  = ranova(rm, 'WithinModel','Time');
            betweenTbl = anova(rm);
        catch ME
            warning('fitrm failed for state %s: %s', st, ME.message);
            rm         = [];
            withinTbl  = table();
            betweenTbl = table();
        end
    end

    OUT.anova.(stField).rm       = rm;
    OUT.anova.(stField).within   = withinTbl;
    OUT.anova.(stField).between  = betweenTbl;

    fprintf('\n===== 2-way RM ANOVA (Time x Genotype) for %s bouts/hour (drugs) =====\n', st);

    % Time and interaction
    try
        if ~isempty(withinTbl) && any(strcmp(withinTbl.Properties.VariableNames,"Term"))
            rowTime = strcmp(withinTbl.Term, 'Time');
            rowInt  = strcmp(withinTbl.Term, 'Time:genotype');

            if any(rowTime) && any(strcmp(withinTbl.Properties.VariableNames,"pValueGG"))
                p_time = withinTbl.pValueGG(rowTime);
                fprintf('Time effect:           p = %.4g (GG-corrected)\n', p_time);
            end
            if any(rowInt) && any(strcmp(withinTbl.Properties.VariableNames,"pValueGG"))
                p_int = withinTbl.pValueGG(rowInt);
                fprintf('Time x Genotype:      p = %.4g (GG-corrected)\n', p_int);
            end
        end
    catch
        fprintf('  (Could not extract Time / Time×Genotype p-values)\n');
    end

    % Genotype main effect
    try
        if ~isempty(betweenTbl) && any(strcmp(betweenTbl.Properties.VariableNames,"Term"))
            rowGen = strcmp(betweenTbl.Term, 'genotype');
            if any(rowGen) && any(strcmp(betweenTbl.Properties.VariableNames,"pValue"))
                p_gen = betweenTbl.pValue(rowGen);
                fprintf('Genotype main effect: p = %.4g\n', p_gen);
            end
        end
    catch
        fprintf('  (Could not extract genotype main effect p-value)\n');
    end

    % ----------- 3) Per-hour means/SEMs + stats -----------
    meanWT  = nan(1, nH); semWT  = nan(1, nH);
    meanAPP = nan(1, nH); semAPP = nan(1, nH);

    p_t      = nan(1, nH);
    p_rs     = nan(1, nH);
    cohen_d  = nan(1, nH);
    nWT_vec  = nan(1, nH);
    nAPP_vec = nan(1, nH);
    max_y    = nan(1, nH);

    for h = 1:nH
        hr = all_hours(h);

        valsWT  = [];
        valsAPP = [];

        if hasWT
            maskWT = (G.genotype == "WT") & (G.hour_idx == hr);
            valsWT = G.bouts_per_h(maskWT);
            valsWT = valsWT(~isnan(valsWT));
            if ~isempty(valsWT)
                meanWT(h) = mean(valsWT);
                semWT(h)  = std(valsWT) / sqrt(numel(valsWT));
            end
        end

        if hasAPP
            maskAPP = (G.genotype == "APP") & (G.hour_idx == hr);
            valsAPP = G.bouts_per_h(maskAPP);
            valsAPP = valsAPP(~isnan(valsAPP));
            if ~isempty(valsAPP)
                meanAPP(h) = mean(valsAPP);
                semAPP(h)  = std(valsAPP) / sqrt(numel(valsAPP));
            end
        end

        if ~isempty(valsWT) && ~isempty(valsAPP)
            nWT_vec(h)  = numel(valsWT);
            nAPP_vec(h) = numel(valsAPP);

            [~, p_t(h)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            p_rs(h) = ranksum(valsWT, valsAPP);

            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
            cohen_d(h) = (m2 - m1) / sp;
        end

        all_vals = [];
        if ~isempty(valsWT),  all_vals = [all_vals; valsWT];  end %#ok<AGROW>
        if ~isempty(valsAPP), all_vals = [all_vals; valsAPP]; end %#ok<AGROW>
        if ~isempty(all_vals)
            max_y(h) = max(all_vals);
        end
    end

    % BH–FDR correction across hours
    p_t_fdr = nan(size(p_t));
    valid   = ~isnan(p_t);
    pvals   = p_t(valid);
    if ~isempty(pvals)
        [sorted_p, sort_idx] = sort(pvals(:));
        m = numel(sorted_p);
        adj = sorted_p .* (m ./ (1:m)');
        for i2 = m-1:-1:1
            adj(i2) = min(adj(i2), adj(i2+1));
        end
        adj(adj>1) = 1;
        p_fdr_vals = nan(size(pvals));
        p_fdr_vals(sort_idx) = adj;
        p_t_fdr(valid) = p_fdr_vals;
    end

    % ----------- 4) Time-course plot for this state -----------
    figure('Color','w'); hold on;

    useBands = any(pointStyle == ["bands","both"]);
    useDots  = any(pointStyle == ["dots","both"]);

    % --- SEM bands (WT, APP) ---
    if useBands
        % WT band
        if hasWT
            validWT = ~isnan(meanWT) & ~isnan(semWT);
            xWT = all_hours(validWT);
            mw  = meanWT(validWT);
            sw  = semWT(validWT);
            if numel(xWT) >= 2 && numel(mw) == numel(sw)
                yu = mw + sw;
                yl = mw - sw;
                xp = [xWT(:); flipud(xWT(:))];
                yp = [yl(:);  flipud(yu(:))];
                if numel(xp) == numel(yp)
                    patch(xp, yp, COL_WT, 'FaceAlpha', bandAlpha, 'EdgeColor','none');
                end
            end
        end
        % APP band
        if hasAPP
            validAPP = ~isnan(meanAPP) & ~isnan(semAPP);
            xAPP = all_hours(validAPP);
            ma   = meanAPP(validAPP);
            sa   = semAPP(validAPP);
            if numel(xAPP) >= 2 && numel(ma) == numel(sa)
                yu = ma + sa;
                yl = ma - sa;
                xp = [xAPP(:); flipud(xAPP(:))];
                yp = [yl(:);   flipud(yu(:))];
                if numel(xp) == numel(yp)
                    patch(xp, yp, COL_APP, 'FaceAlpha', bandAlpha, 'EdgeColor','none');
                end
            end
        end
    end

    % --- Mean ± SEM lines (errorbars) ---
    if hasWT
        errorbar(all_hours, meanWT, semWT, '-o', ...
            'Color', COL_WT, ...
            'MarkerFaceColor', COL_WT, ...
            'MarkerSize', 5, ...
            'LineWidth', 1.2);
    end
    if hasAPP
        errorbar(all_hours, meanAPP, semAPP, '-o', ...
            'Color', COL_APP, ...
            'MarkerFaceColor', COL_APP, ...
            'MarkerSize', 5, ...
            'LineWidth', 1.2);
    end

    % --- Dots (+ optional IDs) ---
    jitterFrac = 0.25;
    y_offset   = 0.5;

    if useDots
        for h = 1:nH
            hr = all_hours(h);

            if hasWT
                maskWT  = (G.genotype == "WT") & (G.hour_idx == hr);
                valsWT  = G.bouts_per_h(maskWT);
                mWT     = G.mouse(maskWT);
                if ~isempty(valsWT)
                    xw = hr - 0.1 + (rand(size(valsWT)) - 0.5) * jitterFrac;
                    plot(xw, valsWT, '.', 'Color', COL_WT_DOT, 'MarkerSize', 10);
                    if showIDs
                        for j = 1:numel(valsWT)
                            thisID = char(mWT(j));
                            text(xw(j), valsWT(j) + y_offset, thisID, ...
                                'Rotation', 45, ...
                                'HorizontalAlignment','left', ...
                                'VerticalAlignment','bottom', ...
                                'FontSize', 8, ...
                                'Color', COL_WT_DOT);
                        end
                    end
                end
            end

            if hasAPP
                maskAPP = (G.genotype == "APP") & (G.hour_idx == hr);
                valsAPP = G.bouts_per_h(maskAPP);
                mAPP    = G.mouse(maskAPP);
                if ~isempty(valsAPP)
                    xa = hr + 0.1 + (rand(size(valsAPP)) - 0.5) * jitterFrac;
                    plot(xa, valsAPP, '.', 'Color', COL_APP_DOT, 'MarkerSize', 10);
                    if showIDs
                        for j = 1:numel(valsAPP)
                            thisID = char(mAPP(j));
                            text(xa(j), valsAPP(j) + y_offset, thisID, ...
                                'Rotation', 45, ...
                                'HorizontalAlignment','left', ...
                                'VerticalAlignment','bottom', ...
                                'FontSize', 8, ...
                                'Color', COL_APP_DOT);
                        end
                    end
                end
            end
        end
    end

    xlabel('Hour from recording start');
    ylabel(sprintf('Bouts per hour (%s)', st));
    title(sprintf('drugs: bouts per hour in %s (WT vs APP)', st));

    if hasWT && hasAPP
        legend({'WT mean \pm SEM','APP mean \pm SEM'}, ...
               'Location','northoutside', ...
               'Orientation','horizontal');
    elseif hasWT
        legend({'WT mean \pm SEM'}, 'Location','northoutside');
    elseif hasAPP
        legend({'APP mean \pm SEM'}, 'Location','northoutside');
    end

    set(gca,'Box','off','FontSize',12);
    xlim([min(all_hours)-0.5, max(all_hours)+0.5]);

    % ----------- 5) Add per-hour stars (optional, exploratory) -----------
    if any(~isnan(max_y))
        global_max_y = max(max_y(~isnan(max_y)));
    else
        global_max_y = max(G.bouts_per_h, [], 'omitnan');
    end
    if isempty(global_max_y) || isnan(global_max_y)
        global_max_y = 1;
    end

    y_top  = global_max_y + 3;
    ylim([0, y_top]);

    if showStars
        for h = 1:nH
            if useFDRforStars
                p_here = p_t_fdr(h);
            else
                p_here = p_t(h);
            end
            if isnan(p_here) || p_here >= 0.05
                continue;
            end
            if nWT_vec(h) < minNperGroupForStats || nAPP_vec(h) < minNperGroupForStats
                continue;
            end

            if p_here < 0.001
                stars = '***';
            elseif p_here < 0.01
                stars = '**';
            else
                stars = '*';
            end

            y_star = max_y(h);
            if isnan(y_star)
                y_star = global_max_y * 0.8;
            end
            y_star = y_star + 1.0;

            text(all_hours(h), y_star, stars, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize', 12, ...
                'FontWeight','bold');
        end
    end

    % ----------- 6) Save time-course figure for this state -----------
    if showIDs
        id_suffix = '_withIDs';
    else
        id_suffix = '_noIDs';
    end

    band_suffix = '';
    switch char(pointStyle)
        case 'dots'
            band_suffix = '_dots';
        case 'bands'
            band_suffix = '_bands';
        case 'both'
            band_suffix = '_bands_dots';
        case 'none'
            band_suffix = '_nolocal';
    end

    fname = sprintf('drugs_bouts_per_hour_%s_APPvsWT%s%s.png', ...
                    lower(st), id_suffix, band_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    OUT.states_plotted = [OUT.states_plotted, st];
    OUT.files.(stField) = out_file;

    fprintf('✅ drugs bouts/hour time-course plot for %s saved to: %s\n', st, out_file);

    % ----------- 7) Store per-hour stats table for this state -----------
    stats_idx = ~isnan(p_t);
    if any(stats_idx)
        stats_tbl = table( ...
            all_hours(stats_idx), ...
            nWT_vec(stats_idx)', ...
            nAPP_vec(stats_idx)', ...
            meanWT(stats_idx)', ...
            meanAPP(stats_idx)', ...
            p_t(stats_idx)', ...
            p_t_fdr(stats_idx)', ...
            p_rs(stats_idx)', ...
            cohen_d(stats_idx)', ...
            'VariableNames', {'Hour','nWT','nAPP','MeanWT','MeanAPP', ...
                              'p_ttest','p_ttest_FDR','p_ranksum','Cohen_d'});
        fprintf('\ndrugs bouts/hour (%s): WT vs APP per-hour stats (exploratory)\n', st);
        disp(stats_tbl);
        OUT.stats.(stField) = stats_tbl;
    else
        fprintf('\n[Stats %s] No hours with both WT and APP present. No per-hour tests.\n', st);
        OUT.stats.(stField) = table();
    end

    % ----------- 8) BETWEEN-MICE variability plot -----------------------
    % For each mouse: mean bouts/hour across *all drugs hours* for this state
    mean_mouse = nan(nM,1);
    for i = 1:nM
        row = dataMat(i,:);
        row = row(~isnan(row));
        if ~isempty(row)
            mean_mouse(i) = mean(row);
        end
    end

    maskWTm  = (genoVec == "WT");
    maskAPPm = (genoVec == "APP");

    WT_means  = mean_mouse(maskWTm);
    APP_means = mean_mouse(maskAPPm);

    WT_means  = WT_means(~isnan(WT_means));
    APP_means = APP_means(~isnan(APP_means));

    nWTm  = numel(WT_means);
    nAPPm = numel(APP_means);

    mean_WT = mean(WT_means,'omitnan');
    mean_APP = mean(APP_means,'omitnan');

    sd_WT = std(WT_means,'omitnan');
    sd_APP = std(APP_means,'omitnan');

    cv_WT = sd_WT / max(mean_WT, eps);
    cv_APP = sd_APP / max(mean_APP, eps);

    % stats on between-mice variability
    p_var   = NaN;   % equality of variance
    p_tmean = NaN;   % difference in means
    d_mean  = NaN;   % Cohen d (means)

    if nWTm >= 2 && nAPPm >= 2
        % F-test for equality of variances
        try
            [~, p_var] = vartest2(WT_means, APP_means);
        catch
            p_var = NaN;
        end

        % difference in mean bouts/hour
        [~, p_tmean] = ttest2(WT_means, APP_means, 'Vartype','unequal');

        m1 = mean_WT;  m2 = mean_APP;
        s1 = sd_WT;    s2 = sd_APP;
        n1 = nWTm;     n2 = nAPPm;
        sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
        d_mean = (m2 - m1) / sp;
    end

    % variability figure (between-mice)
    figure('Color','w'); hold on;
    Xpos = [1 2];
    barWidth = 0.5;

    bar(Xpos(1), mean_WT,  barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
    bar(Xpos(2), mean_APP, barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');

    % error bars = SD across mice (between-mice variability)
    errorbar(Xpos(1), mean_WT,  sd_WT,  'k', 'LineStyle','none', 'LineWidth',1);
    errorbar(Xpos(2), mean_APP, sd_APP, 'k', 'LineStyle','none', 'LineWidth',1);

    % jittered dots per mouse
    jit = 0.12;
    if nWTm > 0
        xw = Xpos(1) + (rand(size(WT_means))-0.5)*2*jit;
        plot(xw, WT_means, '.', 'Color', COL_WT_DOT, 'MarkerSize', 12);
        if showIDs
            mWT_ids = mice(maskWTm);
            mWT_ids = mWT_ids(~isnan(mean_mouse(maskWTm)));
            for j = 1:numel(WT_means)
                text(xw(j), WT_means(j), char(mWT_ids(j)), ...
                     'Rotation', 45, ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize',8, ...
                     'Color',COL_WT_DOT);
            end
        end
    end
    if nAPPm > 0
        xa = Xpos(2) + (rand(size(APP_means))-0.5)*2*jit;
        plot(xa, APP_means, '.', 'Color', COL_APP_DOT, 'MarkerSize', 12);
        if showIDs
            mAPP_ids = mice(maskAPPm);
            mAPP_ids = mAPP_ids(~isnan(mean_mouse(maskAPPm)));
            for j = 1:numel(APP_means)
                text(xa(j), APP_means(j), char(mAPP_ids(j)), ...
                     'Rotation', 45, ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize',8, ...
                     'Color',COL_APP_DOT);
            end
        end
    end

    xlim([0.5 2.5]);
    ylabel(sprintf('Per-mouse mean bouts/hour (%s)', st));
    set(gca,'XTick',Xpos,'XTickLabel',{'WT','APP'},'FontSize',12);
    title(sprintf('drugs: between-mice variability of %s bouts/hour', st));
    set(gca,'Box','off');

    % annotate p-values (optional simple text)
    ymax = max([WT_means;APP_means],[],'omitnan');
    if isempty(ymax) || isnan(ymax)
        ymax = 1;
    end
    ylim([0, ymax*1.4]);

    % --- Text with numeric p-values (variance + mean) ---
    if ~isnan(p_var) || ~isnan(p_tmean)
        txt = sprintf('var-test p=%.3f | mean t-test p=%.3f | d=%.2f', ...
                      p_var, p_tmean, d_mean);
        text(1.5, ymax*1.25, txt, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',10);
    end

    % --- Significance stars for mean difference (WT vs APP) ---
    starStr = '';
    if ~isnan(p_tmean)
        if p_tmean < 0.001
            starStr = '***';
        elseif p_tmean < 0.01
            starStr = '**';
        elseif p_tmean < 0.05
            starStr = '*';
        end
    end

    if ~isempty(starStr)
        % horizontal line above the two bars
        y_star_line = ymax * 1.10;
        line([1 2], [y_star_line y_star_line], 'Color','k', 'LineWidth', 1.2);

        % star label slightly above the line
        text(1.5, ymax * 1.15, starStr, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize', 14, ...
             'FontWeight','bold');
    end

    % save variability figure
    if showIDs
        id_suffix2 = '_withIDs';
    else
        id_suffix2 = '_noIDs';
    end
    fname_var = sprintf('drugs_bouts_per_hour_%s_APPvsWT_betweenmice%s.png', ...
                        lower(st), id_suffix2);
    out_file_var = fullfile(out_dir, fname_var);
    saveas(gcf, out_file_var);
    OUT.var_files.(stField) = out_file_var;
    fprintf('✅ Between-mice variability plot for %s saved to: %s\n', st, out_file_var);

    % store variability stats
    VAR = struct();
    VAR.per_mouse = table( ...
        mice, genoVec, mean_mouse, ...
        'VariableNames', {'mouse','genotype','mean_bouts_per_h'});
    VAR.genotype_summary = table( ...
        ["WT";"APP"], ...
        [nWTm; nAPPm], ...
        [mean_WT; mean_APP], ...
        [sd_WT; sd_APP], ...
        [cv_WT; cv_APP], ...
        'VariableNames', {'Genotype','n','Mean_bouts_per_h','SD_between_mice','CV_between_mice'});
    VAR.p_vartest2_varEqual = p_var;
    VAR.p_ttest_mean        = p_tmean;
    VAR.Cohen_d_mean        = d_mean;

    OUT.varstats.(stField) = VAR;
end

end
