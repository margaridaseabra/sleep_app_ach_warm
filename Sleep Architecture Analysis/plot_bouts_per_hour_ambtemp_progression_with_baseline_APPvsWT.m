function OUT = plot_bouts_per_hour_ambtemp_progression_with_baseline_APPvsWT( ...
                rows_perhr, meta_tbl, states_to_plot, out_dir, varargin)
% plot_bouts_per_hour_ambtemp_progression_with_baseline_APPvsWT
% -------------------------------------------------------------------------
% Progression of bouts/hour around AMBTEMP start, WT vs APP,
% with BASELINE overlay and per-bin (per hour) 2x2 ANOVA:
%   Condition (baseline vs ambtemp) x Genotype (WT vs APP)
% using mixed-effects:
%   bouts_per_h ~ condition * genotype + (1|mouse)
%
% Stars on the plot:
%   black *    : main effect of condition (baseline vs ambtemp), p<0.05
%   magenta *  : condition x genotype interaction, p<0.05
%   (optionally FDR-corrected across hours if useFDRforStars = true)
%
% INPUTS
%   rows_perhr : group_per_hour table from run_group_sleep_architecture
%                Must contain columns:
%                  hour_idx, hour_start_s, bouts_per_h, state,
%                  condition, mouse, genotype
%
%   meta_tbl   : Excel table with columns (VariableNamingRule 'preserve'):
%                  'Mouse' (numeric or "mouseX"),
%                  'Time started (s)'  (ambtemp onset, seconds from rec start)
%
%   states_to_plot : e.g. ["WK","MA","NREM","REM"]
%
%   out_dir    : folder for figures
%
% NAME–VALUE
%   'relHourRange'        : [min max] hours relative to ambtemp start (default [-3 5])
%   'ambWindow_h'         : duration of manipulation in hours (for shaded patch, default 3)
%   'showIDs'             : (not used here, kept for API compatibility)
%   'useFDRforStars'      : BH-FDR over bins for ANOVA stars (default true)
%   'minNperGroupForStats': (currently not enforced, kept for future use)
%
% OUTPUT
%   OUT.success                 : logical
%   OUT.states                  : states actually plotted
%   OUT.fig_files.(STATE)       : PNG paths for each state
%   OUT.stats_perbin.(STATE)    : per-hour ANOVA stats table
% -------------------------------------------------------------------------

OUT = struct('success',false,'states',[], ...
             'fig_files',struct(), ...
             'stats_perbin',struct());

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
addParameter(p,'ambWindow_h',3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'useFDRforStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p,varargin{:});

relRange = p.Results.relHourRange;
ambWindow = p.Results.ambWindow_h;
useFDR    = p.Results.useFDRforStars;
% minN      = p.Results.minNperGroupForStats;   % reserved if you want to enforce n

% ---------- 1) Normalise rows_perhr ----------
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

% keep only baseline + ambtemp
P = P(P.condition=="baseline" | P.condition=="ambtemp",:);
if isempty(P)
    warning('No baseline or ambtemp rows in rows_perhr.');
    return;
end

P_amb  = P(P.condition=="ambtemp",:);
P_base = P(P.condition=="baseline",:); %#ok<NASGU> % (not used directly but kept for clarity)

if isempty(P_amb)
    warning('No ambtemp rows; cannot anchor relative hours.');
    return;
end

% ---------- 2) numeric mouse_id in P_amb ----------
mouse_id = nan(height(P_amb),1);
for i = 1:height(P_amb)
    mstr = char(P_amb.mouse(i));      % e.g. 'mouse8'
    d    = regexp(mstr,'\d+','match','once');
    mouse_id(i) = str2double(d);
end
P_amb.mouse_id = mouse_id;

% ---------- 3) meta table -> mouse_id + t_start_s ----------
if ~ismember("Mouse", meta_tbl.Properties.VariableNames)
    error('meta_tbl must have a "Mouse" column.');
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

% ---------- 4) Join ambtemp with meta & find anchor hour h0 per mouse ----------
J_amb = innerjoin(P_amb, meta_small, 'Keys','mouse_id');

if isempty(J_amb)
    warning('No overlap between ambtemp rows and meta_tbl (by numeric mouse_id).');
    return;
end

J_amb.rel_hour = (J_amb.hour_start_s - J_amb.t_start_s) / 3600;
J_amb.rel_hour_bin = round(J_amb.rel_hour);

% For each mouse, anchor hour_idx where rel_hour_bin == 0 (ambtemp start)
J0 = J_amb(J_amb.rel_hour_bin==0, :);
if isempty(J0)
    warning('Could not find any bins with rel_hour_bin==0; check meta times.');
    return;
end

h0_tbl = groupsummary(J0, 'mouse', 'median', 'hour_idx');
h0_tbl.Properties.VariableNames{end} = 'h0';
% Ensure mouse type matches P.mouse (string)
h0_tbl.mouse = string(h0_tbl.mouse);

% ---------- 5) Attach h0 to BOTH baseline & ambtemp, compute rel_hour_bin ----------
P_all = P(ismember(P.mouse, h0_tbl.mouse), :);
P_all = outerjoin(P_all, h0_tbl(:,{'mouse','h0'}), ...
                  'Keys','mouse', 'MergeKeys',true);

% rel_hour_bin from hour_idx - h0 per mouse
P_all.rel_hour_bin = P_all.hour_idx - P_all.h0;

% Keep bins in chosen range
P_all = P_all(P_all.rel_hour_bin >= relRange(1) & ...
              P_all.rel_hour_bin <= relRange(2), :);
if isempty(P_all)
    warning('No data left after relHourRange filter.');
    return;
end

all_bins = unique(P_all.rel_hour_bin);
all_bins = all_bins(:).';
nBins    = numel(all_bins);

COL_WT       = [0.6 0.6 0.6];
COL_APP      = [0.39 0.58 0.93];
COL_WT_LIGHT = 0.8*[1 1 1];
COL_APP_LIGHT= [0.7 0.8 1.0];

OUT.states = strings(0,1);

% ---------- 6) Loop over states ----------
for s = 1:numel(states_to_plot)
    st = states_to_plot(s);

    Kst = P_all(P_all.state==st,:);
    if isempty(Kst)
        warning('No rows for state %s; skipping.', st);
        continue;
    end

    hasWT  = any(Kst.genotype=="WT");
    hasAPP = any(Kst.genotype=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP data for %s; skipping.', st);
        continue;
    end

    % -------- 6a) Per-bin means / SEMs for baseline and ambtemp ----------
    meanWT_base  = nan(1,nBins); semWT_base  = nan(1,nBins);
    meanAPP_base = nan(1,nBins); semAPP_base = nan(1,nBins);

    meanWT_amb   = nan(1,nBins); semWT_amb   = nan(1,nBins);
    meanAPP_amb  = nan(1,nBins); semAPP_amb  = nan(1,nBins);

    max_y = nan(1,nBins);

    % ANOVA p-values per bin
    pC  = nan(1,nBins);  % Condition
    pGC = nan(1,nBins);  % Genotype x Condition
    pG  = nan(1,nBins);  % Genotype main

    for b = 1:nBins
        rb = all_bins(b);
        Kb = Kst(Kst.rel_hour_bin==rb, :);

        % -- means / SEMs for plotting --
        for_cond = ["baseline","ambtemp"];
        for g = 1:numel(for_cond)
            cond = for_cond(g);

            valsWT  = Kb.bouts_per_h(Kb.genotype=="WT"  & Kb.condition==cond);
            valsAPP = Kb.bouts_per_h(Kb.genotype=="APP" & Kb.condition==cond);

            valsWT  = valsWT(~isnan(valsWT));
            valsAPP = valsAPP(~isnan(valsAPP));

            if cond=="baseline"
                if ~isempty(valsWT)
                    meanWT_base(b) = mean(valsWT);
                    semWT_base(b)  = std(valsWT)/sqrt(numel(valsWT));
                end
                if ~isempty(valsAPP)
                    meanAPP_base(b) = mean(valsAPP);
                    semAPP_base(b)  = std(valsAPP)/sqrt(numel(valsAPP));
                end
            else
                if ~isempty(valsWT)
                    meanWT_amb(b) = mean(valsWT);
                    semWT_amb(b)  = std(valsWT)/sqrt(numel(valsWT));
                end
                if ~isempty(valsAPP)
                    meanAPP_amb(b) = mean(valsAPP);
                    semAPP_amb(b)  = std(valsAPP)/sqrt(numel(valsAPP));
                end
            end

            % track max y per bin for star placement
            all_vals = [valsWT(:); valsAPP(:)];
            if ~isempty(all_vals)
                max_y(b) = max([max_y(b); all_vals(:)], [], 'omitnan');
            end
        end

        % -- per-bin ANOVA: Condition x Genotype via mixed-effects --
        haveBase = any(Kb.condition=="baseline");
        haveAmb  = any(Kb.condition=="ambtemp");
        haveWT   = any(Kb.genotype=="WT");
        haveAPP2 = any(Kb.genotype=="APP");

        if ~(haveBase && haveAmb && haveWT && haveAPP2)
            continue;
        end

        Kb2 = Kb(~isnan(Kb.bouts_per_h),:);
        if height(Kb2) < 4
            continue;
        end

        Kb2.condition = categorical(Kb2.condition, {'baseline','ambtemp'});
        Kb2.genotype  = categorical(Kb2.genotype);
        Kb2.mouse     = categorical(Kb2.mouse);

        try
            lme = fitlme(Kb2, 'bouts_per_h ~ condition * genotype + (1|mouse)', ...
                         'DummyVarCoding','effects');
            a = anova(lme);

            rowCond = find(strcmp(a.Term,'condition'),1);
            rowGen  = find(strcmp(a.Term,'genotype'),1);
            rowInt  = find(strcmp(a.Term,'condition:genotype'),1);

            if ~isempty(rowCond)
                pC(b) = a.pValue(rowCond);
            end
            if ~isempty(rowGen)
                pG(b) = a.pValue(rowGen);
            end
            if ~isempty(rowInt)
                pGC(b) = a.pValue(rowInt);
            end
        catch ME
            warning('Per-bin ANOVA failed for %s, bin %d: %s', st, rb, ME.message);
        end
    end

    % -------- 6b) FDR across bins for Condition & Interaction ----------
    if useFDR
        pC_FDR  = local_fdr(pC);
        pGC_FDR = local_fdr(pGC);
    else
        pC_FDR  = pC;
        pGC_FDR = pGC;
    end

    % -------- 6c) Plot with baseline overlay ----------
    fig = figure('Color','w','Position',[150 150 900 450]); hold on;

    % Baseline (dashed, light)
    if hasWT
        errorbar(all_bins, meanWT_base, semWT_base, '--', ...
            'Color', COL_WT_LIGHT, 'LineWidth', 1);
    end
    if hasAPP
        errorbar(all_bins, meanAPP_base, semAPP_base, '--', ...
            'Color', COL_APP_LIGHT, 'LineWidth', 1);
    end

    % Ambtemp (solid, stronger)
    if hasWT
        errorbar(all_bins, meanWT_amb, semWT_amb, '-o', ...
            'Color', COL_WT, 'MarkerFaceColor', COL_WT, ...
            'LineWidth',1.5,'MarkerSize',5);
    end
    if hasAPP
        errorbar(all_bins, meanAPP_amb, semAPP_amb, '-o', ...
            'Color', COL_APP, 'MarkerFaceColor', COL_APP, ...
            'LineWidth',1.5,'MarkerSize',5);
    end

    % y limits
    y_global = max([ ...
        max(max_y,[],'omitnan'), ...
        max(meanWT_base+semWT_base,[],'omitnan'), ...
        max(meanAPP_base+semAPP_base,[],'omitnan'), ...
        max(meanWT_amb+semWT_amb,[],'omitnan'), ...
        max(meanAPP_amb+semAPP_amb,[],'omitnan') ], [], 'omitnan');
    if isempty(y_global) || isnan(y_global)
        y_global = 1;
    end
    ylim([0, y_global*1.35]);

    % Shaded ambtemp window [0, ambWindow]
    yl = ylim;
    patch([0 ambWindow ambWindow 0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.93 0.9 0.95], 'EdgeColor','none','FaceAlpha',0.25);
    uistack(findobj(gca,'Type','patch'),'bottom');

    % -------- stars for per-bin ANOVA ----------
    for b = 1:nBins
        rb = all_bins(b);
        if isnan(max_y(b)), continue; end

        y_star_base = max_y(b) + 0.05*y_global;

        % Condition main effect (black star)
        pC_use = pC_FDR(b);
        pI_use = pGC_FDR(b);

        if ~isnan(pC_use) && pC_use < 0.05
            stars = local_p2stars(pC_use);
            text(rb, y_star_base, stars, ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',11,'FontWeight','bold', ...
                 'Color','k');
            y_star_base = y_star_base + 0.05*y_global;
        end

        % Interaction (magenta star)
        if ~isnan(pI_use) && pI_use < 0.05
            stars = local_p2stars(pI_use);
            text(rb, y_star_base, stars, ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',11,'FontWeight','bold', ...
                 'Color',[0.7 0 0.7]);
        end
    end

    xlabel('Hours relative to ambtemp start');
    ylabel(sprintf('Bouts per hour (%s)', st));
    title(sprintf('Ambtemp vs baseline: %s bouts/hour (WT vs APP)', st));

    lg = legend({'WT baseline','APP baseline', ...
                 'WT ambtemp','APP ambtemp'}, ...
                'Location','northoutside','Orientation','horizontal');
    legend boxoff;

    % Text explaining what the stars mean
    star_txt = '* black: Condition (baseline vs ambtemp) p<0.05;   ' + ...
               "* magenta: Condition×Genotype p<0.05  (per hour, ";
    if useFDR
        star_txt = star_txt + "BH-FDR corrected)";
    else
        star_txt = star_txt + "uncorrected)";
    end
    annotation('textbox',[0.1 0.88 0.8 0.06], ...
        'String', star_txt, ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',9);

    set(gca,'Box','off','FontSize',12);

    fname = sprintf('ambtemp_progression_withBaseline_%s_APPvsWT.png', lower(st));
    out_file = fullfile(out_dir, fname);
    saveas(fig, out_file);  % leaves figure open

    OUT.fig_files.(matlab.lang.makeValidName(st)) = out_file;
    OUT.states = [OUT.states; st];

    % -------- 6d) per-bin stats table for output ----------
    stats_tbl = table( ...
        all_bins(:), ...
        pG(:), pC(:), pGC(:), ...
        pC_FDR(:), pGC_FDR(:), ...
        'VariableNames', {'Rel_hour', ...
                          'p_genotype','p_condition','p_GxC', ...
                          'p_condition_FDR','p_GxC_FDR'});
    OUT.stats_perbin.(matlab.lang.makeValidName(st)) = stats_tbl;

    % also write CSV
    csv_out = fullfile(out_dir, sprintf('ambtemp_progression_perbin_ANOVA_%s.csv', lower(st)));
    writetable(stats_tbl, csv_out);

    fprintf('[ANOVA per-bin] %s saved: %s\n', st, csv_out);
end

OUT.success = true;
fprintf('✅ Ambtemp progression with baseline overlay + per-bin ANOVA saved in %s\n', out_dir);

end
% -------------------------------------------------------------------------
function p_adj = local_fdr(p)
% BH-FDR along a vector, preserving NaNs
p_adj = nan(size(p));
valid = ~isnan(p);
pv    = p(valid);
if isempty(pv), return; end
[sp,idx] = sort(pv(:));
m = numel(sp);
adj = sp .* (m./(1:m))';
for k = m-1:-1:1
    adj(k) = min(adj(k),adj(k+1));
end
adj(adj>1) = 1;
out = nan(size(pv));
out(idx) = adj;
p_adj(valid) = out;
end
% -------------------------------------------------------------------------
function stars = local_p2stars(p)
if p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
else
    stars = '*';
end
end
