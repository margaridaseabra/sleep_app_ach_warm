function OUT = make_baseline_bouts_earlyLate_APPvsWT(rows_perhr, out_dir, states_to_plot, varargin)
% make_baseline_bouts_earlyLate_APPvsWT
% -------------------------------------------------------------------------
% For BASELINE recordings only:
%   For each requested STATE (WK, MA, NREM, REM):
%
%   1) Compute, per mouse:
%         Early_mean = mean(bouts_per_h) over early_hours (default 0–3)
%         Late_mean  = mean(bouts_per_h) over late_hours  (default ≥4)
%
%   2) Compare WT vs APP separately for each phase:
%         - Welch t-test (APP vs WT)
%         - Mann–Whitney (ranksum)
%         - Cohen's d (APP - WT)
%         - Group SD, CV, SD ratio, F-test on variances (exploratory)
%
%   3) Plot per state:
%         - Bars: mean±SEM (WT vs APP) for Early and Late
%         - Dots: individual mice (optional IDs)
%         - Stars if p < 0.05 for Early / Late comparison
%
% INPUT
%   rows_perhr    : table from run_group_sleep_architecture (group_per_hour)
%   out_dir       : output folder
%   states_to_plot: string/cell array, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'earlyHours'     : vector of hour_idx for early phase (default 0:3)
%   'showIDs'        : true/false, label dots with mouse ID (default: false)
%   'minNperGroup'   : min n per genotype to run stats (default: 3)
%
% OUTPUT
%   OUT.success         : logical
%   OUT.states          : states processed
%   OUT.files.(STATE)   : PNG path for each state plot
%   OUT.stats.(STATE)   : table with early/late WT vs APP stats
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
if nargin < 3 || isempty(states_to_plot)
    states_to_plot = ["WK","MA","NREM","REM"];
end
states_to_plot = string(states_to_plot(:)).';

p = inputParser;
addParameter(p,'earlyHours',0:3,@(v)isnumeric(v));
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroup',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

earlyHours   = p.Results.earlyHours;
showIDs      = p.Results.showIDs;
minNperGroup = p.Results.minNperGroup;

% -------- Baseline only --------
P = rows_perhr;
cond_str    = lower(strtrim(P.condition));
is_baseline = cond_str == "baseline";
P = P(is_baseline, :);

if isempty(P)
    warning('No baseline rows in rows_perhr. Nothing to do.');
    OUT = struct('success',false,'msg','no baseline data');
    return;
end

% Colors
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

OUT = struct();
OUT.success = true;
OUT.states  = [];
OUT.files   = struct();
OUT.stats   = struct();

for s = 1:numel(states_to_plot)
    st = states_to_plot(s);
    Pst = P(P.state == st, :);
    if isempty(Pst)
        warning('No baseline rows for state "%s". Skipping.', st);
        continue;
    end

    G = Pst;

    hasWT  = any(G.genotype=="WT");
    hasAPP = any(G.genotype=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state "%s". Skipping.', st);
        continue;
    end

    % ---------- Build per-mouse early/late means ----------
    mice = unique(G.mouse);
    nM   = numel(mice);

    mouseID   = strings(nM,1);
    genotype  = strings(nM,1);
    earlyMean = nan(nM,1);
    lateMean  = nan(nM,1);

    all_hours = unique(G.hour_idx);
    all_hours = sort(all_hours);
    lateHours = setdiff(all_hours, earlyHours);

    for i = 1:nM
        mID = mice(i);
        maskM = G.mouse == mID;

        g_this = unique(G.genotype(maskM));
        if numel(g_this) ~= 1
            warning('Mouse %s has multiple genotypes? Using first.', string(mID));
            genotype(i) = g_this(1);
        else
            genotype(i) = g_this;
        end
        mouseID(i) = mID;

        % early
        maskEarly = maskM & ismember(G.hour_idx, earlyHours);
        valsE = G.bouts_per_h(maskEarly);
        valsE = valsE(~isnan(valsE));
        if ~isempty(valsE)
            earlyMean(i) = mean(valsE);
        end

        % late
        maskLate = maskM & ismember(G.hour_idx, lateHours);
        valsL = G.bouts_per_h(maskLate);
        valsL = valsL(~isnan(valsL));
        if ~isempty(valsL)
            lateMean(i) = mean(valsL);
        end
    end

    % Keep only mice with at least one phase non-NaN
    keep = ~(isnan(earlyMean) & isnan(lateMean));
    mouseID  = mouseID(keep);
    genotype = genotype(keep);
    earlyMean = earlyMean(keep);
    lateMean  = lateMean(keep);

    % ---------- Compute group stats for early / late ----------
    phases = ["early","late"];
    meanWT  = nan(2,1);
    meanAPP = nan(2,1);
    sdWT    = nan(2,1);
    sdAPP   = nan(2,1);
    cvWT    = nan(2,1);
    cvAPP   = nan(2,1);
    nWT_vec = nan(2,1);
    nAPP_vec= nan(2,1);
    p_t     = nan(2,1);
    p_rs    = nan(2,1);
    p_var   = nan(2,1);
    d_eff   = nan(2,1);

    for k = 1:2
        if k == 1
            valsWT  = earlyMean(genotype=="WT");
            valsAPP = earlyMean(genotype=="APP");
        else
            valsWT  = lateMean(genotype=="WT");
            valsAPP = lateMean(genotype=="APP");
        end

        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        nWT_vec(k)  = numel(valsWT);
        nAPP_vec(k) = numel(valsAPP);

        if ~isempty(valsWT)
            meanWT(k) = mean(valsWT);
            sdWT(k)   = std(valsWT);
            cvWT(k)   = sdWT(k) / meanWT(k);
        end
        if ~isempty(valsAPP)
            meanAPP(k) = mean(valsAPP);
            sdAPP(k)   = std(valsAPP);
            cvAPP(k)   = sdAPP(k) / meanAPP(k);
        end

        if numel(valsWT) >= minNperGroup && numel(valsAPP) >= minNperGroup
            % t-test
            [~, p_t(k)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            % ranksum
            p_rs(k) = ranksum(valsWT, valsAPP);
            % Cohen's d (APP - WT)
            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
            d_eff(k) = (m2 - m1) / sp;
            % F-test on variances (exploratory)
            try
                [~, pF] = vartest2(valsWT, valsAPP, 'Tail','both');
                p_var(k) = pF;
            catch
                p_var(k) = NaN;
            end
        end
    end

    % ---------- Plot ----------
    figure('Color','w'); hold on;
    Xpos = [1 2 4 5]; % early WT, early APP, late WT, late APP
    barWidth = 0.7;

    % bars
    bar(Xpos(1), meanWT(1), barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
    bar(Xpos(2), meanAPP(1),barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
    bar(Xpos(3), meanWT(2), barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
    bar(Xpos(4), meanAPP(2),barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');

    % error bars (SEM)
    semWT = sdWT ./ sqrt(max(1,nWT_vec));
    semAPP = sdAPP ./ sqrt(max(1,nAPP_vec));

    errorbar(Xpos(1), meanWT(1), semWT(1), 'k','LineStyle','none','LineWidth',1);
    errorbar(Xpos(2), meanAPP(1),semAPP(1),'k','LineStyle','none','LineWidth',1);
    errorbar(Xpos(3), meanWT(2), semWT(2), 'k','LineStyle','none','LineWidth',1);
    errorbar(Xpos(4), meanAPP(2),semAPP(2),'k','LineStyle','none','LineWidth',1);

    % dots + optional IDs
    jit = 0.15;
    for i = 1:numel(mouseID)
        mID = mouseID(i);
        g   = genotype(i);

        % x pos for this mouse in early/late
        if g=="WT"
            xE = Xpos(1);
            xL = Xpos(3);
            col = COL_WT_DOT;
        else
            xE = Xpos(2);
            xL = Xpos(4);
            col = COL_APP_DOT;
        end

        if ~isnan(earlyMean(i))
            xe = xE + (rand-0.5)*2*jit;
            plot(xe, earlyMean(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xe, earlyMean(i), char(mID), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end

        if ~isnan(lateMean(i))
            xl = xL + (rand-0.5)*2*jit;
            plot(xl, lateMean(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xl, lateMean(i), char(mID), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end
    end

    % x-axis cosmetics
    xlim([0.5 5.5]);
    set(gca,'XTick',[1.5 4.5], ...
            'XTickLabel',{sprintf('Early (h %s)',join(string(earlyHours),',')), ...
                          'Late (rest)'}, ...
            'FontSize',12);
    ylabel(sprintf('Mean bouts/hour (%s)', st));
    title(sprintf('Baseline: early vs late %s bouts/hour (WT vs APP)', st));
    set(gca,'Box','off');

    % stars for p<0.05 per phase
    yl = ylim;
    yStarBase = yl(2);
    yStarEarly = yStarBase - 0.05*(yl(2)-yl(1));
    yStarLate  = yStarBase - 0.05*(yl(2)-yl(1));

    for k = 1:2
        if p_t(k) < 0.05 && nWT_vec(k)>=minNperGroup && nAPP_vec(k)>=minNperGroup
            if p_t(k) < 0.001
                stars = '***';
            elseif p_t(k) < 0.01
                stars = '**';
            else
                stars = '*';
            end
            if k==1
                x_c = mean(Xpos(1:2));
                line([Xpos(1) Xpos(2)], [yStarEarly yStarEarly],'Color','k','LineWidth',1);
                text(x_c, yStarEarly+(0.02*(yl(2)-yl(1))), stars, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            else
                x_c = mean(Xpos(3:4));
                line([Xpos(3) Xpos(4)], [yStarLate yStarLate],'Color','k','LineWidth',1);
                text(x_c, yStarLate+(0.02*(yl(2)-yl(1))), stars, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            end
        end
    end

    % save figure
    if showIDs
        id_suffix = '_withIDs';
    else
        id_suffix = '_noIDs';
    end
    fname = sprintf('baseline_bouts_earlyLate_%s_APPvsWT%s.png', lower(st), id_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    OUT.states = [OUT.states, st];
    OUT.files.(matlab.lang.makeValidName(st)) = out_file;

    % summary table
    sd_ratio = sdAPP ./ sdWT;

    Tstats = table( ...
        phases', ...
        nWT_vec, nAPP_vec, ...
        meanWT, meanAPP, ...
        sdWT, sdAPP, ...
        cvWT, cvAPP, ...
        sd_ratio, ...
        p_t, p_rs, p_var, d_eff, ...
        'VariableNames', {'Phase','nWT','nAPP', ...
                          'MeanWT','MeanAPP', ...
                          'SD_WT','SD_APP', ...
                          'CV_WT','CV_APP', ...
                          'SD_ratio_APP_vs_WT', ...
                          'p_ttest','p_ranksum','p_var_Ftest','Cohen_d'});

    fprintf('\nEarly vs Late bouts/hour (%s): WT vs APP\n', st);
    disp(Tstats);

    OUT.stats.(matlab.lang.makeValidName(st)) = Tstats;
end
end
