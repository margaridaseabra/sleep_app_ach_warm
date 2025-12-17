function OUT = make_baseline_sleeppercent_earlyLate_APPvsWT(rows_perhr, out_dir, states_to_plot, varargin)
% make_baseline_statepercent_earlyLate_APPvsWT
% -------------------------------------------------------------------------
% For BASELINE recordings only:
%
%   For each requested STATE (e.g. WK, MA, NREM, REM, SLEEP):
%
%     Early phase: hours in earlyHours  (default: 0:3)
%     Late phase : all remaining hours  (baseline hours \ earlyHours)
%
%   For each mouse:
%     - Denominator per phase = total time in [WK, MA, NREM, REM]
%     - Numerator per phase:
%         * if state is WK/MA/NREM/REM: time in that state
%         * if state is SLEEP        : time in NREM + REM
%
%     %state_early = 100 * num_early / den_early
%     %state_late  = 100 * num_late  / den_late
%
%   Then for each state and phase (Early, Late) compare WT vs APP:
%       - Welch t-test (APP vs WT)
%       - Mann–Whitney (ranksum)
%       - Cohen's d (APP - WT)
%       - SD, CV, SD ratio (APP/WT)
%       - F-test on variances (exploratory)
%
%   Plot (per state):
%       X: Early vs Late
%          within each phase: WT vs APP bars (mean ± SEM)
%       Dots: individual mice (optional IDs)
%       Stars: p_ttest < 0.05 (uncorrected, small n, exploratory)
%
% INPUT
%   rows_perhr    : OUT.per_hour table from run_group_sleep_architecture
%                   must contain: mouse, genotype, condition, state,
%                   hour_idx, and a duration column (dur_s or total_dur_s)
%   out_dir       : output folder
%   states_to_plot: e.g. ["WK","MA","NREM","REM","SLEEP"]
%
% NAME–VALUE OPTIONS
%   'earlyHours'   : vector of hour_idx for Early phase (default 0:3)
%   'showIDs'      : true/false, show mouse IDs on dots (default: false)
%   'minNperGroup' : minimum n per genotype to run stats (default: 3)
%
% OUTPUT
%   OUT.success         : logical
%   OUT.states          : states processed
%   OUT.files.(STATE)   : PNG path for each state figure
%   OUT.stats.(STATE)   : table with Early/Late WT vs APP stats
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
if nargin < 3 || isempty(states_to_plot)
    states_to_plot = ["WK","MA","NREM","REM","SLEEP"];
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

% -------- Find duration column (per hour, per state) --------
durVar = '';
if ismember('dur_s', P.Properties.VariableNames)
    durVar = 'dur_s';
elseif ismember('total_dur_s', P.Properties.VariableNames)
    durVar = 'total_dur_s';
else
    error('Could not find duration column (dur_s or total_dur_s) in rows_perhr.');
end

% Base states (used for denominator); sleep components for SLEEP
base_states  = ["WK","MA","NREM","REM"];
sleep_states = ["NREM","REM"];

all_hours = unique(P.hour_idx);
all_hours = sort(all_hours);
lateHours = setdiff(all_hours, earlyHours);

if isempty(lateHours)
    warning('No lateHours found (all hours are in earlyHours).');
end

% Colors
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

% -------- Precompute per-mouse denominators (same for all states) --------
mice = unique(P.mouse);
nM   = numel(mice);

mouseID   = strings(nM,1);
genotype  = strings(nM,1);
denE      = zeros(nM,1);
denL      = zeros(nM,1);

for i = 1:nM
    mID = mice(i);
    maskM = P.mouse == mID;

    g_this = unique(P.genotype(maskM));
    if numel(g_this) ~= 1
        warning('Mouse %s has multiple genotypes? Using first.', string(mID));
        genotype(i) = g_this(1);
    else
        genotype(i) = g_this;
    end
    mouseID(i) = mID;

    % Early denominator (WK+MA+NREM+REM)
    maskEarlyBase = maskM & ismember(P.hour_idx, earlyHours) & ismember(P.state, base_states);
    denE(i) = sum(P.(durVar)(maskEarlyBase), 'omitnan');

    % Late denominator
    maskLateBase  = maskM & ismember(P.hour_idx, lateHours) & ismember(P.state, base_states);
    denL(i) = sum(P.(durVar)(maskLateBase), 'omitnan');
end

OUT = struct();
OUT.success = true;
OUT.states  = [];
OUT.files   = struct();
OUT.stats   = struct();

for s = 1:numel(states_to_plot)
    st = states_to_plot(s);
    st_upper = upper(st);

    % -------- Compute per-mouse numerators for this state --------
    pctE = nan(nM,1);
    pctL = nan(nM,1);

    for i = 1:nM
        mID   = mouseID(i);
        maskM = P.mouse == mID;

        if st_upper == "SLEEP"
            % NREM + REM
            maskEarlyNum = maskM & ismember(P.hour_idx, earlyHours) & ismember(P.state, sleep_states);
            maskLateNum  = maskM & ismember(P.hour_idx, lateHours)  & ismember(P.state, sleep_states);
        else
            maskEarlyNum = maskM & ismember(P.hour_idx, earlyHours) & (P.state == st_upper);
            maskLateNum  = maskM & ismember(P.hour_idx, lateHours)  & (P.state == st_upper);
        end

        numE = sum(P.(durVar)(maskEarlyNum), 'omitnan');
        numL = sum(P.(durVar)(maskLateNum), 'omitnan');

        if denE(i) > 0
            pctE(i) = 100 * numE / denE(i);
        else
            pctE(i) = NaN;
        end
        if denL(i) > 0
            pctL(i) = 100 * numL / denL(i);
        else
            pctL(i) = NaN;
        end
    end

    % Drop mice with no info for this state (both NaN)
    keep = ~(isnan(pctE) & isnan(pctL));
    if ~any(keep)
        warning('No usable data for state "%s". Skipping.', st);
        continue;
    end

    mID      = mouseID(keep);
    geno     = genotype(keep);
    pctE_ks  = pctE(keep);
    pctL_ks  = pctL(keep);

    hasWT  = any(geno=="WT");
    hasAPP = any(geno=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state "%s". Skipping.', st);
        continue;
    end

    % -------- Group-wise stats for Early / Late --------
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
            valsWT  = pctE_ks(geno=="WT");
            valsAPP = pctE_ks(geno=="APP");
        else
            valsWT  = pctL_ks(geno=="WT");
            valsAPP = pctL_ks(geno=="APP");
        end

        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        nWT_vec(k)  = numel(valsWT);
        nAPP_vec(k) = numel(valsAPP);

        if ~isempty(valsWT)
            meanWT(k) = mean(valsWT);
            sdWT(k)   = std(valsWT);
            cvWT(k)   = sdWT(k) / max(eps, meanWT(k));
        end
        if ~isempty(valsAPP)
            meanAPP(k) = mean(valsAPP);
            sdAPP(k)   = std(valsAPP);
            cvAPP(k)   = sdAPP(k) / max(eps, meanAPP(k));
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

    % -------- Plot for this state --------
    figure('Color','w'); hold on;
    Xpos = [1 2 4 5];  % early WT, early APP, late WT, late APP
    barWidth = 0.7;

    % Bars  (store handles!)
    hBarWT  = gobjects(2,1);
    hBarAPP = gobjects(2,1);

    if hasWT
        hBarWT(1) = bar(Xpos(1), meanWT(1), barWidth, ...
                    'FaceColor', COL_WT,  'EdgeColor','none');
        hBarWT(2) = bar(Xpos(3), meanWT(2), barWidth, ...
                    'FaceColor', COL_WT,  'EdgeColor','none');
    end
    if hasAPP
        hBarAPP(1) = bar(Xpos(2), meanAPP(1), barWidth, ...
                        'FaceColor', COL_APP, 'EdgeColor','none');
        hBarAPP(2) = bar(Xpos(4), meanAPP(2), barWidth, ...
                        'FaceColor', COL_APP, 'EdgeColor','none');
    end

    % Dots + optional IDs
    jit = 0.15;
    for i = 1:numel(mID)
        thisID = mID(i);
        g      = geno(i);

        if g=="WT"
            xE = Xpos(1);
            xL = Xpos(3);
            col = COL_WT_DOT;
        else
            xE = Xpos(2);
            xL = Xpos(4);
            col = COL_APP_DOT;
        end

        if ~isnan(pctE_ks(i))
            xe = xE + (rand-0.5)*2*jit;
            plot(xe, pctE_ks(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xe, pctE_ks(i), char(thisID), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end

        if ~isnan(pctL_ks(i))
            xl = xL + (rand-0.5)*2*jit;
            plot(xl, pctL_ks(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xl, pctL_ks(i), char(thisID), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end
    end

    % x-axis cosmetics
    xlim([0.5 5.5]);
    earlyLabel = sprintf('Early (h %s)', strjoin(string(earlyHours),','));
    set(gca,'XTick',[1.5 4.5], ...
            'XTickLabel',{earlyLabel, 'Late (rest)'}, ...
            'FontSize',12);

    ylabel(sprintf('%% time in %s', st_upper));
    title(sprintf('Baseline: early vs late %%time in %s (WT vs APP)', st_upper));
    set(gca,'Box','off');

    legend({'WT mean\pmSEM','APP mean\pmSEM'}, ...
           'Location','northoutside','Orientation','horizontal');

    % Stars for p<0.05 per phase (uncorrected, exploratory)
    yl = ylim;
    ySpan = yl(2) - yl(1);
    yStarEarly = yl(2) - 0.08*ySpan;
    yStarLate  = yl(2) - 0.08*ySpan;

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
                text(x_c, yStarEarly+0.02*ySpan, stars, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            else
                x_c = mean(Xpos(3:4));
                line([Xpos(3) Xpos(4)], [yStarLate yStarLate],'Color','k','LineWidth',1);
                text(x_c, yStarLate+0.02*ySpan, stars, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            end
        end
    end
    % ---- Legend: WT grey, APP blue, black line = significance ----
    legH = [];
    legStr = {};

    if hasWT
        legH(end+1)   = hBarWT(1);          % any WT bar
        legStr{end+1} = 'WT mean\pmSEM';
    end
    if hasAPP
        legH(end+1)   = hBarAPP(1);         % any APP bar
        legStr{end+1} = 'APP mean\pmSEM';
    end

    % dummy line so legend shows a black line for significance
    hSigLeg = plot(NaN,NaN,'k-');
    legH(end+1)   = hSigLeg;
    legStr{end+1} = 'Significant WT vs APP (p<0.05)';

    legend(legH, legStr, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off');
    % Save figure
    if showIDs
        id_suffix = '_withIDs';
    else
        id_suffix = '_noIDs';
    end
    fname    = sprintf('baseline_statepercent_earlyLate_%s_APPvsWT%s.png', lower(st_upper), id_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    % Summary table
    sd_ratio = sdAPP ./ sdWT;
    Tstats   = table( ...
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

    fprintf('\nEarly vs Late %%time in %s: WT vs APP\n', st_upper);
    disp(Tstats);

    OUT.states = [OUT.states, st_upper];
    OUT.files.(matlab.lang.makeValidName(st_upper)) = out_file;
    OUT.stats.(matlab.lang.makeValidName(st_upper)) = Tstats;
end
end
