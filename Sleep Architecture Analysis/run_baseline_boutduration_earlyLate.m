function OUT = run_baseline_boutduration_earlyLate(PERHOUR, out_dir, states_to_use, varargin)
% run_baseline_boutduration_earlyLate
% -------------------------------------------------------------------------
% Baseline-only analysis of mean bout duration per state:
%
%   Window 1 (Early) : first 3 hours  -> hour_idx <= earlyMaxHour (default: 2 => 0–3 h)
%   Window 2 (Rest)  : remaining baseline hours -> hour_idx > earlyMaxHour
%
% For each STATE and each mouse:
%   - Mean bout duration in Early:
%       meanBout_E = sum(dur_s) / sum(bouts_per_h)  over early hours
%   - Mean bout duration in Rest:
%       meanBout_R = sum(dur_s) / sum(bouts_per_h)  over rest hours
%     (set to NaN if no bouts in that window)
%
% Then:
%   1) Between-genotype (WT vs APP) per window (Early/Rest):
%        - Welch t-test
%        - Mann–Whitney (ranksum)
%        - Cohen's d (APP - WT)
%        - SD, CV, SD ratio (APP/WT)
%        - F-test on variances (exploratory)
%
%   2) Within-genotype paired test: Early vs Rest (per mouse):
%        - paired t-test for WT
%        - paired t-test for APP
%
% For each state, we also create a figure:
%   - X: Early vs Rest groups
%     Within each: WT vs APP bars (mean±SEM)
%   - Dots = individual mice (optional IDs)
%   - Stars (p<0.05, uncorrected, exploratory) on WT vs APP comparisons
%
% INPUT
%   PERHOUR       : table from run_group_sleep_architecture (OUT.per_hour)
%                   must contain: hour_idx, dur_s, bouts_per_h,
%                   state, condition, mouse, genotype
%   out_dir       : folder to save CSVs/figures
%   states_to_use : string/cell array of states, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'earlyMaxHour' : max hour_idx included in EARLY (default: 2 -> 0–3 h)
%   'showIDs'      : true/false, show mouse IDs on dots (default: false)
%   'minNperGroup' : minimum n per genotype to run tests (default: 3)
%
% OUTPUT
%   OUT.success                 : logical
%   OUT.states                  : states processed
%   OUT.between.(STATE)         : table WT vs APP (Early/Rest)
%   OUT.within.(STATE)          : table Early vs Rest (WT, APP)
%   OUT.files.(STATE)           : PNG figure paths
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
if nargin < 3 || isempty(states_to_use)
    states_to_use = ["WK","MA","NREM","REM"];
end
states_to_use = string(states_to_use(:)).';

p = inputParser;
addParameter(p,'earlyMaxHour',2,@(x)isscalar(x) && isnumeric(x));
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroup',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

earlyMaxHour = p.Results.earlyMaxHour;   % default: 2 -> 0,1,2 (0–3 h)
showIDs      = p.Results.showIDs;
minNperGroup = p.Results.minNperGroup;

PH = PERHOUR;

% ---------- Baseline only ----------
if ~ismember('condition', PH.Properties.VariableNames)
    error('PERHOUR must contain a "condition" column.');
end
cond_str    = lower(strtrim(PH.condition));
is_baseline = cond_str == "baseline";
PH = PH(is_baseline, :);

if isempty(PH)
    warning('No baseline rows in PERHOUR. Nothing to do.');
    OUT = struct('success',false,'msg','no baseline data');
    return;
end

% ---------- Check required columns ----------
needed = {'hour_idx','dur_s','bouts_per_h','state','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('PERHOUR table missing required columns for run_baseline_boutduration_earlyLate.');
end

% Normalize text columns
PH.state    = string(PH.state);
PH.mouse    = string(PH.mouse);
PH.genotype = string(PH.genotype);

% Collapse any non-WT to APP (if needed)
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% Color scheme
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

OUT = struct();
OUT.success = true;
OUT.states  = [];
OUT.between = struct();
OUT.within  = struct();
OUT.files   = struct();

% ---------- Loop over states ----------
for s = 1:numel(states_to_use)
    st = states_to_use(s);
    st_upper = upper(st);

    PHs = PH(PH.state == st_upper, :);
    if isempty(PHs)
        warning('No baseline rows for state "%s". Skipping.', st_upper);
        continue;
    end

    % --- define early/rest within THIS state subtable ---
    h_s       = double(PHs.hour_idx);
    isEarly_s = h_s <= earlyMaxHour;
    isRest_s  = h_s >  earlyMaxHour;

    if ~any(isEarly_s)
        warning('No early hours for state %s (all rest).', st_upper);
    end

    % ---------- Per-mouse, Early/Rest totals ----------
    mice = unique(PHs.mouse);
    nM   = numel(mice);

    mouseID  = strings(nM,1);
    geno_m   = strings(nM,1);
    meanE    = nan(nM,1);
    meanR    = nan(nM,1);

    for i = 1:nM
        mID    = mice(i);
        maskM  = PHs.mouse == mID;
        mouseID(i) = mID;

        g_this = unique(PHs.geno_group(maskM));
        if numel(g_this) ~= 1
            warning('Mouse %s has multiple genotypes in state %s? Using first.', ...
                string(mID), st_upper);
            geno_m(i) = g_this(1);
        else
            geno_m(i) = g_this;
        end

        % Early window: sum duration and sum bouts (within PHs)
        maskME = maskM & isEarly_s;
        tot_dur_E   = sum(double(PHs.dur_s(maskME)),      'omitnan');
        tot_bouts_E = sum(double(PHs.bouts_per_h(maskME)),'omitnan');
        if tot_bouts_E > 0
            meanE(i) = tot_dur_E / tot_bouts_E;
        else
            meanE(i) = NaN;
        end

        % Rest window
        maskMR = maskM & isRest_s;
        tot_dur_R   = sum(double(PHs.dur_s(maskMR)),      'omitnan');
        tot_bouts_R = sum(double(PHs.bouts_per_h(maskMR)),'omitnan');
        if tot_bouts_R > 0
            meanR(i) = tot_dur_R / tot_bouts_R;
        else
            meanR(i) = NaN;
        end
    end

    % Keep only mice with at least one window non-NaN
    keep = ~(isnan(meanE) & isnan(meanR));
    if ~any(keep)
        warning('All mice have NaN bout durations for state "%s". Skipping.', st_upper);
        continue;
    end
    mouseID = mouseID(keep);
    geno_m  = geno_m(keep);
    meanE   = meanE(keep);
    meanR   = meanR(keep);

    hasWT  = any(geno_m=="WT");
    hasAPP = any(geno_m=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP mice left for state "%s". Skipping.', st_upper);
        continue;
    end

    % ---------- BETWEEN-genotype stats (WT vs APP) per window ----------
    phases  = ["Early","Rest"];
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
        if k==1
            valsWT  = meanE(geno_m=="WT");
            valsAPP = meanE(geno_m=="APP");
        else
            valsWT  = meanR(geno_m=="WT");
            valsAPP = meanR(geno_m=="APP");
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

        if nWT_vec(k) >= minNperGroup && nAPP_vec(k) >= minNperGroup
            % Welch t-test
            [~, p_t(k)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            % Mann–Whitney
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

    sd_ratio = sdAPP ./ sdWT;
    between_tbl = table( ...
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

    fprintf('\n[Between] Baseline bout duration (state=%s): WT vs APP\n', st_upper);
    disp(between_tbl);

    % ---------- WITHIN-genotype stats (Early vs Rest) ----------
    within_rows = {};
    for gname = ["WT","APP"]
        maskG = geno_m == gname;
        if ~any(maskG), continue; end

        vE = meanE(maskG);
        vR = meanR(maskG);

        valid = ~isnan(vE) & ~isnan(vR);
        if sum(valid) < 2
            continue;
        end

        [~, p_paired] = ttest(vE(valid), vR(valid));   % Early vs Rest

        meanE_g = mean(vE(valid),'omitnan');
        meanR_g = mean(vR(valid),'omitnan');
        diff_g  = meanR_g - meanE_g;   % Rest - Early

        within_rows(end+1,:) = {char(gname), sum(valid), ...
                                meanE_g, meanR_g, diff_g, p_paired}; %#ok<AGROW>
    end

    if isempty(within_rows)
        within_tbl = table('Size',[0 6], ...
            'VariableTypes',{'string','double','double','double','double','double'}, ...
            'VariableNames',{'Genotype','n','MeanEarly','MeanRest', ...
                             'Diff_Rest_minus_Early','p_paired'});
    else
        within_tbl = cell2table(within_rows, ...
            'VariableNames',{'Genotype','n','MeanEarly','MeanRest', ...
                             'Diff_Rest_minus_Early','p_paired'});
    end

    fprintf('\n[Within] Baseline bout duration (state=%s): Early vs Rest\n', st_upper);
    disp(within_tbl);

    % ---------- Plot for this state ----------
    figure('Color','w'); hold on;
    Xpos = [1 2 4 5];  % Early WT, Early APP, Rest WT, Rest APP
    barWidth = 0.7;

    % Bars
    if hasWT
        bar(Xpos(1), meanWT(1), barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
        bar(Xpos(3), meanWT(2), barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
    end
    if hasAPP
        bar(Xpos(2), meanAPP(1), barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
        bar(Xpos(4), meanAPP(2), barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
    end

    % Error bars (SEM)
    semWT  = sdWT  ./ sqrt(max(1,nWT_vec));
    semAPP = sdAPP ./ sqrt(max(1,nAPP_vec));

    if hasWT
        errorbar(Xpos(1), meanWT(1),  semWT(1),  'k','LineStyle','none','LineWidth',1);
        errorbar(Xpos(3), meanWT(2),  semWT(2),  'k','LineStyle','none','LineWidth',1);
    end
    if hasAPP
        errorbar(Xpos(2), meanAPP(1), semAPP(1), 'k','LineStyle','none','LineWidth',1);
        errorbar(Xpos(4), meanAPP(2), semAPP(2), 'k','LineStyle','none','LineWidth',1);
    end

    % Dots + optional IDs
    jit = 0.15;
    for i = 1:numel(mouseID)
        id_i = mouseID(i);
        g    = geno_m(i);

        if g=="WT"
            xE = Xpos(1);
            xR = Xpos(3);
            col = COL_WT_DOT;
        else
            xE = Xpos(2);
            xR = Xpos(4);
            col = COL_APP_DOT;
        end

        if ~isnan(meanE(i))
            xe = xE + (rand-0.5)*2*jit;
            plot(xe, meanE(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xe, meanE(i), char(id_i), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end

        if ~isnan(meanR(i))
            xr = xR + (rand-0.5)*2*jit;
            plot(xr, meanR(i), '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xr, meanR(i), char(id_i), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end
    end

    % Axis cosmetics
    xlim([0.5 5.5]);
    earlyLabel = sprintf('Early (h 0–%.0f)', earlyMaxHour+1);  % e.g. 0–3 h if earlyMaxHour=2
    set(gca,'XTick',[1.5 4.5], ...
            'XTickLabel',{earlyLabel, 'Rest (remaining h)'}, ...
            'FontSize',12);

    ylabel(sprintf('Mean bout duration in %s (s)', st_upper));
    title(sprintf('Baseline: Early vs Rest bout duration (%s, WT vs APP)', st_upper));
    set(gca,'Box','off');

    legend({'WT mean\pmSEM','APP mean\pmSEM'}, ...
           'Location','northoutside','Orientation','horizontal');

    % Stars for WT vs APP comparisons (per phase; exploratory)
    yl = ylim;
    ySpan = yl(2) - yl(1);
    yStarEarly = yl(2) - 0.08*ySpan;
    yStarRest  = yl(2) - 0.08*ySpan;

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
                line([Xpos(3) Xpos(4)], [yStarRest yStarRest],'Color','k','LineWidth',1);
                text(x_c, yStarRest+0.02*ySpan, stars, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',12,'FontWeight','bold');
            end
        end
    end

    % Save figure
    if showIDs
        id_suffix = '_withIDs';
    else
        id_suffix = '_noIDs';
    end
    fname    = sprintf('baseline_boutdur_earlyRest_%s_APPvsWT%s.png', lower(st_upper), id_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    % Store in OUT
    OUT.states = [OUT.states, st_upper];
    OUT.between.(matlab.lang.makeValidName(st_upper)) = between_tbl;
    OUT.within.(matlab.lang.makeValidName(st_upper))  = within_tbl;
    OUT.files.(matlab.lang.makeValidName(st_upper))   = out_file;
end
end
