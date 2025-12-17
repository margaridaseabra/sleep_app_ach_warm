function OUT = run_drugs_boutduration_threePhase_APPvsWT_prism(PERHOUR, out_dir, states_to_use, varargin)
% run_drugs_boutduration_threePhase_APPvsWT_prism
% -------------------------------------------------------------------------
% drugs-only analysis of mean bout duration per state, in THREE windows:
%
%   Phase 1 (0–3 h):  hour_idx in [0 1 2]
%   Phase 2 (3–6 h):  hour_idx in [3 4 5]
%   Washout:          all remaining drugs hours
%
% For each STATE and each mouse:
%   - For each phase p:
%       meanBout_p = sum(dur_s) / sum(bouts_per_h) over that phase
%       (set to NaN if no bouts in that phase)
%
% BETWEEN-genotype stats (WT vs APP) per phase:
%   - Welch t-test
%   - Mann–Whitney (ranksum)
%   - Cohen's d (APP - WT)
%   - SD, CV, SD ratio (APP/WT)
%   - F-test on variances (exploratory)
%
% WITHIN-genotype stats (3-phase paired):
%   For each genotype (WT, APP), paired t-tests on matched mice:
%     - 0–3 h vs 3–6 h
%     - 0–3 h vs washout
%     - 3–6 h vs washout
%
% Figure per state:
%   - X: 3 clusters = (0–3 h, 3–6 h, washout)
%     Within each cluster: WT vs APP bars (mean±SEM)
%   - Dots = individual mice (optional IDs)
%   - Stars (p<0.05, uncorrected) on WT vs APP comparisons, per phase
%
% PRISM EXPORT (per state):
%   drugs_boutdur_threePhase_<state>_forPrism.csv
%
%   Layout:
%     Phase, GENO_mouseID, GENO_mouseID2, ...
%     0_3h,  <mean>,      <mean>,        ...
%     3_6h,  ...
%     washout,...
%
% INPUT
%   PERHOUR       : table from run_group_sleep_architecture (OUT.per_hour)
%                   must contain: hour_idx, dur_s, bouts_per_h,
%                   state, condition, mouse, genotype
%   out_dir       : folder to save CSVs/figures
%   states_to_use : string/cell array of states, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'showIDs'      : true/false, show mouse IDs on dots (default: false)
%   'minNperGroup' : minimum n per genotype to run tests (default: 3)
%
% OUTPUT
%   OUT.success                 : logical
%   OUT.states                  : states processed
%   OUT.between.(STATE)         : table WT vs APP (3 phases)
%   OUT.within.(STATE)          : table within-genotype paired stats (3 phases)
%   OUT.files.(STATE)           : PNG figure paths
%   OUT.prism.(STATE).table     : Prism-style phase×mouse table
%   OUT.prism.(STATE).csv       : CSV path
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

% actually used
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroup',3,@(x)isscalar(x)&&x>=1);

% compatibility with other plotting functions (accepted but mostly ignored here)
addParameter(p,'pointStyle','dots',@(s)ischar(s) || isstring(s));
addParameter(p,'showStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'useFDRforStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'bandAlpha',0.2,@(x)isscalar(x)&&x>=0&&x<=1);

parse(p, varargin{:});

showIDs      = p.Results.showIDs;
minNperGroup = p.Results.minNperGroup;

% optional args kept only so inputParser doesn't complain
pointStyle         = string(p.Results.pointStyle);        %#ok<NASGU>
showStarsFlag      = p.Results.showStars;                 %#ok<NASGU>
useFDRforStarsFlag = p.Results.useFDRforStars;            %#ok<NASGU>
minN_forStats      = p.Results.minNperGroupForStats;      %#ok<NASGU>
bandAlpha          = p.Results.bandAlpha;                 %#ok<NASGU>



PH = PERHOUR;

% ---------- drugs only ----------
if ~ismember('condition', PH.Properties.VariableNames)
    error('PERHOUR must contain a "condition" column.');
end
cond_str    = lower(strtrim(PH.condition));
is_drugs = cond_str == "drugs";
PH = PH(is_drugs, :);

if isempty(PH)
    warning('No drugs rows in PERHOUR. Nothing to do.');
    OUT = struct('success',false,'msg','no drugs data');
    return;
end

% ---------- Check required columns ----------
needed = {'hour_idx','dur_s','bouts_per_h','state','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('PERHOUR table missing required columns for run_drugs_boutduration_threePhase_APPvsWT_prism.');
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

% ---------- Overall bout duration across entire drugs (no phases) ----------
allMice = unique(PH.mouse);
nAllM   = numel(allMice);

% genotype per mouse (WT / APP)
geno_all = strings(nAllM,1);
for i = 1:nAllM
    mID   = allMice(i);
    maskM = PH.mouse == mID;
    g_this = unique(PH.geno_group(maskM));
    if numel(g_this) ~= 1
        warning('Mouse %s has multiple genotypes? Using first.', string(mID));
        geno_all(i) = g_this(1);
    else
        geno_all(i) = g_this;
    end
end

nStates = numel(states_to_use);
overallMean = nan(nStates, nAllM);   % rows = states, cols = mice

for si = 1:nStates
    st = states_to_use(si);
    st_upper = upper(st);

    PH_state = PH(PH.state == st_upper, :);
    if isempty(PH_state)
        continue;
    end

    for i = 1:nAllM
        mID   = allMice(i);
        maskM = PH_state.mouse == mID;
        if ~any(maskM)
            continue;
        end

        tot_dur   = sum(double(PH_state.dur_s(maskM)),       'omitnan');
        tot_bouts = sum(double(PH_state.bouts_per_h(maskM)), 'omitnan');

        if tot_bouts > 0
            overallMean(si,i) = tot_dur / tot_bouts;   % mean bout duration (s)
        else
            overallMean(si,i) = NaN;
        end
    end
end

% Remove states with no data at all (optional, but nicer for Prism)
valid_state_rows = ~all(isnan(overallMean),2);
overallMean_valid = overallMean(valid_state_rows, :);
stateNames_valid  = upper(states_to_use(valid_state_rows));

% Build Prism-style table: rows = states, columns = mice
colNames_raw_overall = cell(1, nAllM);
for i = 1:nAllM
    colNames_raw_overall{i} = sprintf('%s_%s', char(geno_all(i)), char(allMice(i)));
end
colNames_valid_overall = matlab.lang.makeValidName(colNames_raw_overall, 'ReplacementStyle','delete');

T_overallPrism = array2table(overallMean_valid, 'VariableNames', colNames_valid_overall);
T_overallPrism = addvars(T_overallPrism, stateNames_valid(:), ...
    'Before',1, 'NewVariableNames','State');

overall_csvPath = fullfile(out_dir, 'drugs_boutdur_wholeRecording_forPrism.csv');
try
    writetable(T_overallPrism, overall_csvPath);
    fprintf('📄 Prism CSV for WHOLE recording bout duration saved to: %s\n', overall_csvPath);
catch ME
    warning('run_drugs_boutduration_threePhase_APPvsWT_prism:OverallPrismWriteFailed', ...
        'Could not write overall Prism CSV (%s): %s', overall_csvPath, ME.message);
    overall_csvPath = '';
end



% Phase definitions (in hour_idx units)
phaseHours   = {0:2, 3:5, []};     % washout will be filled per state
phaseTags    = ["0_3h","3_6h","washout"];
phaseLabels  = ["0–3 h","3–6 h","Washout"];
nPhases      = numel(phaseHours);

OUT = struct();
OUT.success = true;
OUT.states  = [];
OUT.between = struct();
OUT.within  = struct();
OUT.files   = struct();
OUT.prism   = struct();

% ---------- Loop over states ----------
for s = 1:numel(states_to_use)
    st = states_to_use(s);
    st_upper = upper(st);

    PHs = PH(PH.state == st_upper, :);
    if isempty(PHs)
        warning('No drugs rows for state "%s". Skipping.', st_upper);
        continue;
    end

    % Hours present for THIS state
    h_all   = sort(unique(double(PHs.hour_idx)));
    base12  = unique([0:2, 3:5]);
    washH   = setdiff(h_all, base12);   % washout = any other drugs hours
    phaseHours{3} = washH;

    % ---------- Per-mouse, phase totals ----------
    mice = unique(PHs.mouse);
    nM   = numel(mice);

    mouseID  = strings(nM,1);
    geno_m   = strings(nM,1);
    meanP    = nan(nM, nPhases);   % mouse × phase

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

        for ph = 1:nPhases
            hSet = phaseHours{ph};
            if isempty(hSet)
                continue;
            end
            maskPhase = maskM & ismember(double(PHs.hour_idx), hSet);

            tot_dur   = sum(double(PHs.dur_s(maskPhase)),       'omitnan');
            tot_bouts = sum(double(PHs.bouts_per_h(maskPhase)), 'omitnan');

            if tot_bouts > 0
                meanP(i,ph) = tot_dur / tot_bouts;
            else
                meanP(i,ph) = NaN;
            end
        end
    end

    % Keep only mice with at least one non-NaN phase
    keep = ~all(isnan(meanP),2);
    if ~any(keep)
        warning('All mice have NaN bout durations for state "%s". Skipping.', st_upper);
        continue;
    end
    mouseID = mouseID(keep);
    geno_m  = geno_m(keep);
    meanP   = meanP(keep,:);   % nKept × 3
    nKept   = numel(mouseID);

    hasWT  = any(geno_m=="WT");
    hasAPP = any(geno_m=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP mice left for state "%s". Skipping.', st_upper);
        continue;
    end

    % ---------- BETWEEN-genotype stats (WT vs APP) per phase ----------
    meanWT  = nan(nPhases,1);
    meanAPP = nan(nPhases,1);
    sdWT    = nan(nPhases,1);
    sdAPP   = nan(nPhases,1);
    cvWT    = nan(nPhases,1);
    cvAPP   = nan(nPhases,1);
    nWT_vec = nan(nPhases,1);
    nAPP_vec= nan(nPhases,1);
    p_t     = nan(nPhases,1);
    p_rs    = nan(nPhases,1);
    p_var   = nan(nPhases,1);
    d_eff   = nan(nPhases,1);

    for ph = 1:nPhases
        valsWT  = meanP(geno_m=="WT",  ph);
        valsAPP = meanP(geno_m=="APP", ph);

        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        nWT_vec(ph)  = numel(valsWT);
        nAPP_vec(ph) = numel(valsAPP);

        if ~isempty(valsWT)
            meanWT(ph) = mean(valsWT);
            sdWT(ph)   = std(valsWT);
            cvWT(ph)   = sdWT(ph) / max(eps, meanWT(ph));
        end
        if ~isempty(valsAPP)
            meanAPP(ph) = mean(valsAPP);
            sdAPP(ph)   = std(valsAPP);
            cvAPP(ph)   = sdAPP(ph) / max(eps, meanAPP(ph));
        end

        if nWT_vec(ph) >= minNperGroup && nAPP_vec(ph) >= minNperGroup
            % Welch t-test
            [~, p_t(ph)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            % Mann–Whitney
            p_rs(ph) = ranksum(valsWT, valsAPP);
            % Cohen's d (APP - WT)
            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
            d_eff(ph) = (m2 - m1) / sp;
            % F-test on variances (exploratory)
            try
                [~, pF] = vartest2(valsWT, valsAPP, 'Tail','both');
                p_var(ph) = pF;
            catch
                p_var(ph) = NaN;
            end
        end
    end

    sd_ratio = sdAPP ./ sdWT;
    between_tbl = table( ...
        phaseTags', ...
        nWT_vec, nAPP_vec, ...
        meanWT, meanAPP, ...
        sdWT, sdAPP, ...
        cvWT, cvAPP, ...
        sd_ratio, ...
        p_t, p_rs, p_var, d_eff, ...
        'VariableNames', {'PhaseTag','nWT','nAPP', ...
                          'MeanWT','MeanAPP', ...
                          'SD_WT','SD_APP', ...
                          'CV_WT','CV_APP', ...
                          'SD_ratio_APP_vs_WT', ...
                          'p_ttest','p_ranksum','p_var_Ftest','Cohen_d'});

    fprintf('\n[Between] Three-phase bout duration (state=%s): WT vs APP\n', st_upper);
    disp(between_tbl);

    % ---------- WITHIN-genotype stats (paired, 3 phases) ----------
    % Pairwise comparisons: 0_3 vs 3_6, 0_3 vs wash, 3_6 vs wash
    within_rows = {};
    for gname = ["WT","APP"]
        maskG = geno_m == gname;
        if ~any(maskG), continue; end

        V = meanP(maskG,:);   % nG × 3

        % 0_3 vs 3_6
        v1 = V(:,1); v2 = V(:,2);
        valid12 = ~isnan(v1) & ~isnan(v2);
        n12 = sum(valid12);
        if n12 >= 2
            [~, p12] = ttest(v1(valid12), v2(valid12));
        else
            p12 = NaN;
        end

        % 0_3 vs wash
        v3 = V(:,3);
        valid13 = ~isnan(v1) & ~isnan(v3);
        n13 = sum(valid13);
        if n13 >= 2
            [~, p13] = ttest(v1(valid13), v3(valid13));
        else
            p13 = NaN;
        end

        % 3_6 vs wash
        valid23 = ~isnan(v2) & ~isnan(v3);
        n23 = sum(valid23);
        if n23 >= 2
            [~, p23] = ttest(v2(valid23), v3(valid23));
        else
            p23 = NaN;
        end

        mean0 = mean(v1,'omitnan');
        mean1 = mean(v2,'omitnan');
        mean2 = mean(v3,'omitnan');

        within_rows(end+1,:) = { char(gname), ...
                                 size(V,1), ...
                                 mean0, mean1, mean2, ...
                                 n12, p12, ...
                                 n13, p13, ...
                                 n23, p23 }; %#ok<AGROW>
    end

    if isempty(within_rows)
        within_tbl = table('Size',[0 11], ...
            'VariableTypes',{'string','double','double','double','double', ...
                             'double','double','double','double','double','double'}, ...
            'VariableNames',{'Genotype','nMice', ...
                             'Mean_0_3h','Mean_3_6h','Mean_washout', ...
                             'nPairs_0_3_vs_3_6','p_0_3_vs_3_6', ...
                             'nPairs_0_3_vs_wash','p_0_3_vs_wash', ...
                             'nPairs_3_6_vs_wash','p_3_6_vs_wash'});
    else
        within_tbl = cell2table(within_rows, ...
            'VariableNames',{'Genotype','nMice', ...
                             'Mean_0_3h','Mean_3_6h','Mean_washout', ...
                             'nPairs_0_3_vs_3_6','p_0_3_vs_3_6', ...
                             'nPairs_0_3_vs_wash','p_0_3_vs_wash', ...
                             'nPairs_3_6_vs_wash','p_3_6_vs_wash'});
    end

    fprintf('\n[Within] Three-phase bout duration (state=%s): within-genotype paired t-tests\n', st_upper);
    disp(within_tbl);

    % ---------- Plot for this state (3 clusters) ----------
    figure('Color','w'); hold on;

    clusterCenters = [1.5, 4.5, 7.5];
    Xpos = [ ...
        clusterCenters(1)-0.3, clusterCenters(1)+0.3; ...
        clusterCenters(2)-0.3, clusterCenters(2)+0.3; ...
        clusterCenters(3)-0.3, clusterCenters(3)+0.3  ...
    ]; % rows: phase; cols: WT, APP

    barWidth = 0.5;

    for ph = 1:nPhases
        if hasWT && ~isnan(meanWT(ph))
            bar(Xpos(ph,1), meanWT(ph), barWidth, 'FaceColor', COL_WT,  'EdgeColor','none');
        end
        if hasAPP && ~isnan(meanAPP(ph))
            bar(Xpos(ph,2), meanAPP(ph), barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
        end
    end

    % SEM
    semWT  = sdWT  ./ sqrt(max(1,nWT_vec));
    semAPP = sdAPP ./ sqrt(max(1,nAPP_vec));
    for ph = 1:nPhases
        if hasWT && ~isnan(meanWT(ph)) && nWT_vec(ph)>1
            errorbar(Xpos(ph,1), meanWT(ph), semWT(ph), 'k','LineStyle','none','LineWidth',1);
        end
        if hasAPP && ~isnan(meanAPP(ph)) && nAPP_vec(ph)>1
            errorbar(Xpos(ph,2), meanAPP(ph), semAPP(ph), 'k','LineStyle','none','LineWidth',1);
        end
    end

    % Dots + optional IDs
    jit = 0.15;
    for i = 1:nKept
        id_i = mouseID(i);
        g    = geno_m(i);

        for ph = 1:nPhases
            val = meanP(i,ph);
            if isnan(val), continue; end

            if g=="WT"
                xBase = Xpos(ph,1);
                col   = COL_WT_DOT;
            else
                xBase = Xpos(ph,2);
                col   = COL_APP_DOT;
            end

            xp = xBase + (rand-0.5)*2*jit;
            plot(xp, val, '.', 'Color', col, 'MarkerSize',10);
            if showIDs
                text(xp, val, char(id_i), ...
                    'Rotation',45, ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8, 'Color',col);
            end
        end
    end

    xlim([min(clusterCenters)-1, max(clusterCenters)+1]);
    set(gca,'XTick',clusterCenters,'XTickLabel',phaseLabels,'FontSize',12);
    ylabel(sprintf('Mean bout duration in %s (s)', st_upper));
    title(sprintf('drugs: bout duration (%s), three phases, WT vs APP', st_upper));
    set(gca,'Box','off');

    % Stars for WT vs APP per phase
    yl = ylim; ySpan = yl(2) - yl(1);
    yStar = yl(2) - 0.08*ySpan;
    for ph = 1:nPhases
        if p_t(ph) < 0.05 && nWT_vec(ph)>=minNperGroup && nAPP_vec(ph)>=minNperGroup
            if p_t(ph) < 0.001
                stars = '***';
            elseif p_t(ph) < 0.01
                stars = '**';
            else
                stars = '*';
            end
            x_c = mean(Xpos(ph,:));
            line([Xpos(ph,1) Xpos(ph,2)], [yStar yStar],'Color','k','LineWidth',1);
            text(x_c, yStar+0.02*ySpan, stars, ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'FontSize',12,'FontWeight','bold');
        end
    end

    legend({'WT mean\pmSEM','APP mean\pmSEM'}, ...
           'Location','northoutside','Orientation','horizontal');

    % Save figure
    if showIDs
        id_suffix = '_withIDs';
    else
        id_suffix = '_noIDs';
    end
    fname    = sprintf('drugs_boutdur_threePhase_%s_APPvsWT%s.png', lower(st_upper), id_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    % Store in OUT
    OUT.states = [OUT.states, st_upper];
    key = matlab.lang.makeValidName(st_upper);
    OUT.between.(key) = between_tbl;
    OUT.within.(key)  = within_tbl;
    OUT.files.(key)   = out_file;

    % ---------- PRISM EXPORT for this state ----------
    % meanP: nKept × 3  (mouse × phase)
    prismMat = meanP.';              % 3 × nKept  (phase × mouse)
    PhaseCol = phaseTags(:);         % {"0_3h";"3_6h";"washout"}

    % Column names like WT_mouse12, APP_mouse03 ...
    colNames_raw = cell(1, nKept);
    for i = 1:nKept
        colNames_raw{i} = sprintf('%s_%s', char(geno_m(i)), char(mouseID(i)));
    end
    colNames_valid = matlab.lang.makeValidName(colNames_raw, 'ReplacementStyle','delete');

    Tprism = array2table(prismMat, 'VariableNames', colNames_valid);
    Tprism = addvars(Tprism, PhaseCol, 'Before',1, 'NewVariableNames','Phase');

    csvName = sprintf('drugs_boutdur_threePhase_%s_forPrism.csv', lower(st_upper));
    csvPath = fullfile(out_dir, csvName);
    try
        writetable(Tprism, csvPath);
        fprintf('📄 Prism CSV for bout duration (%s) saved to: %s\n', st_upper, csvPath);
    catch ME
        warning('run_drugs_boutduration_threePhase_APPvsWT_prism:PrismWriteFailed', ...
                'Could not write Prism CSV for %s (%s): %s', st_upper, csvPath, ME.message);
        csvPath = '';
    end

    OUT.prism.(key).table = Tprism;
    OUT.prism.(key).csv   = csvPath;
    % Store overall Prism export in OUT
    OUT.overall_prism.table = T_overallPrism;
    OUT.overall_prism.csv   = overall_csvPath;

end
end
