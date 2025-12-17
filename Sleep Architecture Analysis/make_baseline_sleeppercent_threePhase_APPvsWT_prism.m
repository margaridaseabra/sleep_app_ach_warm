function OUT = make_baseline_sleeppercent_threePhase_APPvsWT_prism(rows_perhr, out_dir, states_to_plot, varargin)
% make_baseline_sleeppercent_threePhase_APPvsWT
% -------------------------------------------------------------------------
% For BASELINE recordings only:
%
%   For each requested STATE (e.g. WK, MA, NREM, REM, SLEEP):
%
%     Phase 1 (0–3 h):  hour_idx in [0 1 2]
%     Phase 2 (3–6 h):  hour_idx in [3 4 5]
%     Washout:          all remaining baseline hours (hour_idx not in above)
%
%   For each mouse & phase:
%     - Denominator per phase = total time in [WK, MA, NREM, REM]
%     - Numerator per phase:
%         * if state is WK/MA/NREM/REM: time in that state
%         * if state is SLEEP        : time in NREM + REM
%
%     pct_phase = 100 * num_phase / den_phase
%
%   Then for each state and phase compare WT vs APP:
%       - Welch t-test (APP vs WT)
%       - Mann–Whitney (ranksum)
%       - Cohen's d (APP - WT)
%       - SD, CV, SD ratio (APP/WT)
%       - F-test on variances (exploratory)
%
%   Plot (per state):
%       X: three clusters = 0–3 h, 3–6 h, washout
%          within each cluster: WT vs APP bars (mean ± SEM)
%       Dots: individual mice (optional IDs)
%       Stars: p_ttest < 0.05 (uncorrected, exploratory)
%
%   Prism export:
%       For each phase, write a CSV:
%           baseline_sleeppercent_0to3h_forPrism.csv
%           baseline_sleeppercent_3to6h_forPrism.csv
%           baseline_sleeppercent_washout_forPrism.csv
%
%       Layout:
%           State, mouseXX, mouseYY, ...
%           WK,   <%,   <%, ...
%           MA,   ...
%           NREM, ...
%           REM,  ...
%           SLEEP,...
%
% INPUT
%   rows_perhr    : OUT.per_hour table from run_group_sleep_architecture
%                   must contain: mouse, genotype, condition, state,
%                   hour_idx, and a duration column (dur_s or total_dur_s)
%   out_dir       : output folder
%   states_to_plot: e.g. ["WK","MA","NREM","REM","SLEEP"]
%
% NAME–VALUE OPTIONS
%   'showIDs'      : true/false, show mouse IDs on dots (default: false)
%   'minNperGroup' : minimum n per genotype to run stats (default: 3)
%
% OUTPUT
%   OUT.success           : logical
%   OUT.states            : states processed
%   OUT.files.(STATE)     : PNG path for each state figure
%   OUT.stats.(STATE)     : table with WT vs APP stats for each phase
%
%   OUT.prism.phaseTags   : {'h0_3','h3_6','washout'}
%   OUT.prism.phaseFiles  : struct with CSV file paths for each phaseTag
%   OUT.prism.phaseTables : struct with tables used for CSVs
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
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroup',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

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

% Define phases
all_hours    = sort(unique(P.hour_idx));
phase1Hours  = 0:3;      % 0–3 h (0–1,1–2,2–3)
phase2Hours  = 4:6;      % 3–6 h (3–4,4–5,5–6)
phase12_all  = unique([phase1Hours(:); phase2Hours(:)]).';
washHours    = setdiff(all_hours, phase12_all);   % washout = everything else

phaseHours   = {phase1Hours, phase2Hours, washHours};
phaseTags    = {'h0_3','h3_6','washout'};
phaseLabels  = {'0–3 h','3–6 h','Washout'};
nPhases      = numel(phaseHours);

% Colors
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

% -------- Per mouse --------
mice = unique(P.mouse);
nM   = numel(mice);

mouseID  = strings(nM,1);
genotype = strings(nM,1);
den      = zeros(nM, nPhases);   % denominators per mouse × phase

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

    % Denominator per phase (WK+MA+NREM+REM)
    for ph = 1:nPhases
        hSet = phaseHours{ph};
        maskBase = maskM & ismember(P.hour_idx, hSet) & ismember(P.state, base_states);
        den(i,ph) = sum(P.(durVar)(maskBase), 'omitnan');
    end
end

OUT = struct();
OUT.success = true;
OUT.states  = [];
OUT.files   = struct();
OUT.stats   = struct();
OUT.prism   = struct();
OUT.prism.phaseTags = phaseTags;

nStates = numel(states_to_plot);
% store all percentages: (state × mouse × phase)
pct_all = nan(nStates, nM, nPhases);

for s = 1:nStates
    st = states_to_plot(s);
    st_upper = upper(st);

    % -------- Compute per-mouse numerators for this state, all phases --------
    pct = nan(nM, nPhases);   % mouse × phase

    for i = 1:nM
        mID   = mouseID(i);
        maskM = P.mouse == mID;

        for ph = 1:nPhases
            hSet = phaseHours{ph};

            if st_upper == "SLEEP"
                % NREM + REM
                maskNum = maskM & ismember(P.hour_idx, hSet) & ...
                          ismember(P.state, sleep_states);
            else
                maskNum = maskM & ismember(P.hour_idx, hSet) & ...
                          (P.state == st_upper);
            end

            num_here = sum(P.(durVar)(maskNum), 'omitnan');

            if den(i,ph) > 0
                pct(i,ph) = 100 * num_here / den(i,ph);
            else
                pct(i,ph) = NaN;
            end
        end
    end

    % Store for Prism export (state index s)
    pct_all(s,:,:) = pct;

    % -------- For stats & plotting, drop mice with all-NaN for this state --------
    keep = ~all(isnan(pct),2);
    if ~any(keep)
        warning('No usable data for state "%s". Skipping.', st);
        continue;
    end

    mID_kept   = mouseID(keep);
    geno_kept  = genotype(keep);
    pct_kept   = pct(keep,:);        % size: nKept × nPhases

    hasWT  = any(geno_kept=="WT");
    hasAPP = any(geno_kept=="APP");
    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state "%s". Skipping.', st);
        continue;
    end

    % -------- Group-wise stats for each phase --------
    phaseNames = ["0_3h","3_6h","wash"];
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
        valsWT  = pct_kept(geno_kept=="WT", ph);
        valsAPP = pct_kept(geno_kept=="APP", ph);

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

        if numel(valsWT) >= minNperGroup && numel(valsAPP) >= minNperGroup
            % t-test
            [~, p_t(ph)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            % ranksum
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

    % -------- Plot for this state (3 phases) --------
    figure('Color','w'); hold on;

    % cluster centers for phases
    clusterCenters = [1.5, 4.5, 7.5];   % 3 clusters
    Xpos = [ ...
        clusterCenters(1)-0.3, clusterCenters(1)+0.3; ...
        clusterCenters(2)-0.3, clusterCenters(2)+0.3; ...
        clusterCenters(3)-0.3, clusterCenters(3)+0.3  ...
    ]; % rows: phase; cols: WT, APP

    barWidth = 0.5;
    hBarWT  = gobjects(nPhases,1);
    hBarAPP = gobjects(nPhases,1);

    for ph = 1:nPhases
        if hasWT
            hBarWT(ph) = bar(Xpos(ph,1), meanWT(ph), barWidth, ...
                             'FaceColor', COL_WT,  'EdgeColor','none');
        end
        if hasAPP
            hBarAPP(ph) = bar(Xpos(ph,2), meanAPP(ph), barWidth, ...
                               'FaceColor', COL_APP, 'EdgeColor','none');
        end
    end

    % Error bars (SEM)
    for ph = 1:nPhases
        if hasWT && ~isnan(meanWT(ph)) && ~isnan(sdWT(ph)) && nWT_vec(ph)>1
            semWT = sdWT(ph) / sqrt(nWT_vec(ph));
            errorbar(Xpos(ph,1), meanWT(ph), semWT, 'k', ...
                     'LineStyle','none','LineWidth',1);
        end
        if hasAPP && ~isnan(meanAPP(ph)) && ~isnan(sdAPP(ph)) && nAPP_vec(ph)>1
            semAPP = sdAPP(ph) / sqrt(nAPP_vec(ph));
            errorbar(Xpos(ph,2), meanAPP(ph), semAPP, 'k', ...
                     'LineStyle','none','LineWidth',1);
        end
    end

    % Dots + optional IDs
    jit = 0.15;
    nKept = numel(mID_kept);
    for i = 1:nKept
        thisID = mID_kept(i);
        g      = geno_kept(i);

        for ph = 1:nPhases
            val = pct_kept(i,ph);
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
                text(xp, val, char(thisID), ...
                     'Rotation',45, ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize',8, 'Color',col);
            end
        end
    end

    % x-axis cosmetics
    xlim([min(clusterCenters)-1, max(clusterCenters)+1]);
    set(gca,'XTick',clusterCenters, ...
            'XTickLabel',phaseLabels, ...
            'FontSize',12);

    ylabel(sprintf('%% time in %s', st_upper));
    title(sprintf('Baseline: %%time in %s, 3 phases (WT vs APP)', st_upper));
    set(gca,'Box','off');

    % Stars for p<0.05 per phase (uncorrected, exploratory)
    yl = ylim;
    ySpan = yl(2) - yl(1);
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

    % Legend: WT grey, APP blue, black line = significance
    legH = [];
    legStr = {};
    if hasWT
        legH(end+1)   = hBarWT(find(~isnan(meanWT),1)); %#ok<FNDSB>
        legStr{end+1} = 'WT mean\pmSEM';
    end
    if hasAPP
        legH(end+1)   = hBarAPP(find(~isnan(meanAPP),1)); %#ok<FNDSB>
        legStr{end+1} = 'APP mean\pmSEM';
    end
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
    fname    = sprintf('baseline_sleeppercent_threePhase_%s_APPvsWT%s.png', lower(st_upper), id_suffix);
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    % Summary table
    sd_ratio = sdAPP ./ sdWT;
    Tstats   = table( ...
        (1:nPhases)', string(phaseNames(:)), ...
        nWT_vec, nAPP_vec, ...
        meanWT, meanAPP, ...
        sdWT, sdAPP, ...
        cvWT, cvAPP, ...
        sd_ratio, ...
        p_t, p_rs, p_var, d_eff, ...
        'VariableNames', {'PhaseIndex','PhaseName','nWT','nAPP', ...
                          'MeanWT','MeanAPP', ...
                          'SD_WT','SD_APP', ...
                          'CV_WT','CV_APP', ...
                          'SD_ratio_APP_vs_WT', ...
                          'p_ttest','p_ranksum','p_var_Ftest','Cohen_d'});

    fprintf('\nThree-phase %%time in %s: WT vs APP\n', st_upper);
    disp(Tstats);

    OUT.states = [OUT.states, st_upper];
    OUT.files.(matlab.lang.makeValidName(st_upper)) = out_file;
    OUT.stats.(matlab.lang.makeValidName(st_upper)) = Tstats;
end

% -------------------------------------------------------------------------
%  Prism EXPORT: 3 CSVs, one per phase
% -------------------------------------------------------------------------
OUT.prism.phaseFiles  = struct();
OUT.prism.phaseTables = struct();

mouseVarNames = matlab.lang.makeValidName(cellstr(mouseID));

for ph = 1:nPhases
    thisPct = squeeze(pct_all(:,:,ph));  % (state × mouse)
    StateCol = states_to_plot(:);

    Tphase = array2table(thisPct, 'VariableNames', mouseVarNames);
    Tphase = addvars(Tphase, StateCol, 'Before', 1, 'NewVariableNames','State');

    switch ph
        case 1
            tag  = 'h0_3';
            base = 'baseline_sleeppercent_0to3h_forPrism.csv';
        case 2
            tag  = 'h3_6';
            base = 'baseline_sleeppercent_3to6h_forPrism.csv';
        otherwise
            tag  = 'washout';
            base = 'baseline_sleeppercent_washout_forPrism.csv';
    end

    csvPath = fullfile(out_dir, base);
    try
        writetable(Tphase, csvPath);
    catch ME
        warning('Could not write Prism CSV for phase %s (%s): %s', tag, csvPath, ME.message);
        csvPath = '';
    end

    OUT.prism.phaseFiles.(tag)  = csvPath;
    OUT.prism.phaseTables.(tag) = Tphase;
end

end
