function OUT = run_MA_fragmentation_threePhase_from_perhour_prism(PERHOUR, out_dir, varargin)
% run_MA_fragmentation_threePhase_from_perhour_prism
% -------------------------------------------------------------------------
% MA fragmentation index over THREE phases, using the PER-HOUR table:
%
%   FI_NREM = (# MA bouts) / (hours of NREM sleep)
%
% Phases (in hour_idx):
%   Phase 1 (0–3 h):   0,1,2
%   Phase 2 (3–6 h):   3,4,5
%   Phase 3 (Washout): all other drugs hours
%
% For each mouse & phase:
%   - Numerator = total number of MA bouts in that phase
%   - Denominator = total NREM duration (h) in that phase
%   - FI_NREM_phase = MA_bouts_phase / NREM_hours_phase
%
%   If NREM_hours_phase == 0 => FI = NaN for that mouse/phase.
%
% BETWEEN-genotype (WT vs APP) stats, per phase:
%   - Welch t-test
%   - Mann–Whitney (ranksum)
%   - Cohen's d (APP - WT)
%
% Plot:
%   - 3 clusters: 0–3 h, 3–6 h, washout
%   - Within each cluster: WT vs APP bars (mean±SEM) + dots
%
% PRISM EXPORT:
%   drugs_MAfrag_threePhase_forPrism.csv
%
%   Layout:
%       Phase, WT_mouseXX, WT_mouseYY, APP_mouseZZ, ...
%       0_3h,  FI,         FI,         FI,          ...
%       3_6h,  ...
%       washout,...
%
% INPUT
%   PERHOUR : OUT.per_hour from run_group_sleep_architecture
%             must contain:
%                 condition, state, mouse, genotype,
%                 hour_idx, dur_s,
%                 AND EITHER n_bouts OR bouts_per_h
%   out_dir : output directory
%
% NAME–VALUE OPTIONS
%   'showIDs' : true/false, show mouse IDs next to dots (default: false)
%
% OUTPUT
%   OUT.success             : logical
%   OUT.phaseTags           : ["0_3h","3_6h","washout"]
%   OUT.per_mouse_table     : table (mouse, genotype, FI_0_3h, FI_3_6h, FI_washout)
%   OUT.between_stats       : per-phase WT vs APP stats table
%   OUT.fig_file            : PNG path of the plot
%   OUT.prism.table         : Prism-style table (Phase × mouse)
%   OUT.prism.csv           : CSV path for Prism
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

p = inputParser;
addParameter(p,'showIDs',false,@(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
showIDs = p.Results.showIDs;

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
needed_basic = {'hour_idx','state','mouse','genotype'};
if ~all(ismember(needed_basic, PH.Properties.VariableNames))
    error('PERHOUR missing one of required columns: hour_idx, state, mouse, genotype.');
end

if ~ismember('dur_s', PH.Properties.VariableNames)
    error('PERHOUR must contain dur_s for NREM duration.');
end

has_n_bouts   = ismember('n_bouts', PH.Properties.VariableNames);
has_bouts_per = ismember('bouts_per_h', PH.Properties.VariableNames);

if ~has_n_bouts && ~has_bouts_per
    error('PERHOUR must contain either n_bouts or bouts_per_h to compute MA bouts.');
end

% ---------- Normalize text columns ----------
PH.state     = string(PH.state);
PH.mouse     = string(PH.mouse);
PH.genotype  = string(PH.genotype);
PH.condition = string(PH.condition);

% Collapse any non-WT to APP (if you only care WT vs APP)
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% ---------- Phase definitions ----------
all_hours = sort(unique(double(PH.hour_idx)));
base12    = unique([0 1 2 3 4 5]);  % first 6 hours
washH     = setdiff(all_hours, base12);

phaseHours  = {0:2, 3:5, washH};
phaseTags   = ["0_3h","3_6h","washout"];
phaseLabels = ["0–3 h","3–6 h","Washout"];
nPhases     = numel(phaseHours);

% ---------- Precompute per-row n_bouts (for MA rows) if needed ----------
if has_n_bouts
    n_bouts_row = double(PH.n_bouts);
else
    % Estimate n_bouts from bouts_per_h * (dur_s / 3600)
    hrFrac       = double(PH.dur_s) / 3600;
    n_bouts_row  = double(PH.bouts_per_h) .* hrFrac;
end

% ---------- Per-mouse fragmentation per phase ----------
mice = unique(PH.mouse);
nM   = numel(mice);

mouseID = strings(nM,1);
geno_m  = strings(nM,1);
FI_mat  = nan(nM, nPhases);   % mouse × phase

for i = 1:nM
    mID        = mice(i);
    mouseID(i) = mID;

    maskM = PH.mouse == mID;

    g_this = unique(PH.geno_group(maskM));
    if numel(g_this) ~= 1
        warning('Mouse %s has multiple genotypes? Using first.', string(mID));
        geno_m(i) = g_this(1);
    else
        geno_m(i) = g_this;
    end

    for ph = 1:nPhases
        hSet = phaseHours{ph};
        if isempty(hSet)
            continue;
        end

        % Numerator: MA bouts in this phase
        maskMA = maskM & PH.state=="MA" & ismember(double(PH.hour_idx), hSet);
        n_MA_phase = sum(n_bouts_row(maskMA), 'omitnan');

        % Denominator: NREM duration (in hours) in this phase
        maskNREM   = maskM & PH.state=="NREM" & ismember(double(PH.hour_idx), hSet);
        durNREM_s  = sum(double(PH.dur_s(maskNREM)), 'omitnan');
        hoursNREM  = durNREM_s / 3600;

        if hoursNREM > 0
            FI_mat(i, ph) = n_MA_phase / hoursNREM;
        else
            FI_mat(i, ph) = NaN;
        end
    end
end

% Keep only mice with at least one non-NaN FI
keep = ~all(isnan(FI_mat),2);
if ~any(keep)
    warning('All mice have NaN fragmentation index in all phases. Nothing to do.');
    OUT = struct('success',false,'msg','no FI data');
    return;
end
mouseID = mouseID(keep);
geno_m  = geno_m(keep);
FI_mat  = FI_mat(keep,:);
nKept   = numel(mouseID);

hasWT  = any(geno_m=="WT");
hasAPP = any(geno_m=="APP");

if ~hasWT && ~hasAPP
    warning('No WT or APP mice present after filtering. Nothing to do.');
    OUT = struct('success',false,'msg','no WT/APP data');
    return;
end

% ---------- BETWEEN-genotype stats per phase ----------
meanWT  = nan(nPhases,1);
meanAPP = nan(nPhases,1);
sdWT    = nan(nPhases,1);
sdAPP   = nan(nPhases,1);
semWT   = nan(nPhases,1);
semAPP  = nan(nPhases,1);
nWT_vec = nan(nPhases,1);
nAPP_vec= nan(nPhases,1);
p_t     = nan(nPhases,1);
p_rs    = nan(nPhases,1);
d_eff   = nan(nPhases,1);

for ph = 1:nPhases
    valsWT  = FI_mat(geno_m=="WT",  ph);
    valsAPP = FI_mat(geno_m=="APP", ph);

    valsWT  = valsWT(~isnan(valsWT));
    valsAPP = valsAPP(~isnan(valsAPP));

    nWT_vec(ph)  = numel(valsWT);
    nAPP_vec(ph) = numel(valsAPP);

    if ~isempty(valsWT)
        meanWT(ph) = mean(valsWT);
        sdWT(ph)   = std(valsWT);
        semWT(ph)  = sdWT(ph) / max(1,sqrt(numel(valsWT)));
    end
    if ~isempty(valsAPP)
        meanAPP(ph) = mean(valsAPP);
        sdAPP(ph)   = std(valsAPP);
        semAPP(ph)  = sdAPP(ph) / max(1,sqrt(numel(valsAPP)));
    end

    if nWT_vec(ph) >= 2 && nAPP_vec(ph) >= 2
        [~, p_t(ph)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
        p_rs(ph)     = ranksum(valsWT, valsAPP);

        m1 = mean(valsWT);  m2 = mean(valsAPP);
        s1 = std(valsWT);   s2 = std(valsAPP);
        n1 = numel(valsWT); n2 = numel(valsAPP);
        sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
        d_eff(ph) = (m2 - m1) / sp;
    end
end

between_tbl = table( ...
    phaseTags', ...
    nWT_vec, nAPP_vec, ...
    meanWT, meanAPP, ...
    sdWT, sdAPP, ...
    semWT, semAPP, ...
    p_t, p_rs, d_eff, ...
    'VariableNames', {'PhaseTag','nWT','nAPP', ...
                      'MeanWT','MeanAPP', ...
                      'SD_WT','SD_APP', ...
                      'SEM_WT','SEM_APP', ...
                      'p_ttest','p_ranksum','Cohen_d'});

fprintf('\n=== MA fragmentation (MA per NREM hour), 3 phases (drugs) ===\n');
disp(between_tbl);

% ---------- Plot: 3 clusters WT vs APP ----------
COL_WT      = [0.6 0.6 0.6];
COL_APP     = [0.39 0.58 0.93];
COL_WT_DOT  = [0.3 0.3 0.3];
COL_APP_DOT = [0.1 0.2 0.6];

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
        if nWT_vec(ph) > 1
            errorbar(Xpos(ph,1), meanWT(ph), semWT(ph), 'k','LineStyle','none','LineWidth',1);
        end
    end
    if hasAPP && ~isnan(meanAPP(ph))
        bar(Xpos(ph,2), meanAPP(ph), barWidth, 'FaceColor', COL_APP, 'EdgeColor','none');
        if nAPP_vec(ph) > 1
            errorbar(Xpos(ph,2), meanAPP(ph), semAPP(ph), 'k','LineStyle','none','LineWidth',1);
        end
    end
end

% Dots
jit = 0.15;
for i = 1:nKept
    id_i = mouseID(i);
    g    = geno_m(i);

    for ph = 1:nPhases
        val = FI_mat(i,ph);
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
ylabel('MA bouts per hour of NREM sleep');
title('drugs MA fragmentation (MA per NREM hour), 3 phases, WT vs APP');
set(gca,'Box','off');

% Add stars for WT vs APP per phase (based on p_t)
yl = ylim; ySpan = yl(2) - yl(1);
yStar = yl(2) - 0.08*ySpan;
for ph = 1:nPhases
    if isnan(p_t(ph)) || p_t(ph) >= 0.05
        continue;
    end

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

legend({'WT mean\pmSEM','APP mean\pmSEM'}, ...
       'Location','northoutside','Orientation','horizontal');

fig_file = fullfile(out_dir, 'drugs_MAfrag_threePhase_WTvsAPP.png');
saveas(gcf, fig_file);

% ---------- Per-mouse wide table (for your sanity) ----------
per_mouse_tbl = table(mouseID, geno_m, ...
    FI_mat(:,1), FI_mat(:,2), FI_mat(:,3), ...
    'VariableNames', {'Mouse','Genotype','FI_0_3h','FI_3_6h','FI_washout'});

% ---------- Prism wide table: Phase × mouse ----------
prismMat = FI_mat.';              % 3 × nKept
PhaseCol = phaseTags(:);          % "0_3h","3_6h","washout"

colNames_raw = cell(1, nKept);
for i = 1:nKept
    colNames_raw{i} = sprintf('%s_%s', char(geno_m(i)), char(mouseID(i)));
end
colNames_valid = matlab.lang.makeValidName(colNames_raw, 'ReplacementStyle','delete');

Tprism = array2table(prismMat, 'VariableNames', colNames_valid);
Tprism = addvars(Tprism, PhaseCol, 'Before',1, 'NewVariableNames','Phase');

csvName = 'drugs_MAfrag_threePhase_forPrism.csv';
csvPath = fullfile(out_dir, csvName);
try
    writetable(Tprism, csvPath);
    fprintf('📄 Prism 3-phase MA fragmentation CSV saved to: %s\n', csvPath);
catch ME
    warning('run_MA_fragmentation_threePhase_from_perhour_prism:PrismWriteFailed', ...
        'Could not write Prism CSV (%s): %s', csvPath, ME.message);
    csvPath = '';
end

% ---------- Pack OUT ----------
OUT = struct();
OUT.success         = true;
OUT.phaseTags       = phaseTags;
OUT.per_mouse_table = per_mouse_tbl;
OUT.between_stats   = between_tbl;
OUT.fig_file        = fig_file;
OUT.prism           = struct('table', Tprism, 'csv', csvPath);

end
