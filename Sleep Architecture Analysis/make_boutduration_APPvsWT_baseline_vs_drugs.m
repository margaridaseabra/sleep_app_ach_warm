function OUT = make_boutduration_APPvsWT_baseline_vs_drugs(rows_overall_drugs, out_dir, states_to_use, varargin)
% make_boutduration_APPvsWT_baseline_vs_drugs
% -------------------------------------------------------------------------
% PERHOUR is the OVERALL table from run_group_sleep_architecture_ambtemp_overall:
%   state, n_bouts, total_dur_s, mean_bout_dur_s,
%   file, date, condition, mouse, genotype
%
% Plot:
%   4 bars per state:
%       1) WT baseline
%       2) WT drugs
%       3) APP baseline
%       4) APP drugs
%
% Stats (per state):
%   ***ONLY*** APP baseline vs APP drugs:
%     - Prefer paired t-test (per mouse) if baseline & drugs both present
%     - Otherwise fall back to unpaired t-test
%   Stars shown ONLY above APP baseline vs APP drugs:
%       p < 0.05  -> *
%       p < 0.01  -> **
%       p < 0.001 -> ***
%
% WT bars are shown for context but NOT used in the tests.
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

% ---------- optional args ----------
p = inputParser;
addParameter(p,'showIDs',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

showIDs              = p.Results.showIDs;
minNperGroupForStats = p.Results.minNperGroupForStats;

PH = rows_overall_drugs;

% ---------- 1) Keep only BASELINE and DRUGS ----------
if ~ismember('condition', PH.Properties.VariableNames)
    error('Input table must contain a "condition" column.');
end

PH.condition = string(lower(strtrim(PH.condition)));

% map any drug-like labels to "drugs"
isDrugLike = (PH.condition == "drug")    | ...
             (PH.condition == "drugs")   | ...
             (PH.condition == "washout");
PH.condition(isDrugLike) = "drugs";

% keep only baseline and drugs
keep_cond = (PH.condition == "baseline") | (PH.condition == "drugs");
PH        = PH(keep_cond, :);

if isempty(PH)
    warning('No baseline/drugs rows in input table. Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline/drugs data');
    return;
end

% ---------- 2) Check required columns ----------
needed = {'total_dur_s','n_bouts','state','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('Input table missing required columns for make_boutduration_APPvsWT_baseline_vs_drugs.');
end

% Normalize text columns
PH.state    = string(PH.state);
PH.mouse    = string(PH.mouse);
PH.genotype = string(PH.genotype);

% Collapse any non-WT to APP
geno = PH.genotype;
geno(geno ~= "WT") = "APP";
PH.geno_group = geno;

% ---------- 3) Restrict to desired states ----------
avail_states = unique(PH.state);
state_order  = states_to_use(ismember(states_to_use, avail_states));

if isempty(state_order)
    warning('None of the requested states found in input data.');
    OUT = struct('success', false, 'msg', 'no states found');
    return;
end

% Fixed condition/genotype order for plotting / stats
cond_order = ["baseline","drugs"];
geno_order = ["WT","APP"];

% ---------- 4) Per-mouse mean bout duration per state & condition ----------
[gid, m, gen, cond, st] = findgroups( ...
    PH.mouse, PH.geno_group, PH.condition, PH.state);

tot_dur   = splitapply(@(x) sum(double(x),'omitnan'), PH.total_dur_s, gid);
tot_bouts = splitapply(@(x) sum(double(x),'omitnan'), PH.n_bouts,      gid);

mean_bout = tot_dur ./ max(tot_bouts, 1);
mean_bout(tot_bouts == 0) = NaN;

Tmouse = table(m, gen, cond, st, mean_bout, ...
    'VariableNames', {'mouse','geno','condition','state','mean_bout_dur_s'});

% Keep only states & conditions of interest
Tmouse = Tmouse(ismember(Tmouse.state, state_order) & ...
                ismember(Tmouse.condition, cond_order), :);

if isempty(Tmouse)
    warning('No per-mouse bout durations for requested states/conditions.');
    OUT = struct('success', false, 'msg', 'no per-mouse data');
    return;
end

% ---------- Colors ----------
% 1: WT baseline
% 2: WT drugs
% 3: APP baseline
% 4: APP drugs
COL_BAR = [
    0.6  0.6  0.6 ;  % WT baseline
    0.3  0.3  0.3 ;  % WT drugs (darker)
    0.39 0.58 0.93;  % APP baseline
    0.19 0.36 0.70]; % APP drugs (darker)

COL_DOT = [
    0.2  0.2  0.2 ;  % WT baseline
    0.1  0.1  0.1 ;  % WT drugs
    0.1  0.2  0.6 ;  % APP baseline
    0.05 0.1  0.4]; % APP drugs

% ---------- 5) Compute per-cell stats + APP-only tests ----------
nStates = numel(state_order);
nConds  = numel(cond_order);
nGen    = numel(geno_order);

% mean/SD/n arrays: (state, condition, genotype)
meanVals = nan(nStates, nConds, nGen);
sdVals   = nan(nStates, nConds, nGen);
nVals    = nan(nStates, nConds, nGen);

% For star placement
max_y_state = nan(1, nStates);

% Per-state summary (WT + APP)
n_WT_baseline   = nan(1, nStates);
n_WT_drugs      = nan(1, nStates);
n_APP_baseline  = nan(1, nStates);
n_APP_drugs     = nan(1, nStates);

mean_WT_baseline   = nan(1, nStates);
mean_WT_drugs      = nan(1, nStates);
mean_APP_baseline  = nan(1, nStates);
mean_APP_drugs     = nan(1, nStates);

sd_WT_baseline   = nan(1, nStates);
sd_WT_drugs      = nan(1, nStates);
sd_APP_baseline  = nan(1, nStates);
sd_APP_drugs     = nan(1, nStates);

% APP-only test results
p_APP        = nan(1, nStates);
test_type_APP = strings(1, nStates);  % "paired" or "unpaired"
d_APP        = nan(1, nStates);       % Cohen's d (unpaired) / approx effect

for i = 1:nStates
    st_i = state_order(i);
    Ti   = Tmouse(Tmouse.state == st_i, :);
    if isempty(Ti), continue; end

    % ---------- per cell stats ----------
    for j = 1:nConds
        cond_j = cond_order(j);
        Tij    = Ti(Ti.condition == cond_j, :);

        for k = 1:nGen
            gen_k = geno_order(k);
            y     = Tij.mean_bout_dur_s(Tij.geno == gen_k);
            y     = y(~isnan(y));

            nVals(i,j,k) = numel(y);
            if ~isempty(y)
                meanVals(i,j,k) = mean(y);
                sdVals(i,j,k)   = std(y);
            end
        end
    end

    % max y for plotting (use all mice)
    max_y_state(i) = max(Ti.mean_bout_dur_s, [], 'omitnan');

    % store per-cell stats
    idxB  = 1; idxD = 2; idxWT = 1; idxAPP = 2;

    n_WT_baseline(i)   = nVals(i,idxB,idxWT);
    n_WT_drugs(i)      = nVals(i,idxD,idxWT);
    n_APP_baseline(i)  = nVals(i,idxB,idxAPP);
    n_APP_drugs(i)     = nVals(i,idxD,idxAPP);

    mean_WT_baseline(i)   = meanVals(i,idxB,idxWT);
    mean_WT_drugs(i)      = meanVals(i,idxD,idxWT);
    mean_APP_baseline(i)  = meanVals(i,idxB,idxAPP);
    mean_APP_drugs(i)     = meanVals(i,idxD,idxAPP);

    sd_WT_baseline(i)   = sdVals(i,idxB,idxWT);
    sd_WT_drugs(i)      = sdVals(i,idxD,idxWT);
    sd_APP_baseline(i)  = sdVals(i,idxB,idxAPP);
    sd_APP_drugs(i)     = sdVals(i,idxD,idxAPP);

    % ---------- APP-only stats: baseline vs drugs ----------
    Ti_APP = Ti(Ti.geno=="APP", :);
    if isempty(Ti_APP), continue; end

    % Try paired t-test first (same mouse in both conditions)
    try
        Twide = unstack(Ti_APP(:, {'mouse','condition','mean_bout_dur_s'}), ...
                        'mean_bout_dur_s', 'condition');
    catch
        Twide = table();
    end

    hasBase = ismember("baseline", Twide.Properties.VariableNames);
    hasDrug = ismember("drugs",    Twide.Properties.VariableNames);

    if hasBase && hasDrug
        base = Twide.baseline;
        drug = Twide.drugs;
        ok   = ~isnan(base) & ~isnan(drug);
        base = base(ok);
        drug = drug(ok);

        if numel(base) >= minNperGroupForStats
            [~, pval] = ttest(base, drug);
            p_APP(i)         = pval;
            test_type_APP(i) = "paired";

            % approximate effect size (using SD of difference)
            diffBD = drug - base;
            m1 = mean(base);
            m2 = mean(drug);
            sd_diff = std(diffBD);
            d_APP(i) = (m2 - m1) / max(sd_diff, eps);
        end
    end

    % If paired test not possible, fall back to unpaired t-test
    if isnan(p_APP(i))
        base = Ti_APP.mean_bout_dur_s(Ti_APP.condition=="baseline");
        drug = Ti_APP.mean_bout_dur_s(Ti_APP.condition=="drugs");
        base = base(~isnan(base));
        drug = drug(~isnan(drug));
        if numel(base) >= minNperGroupForStats && ...
           numel(drug) >= minNperGroupForStats
            [~, pval] = ttest2(base, drug, 'Vartype','unequal');
            p_APP(i)         = pval;
            test_type_APP(i) = "unpaired";

            m1 = mean(base); m2 = mean(drug);
            s1 = std(base);  s2 = std(drug);
            n1 = numel(base); n2 = numel(drug);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
            d_APP(i) = (m2 - m1) / max(sp, eps);
        end
    end
end

% ---------- 6) Prepare data for plotting ----------
x = 1:nStates;
idxB  = 1; idxD = 2; idxWT = 1; idxAPP = 2;

Ybar = [
    squeeze(meanVals(:,idxB, idxWT)), ...
    squeeze(meanVals(:,idxD, idxWT)), ...
    squeeze(meanVals(:,idxB, idxAPP)), ...
    squeeze(meanVals(:,idxD, idxAPP))];

Ysd = [
    squeeze(sdVals(:,idxB, idxWT)), ...
    squeeze(sdVals(:,idxD, idxWT)), ...
    squeeze(sdVals(:,idxB, idxAPP)), ...
    squeeze(sdVals(:,idxD, idxAPP))];

Yn  = [
    squeeze(nVals(:,idxB, idxWT)), ...
    squeeze(nVals(:,idxD, idxWT)), ...
    squeeze(nVals(:,idxB, idxAPP)), ...
    squeeze(nVals(:,idxD, idxAPP))];

Ysem = Ysd ./ sqrt(max(Yn, 1));

% ---------- 7) Plot ----------
figure('Color','w'); hold on;

nBars      = 4;
groupWidth = 0.8;
barWidth   = groupWidth / nBars;

barHandles = gobjects(1, nBars);

for b = 1:nBars
    x_b = x + (b - (nBars+1)/2) * barWidth;
    barHandles(b) = bar(x_b, Ybar(:,b), barWidth, ...
        'FaceColor', COL_BAR(b,:), 'EdgeColor','none');
end

% Error bars
for b = 1:nBars
    x_b = x + (b - (nBars+1)/2) * barWidth;
    errorbar(x_b, Ybar(:,b), Ysem(:,b), 'k', ...
        'LineStyle','none','LineWidth',1);
end

% Dots + optional IDs
jitterFrac = 0.4;

for i = 1:nStates
    st_i = state_order(i);
    Ti   = Tmouse(Tmouse.state == st_i, :);
    if isempty(Ti), continue; end

    for j = 1:nConds
        cond_j = cond_order(j);
        for k = 1:nGen
            gen_k = geno_order(k);

            if k == 1 && j == 1      % WT baseline
                b = 1;
            elseif k == 1 && j == 2  % WT drugs
                b = 2;
            elseif k == 2 && j == 1  % APP baseline
                b = 3;
            else                     % APP drugs
                b = 4;
            end

            x_center = x(i) + (b - (nBars+1)/2) * barWidth;

            Tij = Ti(Ti.condition == cond_j & Ti.geno == gen_k, :);
            if isempty(Tij), continue; end

            vals   = Tij.mean_bout_dur_s;
            miceID = Tij.mouse;

            x_dots = x_center + (rand(size(vals)) - 0.5) * barWidth * jitterFrac;

            plot(x_dots, vals, '.', ...
                'Color', COL_DOT(b,:), 'MarkerSize', 10);

            if showIDs
                for q = 1:numel(vals)
                    thisID = char(miceID(q));
                    text(x_dots(q), vals(q), thisID, ...
                        'Rotation', 45, ...
                        'HorizontalAlignment','left', ...
                        'VerticalAlignment','bottom', ...
                        'FontSize',8, ...
                        'Color', COL_DOT(b,:));
                end
            end
        end
    end
end

% Y-limits & star positioning
all_vals = Tmouse.mean_bout_dur_s;
global_max_y = max(all_vals, [], 'omitnan');
if isempty(global_max_y) || isnan(global_max_y)
    global_max_y = 1;
end
y_top = global_max_y * 1.3;
ylim([0, y_top]);

y_offset_abs = 0.03 * y_top;  % small offset for text

% Stars ONLY for APP baseline vs APP drugs
for i = 1:nStates
    p_here = p_APP(i);
    if isnan(p_here) || p_here >= 0.05
        continue;
    end

    if p_here < 0.001
        stars = '***';
    elseif p_here < 0.01
        stars = '**';
    else
        stars = '*';
    end

    y_star = max_y_state(i);
    if isnan(y_star)
        y_star = global_max_y * 0.8;
    end
    y_star = y_star + y_offset_abs;

    % APP baseline (b=3) and APP drugs (b=4)
    x3 = x(i) + (3 - (nBars+1)/2) * barWidth;
    x4 = x(i) + (4 - (nBars+1)/2) * barWidth;

    line([x3, x4], [y_star, y_star], ...
         'Color','k','LineWidth',1.0);
    text(mean([x3 x4]), y_star + y_offset_abs, stars, ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',14, ...
         'FontWeight','bold');
end

% Axes & labels
xticks(x);
xticklabels(state_order);
xlabel('State');
ylabel('Mean bout duration (s)');
title('Mean bout duration per state: WT vs APP, baseline vs drugs');
set(gca,'Box','off','FontSize',12);

legend(barHandles, ...
    {'WT baseline','WT drugs','APP baseline','APP drugs'}, ...
    'Location','northoutside','Orientation','horizontal');

% ---------- 8) Stats table ----------
stats_tbl = table( ...
    state_order', ...
    n_WT_baseline',   n_WT_drugs',   n_APP_baseline',   n_APP_drugs', ...
    mean_WT_baseline',mean_WT_drugs',mean_APP_baseline',mean_APP_drugs', ...
    sd_WT_baseline',  sd_WT_drugs',  sd_APP_baseline',  sd_APP_drugs', ...
    p_APP', test_type_APP', d_APP', ...
    'VariableNames', { ...
        'State', ...
        'n_WT_baseline','n_WT_drugs','n_APP_baseline','n_APP_drugs', ...
        'Mean_WT_baseline','Mean_WT_drugs', ...
        'Mean_APP_baseline','Mean_APP_drugs', ...
        'SD_WT_baseline','SD_WT_drugs', ...
        'SD_APP_baseline','SD_APP_drugs', ...
        'p_APP_baseline_vs_drugs','test_type_APP','Cohen_d_APP'});

fprintf('\nAPP baseline vs APP drugs: per-state stats\n');
disp(stats_tbl);

% ---------- 9) Save figure ----------
if showIDs
    id_suffix = '_withIDs';
else
    id_suffix = '_noIDs';
end
fname    = sprintf('boutduration_APPvsWT_baseline_vs_drugs%s.png', id_suffix);
out_file = fullfile(out_dir, fname);
saveas(gcf, out_file);

% ---------- OUT ----------
OUT = struct();
OUT.success         = true;
OUT.out_file        = out_file;
OUT.state_order     = state_order;
OUT.condition_order = cond_order;
OUT.stats           = stats_tbl;

fprintf('✅ Bout duration (baseline vs drugs; APP stats only) plot saved to: %s\n', out_file);
end
