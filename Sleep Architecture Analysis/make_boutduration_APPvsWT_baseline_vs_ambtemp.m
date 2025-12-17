function OUT = make_boutduration_APPvsWT_baseline_vs_ambtemp(PERHOUR, out_dir, states_to_use, varargin)
% make_boutduration_APPvsWT_baseline_vs_ambtemp
% -------------------------------------------------------------------------
% Comparison of:
%   - Genotype: WT vs APP (all non-WT collapsed to APP)
%   - Condition: baseline vs ambtemp
%
% For each mouse, condition, genotype and state:
%   mean_bout_dur = total_dur_s_across_hours / total_number_of_bouts
%   (number of bouts approximated as sum(bouts_per_h) across hours)
%
% Plot (single figure):
%   - X: states (e.g., WK, MA, NREM, REM)
%   - Bars: 4 per state
%       1) WT baseline
%       2) WT ambtemp
%       3) APP baseline
%       4) APP ambtemp
%     with mean ± SEM
%   - Dots: individual mice (optional IDs)
%
% Stats per state:
%   - n and mean/SD for each of the 4 cells:
%       WT_baseline, WT_ambtemp, APP_baseline, APP_ambtemp
%   - 2-way ANOVA (per state, using per-mouse means):
%       factors: Genotype (WT vs APP), Condition (baseline vs ambtemp)
%       outputs: p_genotype, p_condition, p_interaction
%
% INPUT
%   PERHOUR       : OUT.per_hour table from run_group_sleep_architecture
%   out_dir       : folder to save the figure (default: pwd)
%   states_to_use : string/cell array, e.g. ["WK","MA","NREM","REM"]
%
% NAME–VALUE OPTIONS
%   'showIDs'             : true/false, show mouse IDs next to dots (default: true)
%   'minNperGroupForStats': not critical here, kept for compatibility (default: 3)
%
% OUTPUT
%   OUT.success        : logical
%   OUT.out_file       : PNG filename
%   OUT.state_order    : states actually plotted
%   OUT.condition_order: ["baseline","ambtemp"]
%   OUT.stats          : stats table (one row per state)
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

PH = PERHOUR;

% ---------- 1) Keep only BASELINE and AMBTEMP ----------
if ~ismember('condition', PH.Properties.VariableNames)
    error('PERHOUR must contain a "condition" column.');
end

PH.condition = string(lower(strtrim(PH.condition)));
keep_cond    = PH.condition == "baseline" | PH.condition == "ambtemp";
PH           = PH(keep_cond, :);

if isempty(PH)
    warning('No baseline/ambtemp rows in PERHOUR. Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline/ambtemp data');
    return;
end

% ---------- 2) Check required columns ----------
needed = {'dur_s','bouts_per_h','state','mouse','genotype'};
if ~all(ismember(needed, PH.Properties.VariableNames))
    error('PERHOUR table missing required columns for make_boutduration_APPvsWT_baseline_vs_ambtemp.');
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
    warning('None of the requested states found in PERHOUR data.');
    OUT = struct('success', false, 'msg', 'no states found');
    return;
end

% Fixed condition/genotype order for plotting / stats
cond_order = ["baseline","ambtemp"];
geno_order = ["WT","APP"];

% ---------- 4) Per-mouse mean bout duration per state & condition ----------
[gid, m, gen, cond, st] = findgroups( ...
    PH.mouse, PH.geno_group, PH.condition, PH.state);

tot_dur   = splitapply(@(x) sum(double(x),'omitnan'), PH.dur_s,       gid);
tot_bouts = splitapply(@(x) sum(double(x),'omitnan'), PH.bouts_per_h, gid);

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
% 2: WT ambtemp
% 3: APP baseline
% 4: APP ambtemp
COL_BAR = [
    0.6  0.6  0.6 ;  % WT baseline
    0.3  0.3  0.3 ;  % WT ambtemp (darker)
    0.39 0.58 0.93;  % APP baseline
    0.19 0.36 0.70]; % APP ambtemp (darker)

COL_DOT = [
    0.2  0.2  0.2 ;  % WT baseline
    0.1  0.1  0.1 ;  % WT ambtemp
    0.1  0.2  0.6 ;  % APP baseline
    0.05 0.1  0.4]; % APP ambtemp

% ---------- 5) Compute per-cell stats + ANOVA per state ----------
nStates = numel(state_order);
nConds  = numel(cond_order);
nGen    = numel(geno_order);

% mean/SD/n arrays: (state, condition, genotype)
meanVals = nan(nStates, nConds, nGen);
sdVals   = nan(nStates, nConds, nGen);
nVals    = nan(nStates, nConds, nGen);

% Per-state ANOVA p-values
p_genotype   = nan(1, nStates);
p_condition  = nan(1, nStates);
p_interact   = nan(1, nStates);

% For star placement
max_y_state = nan(1, nStates);

% For output table, we keep per state / cell
n_WT_baseline   = nan(1, nStates);
n_WT_ambtemp    = nan(1, nStates);
n_APP_baseline  = nan(1, nStates);
n_APP_ambtemp   = nan(1, nStates);

mean_WT_baseline   = nan(1, nStates);
mean_WT_ambtemp    = nan(1, nStates);
mean_APP_baseline  = nan(1, nStates);
mean_APP_ambtemp   = nan(1, nStates);

sd_WT_baseline   = nan(1, nStates);
sd_WT_ambtemp    = nan(1, nStates);
sd_APP_baseline  = nan(1, nStates);
sd_APP_ambtemp   = nan(1, nStates);

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

    % Flatten for ANOVA
    y_all = Ti.mean_bout_dur_s;
    valid = ~isnan(y_all);
    Ti    = Ti(valid,:);
    y_all = y_all(valid);

    if numel(y_all) >= 3 && ...
       numel(unique(Ti.geno)) >= 2 && ...
       numel(unique(Ti.condition)) >= 2

        try
            % Two-way ANOVA (Genotype × Condition)
            [p, ~, ~] = anovan(y_all, ...
                {Ti.geno, Ti.condition}, ...
                'model','interaction', ...
                'display','off', ...
                'varnames',{'Genotype','Condition'});
            % p(1): Genotype, p(2): Condition, p(3): Interaction
            p_genotype(i)  = p(1);
            p_condition(i) = p(2);
            p_interact(i)  = p(3);
        catch
            warning('ANOVA failed for state %s.', st_i);
        end
    end

    % max y for plotting
    max_y_state(i) = max(Ti.mean_bout_dur_s, [], 'omitnan');

    % store per-cell stats for output
    % indices: cond_order = ["baseline","ambtemp"]; geno_order = ["WT","APP"];
    idxB  = 1; idxA = 2; idxWT = 1; idxAPP = 2;

    n_WT_baseline(i)   = nVals(i,idxB,idxWT);
    n_WT_ambtemp(i)    = nVals(i,idxA,idxWT);
    n_APP_baseline(i)  = nVals(i,idxB,idxAPP);
    n_APP_ambtemp(i)   = nVals(i,idxA,idxAPP);

    mean_WT_baseline(i)   = meanVals(i,idxB,idxWT);
    mean_WT_ambtemp(i)    = meanVals(i,idxA,idxWT);
    mean_APP_baseline(i)  = meanVals(i,idxB,idxAPP);
    mean_APP_ambtemp(i)   = meanVals(i,idxA,idxAPP);

    sd_WT_baseline(i)   = sdVals(i,idxB,idxWT);
    sd_WT_ambtemp(i)    = sdVals(i,idxA,idxWT);
    sd_APP_baseline(i)  = sdVals(i,idxB,idxAPP);
    sd_APP_ambtemp(i)   = sdVals(i,idxA,idxAPP);
end

% ---------- 6) Prepare data for plotting ----------
x = 1:nStates;

% Map into 4 columns: [WT_baseline, WT_ambtemp, APP_baseline, APP_ambtemp]
idxB  = 1; idxA = 2; idxWT = 1; idxAPP = 2;

Ybar = [
    squeeze(meanVals(:,idxB, idxWT)), ...
    squeeze(meanVals(:,idxA, idxWT)), ...
    squeeze(meanVals(:,idxB, idxAPP)), ...
    squeeze(meanVals(:,idxA, idxAPP))];

Ysd = [
    squeeze(sdVals(:,idxB, idxWT)), ...
    squeeze(sdVals(:,idxA, idxWT)), ...
    squeeze(sdVals(:,idxB, idxAPP)), ...
    squeeze(sdVals(:,idxA, idxAPP))];

Yn  = [
    squeeze(nVals(:,idxB, idxWT)), ...
    squeeze(nVals(:,idxA, idxWT)), ...
    squeeze(nVals(:,idxB, idxAPP)), ...
    squeeze(nVals(:,idxA, idxAPP))];

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

            % determine which bar index this cell corresponds to
            if k == 1 && j == 1      % WT baseline
                b = 1;
            elseif k == 1 && j == 2  % WT ambtemp
                b = 2;
            elseif k == 2 && j == 1  % APP baseline
                b = 3;
            else                     % APP ambtemp
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

% Y-limits & star positioning (use interaction p for stars)
all_vals = Tmouse.mean_bout_dur_s;
global_max_y = max(all_vals, [], 'omitnan');
if isempty(global_max_y) || isnan(global_max_y)
    global_max_y = 1;
end
y_top = global_max_y * 1.3;
ylim([0, y_top]);

y_offset_abs = 0.03 * y_top;  % small offset for text

for i = 1:nStates
    p_here = p_interact(i);   % stars reflect Genotype×Condition interaction
    if isnan(p_here) || p_here >= 0.05
        continue;
    end

    % decide number of stars
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

    % horizontal line over the 4 bars
    x_left  = x(i) - groupWidth/2;
    x_right = x(i) + groupWidth/2;
    line([x_left, x_right], [y_star, y_star], ...
         'Color','k','LineWidth',1.0);
    text(x(i), y_star + y_offset_abs, stars, ...
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
title('Mean bout duration per state: WT vs APP, baseline vs ambtemp');
set(gca,'Box','off','FontSize',12);

legend(barHandles, ...
    {'WT baseline','WT ambtemp','APP baseline','APP ambtemp'}, ...
    'Location','northoutside','Orientation','horizontal');

% ---------- 8) Stats table ----------
stats_tbl = table( ...
    state_order', ...
    n_WT_baseline',   n_WT_ambtemp',   n_APP_baseline',   n_APP_ambtemp', ...
    mean_WT_baseline',mean_WT_ambtemp',mean_APP_baseline',mean_APP_ambtemp', ...
    sd_WT_baseline',  sd_WT_ambtemp',  sd_APP_baseline',  sd_APP_ambtemp', ...
    p_genotype', p_condition', p_interact', ...
    'VariableNames', { ...
        'State', ...
        'n_WT_baseline','n_WT_ambtemp','n_APP_baseline','n_APP_ambtemp', ...
        'Mean_WT_baseline','Mean_WT_ambtemp', ...
        'Mean_APP_baseline','Mean_APP_ambtemp', ...
        'SD_WT_baseline','SD_WT_ambtemp', ...
        'SD_APP_baseline','SD_APP_ambtemp', ...
        'p_Anova_Genotype','p_Anova_Condition','p_Anova_Interaction'});

fprintf('\nMean bout duration: WT vs APP, baseline vs ambtemp, per-state ANOVA stats\n');
disp(stats_tbl);

% ---------- 9) Save figure ----------
if showIDs
    id_suffix = '_withIDs';
else
    id_suffix = '_noIDs';
end
fname    = sprintf('boutduration_APPvsWT_baseline_vs_ambtemp%s.png', id_suffix);
out_file = fullfile(out_dir, fname);
saveas(gcf, out_file);

% ---------- OUT ----------
OUT = struct();
OUT.success        = true;
OUT.out_file       = out_file;
OUT.state_order    = state_order;
OUT.condition_order = cond_order;
OUT.stats          = stats_tbl;

fprintf('✅ Bout duration (baseline vs ambtemp, WT vs APP) plot saved to: %s\n', out_file);
end
