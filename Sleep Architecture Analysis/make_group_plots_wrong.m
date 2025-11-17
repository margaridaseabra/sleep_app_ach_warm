function make_group_plots_wrong(rows_overall, rows_perhr, out_dir)
% MAKE_GROUP_PLOTS
% -------------------------------------------------------------------------
% 1) One big 3x3 grid figure summarising:
%       - total time per state & condition
%       - mean bout duration per state & condition
%       - bouts per hour per state & condition
%    Colours: Wake=purple, NREM=red, REM=green.
%
% 2) Three genotype comparison figures (WT vs APP):
%       - Fig10: total time
%       - Fig11: mean bout duration
%       - Fig12: bouts per hour
%    WT = grey, APP = orange.
%
% INPUTS
%   rows_overall : table from run_group_sleep_architecture.overall
%   rows_perhr   : table from run_group_sleep_architecture.per_hour
%   out_dir      : folder where figures should be saved
% -------------------------------------------------------------------------

if nargin < 3 || isempty(out_dir), out_dir = pwd; end
if ~isfolder(out_dir), mkdir(out_dir); end

% ---- constants ---------------------------------------------------------
stateOrder  = ["WK","NREM","REM"];
condOrder   = unique(rows_overall.condition, 'stable');  % keep file order

COL_WAKE = [0.55 0.10 0.85];   % purple
COL_NREM = [0.85 0.25 0.25];   % red
COL_REM  = [0.20 0.70 0.30];   % green

colMap = containers.Map( ...
    {'WK','NREM','REM'}, ...
    {COL_WAKE, COL_NREM, COL_REM});

COL_WT  = [0.60 0.60 0.60];    % grey
COL_APP = [1.00 0.55 0.00];    % orange

%% ================= 1) 3x3 GRID: STATES x METRICS ======================
fig_grid = figure('Name','Sleep architecture summary (3x3 grid)', ...
                  'Color','w', 'Units','normalized','Position',[0.05 0.05 0.9 0.85]);

metricNames = {'Total time (h)','Mean bout duration (s)','Bouts per hour'};

for si = 1:numel(stateOrder)
    st = stateOrder(si);
    col = colMap(char(st));

    % ---------- metric 1: total time (h) ----------
    [mu, se, condLabels] = compute_total_time(rows_overall, st, condOrder);
    subplot(3,3,si); hold on;
    plot_metric_bar(condLabels, mu, se, col);
    if si == 1
        ylabel(metricNames{1});
    end
    title(sprintf('%s – total time', state_label(st)));

    % ---------- metric 2: mean bout duration (s) ----------
    [mu, se, condLabels] = compute_mean_bout_dur(rows_overall, st, condOrder);
    subplot(3,3,3+si); hold on;
    plot_metric_bar(condLabels, mu, se, col);
    if si == 1
        ylabel(metricNames{2});
    end
    title(sprintf('%s – mean bout duration', state_label(st)));

    % ---------- metric 3: bouts per hour ----------
    [mu, se, condLabels] = compute_bouts_per_hour(rows_perhr, st, condOrder);
    subplot(3,3,6+si); hold on;
    plot_metric_bar(condLabels, mu, se, col);
    if si == 1
        ylabel(metricNames{3});
    end
    title(sprintf('%s – bouts per hour', state_label(st)));
end

sgtitle('Group sleep architecture – all mice / both genotypes','FontWeight','bold');

grid_file = fullfile(out_dir, 'Fig1_9_sleep_architecture_grid.png');
saveas(fig_grid, grid_file);

%% ============== 2) GENOTYPE COMPARISONS (WT vs APP) ====================

% Normalise genotype labels: WT vs APP (everything that is not WT)
rows_overall.geno_group = rows_overall.genotype;
rows_overall.geno_group(rows_overall.genotype ~= "WT") = "APP";
rows_perhr.geno_group = rows_perhr.genotype;
rows_perhr.geno_group(rows_perhr.genotype ~= "WT") = "APP";

% ---- helper for genotype bar plotting ----
plot_genotype_figure(@compute_total_time_by_genotype, ...
    rows_overall, rows_perhr, condOrder, stateOrder, ...
    COL_WT, COL_APP, out_dir, ...
    'Fig10_total_time_genotype.png', ...
    'Total sleep time – WT vs APP (per condition)');

plot_genotype_figure(@compute_mean_bout_dur_by_genotype, ...
    rows_overall, rows_perhr, condOrder, stateOrder, ...
    COL_WT, COL_APP, out_dir, ...
    'Fig11_mean_bout_duration_genotype.png', ...
    'Mean bout duration – WT vs APP (per condition)');

plot_genotype_figure(@compute_bouts_per_hour_by_genotype, ...
    rows_overall, rows_perhr, condOrder, stateOrder, ...
    COL_WT, COL_APP, out_dir, ...
    'Fig12_bouts_per_hour_genotype.png', ...
    'Bouts per hour – WT vs APP (per condition)');

end

% =====================================================================
% ====================== METRIC COMPUTATION ===========================
% =====================================================================

function [mu, se, condLabels] = compute_total_time(overall, state, condOrder)
mask = overall.state == state;
T = overall(mask,:);
if isempty(T)
    mu = []; se = []; condLabels = [];
    return;
end
vals_h = T.total_dur_s / 3600;  % seconds -> hours
[G, condLabels] = findgroups(T.condition);
mu = splitapply(@mean, vals_h, G);
se = splitapply(@(x) std(x)/sqrt(numel(x)), vals_h, G);

% reorder to condOrder
[condLabels, idx] = reorder_by_list(condLabels, condOrder);
mu = mu(idx);
se = se(idx);
end

function [mu, se, condLabels] = compute_mean_bout_dur(overall, state, condOrder)
mask = overall.state == state;
T = overall(mask,:);
if isempty(T)
    mu = []; se = []; condLabels = [];
    return;
end
vals_s = T.mean_bout_dur_s;
[G, condLabels] = findgroups(T.condition);
mu = splitapply(@mean, vals_s, G);
se = splitapply(@(x) std(x)/sqrt(numel(x)), vals_s, G);

[condLabels, idx] = reorder_by_list(condLabels, condOrder);
mu = mu(idx);
se = se(idx);
end

function [mu, se, condLabels] = compute_bouts_per_hour(perhr, state, condOrder)
mask = perhr.state == state;
T = perhr(mask,:);
if isempty(T)
    mu = []; se = []; condLabels = [];
    return;
end

% first compute mean bouts_per_h for each file (across hours)
[Gfile, fileNames] = findgroups(T.file);
bph_file = splitapply(@mean, T.bouts_per_h, Gfile);
cond_file = splitapply(@(x) x(1), T.condition, Gfile);  % one condition per file

% now summarise per condition
[Gc, condLabels] = findgroups(cond_file);
mu = splitapply(@mean, bph_file, Gc);
se = splitapply(@(x) std(x)/sqrt(numel(x)), bph_file, Gc);

[condLabels, idx] = reorder_by_list(condLabels, condOrder);
mu = mu(idx);
se = se(idx);
end

% ---------- genotype versions ----------

function [muWT, seWT, muAPP, seAPP, condLabels] = ...
    compute_total_time_by_genotype(overall, perhr, state, condOrder) %#ok<INUSD>
mask = overall.state == state;
T = overall(mask,:);
if isempty(T)
    muWT = []; seWT = []; muAPP = []; seAPP = []; condLabels = [];
    return;
end

vals_h = T.total_dur_s/3600;

[~, condLabels, genoLabs] = findgroups(T.condition, T.geno_group);

% helper: function to pull stats for one geno
[muWT, seWT, condLabels]  = stats_for_geno(vals_h, T.condition, T.geno_group, "WT",  condOrder);
[muAPP,seAPP,condLabels]  = stats_for_geno(vals_h, T.condition, T.geno_group, "APP", condOrder);
end

function [muWT, seWT, muAPP, seAPP, condLabels] = ...
    compute_mean_bout_dur_by_genotype(overall, perhr, state, condOrder) %#ok<INUSD>
mask = overall.state == state;
T = overall(mask,:);
if isempty(T)
    muWT = []; seWT = []; muAPP = []; seAPP = []; condLabels = [];
    return;
end

vals = T.mean_bout_dur_s;
[muWT, seWT, condLabels]  = stats_for_geno(vals, T.condition, T.geno_group, "WT",  condOrder);
[muAPP,seAPP,condLabels]  = stats_for_geno(vals, T.condition, T.geno_group, "APP", condOrder);
end

function [muWT, seWT, muAPP, seAPP, condLabels] = ...
    compute_bouts_per_hour_by_genotype(overall, perhr, state, condOrder)
mask = perhr.state == state;
T = perhr(mask,:);
if isempty(T)
    muWT = []; seWT = []; muAPP = []; seAPP = []; condLabels = [];
    return;
end

% per file first
[Gfile, ~] = findgroups(T.file);
bph_file  = splitapply(@mean, T.bouts_per_h, Gfile);
cond_file = splitapply(@(x)x(1), T.condition,  Gfile);
geno_file = splitapply(@(x)x(1), T.geno_group, Gfile);

[muWT, seWT, condLabels]  = stats_for_geno(bph_file, cond_file, geno_file, "WT",  condOrder);
[muAPP,seAPP,condLabels]  = stats_for_geno(bph_file, cond_file, geno_file, "APP", condOrder);
end

function [mu, se, condLabelsOut] = stats_for_geno(vals, cond, geno, targetGeno, condOrder)
mask = (geno == targetGeno);
if ~any(mask)
    mu = []; se = []; condLabelsOut = condOrder;
    mu = nan(size(condOrder)); se = nan(size(condOrder));
    return;
end
valsG  = vals(mask);
condG  = cond(mask);
[Gc, condLabels] = findgroups(condG);
mu = splitapply(@mean, valsG, Gc);
se = splitapply(@(x) std(x)/sqrt(numel(x)), valsG, Gc);

[condLabelsOut, idx] = reorder_by_list(condLabels, condOrder);
mu2 = nan(numel(condOrder),1);
se2 = nan(numel(condOrder),1);
% put into correct positions
[~,loc] = ismember(condLabels, condLabelsOut);
for i = 1:numel(condLabels)
    mu2(loc(i)) = mu(i);
    se2(loc(i)) = se(i);
end
mu = mu2;
se = se2;
end

% =====================================================================
% ======================== PLOTTING HELPERS ===========================
% =====================================================================

function plot_metric_bar(condLabels, mu, se, col)
if isempty(mu)
    text(0.5,0.5,'no data','HorizontalAlignment','center');
    axis off;
    return;
end
x = 1:numel(mu);
b = bar(x, mu, 'FaceColor', col, 'EdgeColor','none'); %#ok<NASGU>
hold on;
errorbar(x, mu, se, 'k','LineStyle','none','LineWidth',1,'CapSize',8);
set(gca,'XTick',x,'XTickLabel',cellstr(condLabels));
xtickangle(30);
box off;
end

function [labelsOut, idx] = reorder_by_list(labelsIn, desiredOrder)
% labelsIn, desiredOrder are string arrays
[~, idx] = ismember(labelsIn, desiredOrder);
[~, sortIdx] = sort(idx);
labelsOut = labelsIn(sortIdx);
idx = sortIdx;
end

function txt = state_label(st)
switch char(st)
    case 'WK',   txt = 'Wake';
    case 'NREM', txt = 'NREM';
    case 'REM',  txt = 'REM';
    otherwise,   txt = char(st);
end
end

function plot_genotype_figure(metricFun, rows_overall, rows_perhr, condOrder, stateOrder, ...
                              COL_WT, COL_APP, out_dir, filename, mainTitle)

% Create a 1x3 figure (Wake, NREM, REM) comparing WT vs APP.

fig = figure('Name',mainTitle,'Color','w', ...
             'Units','normalized','Position',[0.1 0.1 0.8 0.35]);

for si = 1:numel(stateOrder)
    st = stateOrder(si);

    [muWT, seWT, muAPP, seAPP, condLabels] = ...
        metricFun(rows_overall, rows_perhr, st, condOrder);

    subplot(1,3,si); hold on;
    if isempty(muWT) && isempty(muAPP)
        text(0.5,0.5,'no data','HorizontalAlignment','center');
        axis off; continue;
    end

    x = 1:numel(condLabels);
    width = 0.35;

    % WT left bar
    bar(x - width/2, muWT, width, 'FaceColor', COL_WT,  'EdgeColor','none');
    % APP right bar
    bar(x + width/2, muAPP, width, 'FaceColor', COL_APP,'EdgeColor','none');

    % errorbars
    errorbar(x - width/2, muWT, seWT, 'k','LineStyle','none','LineWidth',1,'CapSize',8);
    errorbar(x + width/2, muAPP, seAPP, 'k','LineStyle','none','LineWidth',1,'CapSize',8);

    set(gca,'XTick',x,'XTickLabel',cellstr(condLabels));
    xtickangle(30);
    title(state_label(st));
    box off;
    if si == 1
        legend({'WT','APP'},'Location','best');
    end
end

sgtitle(mainTitle,'FontWeight','bold');

saveas(fig, fullfile(out_dir, filename));
end
