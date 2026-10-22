%% RUN_WARMING_MOUSE_DIAGNOSTICS
%
% Individual-mouse diagnostics for the ambient-warming REM result.
%
% PURPOSE
% -------
% Identify which mice:
%   1) show / do not show a pre -> warming increase in REM%;
%   2) show / do not show a pre -> warming increase in REM bout duration;
%   3) show / do not show a positive difference-in-differences (DiD);
%   4) contribute most strongly to the APP-vs-WT response contrast;
%   5) look unusual because of baseline values rather than warming response.
%
% IMPORTANT
% ---------
% This script is for QC / sensitivity / interpretation.
% It DOES NOT automatically exclude animals.
% Outcome direction alone is NOT a valid exclusion criterion.
%
% Run after:
%   run('pipelines/run_warming_complete.m')
%
% Outputs:
%   outputs/sleep_architecture/warming_complete/mouse_diagnostics/
%
% -------------------------------------------------------------------------

clear;
clc;

pipelineDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(pipelineDir);

addpath(fullfile(repoDir,'config'));
P = paths_config();

inDir = fullfile(P.sleep_output,'warming_complete');
outDir = fullfile(inDir,'mouse_diagnostics');
figDir = fullfile(outDir,'figures');

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end

deltaFile = fullfile(inDir,'warming_delta_from_pre.csv');
didFile   = fullfile(inDir,'warming_difference_in_differences.csv');
rawFile   = fullfile(inDir,'warming_raw_metrics_long.csv');
survFile  = fullfile(inDir,'survival','warming_survival_dataset.csv');

assert(isfile(deltaFile),'Missing %s. Run run_warming_complete.m first.',deltaFile);
assert(isfile(didFile),'Missing %s. Run run_warming_complete.m first.',didFile);
assert(isfile(rawFile),'Missing %s. Run run_warming_complete.m first.',rawFile);
assert(isfile(survFile),'Missing %s. Run run_warming_complete.m first.',survFile);

DELTA = readtable(deltaFile,'TextType','string','VariableNamingRule','preserve');
DID   = readtable(didFile,'TextType','string','VariableNamingRule','preserve');
RAW   = readtable(rawFile,'TextType','string','VariableNamingRule','preserve');
SURV  = readtable(survFile,'TextType','string','VariableNamingRule','preserve');
MICE  = readtable(P.mice,'TextType','string','VariableNamingRule','preserve');

mice = unique(MICE.mouse_id);

%% ------------------------------------------------------------------------
% Build one-row-per-mouse REM diagnostic table
% -------------------------------------------------------------------------

D = MICE(:,["mouse_id","genotype","sex","cohort"]);
D.cohort_display = D.cohort;
D.cohort_display(D.cohort_display=="thesis") = "original";

numericVars = [ ...
    "rem_pct_warming_pre", ...
    "rem_pct_warming_warm", ...
    "rem_pct_delta", ...
    "rem_pct_baseline_pre", ...
    "rem_pct_baseline_warmZT", ...
    "rem_pct_baseline_change", ...
    "rem_pct_DiD", ...
    "rem_dur_warming_pre_s", ...
    "rem_dur_warming_warm_s", ...
    "rem_dur_delta_s", ...
    "rem_dur_baseline_pre_s", ...
    "rem_dur_baseline_warmZT_s", ...
    "rem_dur_baseline_change_s", ...
    "rem_dur_DiD_s", ...
    "rem_rate_warming_pre_per_h", ...
    "rem_rate_warming_warm_per_h", ...
    "rem_rate_delta_per_h", ...
    "rem_rate_DiD_per_h", ...
    "rem_surv_baseline_median_s", ...
    "rem_surv_warming_median_s", ...
    "rem_surv_change_s", ...
    "rem_surv_baseline_n_bouts", ...
    "rem_surv_warming_n_bouts" ...
    ];

for v = numericVars
    D.(v) = nan(height(D),1);
end

for i = 1:height(D)

    m = D.mouse_id(i);

    % REM percentage
    q = DELTA( ...
        DELTA.mouse_id==m & ...
        DELTA.domain=="state_percent" & ...
        DELTA.state=="REM" & ...
        DELTA.metric=="percent_time",:);

    D.rem_pct_warming_pre(i)  = pick(q,"ZT3_6","value");
    D.rem_pct_warming_warm(i) = pick(q,"ZT6_9","value");
    D.rem_pct_delta(i)        = pick(q,"ZT6_9","delta_from_pre");

    qdid = DID( ...
        DID.mouse_id==m & ...
        DID.domain=="state_percent" & ...
        DID.state=="REM" & ...
        DID.metric=="percent_time" & ...
        DID.time_bin=="ZT6_9",:);

    if height(qdid)==1
        D.rem_pct_baseline_pre(i) = double(qdid.baseline_pre);
        D.rem_pct_baseline_warmZT(i) = double(qdid.baseline_target);
        D.rem_pct_baseline_change(i) = double(qdid.baseline_change);
        D.rem_pct_DiD(i) = double(qdid.difference_in_differences);
    end

    % REM mean bout duration
    q = DELTA( ...
        DELTA.mouse_id==m & ...
        DELTA.domain=="bout_metric" & ...
        DELTA.state=="REM" & ...
        DELTA.metric=="mean_bout_duration_s",:);

    D.rem_dur_warming_pre_s(i)  = pick(q,"ZT3_6","value");
    D.rem_dur_warming_warm_s(i) = pick(q,"ZT6_9","value");
    D.rem_dur_delta_s(i)        = pick(q,"ZT6_9","delta_from_pre");

    qdid = DID( ...
        DID.mouse_id==m & ...
        DID.domain=="bout_metric" & ...
        DID.state=="REM" & ...
        DID.metric=="mean_bout_duration_s" & ...
        DID.time_bin=="ZT6_9",:);

    if height(qdid)==1
        D.rem_dur_baseline_pre_s(i) = double(qdid.baseline_pre);
        D.rem_dur_baseline_warmZT_s(i) = double(qdid.baseline_target);
        D.rem_dur_baseline_change_s(i) = double(qdid.baseline_change);
        D.rem_dur_DiD_s(i) = double(qdid.difference_in_differences);
    end

    % REM bout rate
    q = DELTA( ...
        DELTA.mouse_id==m & ...
        DELTA.domain=="bout_metric" & ...
        DELTA.state=="REM" & ...
        DELTA.metric=="bouts_per_h",:);

    D.rem_rate_warming_pre_per_h(i)  = pick(q,"ZT3_6","value");
    D.rem_rate_warming_warm_per_h(i) = pick(q,"ZT6_9","value");
    D.rem_rate_delta_per_h(i)        = pick(q,"ZT6_9","delta_from_pre");

    qdid = DID( ...
        DID.mouse_id==m & ...
        DID.domain=="bout_metric" & ...
        DID.state=="REM" & ...
        DID.metric=="bouts_per_h" & ...
        DID.time_bin=="ZT6_9",:);

    if height(qdid)==1
        D.rem_rate_DiD_per_h(i) = double(qdid.difference_in_differences);
    end

    % Per-mouse REM survival medians, ZT6_9, baseline vs warming
    for condition = ["baseline","warming"]

        qs = SURV( ...
            SURV.mouse_id==m & ...
            SURV.condition==condition & ...
            SURV.time_bin=="ZT6_9" & ...
            SURV.state=="REM",:);

        if isempty(qs)
            continue;
        end

        med = km_median(double(qs.time_s),logical(qs.event));

        if condition=="baseline"
            D.rem_surv_baseline_median_s(i) = med;
            D.rem_surv_baseline_n_bouts(i) = height(qs);
        else
            D.rem_surv_warming_median_s(i) = med;
            D.rem_surv_warming_n_bouts(i) = height(qs);
        end
    end

    D.rem_surv_change_s(i) = ...
        D.rem_surv_warming_median_s(i)-D.rem_surv_baseline_median_s(i);
end

%% Direction flags (descriptive only)

D.rem_pct_increased = D.rem_pct_delta > 0;
D.rem_pct_DiD_positive = D.rem_pct_DiD > 0;

D.rem_dur_increased = D.rem_dur_delta_s > 0;
D.rem_dur_DiD_positive = D.rem_dur_DiD_s > 0;

D.rem_survival_median_increased = D.rem_surv_change_s > 0;

%% ------------------------------------------------------------------------
% Leave-one-mouse-out influence on the APP-vs-WT rescue contrast
%
% Positive rescue contrast:
%       mean(APP response) - mean(WT response) > 0
%
% influence_on_rescue =
%       full contrast - contrast without that mouse
%
% > 0 : mouse supports / increases the APP-selective rescue contrast
% < 0 : mouse suppresses / decreases the APP-selective rescue contrast
% -------------------------------------------------------------------------

metrics = { ...
    'rem_pct_delta', ...
    'rem_pct_DiD', ...
    'rem_dur_delta_s', ...
    'rem_dur_DiD_s', ...
    'rem_surv_change_s' ...
    };

for k = 1:numel(metrics)

    var = metrics{k};

    [fullContrast,withoutMouse,influence] = ...
        loo_influence(D,var);

    D.([var '_full_APPminusWT']) = repmat(fullContrast,height(D),1);
    D.([var '_APPminusWT_without_mouse']) = withoutMouse;
    D.([var '_influence_on_rescue']) = influence;
end

%% Composite descriptive inspection score
%
% This is NOT an exclusion score.
% Negative values mean the mouse tends to pull several REM metrics
% away from an APP-selective rescue pattern.

influenceVars = [ ...
    "rem_pct_delta_influence_on_rescue", ...
    "rem_pct_DiD_influence_on_rescue", ...
    "rem_dur_delta_s_influence_on_rescue", ...
    "rem_dur_DiD_s_influence_on_rescue", ...
    "rem_surv_change_s_influence_on_rescue" ...
    ];

Z = nan(height(D),numel(influenceVars));

for j = 1:numel(influenceVars)
    x = double(D.(influenceVars(j)));
    Z(:,j) = robust_z(x);
end

D.rem_rescue_influence_composite = mean(Z,2,'omitnan');

%% Rank mice

D = sortrows(D,'rem_rescue_influence_composite','ascend');

writetable(D,fullfile(outDir,'REM_mouse_diagnostics.csv'));

APP = D(D.genotype=="APP",:);
WT  = D(D.genotype=="WT",:);

writetable(APP,fullfile(outDir,'REM_APP_mouse_diagnostics.csv'));
writetable(WT, fullfile(outDir,'REM_WT_mouse_diagnostics.csv'));

%% Compact review table

reviewVars = [ ...
    "mouse_id","genotype","cohort_display","sex", ...
    "rem_pct_warming_pre","rem_pct_warming_warm","rem_pct_delta","rem_pct_DiD", ...
    "rem_dur_warming_pre_s","rem_dur_warming_warm_s","rem_dur_delta_s","rem_dur_DiD_s", ...
    "rem_surv_baseline_median_s","rem_surv_warming_median_s","rem_surv_change_s", ...
    "rem_rescue_influence_composite" ...
    ];

REVIEW = D(:,reviewVars);
writetable(REVIEW,fullfile(outDir,'REM_mouse_review_compact.csv'));

%% ------------------------------------------------------------------------
% Figures
% -------------------------------------------------------------------------

% 1. REM percentage pre -> warming, individual mice
plot_pre_warm(D, ...
    'rem_pct_warming_pre','rem_pct_warming_warm', ...
    'REM time (%)','REM percentage: pre -> warming', ...
    fullfile(figDir,'REM_percent_pre_to_warming.png'));

% 2. REM mean bout duration pre -> warming
plot_pre_warm(D, ...
    'rem_dur_warming_pre_s','rem_dur_warming_warm_s', ...
    'REM mean bout duration (s)','REM bout duration: pre -> warming', ...
    fullfile(figDir,'REM_duration_pre_to_warming.png'));

% 3. DiD per mouse: REM%
plot_mouse_response(D,'rem_pct_DiD', ...
    'REM % difference-in-differences', ...
    'REM percent: warming response beyond baseline-day change', ...
    fullfile(figDir,'REM_percent_DiD_by_mouse.png'));

% 4. DiD per mouse: REM duration
plot_mouse_response(D,'rem_dur_DiD_s', ...
    'REM duration DiD (s)', ...
    'REM bout duration: warming response beyond baseline-day change', ...
    fullfile(figDir,'REM_duration_DiD_by_mouse.png'));

% 5. Per-mouse REM survival median change
plot_mouse_response(D,'rem_surv_change_s', ...
    'Change in REM KM median (s)', ...
    'REM survival median: baseline -> warming', ...
    fullfile(figDir,'REM_survival_median_change_by_mouse.png'));

% 6. Influence ranking
fig = figure('Color','w','Position',[80 80 1100 520]);
ax = axes(fig);
hold(ax,'on');

x = 1:height(D);
y = D.rem_rescue_influence_composite;

scatter(ax,x,y,55,'filled');
yline(ax,0,'k:');

xticks(ax,x);
xticklabels(ax,D.mouse_id);
xtickangle(ax,45);
ylabel(ax,'Composite influence on APP-selective REM rescue');
title(ax,['Individual-mouse influence diagnostic (descriptive only)' newline ...
          'negative = pulls against rescue; positive = supports rescue']);
box(ax,'off');

for i = 1:height(D)
    text(ax,x(i),y(i),sprintf(' %s',D.genotype(i)), ...
        'FontSize',8,'VerticalAlignment','bottom');
end

exportgraphics(fig,fullfile(figDir,'REM_rescue_influence_ranking.png'), ...
    'Resolution',300);

%% ------------------------------------------------------------------------
% Console output
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('REM INDIVIDUAL-MOUSE DIAGNOSTICS\n');
fprintf('============================================================\n\n');

fprintf('APP mice ranked from strongest pull AGAINST expected rescue\n');
fprintf('to strongest support FOR expected rescue:\n\n');

disp(APP(:,[ ...
    "mouse_id","cohort_display","sex", ...
    "rem_pct_delta","rem_pct_DiD", ...
    "rem_dur_delta_s","rem_dur_DiD_s", ...
    "rem_surv_change_s", ...
    "rem_rescue_influence_composite" ...
    ]));

fprintf('\nWT mice:\n\n');

disp(WT(:,[ ...
    "mouse_id","cohort_display","sex", ...
    "rem_pct_delta","rem_pct_DiD", ...
    "rem_dur_delta_s","rem_dur_DiD_s", ...
    "rem_surv_change_s", ...
    "rem_rescue_influence_composite" ...
    ]));

fprintf('\nAPP direction counts:\n');

fprintf('  REM%% increased pre->warm:             %d/%d\n', ...
    sum(APP.rem_pct_increased),height(APP));

fprintf('  REM%% positive DiD:                   %d/%d\n', ...
    sum(APP.rem_pct_DiD_positive),height(APP));

fprintf('  REM duration increased pre->warm:     %d/%d\n', ...
    sum(APP.rem_dur_increased),height(APP));

fprintf('  REM duration positive DiD:            %d/%d\n', ...
    sum(APP.rem_dur_DiD_positive),height(APP));

fprintf('  REM survival median increased:        %d/%d\n', ...
    sum(APP.rem_survival_median_increased),height(APP));

fprintf('\nOutputs:\n%s\n\n',outDir);

%% ========================================================================
% Local functions
% =========================================================================

function value = pick(T,timeBin,varName)

value = NaN;

if isempty(T)
    return;
end

x = T.(varName)(T.time_bin==timeBin);

if numel(x)==1
    value = double(x);
end
end


function med = km_median(time,event)

med = NaN;

time = double(time(:));
event = logical(event(:));

valid = isfinite(time) & time>0;
time = time(valid);
event = event(valid);

if isempty(time)
    return;
end

eventTimes = sort(unique(time(event)));
S = 1;

for i = 1:numel(eventTimes)

    t = eventTimes(i);
    nRisk = sum(time>=t);
    nEvent = sum(time==t & event);

    if nRisk>0
        S = S*(1-nEvent/nRisk);
    end

    if S<=0.5
        med = t;
        return;
    end
end
end


function [fullContrast,withoutMouse,influence] = loo_influence(D,varName)

x = double(D.(varName));

valid = isfinite(x);

fullContrast = ...
    mean(x(valid & D.genotype=="APP"),'omitnan') - ...
    mean(x(valid & D.genotype=="WT"),'omitnan');

withoutMouse = nan(height(D),1);
influence = nan(height(D),1);

for i = 1:height(D)

    keep = valid;
    keep(i) = false;

    a = x(keep & D.genotype=="APP");
    w = x(keep & D.genotype=="WT");

    if isempty(a) || isempty(w)
        continue;
    end

    withoutMouse(i) = mean(a,'omitnan')-mean(w,'omitnan');
    influence(i) = fullContrast-withoutMouse(i);
end
end


function z = robust_z(x)

x = double(x(:));
z = nan(size(x));

valid = isfinite(x);

if sum(valid)<2
    return;
end

med = median(x(valid));
madv = median(abs(x(valid)-med));

if madv<=0 || ~isfinite(madv)
    s = std(x(valid));
    if s<=0 || ~isfinite(s)
        z(valid)=0;
    else
        z(valid)=(x(valid)-mean(x(valid)))/s;
    end
else
    z(valid)=0.67448975*(x(valid)-med)/madv;
end
end


function plot_pre_warm(D,preVar,warmVar,yLabelText,titleText,outFile)

fig = figure('Color','w','Position',[80 80 1050 500]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for g = 1:2

    ax = nexttile;
    hold(ax,'on');

    if g==1
        G = D(D.genotype=="WT",:);
        groupName = "WT";
    else
        G = D(D.genotype=="APP",:);
        groupName = "APP";
    end

    for i = 1:height(G)

        y = [double(G.(preVar)(i)), double(G.(warmVar)(i))];

        if any(~isfinite(y))
            continue;
        end

        plot(ax,[1 2],y,'-o','LineWidth',1.2);

        text(ax,2.03,y(2),G.mouse_id(i), ...
            'FontSize',8,'VerticalAlignment','middle');
    end

    xlim(ax,[0.8 2.45]);
    xticks(ax,[1 2]);
    xticklabels(ax,{'Pre','Warming'});
    ylabel(ax,yLabelText);
    title(ax,groupName);
    box(ax,'off');
end

title(tl,titleText);

exportgraphics(fig,outFile,'Resolution',300);
end


function plot_mouse_response(D,varName,yLabelText,titleText,outFile)

fig = figure('Color','w','Position',[80 80 1100 500]);
ax = axes(fig);
hold(ax,'on');

% Keep display ordered by genotype, cohort, mouse ID.
T = sortrows(D,{'genotype','cohort_display','mouse_id'});

x = 1:height(T);
y = double(T.(varName));

scatter(ax,x,y,55,'filled');
yline(ax,0,'k:');

xticks(ax,x);
xticklabels(ax,T.mouse_id);
xtickangle(ax,45);
ylabel(ax,yLabelText);
title(ax,titleText);
box(ax,'off');

for i = 1:height(T)
    text(ax,x(i),y(i),sprintf(' %s/%s',T.genotype(i),T.cohort_display(i)), ...
        'FontSize',7,'Rotation',45,'VerticalAlignment','bottom');
end

exportgraphics(fig,outFile,'Resolution',300);
end
