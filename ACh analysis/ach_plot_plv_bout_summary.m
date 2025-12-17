function STATS = ach_plot_plv_bout_summary(PLV_list, cond_labels, title_str)
% ACH_PLOT_PLV_BOUT_SUMMARY
% -------------------------------------------------------------
% Summarise NREM bout-wise PLV across conditions.
%
% INPUTS
%   PLV_list    : struct array, one per condition, with fields:
%                   .plv_all  - [nBouts×1] PLV for each bout
%                   .dur_all  - [nBouts×1] bout duration (s)
%   cond_labels : cellstr, same length as PLV_list (e.g. {'baseline','ambtemp','drugs'})
%   title_str   : figure title (optional)
%
% OUTPUT
%   STATS.plv.p_anova  - p from one-way ANOVA on PLV across conditions
%        .corr.r        - correlation (Pearson) PLV vs duration
%        .corr.p        - p-value for correlation

if nargin < 3 || isempty(title_str)
    title_str = 'Bout-wise PLV summary';
end

nCond = numel(PLV_list);
if numel(cond_labels) ~= nCond
    error('cond_labels must have one label per PLV struct.');
end

% ---- pool data across conditions ----
plv_all   = [];
dur_all   = [];
groups    = [];
condNames = cond_labels;

for c = 1:nCond
    if ~isfield(PLV_list(c),'plv_all') || isempty(PLV_list(c).plv_all)
        continue;
    end
    p = PLV_list(c).plv_all(:);
    d = PLV_list(c).dur_all(:);

    % keep only finite values
    ok = isfinite(p) & isfinite(d);
    p  = p(ok);
    d  = d(ok);

    plv_all = [plv_all; p];
    dur_all = [dur_all; d];
    groups  = [groups;  c*ones(numel(p),1)];
end

STATS = struct();
if isempty(plv_all)
    warning('No PLV data to summarise.');
    STATS.plv.p_anova = NaN;
    STATS.corr.r      = NaN;
    STATS.corr.p      = NaN;
    return;
end

% ---- ANOVA on PLV across conditions ----
if numel(unique(groups)) > 1
    p_anova = anova1(plv_all, groups, 'off');
else
    p_anova = NaN;
end
STATS.plv.p_anova = p_anova;

% ---- correlation PLV vs duration ----
[r_corr, p_corr] = corr(plv_all, dur_all, 'type','Pearson','rows','complete');
STATS.corr.r = r_corr;
STATS.corr.p = p_corr;

% ---- plotting ----
figure('Color','w','Position',[100 100 1000 400]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (1) PLV per condition: bar + scatter
nexttile; hold on;
mean_plv = nan(1,nCond);
for c = 1:nCond
    if ~isfield(PLV_list(c),'plv_all') || isempty(PLV_list(c).plv_all)
        mean_plv(c) = NaN;
    else
        mean_plv(c) = mean(PLV_list(c).plv_all,'omitnan');
    end
end

x = 1:nCond;
bar(x, mean_plv, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none');

for c = 1:nCond
    if ~isfield(PLV_list(c),'plv_all') || isempty(PLV_list(c).plv_all)
        continue;
    end
    p = PLV_list(c).plv_all(:);
    jitter = (rand(size(p))*0.5 - 0.25);
    plot(c + jitter, p, 'k.', 'MarkerSize',8);
end

set(gca,'XTick',x,'XTickLabel',condNames);
ylabel('PLV (delta–ACh)');
title(sprintf('Bout-wise PLV (ANOVA p = %.3f)', p_anova));
box on; grid on;

% (2) PLV vs duration scatter (all bouts)
nexttile; hold on;
cols = lines(nCond);

for c = 1:nCond
    if ~isfield(PLV_list(c),'plv_all') || isempty(PLV_list(c).plv_all)
        continue;
    end
    p = PLV_list(c).plv_all(:);
    d = PLV_list(c).dur_all(:);
    ok = isfinite(p) & isfinite(d);
    scatter(d(ok), p(ok), 25, cols(c,:), 'filled', ...
        'MarkerFaceAlpha',0.7);
end

xlabel('Bout duration (s)');
ylabel('PLV (delta–ACh)');
title(sprintf('PLV vs duration (r = %.2f, p = %.3f)', r_corr, p_corr));
legend(condNames,'Location','best','Box','off');
box on; grid on;

sgtitle(title_str, 'FontWeight','bold');
end
