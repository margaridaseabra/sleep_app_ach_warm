function ach_plv_heatmap_per_condition(PLV_list, cond_labels, title_str)
% ACH_PLV_HEATMAP_PER_CONDITION
%   For a given transition (e.g. NREM→REM theta–ACh),
%   show PLV event heatmaps + mean±SEM per condition.
%
% INPUTS
%   PLV_list    : 1×N struct array, each with fields:
%                 .t_rel (1×T), .plv_mat (nEvents×T)
%   cond_labels : 1×N cellstr, condition per element of PLV_list
%   title_str   : title string

if ~iscell(cond_labels)
    cond_labels = cellstr(cond_labels);
end

[t_unique,~,~] = unique(cond_labels,'stable');
condNames = t_unique;
nCond = numel(condNames);

% Assume same time axis for all
t_rel = PLV_list(1).t_rel(:)';

figure('Color','w','Position',[100 100 1200 400*nCond]);
tiledlayout(nCond,2,'TileSpacing','compact','Padding','compact');

cols = lines(nCond);

for c = 1:nCond
    this_cond = condNames{c};
    idx = find(strcmp(cond_labels,this_cond));

    % collect all events for this condition
    allMat = [];
    for k = idx
        M = PLV_list(k).plv_mat;    % nEv × nTime
        if isempty(M), continue; end
        allMat = [allMat; M]; %#ok<AGROW>
    end
    if isempty(allMat)
        continue;
    end

    % ---- left: heatmap of all events ----
    ax1 = nexttile((c-1)*2 + 1);
    imagesc(ax1, t_rel, 1:size(allMat,1), allMat);
    axis(ax1,'xy');
    colormap(ax1, turbo);
    caxis(ax1,[0 1]); % PLV range
    ylabel(ax1, sprintf('%s events', this_cond));
    if c == nCond
        xlabel(ax1,'Time from transition (s)');
    end
    hold(ax1,'on');
    plot(ax1,[0 0],[0 size(allMat,1)+1],'k--','LineWidth',1);

    % ---- right: mean ± SEM ----
    mu  = mean(allMat,1,'omitnan');
    se  = std(allMat,[],1,'omitnan')/sqrt(size(allMat,1));

    ax2 = nexttile((c-1)*2 + 2);
    hold(ax2,'on');
    patch(ax2,[t_rel fliplr(t_rel)], [mu-se fliplr(mu+se)], ...
          cols(c,:), 'FaceAlpha',0.15,'EdgeColor','none');
    plot(ax2, t_rel, mu, 'Color',cols(c,:),'LineWidth',2);
    plot(ax2, [0 0], ylim(ax2),'k--');
    if c == nCond
        xlabel(ax2,'Time from transition (s)');
    end
    ylabel(ax2,'PLV');
    title(ax2, sprintf('%s (nEvents = %d)', ...
           this_cond, size(allMat,1)));
    grid(ax2,'on');
end

sgtitle(title_str,'FontWeight','bold');
end
