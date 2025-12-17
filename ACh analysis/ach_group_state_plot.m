function GROUP_STATE = ach_group_state_plot(groupTransMat, stateName, t_pre, t_post, dt)
% ACH_GROUP_STATE_PLOT
%   Group-level peri-onset ACh plot for one sleep state
%   (Wake, NREM, or REM) using ACh_periTransitions_all.mat
%
%   GROUP_STATE = ach_group_state_plot('ACh_periTransitions_all.mat','NREM');

% ----------- defaults -----------
if nargin < 1 || isempty(groupTransMat)
    [fname, fpath] = uigetfile('*.mat', 'Select ACh_periTransitions_all.mat');
    if isequal(fname,0), error('No file selected.'); end
    groupTransMat = fullfile(fpath, fname);
end
if nargin < 2 || isempty(stateName), stateName = 'NREM'; end
if nargin < 3 || isempty(t_pre),     t_pre  = 30;  end
if nargin < 4 || isempty(t_post),    t_post = 60;  end
if nargin < 5 || isempty(dt),        dt     = 0.1; end

t_rel_common = -t_pre:dt:t_post;

% ----------- load GROUP_TRANS -----------
S = load(groupTransMat, 'GROUP_TRANS');
GROUP_TRANS = S.GROUP_TRANS;

sessions = GROUP_TRANS.sessions;
OUT_all  = GROUP_TRANS.out;
nSess    = numel(sessions);

genos = cell(1,nSess);
conds = cell(1,nSess);
for k = 1:nSess
    genos{k} = sessions(k).geno;
    conds{k} = sessions(k).cond;
end
genoList = unique(genos);
condList = unique(conds);

% ----------- per-session mean traces -----------
sessMean = nan(nSess, numel(t_rel_common));

for k = 1:nSess
    OUT_T = OUT_all{k};
    if isempty(OUT_T) || ~isfield(OUT_T, stateName), continue; end
    M = OUT_T.(stateName).M;
    if isempty(M), continue; end

    t_sess = OUT_T.t_rel(:).';
    mu_evt = mean(M, 1, 'omitnan');

    if numel(t_sess) ~= numel(mu_evt)
        warning('Session %d: t_rel (%d) ~= trace (%d). Skipping.', ...
            k, numel(t_sess), numel(mu_evt));
        continue;
    end

    mu_interp = interp1(t_sess, mu_evt, t_rel_common, 'linear', NaN);
    sessMean(k,:) = mu_interp;
end

hasData = any(~isnan(sessMean), 2);   % session contributes if ANY non-NaN

% ----------- average per geno x cond -----------
nG = numel(genoList);
nC = numel(condList);

GROUP_STATE = struct();
GROUP_STATE.t_rel      = t_rel_common;
GROUP_STATE.genotypes  = genoList;
GROUP_STATE.conditions = condList;
GROUP_STATE.mean       = cell(nG,nC);
GROUP_STATE.sem        = cell(nG,nC);
GROUP_STATE.nSessions  = zeros(nG,nC);

for g = 1:nG
    for c = 1:nC
        mask = strcmp(genos, genoList{g}) & ...
               strcmp(conds, condList{c}) & ...
               hasData';

        idx = find(mask);
        if isempty(idx)
            continue;
        end

        data = sessMean(idx,:);
        GROUP_STATE.mean{g,c} = mean(data, 1, 'omitnan');
        GROUP_STATE.sem{g,c}  = std(data, [], 1, 'omitnan') ./ sqrt(size(data,1));
        GROUP_STATE.nSessions(g,c) = size(data,1);
    end
end

% ----------- plot with correct legend -----------
figure('Color','w','Position',[100 100 800 400]); hold on;

% One colour per genotype, one line style per condition
baseColors = lines(nG);                  % e.g. APP = colour 1, WT = colour 2
linestyles = {'-','--',':','-.'};        % baseline / ambtemp / drugs / extra

hLines = []; legendEntries = {};

for g = 1:nG
    for c = 1:nC
        mu  = GROUP_STATE.mean{g,c};
        sem = GROUP_STATE.sem{g,c};
        nS  = GROUP_STATE.nSessions(g,c);
        if isempty(mu) || nS == 0, continue; end

        col = baseColors(g,:);
        ls  = linestyles{mod(c-1,numel(linestyles))+1};

        % SEM patch – hidden from legend
        patch([t_rel_common fliplr(t_rel_common)], ...
              [mu - sem, fliplr(mu + sem)], ...
              col, 'FaceAlpha',0.15, 'EdgeColor','none', ...
              'HandleVisibility','off');

        % mean line – this goes in legend
        h = plot(t_rel_common, mu, 'Color', col, ...
                 'LineStyle', ls, 'LineWidth', 2);

        hLines(end+1) = h; %#ok<AGROW>
        legendEntries{end+1} = sprintf('%s, %s (n=%d)', ...
            genoList{g}, condList{c}, nS);
    end
end

plot([0 0], ylim, 'k--','HandleVisibility','off');   % onset

xlabel('Time from onset (s)');
ylabel('\DeltaF/F ACh');
title(sprintf('%s onset: group mean \x00b1 SEM (per session)', stateName));
grid on;

if ~isempty(hLines)
    legend(hLines, legendEntries, 'Location','bestoutside');
end

end
