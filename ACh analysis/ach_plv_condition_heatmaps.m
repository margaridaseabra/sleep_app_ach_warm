function ach_plv_condition_heatmaps(PLV_list, cond_labels, titleStr)
% ACH_PLV_CONDITION_HEATMAPS
% -------------------------------------------------------------
% Visualise PLV across episodes and conditions as heatmaps.
%
% INPUTS
%   PLV_list    : 1×N struct array from ach_eeg_plv_periTransitions_segment
%                 MUST contain:
%                     .t_rel   [1×T] time axis (s)
%                 and SOME numeric vector field of length T that is the
%                 PLV time course (the function will auto-detect it).
%
%   cond_labels : 1×N cellstr, one label per PLV_list element
%   titleStr    : title string for the figure
%
% OUTPUT
%   (figure only)

if isempty(PLV_list)
    warning('PLV_list is empty, nothing to plot.');
    return;
end

if ~iscell(cond_labels)
    cond_labels = cellstr(cond_labels);
end

% ---------- common time axis ----------
if isfield(PLV_list(1),'t_rel')
    t_ref = PLV_list(1).t_rel(:)';   % row vector
else
    error('PLV_list(1) has no field t_rel.');
end
nT = numel(t_ref);

% ---------- conditions ----------
[condNames, ~, condIdx_all] = unique(cond_labels, 'stable');
nCond = numel(condNames);

% prettier names for labels
niceNames = condNames;
for c = 1:nCond
    name = lower(condNames{c});
    switch name
        case 'baseline'
            niceNames{c} = 'Baseline';
        case 'ambtemp'
            niceNames{c} = 'Amb. temp';
        case 'drugs'
            niceNames{c} = 'Drugs';
        otherwise
            niceNames{c} = regexprep(condNames{c}, '(^.)','${upper($1)}');
    end
end

% colours for scalar bar plot
cols_cond = [ ...
    0.85 0.85 0.90;   % baseline
    0.20 0.45 0.85;   % ambtemp
    0.96 0.80 0.60;   % drugs
    ];
if nCond > size(cols_cond,1)
    cols_cond = lines(nCond);
end

% ---------- figure layout ----------
figure('Color','w','Position',[80 80 340*nCond 480]);
tiledlayout(2, nCond, 'TileSpacing','compact', 'Padding','compact');

for c = 1:nCond
    idx_c = find(condIdx_all == c);
    if isempty(idx_c)
        continue;
    end

    H        = nan(numel(idx_c), nT);   % episodes × time PLV
    scalarEp = nan(numel(idx_c), 1);    % mean PLV per episode

    for kk = 1:numel(idx_c)
        P = PLV_list(idx_c(kk));

        % ---------- auto-detect PLV timecourse field ----------
        v = [];
        flds = fieldnames(P);

        % fields we definitely do NOT want to use as PLV timecourse
        bad = {'t_rel','event_t','event_times','scalar_events', ...
               'scalar_mean','nEvents','dt','fs_eeg','fs_ach'};

        for ii = 1:numel(flds)
            fname = flds{ii};
            if any(strcmp(fname,bad)), continue; end

            val = P.(fname);
            if isnumeric(val) && isvector(val) && numel(val) == nT
                v = val(:)';   % use this as PLV(t)
                break;
            end
        end

        if isempty(v)
            warning('Could not find PLV timecourse field for episode %d, skipping in heatmap.', idx_c(kk));
            continue;
        end

        % store in matrix
        H(kk,:) = v;

        % scalar = mean PLV over peri window (simple summary)
        scalarEp(kk) = mean(v, 'omitnan');
    end

    % ===== 1st row: heatmap =====
    nexttile(c); hold on;
    if all(isnan(H(:)))
        text(0.5,0.5,'No PLV data','Units','normalized',...
             'HorizontalAlignment','center');
        axis off;
    else
        imagesc(t_ref, 1:size(H,1), H);
        axis xy;
        xlabel('Time from transition (s)');
        ylabel('Episode #');
        title(niceNames{c});
        colormap(turbo);
        mx = max(H(:),[],'omitnan');
        if ~isempty(mx) && mx > 0
            caxis([0 mx]);
        end
        colorbar;
    end

    % ===== 2nd row: scalar PLV per episode =====
    nexttile(c + nCond); hold on;
    if all(isnan(scalarEp))
        axis off;
    else
        bar(1:numel(idx_c), scalarEp, ...
            'FaceColor', cols_cond(c,:), 'EdgeColor','none');
        plot(1:numel(idx_c), scalarEp, 'k.', 'MarkerSize', 10);
        xlim([0.5 numel(idx_c)+0.5]);
        xlabel('Episode');
        ylabel('Mean PLV');
        title(sprintf('%s – mean PLV', niceNames{c}));
        box on;
    end
end

sgtitle(titleStr, 'FontWeight','bold','Interpreter','none');
end
