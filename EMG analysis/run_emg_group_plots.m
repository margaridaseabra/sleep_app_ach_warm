function run_emg_group_plots(EMG_GROUP)
% RUN_EMG_GROUP_PLOTS
% -------------------------------------------------------------------------
% Group-level APP vs WT figures for EMG bursts:
%   - Bursts per min (Wake/NREM/REM)
%   - Burst peak amplitude (Wake/NREM/REM)

if isfield(EMG_GROUP,'metrics_tbl')
    M = EMG_GROUP.metrics_tbl;
else
    M = struct2table(EMG_GROUP.metrics);
end

M.cond = string(M.cond);
M.geno = string(M.geno);

conds_all = unique(M.cond,'stable');
genos_all = unique(M.geno,'stable');

% order genotypes: WT then APP
if any(genos_all=="WT")
    genoWT  = "WT";
    genoMut = genos_all(genos_all~="WT");
else
    genoWT  = genos_all(1);
    genoMut = genos_all(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes (WT + APP).');
end
genoMut    = genoMut(1);
geno_order = [genoMut genoWT];

COL_APP = [0.6 0.6 0.6];
COL_WT  = [0.392 0.584 0.929];
color_for = @(g) (strcmpi(g,genoWT)*COL_WT) + (~strcmpi(g,genoWT)*COL_APP);

% ------------- burts per minute -----------------------------------------
figure('Name','EMG bursts per minute','Color','w');
states = {'Wake','NREM','REM'};
for sIdx = 1:numel(states)
    subplot(1,3,sIdx); hold on;
    fieldName = ['bursts_per_min_' states{sIdx}];
    [x_all, y_all] = deal([]);

    for gIdx = 1:numel(geno_order)
        g = geno_order(gIdx);
        col = color_for(g);
        for cIdx = 1:numel(conds_all)
            c = conds_all(cIdx);
            mask = (M.cond==c) & (M.geno==g);
            vals = M.(fieldName)(mask);

            mu = mean(vals,'omitnan');
            se = std(vals,'omitnan') / max(1,sqrt(sum(~isnan(vals))));

            x = cIdx + (gIdx - (numel(geno_order)+1)/2)*0.25;
            bar(x, mu, 0.22, 'FaceColor', col);
            errorbar(x, mu, se, 'k','LineStyle','none','CapSize',8);

            x_all = [x_all; repmat(x,size(vals))]; %#ok<AGROW>
            y_all = [y_all; vals(:)];             %#ok<AGROW>
        end
    end

    scatter(x_all, y_all, 15, [0.2 0.2 0.2], 'filled','MarkerFaceAlpha',0.5);

    xlim([0.5 numel(conds_all)+0.5]);
    set(gca,'XTick',1:numel(conds_all),'XTickLabel',conds_all);
    ylabel('Bursts per min');
    title(['EMG bursts – ' states{sIdx}]);
    box off;
end

% legend
axes(subplot(1,3,1)); hold on;
hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
legend([hAPP hWT], {char(genoMut), char(genoWT)}, ...
       'Location','northoutside','Orientation','horizontal','Box','off');
sgtitle('EMG bursts per minute – APP vs WT across conditions');

% ------------- burst amplitudes ----------------------------------------
figure('Name','EMG burst peak amplitude','Color','w');
for sIdx = 1:numel(states)
    subplot(1,3,sIdx); hold on;
    fieldName = ['peak_amp_' states{sIdx}];
    [x_all, y_all] = deal([]);

    for gIdx = 1:numel(geno_order)
        g = geno_order(gIdx);
        col = color_for(g);
        for cIdx = 1:numel(conds_all)
            c = conds_all(cIdx);
            mask = (M.cond==c) & (M.geno==g);
            vals = M.(fieldName)(mask);

            mu = mean(vals,'omitnan');
            se = std(vals,'omitnan') / max(1,sqrt(sum(~isnan(vals))));

            x = cIdx + (gIdx - (numel(geno_order)+1)/2)*0.25;
            bar(x, mu, 0.22, 'FaceColor', col);
            errorbar(x, mu, se, 'k','LineStyle','none','CapSize',8);

            x_all = [x_all; repmat(x,size(vals))]; %#ok<AGROW>
            y_all = [y_all; vals(:)];             %#ok<AGROW>
        end
    end

    scatter(x_all, y_all, 15, [0.2 0.2 0.2], 'filled','MarkerFaceAlpha',0.5);

    xlim([0.5 numel(conds_all)+0.5]);
    set(gca,'XTick',1:numel(conds_all),'XTickLabel',conds_all);
    ylabel('Peak EMG envelope (arb. units)');
    title(['Burst amplitude – ' states{sIdx}]);
    box off;
end

axes(subplot(1,3,1)); hold on;
hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
legend([hAPP hWT], {char(genoMut), char(genoWT)}, ...
       'Location','northoutside','Orientation','horizontal','Box','off');
sgtitle('EMG burst amplitude – APP vs WT across conditions');

end
