function run_ach_group_plots(GROUP)
% run_ach_group_plots
% -------------------------------------------------------------------------
% Make group-level figures comparing APP vs WT across conditions for:
%   - NREM PSD metrics (power, peak Hz, peak amp)
%   - Wake-onset peaks & slopes
%   - State-wise slopes (Wake / NREM / REM)
%
% INPUT
%   GROUP : struct returned by run_ach_batch_auto, containing:
%           GROUP.metrics (struct array) or GROUP.metrics_tbl (table)

% --------- get metrics as a table ---------------------------------------
if isfield(GROUP,'metrics_tbl')
    M = GROUP.metrics_tbl;
else
    M = struct2table(GROUP.metrics);
end

M.cond = string(M.cond);
M.geno = string(M.geno);

conds_all = unique(M.cond,'stable');   % baseline / ambtemp / drugs / ...
genos_all = unique(M.geno,'stable');

% Order genotypes as WT then APP-like
if any(genos_all=="WT")
    genoWT  = "WT";
    genoMut = genos_all(genos_all ~= "WT");
else
    genoWT  = genos_all(1);
    genoMut = genos_all(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes (WT + APP).');
end
genoMut = genoMut(1);          % assume a single mutant group (APP/PS1)
geno_order = [genoMut genoWT]; % for bars: APP first, WT second

% Colours (match paper style)
COL_APP = [1.0 0.5 0.0];       % orange
COL_WT  = [0.4 0.4 0.4];       % grey

% simple function to pick colour by genotype name
color_for = @(g) (strcmpi(g,genoWT) * COL_WT) + ...
                 (~strcmpi(g,genoWT) * COL_APP);

%% ===================== 1) NREM PSD metrics ============================
fig_psd = figure('Name','Group NREM ACh PSD metrics','Color','w');
ax_psd  = gobjects(1,3);

metrics = {'NREM_power','NREM_peakHz','NREM_peakAmp'};
titles  = {'NREM ACh power','NREM ACh peak freq','NREM ACh peak amp'};
ylabels = {'Power (arb. units)','Frequency (Hz)','Amplitude (arb. units)'};

for mIdx = 1:numel(metrics)
    ax_psd(mIdx) = subplot(1,3,mIdx); hold on;
    metricName = metrics{mIdx};

    [x_all, y_all] = deal([]);

    for gIdx = 1:numel(geno_order)
        g = geno_order(gIdx);
        for cIdx = 1:numel(conds_all)
            c = conds_all(cIdx);

            mask = (M.cond==c) & (M.geno==g);
            vals = M.(metricName)(mask);

            mu = mean(vals,'omitnan');
            se = std(vals,'omitnan') / max(1,sqrt(sum(~isnan(vals))));

            % cluster position
            x = cIdx + (gIdx - (numel(geno_order)+1)/2)*0.25;
            bar(x, mu, 0.22, 'FaceColor', color_for(g));
            errorbar(x, mu, se, 'k','LineStyle','none','CapSize',8);

            x_all = [x_all; repmat(x,size(vals))]; %#ok<AGROW>
            y_all = [y_all; vals(:)];             %#ok<AGROW>
        end
    end

    scatter(x_all, y_all, 15, [0.2 0.2 0.2], 'filled','MarkerFaceAlpha',0.5);

    xlim([0.5 numel(conds_all)+0.5]);
    set(gca,'XTick',1:numel(conds_all),'XTickLabel',conds_all);
    ylabel(ylabels{mIdx});
    title(titles{mIdx});
    box off;
end

% clean legend: dummy bars with APP orange + WT grey
axes(ax_psd(1)); hold on;
hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
legend([hAPP hWT], {char(genoMut), char(genoWT)}, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off');

sgtitle('NREM ACh PSD metrics – APP vs WT across conditions');

%% ===================== 2) Wake-onset peaks & slopes ===================
fig_wake = figure('Name','Wake-onset ACh peaks & slopes','Color','w');
ax_wake  = gobjects(1,2);

% Panel 1: Wake-onset peak dF/F
ax_wake(1) = subplot(1,2,1); hold on;
plot_group_bars(M, conds_all, geno_order, genoWT, ...
                COL_APP, COL_WT, ...
                'WakeOn_peak_mean', ...
                'Wake-onset peak dF/F', ...
                'Peak dF/F');

% Panel 2: Wake-onset slope
ax_wake(2) = subplot(1,2,2); hold on;
plot_group_bars(M, conds_all, geno_order, genoWT, ...
                COL_APP, COL_WT, ...
                'WakeOn_slope_mean', ...
                'Wake-onset ACh slope', ...
                'Slope (dF/F per s)');

% clean legend on first subplot
axes(ax_wake(1)); hold on;
hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
legend([hAPP hWT], {char(genoMut), char(genoWT)}, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off');

sgtitle('ACh at wake onset – APP vs WT across conditions');

%% ===================== 3) State-wise ACh slopes =======================
fig_state = figure('Name','State-wise ACh slopes','Color','w');
ax_state  = gobjects(1,3);

state_fields = {'slope_Wake','slope_NREM','slope_REM'};
state_labels = {'Wake','NREM','REM'};

for sIdx = 1:numel(state_fields)
    ax_state(sIdx) = subplot(1,3,sIdx); hold on;
    plot_group_bars(M, conds_all, geno_order, genoWT, ...
                    COL_APP, COL_WT, ...
                    state_fields{sIdx}, ...
                    ['ACh slope – ' state_labels{sIdx}], ...
                    'Slope (dF/F per s)');
end

% legend on first subplot
axes(ax_state(1)); hold on;
hAPP = bar(nan,nan,0.22,'FaceColor',COL_APP);
hWT  = bar(nan,nan,0.22,'FaceColor',COL_WT);
legend([hAPP hWT], {char(genoMut), char(genoWT)}, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off');

sgtitle('State-wise ACh slopes – APP vs WT across conditions');

end

% ======================================================================
% ----------------------- Helper: plot_group_bars -----------------------
% ======================================================================
function plot_group_bars(M, conds, geno_order, genoWT, COL_APP, COL_WT, ...
                         fieldName, titleStr, yLabelStr)
% Generic grouped bar + scatter for one metric field in M

color_for = @(g) (strcmpi(g,genoWT) * COL_WT) + ...
                 (~strcmpi(g,genoWT) * COL_APP);

x_all = [];
y_all = [];

for gIdx = 1:numel(geno_order)
    g = geno_order(gIdx);
    col = color_for(g);

    for cIdx = 1:numel(conds)
        c = conds(cIdx);
        mask = (M.cond==c) & (M.geno==g);
        vals = M.(fieldName)(mask);

        mu = mean(vals,'omitnan');
        se = std(vals,'omitnan') / max(1,sqrt(sum(~isnan(vals))));

        % cluster position: APP left, WT right within each condition
        x = cIdx + (gIdx - (numel(geno_order)+1)/2)*0.25;
        bar(x, mu, 0.22, 'FaceColor', col);
        errorbar(x, mu, se, 'k','LineStyle','none','CapSize',8);

        x_all = [x_all; repmat(x,size(vals))]; %#ok<AGROW>
        y_all = [y_all; vals(:)];             %#ok<AGROW>
    end
end

scatter(x_all, y_all, 15, [0.2 0.2 0.2], 'filled','MarkerFaceAlpha',0.5);

xlim([0.5 numel(conds)+0.5]);
set(gca,'XTick',1:numel(conds),'XTickLabel',conds);
ylabel(yLabelStr);
title(titleStr);
box off;
end
