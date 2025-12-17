function ach_plot_nrem_psd_baseline(GROUP, out_dir, show_mouse_ids)
% Plot baseline NREM ACh PSD (WT vs APP) and bar plots
%
% INPUTS:
%   GROUP          - structure with sessions and metrics
%   out_dir        - directory to save figures (optional)
%   show_mouse_ids - logical, if true show mouse IDs on scatter plots (default: false)

if nargin < 2
    out_dir = [];
end
if nargin < 3 || isempty(show_mouse_ids)
    show_mouse_ids = false;
end

sessions = GROUP.sessions;

% figure out genotypes
geno_vals = unique(string({sessions.geno}),'stable');
if any(geno_vals=="WT")
    genoWT  = "WT";
    genoMut = geno_vals(geno_vals~="WT");
else
    genoWT  = geno_vals(1);
    genoMut = geno_vals(2:end);
end
genoMut = genoMut(1);

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.392 0.584 0.929];

% ---------- collect PSD curves ----------
f_ref    = [];
WT_psd   = [];
APP_psd  = [];

for k = 1:numel(sessions)
    OUT = GROUP.out{k};
    if ~isfield(OUT,'psd') || ~isfield(OUT.psd,'NREM'); continue; end
    P = [];
    if isfield(OUT.psd.NREM,'Pxx')
        P = OUT.psd.NREM.Pxx(:);
    elseif isfield(OUT.psd.NREM,'psd')
        P = OUT.psd.NREM.psd(:);
    elseif isfield(OUT.psd.NREM,'power')
        P = OUT.psd.NREM.power(:);
    end
    if isempty(P) || ~isfield(OUT.psd.NREM,'f'); continue; end

    f_this = OUT.psd.NREM.f(:);

    % reference frequency grid
    if isempty(f_ref)
        f_ref = f_this;
    elseif numel(f_this) ~= numel(f_ref) || any(abs(f_this - f_ref) > 1e-6)
        P = interp1(f_this, P, f_ref, 'linear', 'extrap');
    end

    g = string(sessions(k).geno);
    if g == genoWT
        WT_psd = [WT_psd P]; %#ok<AGROW>
    elseif g == genoMut
        APP_psd = [APP_psd P]; %#ok<AGROW>
    end
end

% ---------- metrics from GROUP.metrics ----------
if isfield(GROUP,'metrics_tbl')
    M = GROUP.metrics_tbl;
else
    M = struct2table(GROUP.metrics);
end
M.geno = string(M.geno);

y_WT_power    = M.NREM_power(M.geno==genoWT);
y_APP_power   = M.NREM_power(M.geno==genoMut);
y_WT_peakHz   = M.NREM_peakHz(M.geno==genoWT);
y_APP_peakHz  = M.NREM_peakHz(M.geno==genoMut);
y_WT_peakAmp  = M.NREM_peakAmp(M.geno==genoWT);
y_APP_peakAmp = M.NREM_peakAmp(M.geno==genoMut);

% ACh cycle frequency (if available)
if ismember('NREM_cycleHz', M.Properties.VariableNames)
    y_WT_cycleHz  = M.NREM_cycleHz(M.geno==genoWT);
    y_APP_cycleHz = M.NREM_cycleHz(M.geno==genoMut);
else
    y_WT_cycleHz  = [];
    y_APP_cycleHz = [];
end

% NREM slope (if available)
if ismember('slope_NREM', M.Properties.VariableNames)
    y_WT_slope  = M.slope_NREM(M.geno==genoWT);
    y_APP_slope = M.slope_NREM(M.geno==genoMut);
else
    y_WT_slope  = [];
    y_APP_slope = [];
end

% ---------- figure ----------
fig = figure('Name','Baseline NREM ACh PSD','Color','w','Position',[100 100 1680 300]);

% PSD curves
subplot(1,6,1); hold on;
if ~isempty(WT_psd)
    m  = mean(WT_psd,2);
    se = std(WT_psd,0,2)/sqrt(size(WT_psd,2));
    fill_between(f_ref, m-se, m+se, COL_WT, 0.2);
    plot(f_ref, m, 'Color', COL_WT, 'LineWidth', 1.8);
end
if ~isempty(APP_psd)
    m  = mean(APP_psd,2);
    se = std(APP_psd,0,2)/sqrt(size(APP_psd,2));
    fill_between(f_ref, m-se, m+se, COL_APP, 0.2);
    plot(f_ref, m, 'Color', COL_APP, 'LineWidth', 1.8);
end
xlabel('Frequency (Hz)');
ylabel('Power (A.U.)');
title('NREM ACh PSD');
box off;
xlim([min(f_ref) max(f_ref)]);
if ~isempty(WT_psd) || ~isempty(APP_psd)
    legend({'WT','APP'},'Box','off','Location','northeast');
end

% Bar plots
subplot(1,6,2);
plot_two_group_bar(y_WT_power, y_APP_power, genoWT, genoMut, ...
                   COL_WT, COL_APP, 'NREM power (A.U.)', ...
                   show_mouse_ids, M(M.geno==genoWT,:), M(M.geno==genoMut,:));

subplot(1,6,3);
plot_two_group_bar(y_WT_peakHz, y_APP_peakHz, genoWT, genoMut, ...
                   COL_WT, COL_APP, 'NREM peak freq (Hz)', ...
                   show_mouse_ids, M(M.geno==genoWT,:), M(M.geno==genoMut,:));

subplot(1,6,4);
plot_two_group_bar(y_WT_peakAmp, y_APP_peakAmp, genoWT, genoMut, ...
                   COL_WT, COL_APP, 'NREM peak amp (A.U.)', ...
                   show_mouse_ids, M(M.geno==genoWT,:), M(M.geno==genoMut,:));

subplot(1,6,5);
if ~isempty(y_WT_cycleHz) || ~isempty(y_APP_cycleHz)
    plot_two_group_bar(y_WT_cycleHz, y_APP_cycleHz, genoWT, genoMut, ...
                       COL_WT, COL_APP, 'ACh cycles (Hz)', ...
                       show_mouse_ids, M(M.geno==genoWT,:), M(M.geno==genoMut,:));
else
    text(0.5, 0.5, 'No cycle data', 'HorizontalAlignment', 'center');
    axis off;
end

subplot(1,6,6);
if ~isempty(y_WT_slope) || ~isempty(y_APP_slope)
    plot_two_group_bar(y_WT_slope, y_APP_slope, genoWT, genoMut, ...
                       COL_WT, COL_APP, 'NREM slope (ΔF/F/s)', ...
                       show_mouse_ids, M(M.geno==genoWT,:), M(M.geno==genoMut,:));
else
    text(0.5, 0.5, 'No slope data', 'HorizontalAlignment', 'center');
    axis off;
end

sgtitle('Baseline NREM ACh – PSD & metrics');

if ~isempty(out_dir)
    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    saveas(fig, fullfile(out_dir, 'Baseline_NREM_ACh_PSD.png'));
end
end

% ----------------------------------------------------------------------
function plot_two_group_bar(yWT, yAPP, genoWT, genoMut, COL_WT, COL_APP, yLabel, show_ids, tbl_WT, tbl_APP)
hold on;

bar(1, mean(yWT,'omitnan'),  0.6, 'FaceColor',COL_WT);
bar(2, mean(yAPP,'omitnan'), 0.6, 'FaceColor',COL_APP);

% scatter points
if ~isempty(yWT)
    xpos = 1 + (rand(size(yWT))-0.5)*0.1;
    scatter(xpos, yWT, 25, [0.2 0.2 0.2], 'filled');
    
    if show_ids && nargin >= 9 && ~isempty(tbl_WT)
        for i = 1:numel(yWT)
            text(xpos(i), yWT(i), sprintf(' %s', tbl_WT.mouse{i}), ...
                 'FontSize', 7, 'HorizontalAlignment', 'left');
        end
    end
end
if ~isempty(yAPP)
    xpos = 2 + (rand(size(yAPP))-0.5)*0.1;
    scatter(xpos, yAPP, 25, [0.2 0.2 0.2], 'filled');
    
    if show_ids && nargin >= 10 && ~isempty(tbl_APP)
        for i = 1:numel(yAPP)
            text(xpos(i), yAPP(i), sprintf(' %s', tbl_APP.mouse{i}), ...
                 'FontSize', 7, 'HorizontalAlignment', 'left');
        end
    end
end

xlim([0.5 2.5]);
set(gca,'XTick',[1 2],'XTickLabel',{char(genoWT), char(genoMut)});
ylabel(yLabel);
box off;

if numel(yWT) > 1 && numel(yAPP) > 1
    [~,p] = ttest2(yWT, yAPP);
    yMax   = max([yWT; yAPP]);
    yMin   = min([yWT; yAPP]);
    yrange = max(yMax - yMin, eps);

    yStar = yMax + 0.15*yrange;
    plot([1 2], [yStar yStar], 'k-', 'LineWidth', 1);
    text(1.5, yStar + 0.03*yrange, p_to_star(p), ...
         'HorizontalAlignment','center', 'FontSize',10);
    text(1.5, yStar + 0.15*yrange, sprintf('p = %.3g', p), ...
         'HorizontalAlignment','center', 'FontSize',8);
end
end

% ----------------------------------------------------------------------
function fill_between(x,y1,y2,col,alphaVal)
x   = x(:);
y1  = y1(:);
y2  = y2(:);
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
     'EdgeColor','none','FaceAlpha',alphaVal);
end

% ----------------------------------------------------------------------
function str = p_to_star(p)
% Convert p-value to star notation for significance
if p < 0.001
    str = '***';
elseif p < 0.01
    str = '**';
elseif p < 0.05
    str = '*';
else
    str = 'ns';
end
end
