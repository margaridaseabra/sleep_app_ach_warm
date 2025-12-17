function OUT = plot_MA_before_REM_by_type(REM_TRANS, out_dir)
% Baseline only.
% Y: MA seconds in the pre-REM window (e.g. last 120 s)
% X: REM type (Short/Medium/Long)
% Bars: WT vs APP, with mouse-level dots + IDs.

if nargin < 2 || isempty(out_dir), out_dir = pwd; end
if ~isfolder(out_dir), mkdir(out_dir); end

T = REM_TRANS;
is_baseline = lower(strtrim(T.condition))=="baseline";
T = T(is_baseline, :);

if isempty(T)
    warning('No baseline REM bouts found.'); OUT.success=false; return;
end

% per-mouse mean per rem_type
G = groupsummary(T, {'mouse','genotype','rem_type'}, 'mean', 'ma_prev_window_sec');
G.Properties.VariableNames{end} = 'mean_ma_prev_sec';

rem_types_all = ["Short","Medium","Long"];
rem_types = rem_types_all(ismember(rem_types_all, unique(G.rem_type)));
if isempty(rem_types)
    warning('No Short/Medium/Long rem_type labels present.');
    OUT.success = false; return;
end

hasWT  = any(G.genotype=="WT");
hasAPP = any(G.genotype=="APP");
if ~hasWT && ~hasAPP
    warning('No WT or APP genotypes found.'); OUT.success=false; return;
end

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

nTypes = numel(rem_types);
meanWT  = nan(1,nTypes); semWT  = nan(1,nTypes);
meanAPP = nan(1,nTypes); semAPP = nan(1,nTypes);

for i = 1:nTypes
    rt = rem_types(i);
    if hasWT
        mWT = (G.genotype=="WT") & (G.rem_type==rt);
        vWT = G.mean_ma_prev_sec(mWT);
        vWT = vWT(~isnan(vWT));
        if ~isempty(vWT)
            meanWT(i) = mean(vWT);
            semWT(i)  = std(vWT)/sqrt(numel(vWT));
        end
    end
    if hasAPP
        mAPP = (G.genotype=="APP") & (G.rem_type==rt);
        vAPP = G.mean_ma_prev_sec(mAPP);
        vAPP = vAPP(~isnan(vAPP));
        if ~isempty(vAPP)
            meanAPP(i) = mean(vAPP);
            semAPP(i)  = std(vAPP)/sqrt(numel(vAPP));
        end
    end
end

figure('Color','w'); hold on;
x = 1:nTypes; barWidth = 0.4;
jitterFrac = 0.25; y_offset = 2;  % seconds

if hasWT
    bar(x - barWidth/2, meanWT, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
    errorbar(x - barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
end
if hasAPP
    bar(x + barWidth/2, meanAPP, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
    errorbar(x + barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
end

for i = 1:nTypes
    rt = rem_types(i);
    if hasWT
        mWT = (G.genotype=="WT") & (G.rem_type==rt);
        v = G.mean_ma_prev_sec(mWT); mids = G.mouse(mWT);
        xw = (x(i)-barWidth/2) + (rand(size(v))-0.5)*barWidth*jitterFrac;
        plot(xw, v, '.', 'Color',[0.3 0.3 0.3]);
        for j = 1:numel(v)
            text(xw(j), v(j)+y_offset, char(mids(j)), ...
                 'Rotation',45,'FontSize',8, ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color',[0.2 0.2 0.2]);
        end
    end
    if hasAPP
        mAPP = (G.genotype=="APP") & (G.rem_type==rt);
        v = G.mean_ma_prev_sec(mAPP); mids = G.mouse(mAPP);
        xa = (x(i)+barWidth/2) + (rand(size(v))-0.5)*barWidth*jitterFrac;
        plot(xa, v, '.', 'Color',[0.1 0.2 0.6]);
        for j = 1:numel(v)
            text(xa(j), v(j)+y_offset, char(mids(j)), ...
                 'Rotation',45,'FontSize',8, ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color',[0.1 0.2 0.6]);
        end
    end
end

xticks(x); xticklabels(rem_types);
xlabel('REM type');
ylabel('MA in pre-REM window (s)');
title('Baseline: MA before REM (WT vs APP, by REM type)');
if hasWT && hasAPP
    legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
end
set(gca,'Box','off','FontSize',12);

out_file = fullfile(out_dir,'baseline_MA_before_REM_byType_APPvsWT.png');
saveas(gcf,out_file);

OUT = struct('success',true,'file',out_file,'rem_types',rem_types);
fprintf('✅ Saved MA-before-REM plot to %s\n', out_file);
end
