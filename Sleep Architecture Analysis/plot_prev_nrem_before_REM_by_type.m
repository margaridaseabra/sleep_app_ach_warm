function OUT = plot_prev_nrem_before_REM_by_type(REM_TRANS, out_dir)
% Baseline only.
% Y: mean NREM duration *immediately before REM* (minutes)
% X: REM type (Short / Medium / Long)
% Bars: WT vs APP
% Dots: each mouse (labelled)

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

T = REM_TRANS;

% ---- baseline + only REM that have NREM directly before ----
is_baseline = lower(strtrim(T.condition)) == "baseline";
T = T(is_baseline & T.prev_state=="NREM", :);

if isempty(T)
    warning('No baseline REM bouts with preceding NREM found.');
    OUT.success = false; return;
end

% convert to minutes
T.prev_nrem_min = T.prev_nrem_dur_s / 60;

% ---- per-mouse mean per rem_type ----
G = groupsummary(T, {'mouse','genotype','rem_type'}, 'mean', 'prev_nrem_min');
G.Properties.VariableNames{end} = 'mean_prev_nrem_min';

rem_types_all = ["Short","Medium","Long"];
rem_types = rem_types_all(ismember(rem_types_all, unique(G.rem_type)));

if isempty(rem_types)
    warning('No Short/Medium/Long rem_type labels present.');
    OUT.success = false; return;
end

% ---- WT vs APP stats ----
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
        valsWT = G.mean_prev_nrem_min(mWT);
        valsWT = valsWT(~isnan(valsWT));
        if ~isempty(valsWT)
            meanWT(i) = mean(valsWT);
            semWT(i)  = std(valsWT)/sqrt(numel(valsWT));
        end
    end
    if hasAPP
        mAPP = (G.genotype=="APP") & (G.rem_type==rt);
        valsAPP = G.mean_prev_nrem_min(mAPP);
        valsAPP = valsAPP(~isnan(valsAPP));
        if ~isempty(valsAPP)
            meanAPP(i) = mean(valsAPP);
            semAPP(i)  = std(valsAPP)/sqrt(numel(valsAPP));
        end
    end
end

% ---- plot ----
figure('Color','w'); hold on;
x = 1:nTypes; barWidth = 0.4;
jitterFrac = 0.25; y_offset = 0.05;

if hasWT
    bar(x - barWidth/2, meanWT, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
    errorbar(x - barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
end
if hasAPP
    bar(x + barWidth/2, meanAPP, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
    errorbar(x + barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
end

% dots + mouse IDs
for i = 1:nTypes
    rt = rem_types(i);

    if hasWT
        mWT  = (G.genotype=="WT") & (G.rem_type==rt);
        vals = G.mean_prev_nrem_min(mWT);
        mids = G.mouse(mWT);
        xw   = (x(i)-barWidth/2) + (rand(size(vals))-0.5)*barWidth*jitterFrac;
        plot(xw, vals, '.', 'Color',[0.3 0.3 0.3]);
        for j = 1:numel(vals)
            text(xw(j), vals(j)+y_offset, char(mids(j)), ...
                 'Rotation',45,'FontSize',8, ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color',[0.2 0.2 0.2]);
        end
    end

    if hasAPP
        mAPP = (G.genotype=="APP") & (G.rem_type==rt);
        vals = G.mean_prev_nrem_min(mAPP);
        mids = G.mouse(mAPP);
        xa   = (x(i)+barWidth/2) + (rand(size(vals))-0.5)*barWidth*jitterFrac;
        plot(xa, vals, '.', 'Color',[0.1 0.2 0.6]);
        for j = 1:numel(vals)
            text(xa(j), vals(j)+y_offset, char(mids(j)), ...
                 'Rotation',45,'FontSize',8, ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color',[0.1 0.2 0.6]);
        end
    end
end

xticks(x); xticklabels(rem_types);
xlabel('REM type');
ylabel('NREM duration before REM (min)');
title('Baseline: NREM before REM (WT vs APP, by REM type)');
if hasWT && hasAPP
    legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
elseif hasWT
    legend({'WT'});
elseif hasAPP
    legend({'APP'});
end
set(gca,'Box','off','FontSize',12);

out_file = fullfile(out_dir,'baseline_prevNREM_before_REM_byType_APPvsWT.png');
saveas(gcf,out_file);

OUT = struct('success',true,'file',out_file,'rem_types',rem_types);
fprintf('✅ Saved NREM-before-REM plot to %s\n', out_file);
end
