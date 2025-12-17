function OUT = plot_NREM_sequence_patterns(NREM_TRANS, out_dir)
% Uses NREM_TRANS from run_REM_NREM_MA_transition_analysis.
% Figure 1: per-mouse probability of each sequence type:
%   NREM_REM_WAKE vs NREM_WAKE  (baseline only, WT vs APP)
% Figure 2: within NREM_REM_WAKE, fraction of Short/Medium/Long REM.

if nargin < 2 || isempty(out_dir), out_dir = pwd; end
if ~isfolder(out_dir), mkdir(out_dir); end

T = NREM_TRANS;
is_baseline = lower(strtrim(T.condition))=="baseline";
T = T(is_baseline, :);

if isempty(T)
    warning('No baseline NREM sequences found.');
    OUT.success = false; return;
end

% ---------------- Fig 1: NREM_REM_WAKE vs NREM_WAKE  ----------------
mice = unique(T.mouse);
SEQ = table();

for m = 1:numel(mice)
    mid = mice(m);
    maskM = T.mouse==mid;
    TM = T(maskM,:);

    geno = unique(TM.genotype); geno = geno(1);
    nAll = height(TM);

    p_nrw   = mean(TM.seq_type=="NREM_REM_WAKE");
    p_nw    = mean(TM.seq_type=="NREM_WAKE");

    SEQ = [SEQ; table(mid, geno, p_nrw, p_nw, ...
        'VariableNames',{'mouse','genotype','p_NREM_REM_WAKE','p_NREM_WAKE'})]; %#ok<AGROW>
end

hasWT  = any(SEQ.genotype=="WT");
hasAPP = any(SEQ.genotype=="APP");

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

figure('Color','w'); hold on;
x = 1:2; barWidth = 0.35;
jitterFrac = 0.25; y_offset = 0.03;

labels = {'NREM→REM→Wake','NREM→Wake'};

% genotype means
meanWT  = nan(1,2); semWT  = nan(1,2);
meanAPP = nan(1,2); semAPP = nan(1,2);

if hasWT
    mWT = SEQ.genotype=="WT";
    v1 = SEQ.p_NREM_REM_WAKE(mWT);
    v2 = SEQ.p_NREM_WAKE(mWT);
    meanWT(1) = mean(v1); semWT(1) = std(v1)/sqrt(numel(v1));
    meanWT(2) = mean(v2); semWT(2) = std(v2)/sqrt(numel(v2));
end
if hasAPP
    mAPP = SEQ.genotype=="APP";
    v1 = SEQ.p_NREM_REM_WAKE(mAPP);
    v2 = SEQ.p_NREM_WAKE(mAPP);
    meanAPP(1) = mean(v1); semAPP(1) = std(v1)/sqrt(numel(v1));
    meanAPP(2) = mean(v2); semAPP(2) = std(v2)/sqrt(numel(v2));
end

if hasWT
    bar(x-barWidth/2, meanWT, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
    errorbar(x-barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
end
if hasAPP
    bar(x+barWidth/2, meanAPP, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
    errorbar(x+barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
end

% mouse dots + IDs
for i = 1:2
    if hasWT
        mWT = SEQ.genotype=="WT";
        if i==1, v = SEQ.p_NREM_REM_WAKE(mWT); else, v = SEQ.p_NREM_WAKE(mWT); end
        mids = SEQ.mouse(mWT);
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
        mAPP = SEQ.genotype=="APP";
        if i==1, v = SEQ.p_NREM_REM_WAKE(mAPP); else, v = SEQ.p_NREM_WAKE(mAPP); end
        mids = SEQ.mouse(mAPP);
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

xticks(x); xticklabels(labels);
ylabel('Proportion of NREM bouts');
title('Baseline: NREM sequence patterns (WT vs APP)');
if hasWT && hasAPP
    legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
end
set(gca,'Box','off','FontSize',12);

out1 = fullfile(out_dir,'baseline_NREM_sequence_probs_APPvsWT.png');
saveas(gcf,out1);

% ---------------- Fig 2: REM size in NREM_REM_WAKE  ----------------
Ttriad = T(T.seq_type=="NREM_REM_WAKE" & NREM_TRANS.seq_rem_dur_s>0, :);
if isempty(Ttriad)
    warning('No NREM_REM_WAKE triads found.'); 
    OUT.success=true; OUT.files = struct('seq_probs',out1,'rem_type_triads',''); 
    return;
end

rem_types_all = ["Short","Medium","Long"];
rem_types = rem_types_all(ismember(rem_types_all, unique(Ttriad.seq_rem_type)));

% per-mouse proportions of each rem_type among triads
TR = table();
mice = unique(Ttriad.mouse);
for m = 1:numel(mice)
    mid = mice(m);
    maskM = Ttriad.mouse==mid;
    TM = Ttriad(maskM,:);
    geno = unique(TM.genotype); geno = geno(1);

    nTriad = height(TM);
    pShort  = mean(TM.seq_rem_type=="Short");
    pMedium = mean(TM.seq_rem_type=="Medium");
    pLong   = mean(TM.seq_rem_type=="Long");

    TR = [TR; table(mid, geno, nTriad, pShort, pMedium, pLong, ...
        'VariableNames',{'mouse','genotype','nTriad','pShort','pMedium','pLong'})]; %#ok<AGROW>
end

hasWT  = any(TR.genotype=="WT");
hasAPP = any(TR.genotype=="APP");

figure('Color','w'); hold on;
x = 1:numel(rem_types); barWidth = 0.35;
jitterFrac = 0.25; y_offset = 0.03;

meanWT  = nan(1,numel(rem_types)); semWT  = nan(1,numel(rem_types));
meanAPP = nan(1,numel(rem_types)); semAPP = nan(1,numel(rem_types));

for i = 1:numel(rem_types)
    rt = rem_types(i);
    pname = ['p' char(rt)];  % pShort, pMedium, pLong

    if hasWT
        mWT = TR.genotype=="WT";
        vWT = TR.(pname)(mWT); vWT = vWT(~isnan(vWT));
        if ~isempty(vWT)
            meanWT(i) = mean(vWT); semWT(i) = std(vWT)/sqrt(numel(vWT));
        end
    end
    if hasAPP
        mAPP = TR.genotype=="APP";
        vAPP = TR.(pname)(mAPP); vAPP = vAPP(~isnan(vAPP));
        if ~isempty(vAPP)
            meanAPP(i) = mean(vAPP); semAPP(i) = std(vAPP)/sqrt(numel(vAPP));
        end
    end
end

if hasWT
    bar(x-barWidth/2, meanWT, barWidth, 'FaceColor',COL_WT,'EdgeColor','none');
    errorbar(x-barWidth/2, meanWT, semWT, 'k','LineStyle','none','LineWidth',1);
end
if hasAPP
    bar(x+barWidth/2, meanAPP, barWidth, 'FaceColor',COL_APP,'EdgeColor','none');
    errorbar(x+barWidth/2, meanAPP, semAPP, 'k','LineStyle','none','LineWidth',1);
end

for i = 1:numel(rem_types)
    rt = rem_types(i);
    pname = ['p' char(rt)];
    if hasWT
        mWT = TR.genotype=="WT";
        v = TR.(pname)(mWT); mids = TR.mouse(mWT);
        xw = (x(i)-barWidth/2)+(rand(size(v))-0.5)*barWidth*jitterFrac;
        plot(xw, v, '.', 'Color',[0.3 0.3 0.3]);
        for j = 1:numel(v)
            text(xw(j), v(j)+y_offset, char(mids(j)), ...
                 'Rotation',45,'FontSize',8, ...
                 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                 'Color',[0.2 0.2 0.2]);
        end
    end
    if hasAPP
        mAPP = TR.genotype=="APP";
        v = TR.(pname)(mAPP); mids = TR.mouse(mAPP);
        xa = (x(i)+barWidth/2)+(rand(size(v))-0.5)*barWidth*jitterFrac;
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
ylabel('Proportion of NREM→REM→Wake triads');
xlabel('REM type');
title('REM size in NREM→REM→Wake sequences (baseline)');
if hasWT && hasAPP
    legend({'WT','APP'},'Location','northoutside','Orientation','horizontal');
end
set(gca,'Box','off','FontSize',12);

out2 = fullfile(out_dir,'baseline_NREM_REM_WAKE_REMsize_APPvsWT.png');
saveas(gcf,out2);

OUT.success = true;
OUT.files = struct('seq_probs',out1,'rem_type_triads',out2);
fprintf('✅ Saved NREM sequence plots to:\n  %s\n  %s\n', out1, out2);
end
