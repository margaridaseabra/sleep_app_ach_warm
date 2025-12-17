function STATS = compare_ach_segment_features(OUT_list, labels)
% COMPARE_ACH_SEGMENT_FEATURES
%   Compare ACh onset features across multiple segments (episodes).
%
%   OUT_list : array or cell array of OUT structs returned by
%              plot_ach_eeg_segment.
%   labels   : cellstr with a label for each segment (optional).
%
% For each state (Wake, NREM, REM) it:
%   - builds boxplots of ΔF/F and slope across segments
%   - performs:
%       * ttest2 if there are 2 segments
%       * ANOVA if more than 2
%
% STATS.(state).deltaF / .slope contain p-values etc.

if ~iscell(OUT_list)
    OUT_list = num2cell(OUT_list);
end
nEp = numel(OUT_list);

if nargin < 2 || isempty(labels)
    labels = arrayfun(@(k) sprintf('Ep%d',k), 1:nEp, 'uni',0);
end

states = {'Wake','NREM','REM'};
STATS = struct();

for s = 1:numel(states)
    stName = states{s};

    all_dF = [];
    all_slope = [];
    grp_dF = [];
    grp_slope = [];

    for e = 1:nEp
        F = OUT_list{e}.features.(stName);
        if isempty(F) || F.n == 0, continue; end

        n1 = numel(F.deltaF_all);
        all_dF   = [all_dF;   F.deltaF_all(:)];
        grp_dF   = [grp_dF;   repmat(e,n1,1)];

        n2 = numel(F.slope_all);
        all_slope = [all_slope; F.slope_all(:)];
        grp_slope = [grp_slope; repmat(e,n2,1)];
    end

    % ---- boxplots ----
    if ~isempty(all_dF)
        figure('Color','w','Position',[300 300 800 350]);
        subplot(1,2,1);
        boxplot(all_dF, grp_dF, 'Labels', labels);
        ylabel('\DeltaF/F at onset');
        title(sprintf('%s: \\DeltaF/F per episode', stName));

        subplot(1,2,2);
        boxplot(all_slope, grp_slope, 'Labels', labels);
        ylabel('Slope (\DeltaF/F per s)');
        title(sprintf('%s: slopes per episode', stName));

        sgtitle(sprintf('%s comparison across %d episodes', stName, nEp), ...
                'FontWeight','bold');
    end

    ST = struct();

    % ---- stats: ΔF/F ----
    if numel(unique(grp_dF)) > 1
        if nEp == 2
            g1 = all_dF(grp_dF==1);
            g2 = all_dF(grp_dF==2);
            [~,p,~,stats] = ttest2(g1, g2, 'Vartype','unequal');
            ST.deltaF.test  = 'ttest2';
            ST.deltaF.p     = p;
            ST.deltaF.stats = stats;
        else
            [p, tbl, astats] = anova1(all_dF, grp_dF, 'off');
            ST.deltaF.test  = 'anova1';
            ST.deltaF.p     = p;
            ST.deltaF.tbl   = tbl;
            ST.deltaF.stats = astats;
        end
    else
        ST.deltaF.test  = 'none';
        ST.deltaF.p     = NaN;
        ST.deltaF.stats = [];
    end

    % ---- stats: slope ----
    if numel(unique(grp_slope)) > 1
        if nEp == 2
            g1 = all_slope(grp_slope==1);
            g2 = all_slope(grp_slope==2);
            [~,p,~,stats] = ttest2(g1, g2, 'Vartype','unequal');
            ST.slope.test  = 'ttest2';
            ST.slope.p     = p;
            ST.slope.stats = stats;
        else
            [p, tbl, astats] = anova1(all_slope, grp_slope, 'off');
            ST.slope.test  = 'anova1';
            ST.slope.p     = p;
            ST.slope.tbl   = tbl;
            ST.slope.stats = astats;
        end
    else
        ST.slope.test  = 'none';
        ST.slope.p     = NaN;
        ST.slope.stats = [];
    end

    STATS.(stName) = ST;
end
end
