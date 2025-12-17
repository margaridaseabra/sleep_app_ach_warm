function plot_bandpower_cond_stats_ambtemp(csvFile, name, outPng, compareMode)
% plot_bandpower_cond_stats_ambtemp(csvFile, name, outPng, compareMode)
%
% Modes:
%   1) Genotype comparison (default, same as your original behaviour):
%        compareMode = 'genotype'
%        name        = conditionName  (e.g. 'baseline' or 'ambtemp')
%
%      -> For the given CONDITION, compares WT vs APP for each band & state.
%
%      Example:
%        plot_bandpower_cond_stats_ambtemp( ...
%            'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
%            'ambtemp', []);
%
%   2) Condition comparison within one genotype:
%        compareMode = 'condition'
%        name        = genotypeName   (e.g. 'WT' or 'APP')
%
%      -> For the given GENOTYPE, compares CONDITIONS (e.g. baseline vs ambtemp)
%         for each band & state.
%
%      Example:
%        plot_bandpower_cond_stats_ambtemp( ...
%            'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
%            'APP', [], 'condition');
%
%   outPng:
%      If empty, a default filename is created in the same folder as csvFile.

    if nargin < 4 || isempty(compareMode)
        compareMode = 'genotype';   % backwards compatible default
    end

    % ---- default output filename if not provided ----
    if nargin < 3 || isempty(outPng)
        [csvDir, ~, ~] = fileparts(csvFile);
        if isempty(csvDir), csvDir = '.'; end

        switch lower(compareMode)
            case 'genotype'
                outPng = fullfile(csvDir, sprintf('GroupBandPower_%s_stats.png', char(name)));
            case 'condition'
                outPng = fullfile(csvDir, sprintf('GroupBandPower_%s_byCondition_stats.png', char(name)));
            otherwise
                error('Unknown compareMode "%s". Use "genotype" or "condition".', compareMode);
        end
    end

    % ---- load table ----
    tbl = readtable(csvFile);

    % enforce string for key columns
    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);
    tbl.State     = string(tbl.State);
    tbl.Band      = string(tbl.Band);

    % basic sets
    genotypes = ["WT","APP"];   % adjust if needed
    bandsUse  = ["Delta","Theta","Sigma","Beta","lGamma1","lGamma2","hGamma"];
    states    = ["Wake","NREM","REM"];

    % common colours for genotype comparison
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.392 0.584 0.929];

    % ============================================================
    %  MODE 1: Compare GENOTYPES within one CONDITION
    % ============================================================
    if strcmpi(compareMode, 'genotype')
        conditionName = string(name);

        % filter condition (case-insensitive)
        maskCond = strcmpi(tbl.Condition, conditionName);
        tblCond  = tbl(maskCond, :);

        if isempty(tblCond)
            error('No rows with Condition == "%s" found.', conditionName);
        end

        nSt   = numel(states);
        nb    = numel(bandsUse);
        ng    = numel(genotypes);

        figure('Color','w','Position',[100 100 1400 400]);

        for si = 1:nSt
            st = states(si);
            subplot(1, nSt, si); hold on;

            maskSt = (tblCond.State == st) & ismember(tblCond.Band, bandsUse);
            sub    = tblCond(maskSt, :);

            if isempty(sub)
                text(0.5,0.5,'No data','HorizontalAlignment','center');
                axis off;
                continue;
            end

            % Mean + SEM: bands x genotypes
            M = nan(nb, ng);
            E = nan(nb, ng);
            P = nan(nb, 1);   % p-values WT vs APP per band

            for bi = 1:nb
                bname = bandsUse(bi);
                valsWT  = sub.Power_dB(sub.Band == bname & sub.Genotype == genotypes(1));
                valsAPP = sub.Power_dB(sub.Band == bname & sub.Genotype == genotypes(2));

                valsWT  = valsWT(~isnan(valsWT));
                valsAPP = valsAPP(~isnan(valsAPP));

                if ~isempty(valsWT)
                    M(bi,1) = mean(valsWT);
                    E(bi,1) = std(valsWT) ./ max(1,sqrt(numel(valsWT)));
                end
                if ~isempty(valsAPP)
                    M(bi,2) = mean(valsAPP);
                    E(bi,2) = std(valsAPP) ./ max(1,sqrt(numel(valsAPP)));
                end

                if numel(valsWT) >= 2 && numel(valsAPP) >= 2
                    [~, p] = ttest2(valsWT, valsAPP);
                    P(bi)  = p;
                else
                    P(bi) = NaN;
                end
            end

            % bar plot
            bh = bar(M, 'grouped');
            bh(1).FaceColor = COL_WT;
            bh(2).FaceColor = COL_APP;

            % Error bars
            hold on;
            for gi = 1:ng
                x = bh(gi).XEndPoints;
                errorbar(x, M(:,gi), E(:,gi), 'k', 'LineStyle','none','LineWidth',1);
            end

            % Stars for WT vs APP per band
            xWT  = bh(1).XEndPoints;
            xAPP = bh(2).XEndPoints;

            allMax = max(M(:) + E(:), [], 'omitnan');
            if isempty(allMax) || isnan(allMax), allMax = 0; end
            allMin = min(M(:) - E(:), [], 'omitnan');
            if isempty(allMin) || isnan(allMin), allMin = 0; end
            yRange = max(allMax - allMin, 1);
            yBase  = allMax;

            for bi = 1:nb
                p = P(bi);
                if isnan(p), continue; end

                if     p < 0.001, stars = '***';
                elseif p < 0.01,  stars = '**';
                elseif p < 0.05,  stars = '*';
                else,             stars = '';
                end
                if isempty(stars), continue; end

                xMid  = mean([xWT(bi), xAPP(bi)]);
                yStar = yBase + 0.05*yRange;

                plot([xWT(bi) xAPP(bi)], [yStar yStar], 'k-', 'LineWidth',1);
                text(xMid, yStar + 0.02*yRange, stars, ...
                     'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                     'FontSize',8, 'FontWeight','bold');
            end

            set(gca, 'XTick', 1:nb, ...
                     'XTickLabel', bandsUse, ...
                     'XTickLabelRotation', 45, ...
                     'FontSize', 8);

            if si == 1
                ylabel('Power (dB)');
            end

            title(st);
            if si == nSt
                legend(genotypes, 'Location','southoutside', 'Orientation','horizontal');
            end

            box off;
            grid on;
        end

        sgtitle(sprintf('EEG band power – %s (WT vs APP)', conditionName));

    % ============================================================
    %  MODE 2: Compare CONDITIONS within one GENOTYPE
    % ============================================================
    elseif strcmpi(compareMode, 'condition')
        genotypeName = string(name);

        % filter genotype (case-insensitive)
        maskG = strcmpi(tbl.Genotype, genotypeName);
        tblG  = tbl(maskG, :);

        if isempty(tblG)
            error('No rows with Genotype == "%s" found.', genotypeName);
        end

        % conditions present for this genotype
        condList = unique(tblG.Condition, 'stable');   % keep original order

        % (Optional) gently prioritise 'baseline' and 'ambtemp' first, if they exist
        priority = ["baseline","ambtemp"];
        ordered  = strings(0,1);
        for pr = priority
            idx = find(strcmpi(condList, pr), 1);
            if ~isempty(idx)
                ordered(end+1,1) = condList(idx); %#ok<AGROW>
            end
        end
        % append any remaining conditions
        rest = condList(~ismember(lower(condList), lower(ordered)));
        condList = [ordered; rest];

        nCond = numel(condList);
        if nCond < 1
            error('No conditions found for genotype %s.', genotypeName);
        end

        nSt = numel(states);
        nb  = numel(bandsUse);

        % colour per condition
        condColors = lines(nCond);

        figure('Color','w','Position',[100 100 1400 400]);

        for si = 1:nSt
            st = states(si);
            subplot(1, nSt, si); hold on;

            maskSt = (tblG.State == st) & ismember(tblG.Band, bandsUse);
            sub    = tblG(maskSt, :);

            if isempty(sub)
                text(0.5,0.5,'No data','HorizontalAlignment','center');
                axis off;
                continue;
            end

            % Mean + SEM: bands x conditions
            M = nan(nb, nCond);
            E = nan(nb, nCond);

            % p-values between *first two* conditions (e.g. baseline vs ambtemp)
            P = nan(nb, 1);

            for bi = 1:nb
                bname = bandsUse(bi);

                for ci = 1:nCond
                    cname = condList(ci);
                    vals = sub.Power_dB(sub.Band == bname & sub.Condition == cname);
                    vals = vals(~isnan(vals));

                    if ~isempty(vals)
                        M(bi,ci) = mean(vals);
                        E(bi,ci) = std(vals) ./ max(1,sqrt(numel(vals)));
                    end
                end

                % only compute stats if we have at least 2 conditions
                if nCond >= 2
                    c1 = condList(1);
                    c2 = condList(2);
                    v1 = sub.Power_dB(sub.Band == bname & sub.Condition == c1);
                    v2 = sub.Power_dB(sub.Band == bname & sub.Condition == c2);
                    v1 = v1(~isnan(v1));
                    v2 = v2(~isnan(v2));

                    if numel(v1) >= 2 && numel(v2) >= 2
                        [~, p] = ttest2(v1, v2);
                        P(bi)  = p;
                    end
                end
            end

            % bar plot: bands x conditions
            bh = bar(M, 'grouped');
            for ci = 1:nCond
                bh(ci).FaceColor = condColors(ci,:);
            end

            % Error bars
            hold on;
            for ci = 1:nCond
                x = bh(ci).XEndPoints;
                errorbar(x, M(:,ci), E(:,ci), 'k', 'LineStyle','none','LineWidth',1);
            end

            % Stars for first two conditions (if present), e.g. baseline vs ambtemp
            if nCond >= 2
                x1 = bh(1).XEndPoints;
                x2 = bh(2).XEndPoints;

                allMax = max(M(:) + E(:), [], 'omitnan');
                if isempty(allMax) || isnan(allMax), allMax = 0; end
                allMin = min(M(:) - E(:), [], 'omitnan');
                if isempty(allMin) || isnan(allMin), allMin = 0; end
                yRange = max(allMax - allMin, 1);
                yBase  = allMax;

                for bi = 1:nb
                    p = P(bi);
                    if isnan(p), continue; end

                    if     p < 0.001, stars = '***';
                    elseif p < 0.01,  stars = '**';
                    elseif p < 0.05,  stars = '*';
                    else,             stars = '';
                    end
                    if isempty(stars), continue; end

                    xMid  = mean([x1(bi), x2(bi)]);
                    yStar = yBase + 0.05*yRange;

                    plot([x1(bi) x2(bi)], [yStar yStar], 'k-', 'LineWidth',1);
                    text(xMid, yStar + 0.02*yRange, stars, ...
                         'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                         'FontSize',8, 'FontWeight','bold');
                end
            end

            set(gca, 'XTick', 1:nb, ...
                     'XTickLabel', bandsUse, ...
                     'XTickLabelRotation', 45, ...
                     'FontSize', 8);

            if si == 1
                ylabel('Power (dB)');
            end

            title(st);

            if si == nSt
                legend(cellstr(condList), 'Location','southoutside', 'Orientation','horizontal');
            end

            box off;
            grid on;
        end

        sgtitle(sprintf('EEG band power – genotype %s (conditions compared)', genotypeName));

    else
        error('Unknown compareMode "%s". Use "genotype" or "condition".', compareMode);
    end

    % ---- save figure ----
    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    saveas(gcf, outPng);
    fprintf('Saved band power stats figure (%s mode, name=%s): %s\n', ...
        compareMode, char(name), outPng);
end
