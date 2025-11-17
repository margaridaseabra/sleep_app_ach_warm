function plot_bandpower_grid(csvFile, outPng)
% plot_bandpower_grid(csvFile, outPng)
%
% Uses EEG_band_power_allmice.csv and makes a single figure:
%   rows = Conditions (in the order they appear)
%   cols = States: Wake, NREM, REM
% Each subplot: grouped bars WT vs APP for each band (delta, theta, ...).

    if nargin < 2 || isempty(outPng)
        outPng = 'EEG_PSD_AllMice/GroupBandPower_grid.png';
    end

    tbl = readtable(csvFile);

    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);
    tbl.State     = string(tbl.State);
    tbl.Band      = string(tbl.Band);

    genotypes = ["WT","APP"];   % adjust if needed
    bandsUse  = ["Delta","Theta","Sigma","Beta","lGamma1","lGamma2","hGamma"];
    states    = ["Wake","NREM","REM"];

    % Conditions in the order they appear in the table
    [~, idxFirst] = unique(tbl.Condition, 'stable');
    conditions    = tbl.Condition(sort(idxFirst))';
    nCond = numel(conditions);
    nSt   = numel(states);
    nb    = numel(bandsUse);
    ng    = numel(genotypes);

    figure('Color','w','Position',[100 100 1400 900]);

    for ci = 1:nCond
        cond = conditions(ci);

        for si = 1:nSt
            st = states(si);
            subplot(nCond, nSt, (ci-1)*nSt + si); hold on;

            mask = tbl.Condition == cond & tbl.State == st & ismember(tbl.Band, bandsUse);
            sub  = tbl(mask, :);

            if isempty(sub)
                text(0.5,0.5,'No data','HorizontalAlignment','center');
                axis off;
                continue;
            end

            % Mean + SEM matrices: bands x genotypes
            M = nan(nb, ng);
            E = nan(nb, ng);

            for bi = 1:nb
                bname = bandsUse(bi);
                for gi = 1:ng
                    gname = genotypes(gi);

                    mBG  = sub.Band == bname & sub.Genotype == gname;
                    vals = sub.Power_dB(mBG);

                    if ~isempty(vals)
                        M(bi,gi) = mean(vals);
                        E(bi,gi) = std(vals) ./ sqrt(numel(vals));
                    end
                end
            end
            % Fixed genotype colours
            COL_WT  = [0.6 0.6 0.6];           % grey
            COL_APP = [0.392 0.584 0.929];     % cornflower blue-ish

            bh = bar(M, 'grouped');  % WT vs APP per band

            % Apply colours: assume genotypes = ["WT","APP"];
            bh(1).FaceColor = COL_WT;   % WT
            bh(2).FaceColor = COL_APP;  % APP
            
            % Add SEM error bars
            for gi = 1:ng
                x = bh(gi).XEndPoints;
                errorbar(x, M(:,gi), E(:,gi), 'k', ...
                         'LineStyle','none','LineWidth',1);
            end

            set(gca, 'XTick', 1:nb, ...
                     'XTickLabel', bandsUse, ...
                     'XTickLabelRotation', 45, ...
                     'FontSize', 8);

            if si == 1
                ylabel('Power (dB)');
            end

            if ci == 1
                title(st);
            end

            if ci == nCond && si == nSt
                legend(genotypes, 'Location','southoutside', 'Orientation','horizontal');
            end

            box off;
            grid on;
        end
    end

    sgtitle('EEG Power Spectrum – Band Power (WT vs APP)');

    % Ensure output directory exists
    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    saveas(gcf, outPng);
    fprintf('Saved band power grid figure: %s\n', outPng);
end
