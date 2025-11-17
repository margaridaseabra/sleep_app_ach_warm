function plot_bandpower_per_condition(csvFile, outDir)
% plot_bandpower_per_condition(csvFile, outDir)
%
% Creates one figure PER CONDITION:
%   columns = Wake, NREM, REM
%   bars    = WT vs APP band power (dB) for each band.
%
% INPUTS:
%   csvFile : path to EEG_band_power_allmice.csv
%   outDir  : output directory for PNGs

    if nargin < 2 || isempty(outDir)
        outDir = 'EEG_PSD_AllMice/BandPower_ByCondition';
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---------- Load table ----------
    tbl = readtable(csvFile);
    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);
    tbl.State     = string(tbl.State);
    tbl.Band      = string(tbl.Band);

    % Genotypes (order in legend/bars)
    genotypes = ["WT","APP"];      % change if your labels differ

    % Bands (x-axis order)
    bandsUse = ["Delta","Theta","Sigma","Beta","lGamma1","lGamma2","hGamma"];

    % States (columns)
    states = ["Wake","NREM","REM"];

    % Conditions (rows) – keep original order
    [~, idxFirst] = unique(tbl.Condition, 'stable');
    conditions    = tbl.Condition(sort(idxFirst))';

    % ---------- Loop over conditions ----------
    for ci = 1:numel(conditions)
        cond = conditions(ci);

        fig = figure('Color','w','Name',cond, ...
                     'Position',[100 100 1400 400]);

        for si = 1:numel(states)
            st = states(si);
            subplot(1, numel(states), si); hold on;

            % Subset for this condition + state
            mask = tbl.Condition == cond & tbl.State == st & ...
                   ismember(tbl.Band, bandsUse);
            sub  = tbl(mask, :);

            if isempty(sub)
                text(0.5,0.5,'No data', ...
                     'HorizontalAlignment','center');
                axis off;
                continue;
            end

            nb = numel(bandsUse);
            ng = numel(genotypes);

            M = nan(nb, ng);   % mean power
            E = nan(nb, ng);   % SEM

            for bi = 1:nb
                bname = bandsUse(bi);
                for gi = 1:ng
                    gname = genotypes(gi);

                    ix = sub.Band == bname & sub.Genotype == gname;
                    vals = sub.Power_dB(ix);

                    if ~isempty(vals)
                        M(bi, gi) = mean(vals);
                        E(bi, gi) = std(vals) ./ sqrt(numel(vals));
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
            % Error bars
            for gi = 1:ng
                x = bh(gi).XEndPoints;
                errorbar(x, M(:,gi), E(:,gi), 'k', ...
                         'LineStyle','none','LineWidth',1);
            end

            set(gca, 'XTick', 1:nb, ...
                     'XTickLabel', bandsUse, ...
                     'XTickLabelRotation', 45, ...
                     'FontSize', 8);

            ylim([min(M(:)-E(:))-5, max(M(:)+E(:))+5]); % quick sensible y-lims
            grid on; box off;

            if si == 1
                ylabel('Power (dB)');
            end
            title(st);

            if si == numel(states)
                legend(genotypes, 'Location','southoutside', ...
                       'Orientation','horizontal');
            end
        end

        sgtitle(sprintf('EEG Power Spectrum – Band Power (%s)', cond));

        % Save figure
        safeCond = regexprep(cond, '[^a-zA-Z0-9]', '');
        outPng   = fullfile(outDir, sprintf('BandPower_%s.png', safeCond));
        saveas(fig, outPng);
        fprintf('Saved band-power figure for %s: %s\n', cond, outPng);
    end
end
