function plot_bandpower_baseline_stats(csvFile, outPng)
% plot_bandpower_baseline_stats(csvFile, outPng)
%
% Focuses ONLY on the "baseline" condition and compares WT vs APP.
% For each state (Wake, NREM, REM) and each band, it:
%   - computes group means ± SEM in dB
%   - runs an unpaired t-test (APP vs WT)
%   - adds significance stars above the bars (*, **, ***)
%
% INPUTS
%   csvFile : path to EEG_band_power_allmice.csv
%   outPng  : output PNG filename
%
% EXAMPLE
%   plot_bandpower_baseline_stats( ...
%       'EEG_PSD_AllMice/EEG_band_power_allmice.csv', []);
%

    if nargin < 1 || isempty(csvFile)
        csvFile = 'EEG_PSD_AllMice/EEG_band_power_allmice.csv';
    end
    if nargin < 2 || isempty(outPng)
        outPng = 'EEG_PSD_AllMice/GroupBandPower_baseline_stats.png';
    end

    % ---------- Load table ----------
    tbl = readtable(csvFile);

    % Make sure key columns are string (for easy filtering)
    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);
    tbl.State     = string(tbl.State);
    tbl.Band      = string(tbl.Band);

    % ---------- Focus on baseline only ----------
    condLower = lower(tbl.Condition);
    maskBase  = condLower == "baseline";

    if ~any(maskBase)
        error('No rows with Condition == "baseline" (case insensitive) were found.');
    end

    tblBase = tbl(maskBase, :);
    disp('Using only rows with Condition == "baseline".');

    % ---------- Settings ----------
    genotypes = ["WT","APP"];   % assuming these two
    bandsUse  = ["Delta","Theta","Sigma","Beta","lGamma1","lGamma2","hGamma"];
    states    = ["Wake","NREM","REM"];

    nb = numel(bandsUse);
    ng = numel(genotypes);
    ns = numel(states);

    % Colors (same as your other functions)
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.392 0.584 0.929];

    % ---------- Figure ----------
    figure('Color','w','Position',[100 100 1400 400]);

    % For printing stats summary later
    statsSummary = table();

    for si = 1:ns
        st = states(si);
        subplot(1, ns, si); hold on;

        % Subset baseline rows for this state
        maskState = tblBase.State == st;
        sub       = tblBase(maskState, :);

        if isempty(sub)
            text(0.5,0.5,'No data','HorizontalAlignment','center');
            axis off;
            continue;
        end

        % Mean + SEM matrices: bands x genotypes
        M = nan(nb, ng);
        E = nan(nb, ng);
        P = nan(nb, 1);   % p-values per band

        % First collect data, compute stats
        for bi = 1:nb
            bname = bandsUse(bi);

            valsWT  = sub.Power_dB(sub.Band == bname & sub.Genotype == genotypes(1));
            valsAPP = sub.Power_dB(sub.Band == bname & sub.Genotype == genotypes(2));

            valsWT  = valsWT(~isnan(valsWT));
            valsAPP = valsAPP(~isnan(valsAPP));

            % Means / SEM
            if ~isempty(valsWT)
                M(bi,1) = mean(valsWT);
                E(bi,1) = std(valsWT) ./ max(1, sqrt(numel(valsWT)));
            end
            if ~isempty(valsAPP)
                M(bi,2) = mean(valsAPP);
                E(bi,2) = std(valsAPP) ./ max(1, sqrt(numel(valsAPP)));
            end

            % Only run t-test if both groups have at least 2 mice
            if numel(valsWT) >= 2 && numel(valsAPP) >= 2
                [~, p] = ttest2(valsWT, valsAPP);
                P(bi)  = p;
            else
                P(bi) = NaN;  % not enough data
            end

            % Add to stats summary table (one row per band per state)
            statsSummary = [statsSummary; table( ...
                st, bname, ...
                numel(valsWT), numel(valsAPP), ...
                meanOrNaN(valsWT), meanOrNaN(valsAPP), ...
                P(bi), ...
                'VariableNames', { ...
                    'State','Band', ...
                    'n_WT','n_APP', ...
                    'Mean_WT_dB','Mean_APP_dB', ...
                    'p_value'})]; %#ok<AGROW>
        end

        % ---------- Plot grouped bars ----------
        bh = bar(M, 'grouped');
        bh(1).FaceColor = COL_WT;   % WT
        bh(2).FaceColor = COL_APP;  % APP

        % Error bars & stars
        hold on;
        for gi = 1:ng
            x = bh(gi).XEndPoints;
            errorbar(x, M(:,gi), E(:,gi), 'k', ...
                     'LineStyle','none','LineWidth',1);
        end

        % Add significance stars for each band
        xWT  = bh(1).XEndPoints;
        xAPP = bh(2).XEndPoints;

        yMaxAll = max(M(:) + E(:), [], 'omitnan');
        if isempty(yMaxAll) || isnan(yMaxAll)
            yMaxAll = 1;
        end
        % Space above bars for stars
        yMargin = 0.05 * abs(yMaxAll);
        if yMargin == 0
            yMargin = 0.5;
        end

        for bi = 1:nb
            p = P(bi);
            if isnan(p)
                continue; % no test
            end

            if p < 0.001
                stars = '***';
            elseif p < 0.01
                stars = '**';
            elseif p < 0.05
                stars = '*';
            else
                stars = '';
            end

            if ~isempty(stars)
                % x position midway between WT and APP bars for this band
                xMid = mean([xWT(bi), xAPP(bi)]);

                % y position a bit above the taller bar
                yBandMax = max(M(bi,:) + E(bi,:), [], 'omitnan');
                if isnan(yBandMax)
                    yBandMax = 0;
                end
                yStar = yBandMax + yMargin;

                % Draw a small horizontal line between the two bars (optional)
                plot([xWT(bi) xAPP(bi)], [yStar yStar], 'k-', 'LineWidth', 1);

                % Put the stars above
                text(xMid, yStar + 0.02*abs(yMaxAll), stars, ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','bottom', ...
                     'FontSize', 10, 'FontWeight','bold');
            end
        end

        set(gca, 'XTick', 1:nb, ...
                 'XTickLabel', bandsUse, ...
                 'XTickLabelRotation', 45, ...
                 'FontSize', 9);

        ylabel('Power (dB)');
        title(st);

        if si == ns
            legend(genotypes, 'Location','southoutside', ...
                   'Orientation','horizontal');
        end

        box off;
        grid on;
    end

    sgtitle('Baseline EEG Band Power (WT vs APP) – with Statistics');

    % Ensure output directory exists
    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    saveas(gcf, outPng);
    fprintf('Saved baseline stats band power figure: %s\n', outPng);

    % Print stats summary
    disp('--- Baseline WT vs APP t-test summary (per State x Band) ---');
    disp(statsSummary);

end

% --------- Helper: meanOrNaN ----------
function m = meanOrNaN(x)
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end
