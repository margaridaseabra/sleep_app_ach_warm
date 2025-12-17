function plot_baseline_bandpower_focus(csvFile, outPng)
% plot_baseline_bandpower_focus(csvFile, outPng)
%
% Baseline-only comparison of WT vs APP in the bands that look most
% interesting from your stats:
%   NREM: Sigma, Beta, lGamma1
%   REM : Sigma, Beta, lGamma1, lGamma2
%
% Creates a figure with:
%   Top row  : NREM bands
%   Bottom row: REM bands
% Each subplot: per-mouse values (scatter + box) with t-test stars.
%
% INPUTS
%   csvFile : path to EEG_band_power_allmice.csv
%   outPng  : output PNG (default:
%             EEG_PSD_AllMice/Bandpower_baseline_focus.png)

    if nargin < 1 || isempty(csvFile)
        csvFile = 'EEG_PSD_AllMice/EEG_band_power_allmice.csv';
    end
    if nargin < 2 || isempty(outPng)
        outPng = 'EEG_PSD_AllMice/Bandpower_baseline_focus.png';
    end

    T = readtable(csvFile);

    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = string(T.Condition);
    T.State     = string(T.State);
    T.Band      = string(T.Band);

    % ---- baseline only ----
    maskBase = lower(T.Condition) == "baseline";
    Tb = T(maskBase, :);
    if isempty(Tb)
        error('No baseline rows found in table.');
    end

    % Bands of interest
    bandsNREM = ["Sigma","Beta","lGamma1"];
    bandsREM  = ["Sigma","Beta","lGamma1","lGamma2"];

    figure('Color','w','Position',[100 100 1200 600]);

    % ---------- NREM row ----------
    for bi = 1:numel(bandsNREM)
        subplot(2, max(numel(bandsNREM), numel(bandsREM)), bi); hold on;

        band = bandsNREM(bi);
        mask = Tb.State == "NREM" & Tb.Band == band;

        if ~any(mask)
            title(sprintf('NREM %s (no data)', band));
            axis off;
            continue;
        end

        sub = Tb(mask, :);
        title(sprintf('NREM %s', band));
        ylabel('Power (dB)');
        [p, d] = box_with_scatter_and_stats(sub);

        text(0.5, 0.9, sprintf('p = %.3g, d = %.2f', p, d), ...
            'Units','normalized', 'HorizontalAlignment','center', ...
            'FontSize', 8);
    end

    % ---------- REM row ----------
    offset = max(numel(bandsNREM), numel(bandsREM));
    for bi = 1:numel(bandsREM)
        subplot(2, offset, offset + bi); hold on;

        band = bandsREM(bi);
        mask = Tb.State == "REM" & Tb.Band == band;

        if ~any(mask)
            title(sprintf('REM %s (no data)', band));
            axis off;
            continue;
        end

        sub = Tb(mask, :);
        title(sprintf('REM %s', band));
        ylabel('Power (dB)');
        [p, d] = box_with_scatter_and_stats(sub);

        text(0.5, 0.9, sprintf('p = %.3g, d = %.2f', p, d), ...
            'Units','normalized', 'HorizontalAlignment','center', ...
            'FontSize', 8);
    end

    sgtitle('Baseline band power (WT vs APP) – key AD-relevant bands');

    % Save
    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    saveas(gcf, outPng);
    fprintf('Saved baseline focus bandpower figure: %s\n', outPng);
end

% ---------- Helper: WT vs APP box + scatter + stats ----------
function [p, d] = box_with_scatter_and_stats(sub)

    genotypes = unique(sub.Genotype);
    if numel(genotypes) ~= 2
        error('Expected exactly 2 genotypes for this band/state.');
    end

    g1 = genotypes(1);
    g2 = genotypes(2);

    vals1 = sub.Power_dB(sub.Genotype == g1);
    vals2 = sub.Power_dB(sub.Genotype == g2);

    COL_1 = [0.6 0.6 0.6];        % WT-ish
    COL_2 = [0.392 0.584 0.929];  % APP-ish

    % Boxplot
    boxplot([vals1; vals2], ...
            [repmat({char(g1)}, numel(vals1), 1); ...
             repmat({char(g2)}, numel(vals2), 1)]);
    hold on;

    % Scatter with jitter
    jitter = 0.05;
    x1 = ones(size(vals1));
    x2 = ones(size(vals2)) * 2;

    scatter(x1 + (rand(size(x1))-0.5)*jitter, vals1, 40, COL_1, 'filled');
    scatter(x2 + (rand(size(x2))-0.5)*jitter, vals2, 40, COL_2, 'filled');

    set(gca, 'XTick', [1 2], 'XTickLabel', {char(g1), char(g2)});
    grid on;

    % t-test & Cohen's d
    if numel(vals1) >= 2 && numel(vals2) >= 2
        [~, p] = ttest2(vals1, vals2);
    else
        p = NaN;
    end

    m1 = mean(vals1);
    m2 = mean(vals2);
    s1 = std(vals1);
    s2 = std(vals2);
    n1 = numel(vals1);
    n2 = numel(vals2);
    sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d  = (m2 - m1) / sp;

    % Stars
    if isnan(p)
        stars = 'n.s.';
    elseif p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end

    allVals = [vals1; vals2];
    yMax = max(allVals);
    yMin = min(allVals);
    yRange = yMax - yMin;
    if yRange == 0
        yRange = 1;
    end
    yStar = yMax + 0.1*yRange;

    plot([1 2],[yStar yStar],'k-','LineWidth',1);
    text(1.5, yStar + 0.02*yRange, stars, ...
         'HorizontalAlignment','center', 'FontWeight','bold');
end
