function plot_bandpower_change_by_genotype(csvFile, state, bandsUse, outPng)
% plot_bandpower_change_by_genotype
%
% For a given STATE (e.g. "NREM" or "REM") and list of BANDS, this:
%   - uses baseline + ambtemp conditions from bandpower CSV
%   - for each mouse & band, computes:
%         Delta = Power_dB(ambtemp) - Power_dB(baseline)
%   - makes a figure with one subplot per band:
%         2 groups: WT vs APP (Delta values)
%         dots = mice, mean±SEM markers
%         stars:
%             WT:  Δ vs 0          (effect of ambtemp in WT)
%             APP: Δ vs 0          (effect of ambtemp in APP)
%             WT vs APP on Δ       (difference in change)
%
% EXAMPLE:
%   plot_bandpower_change_by_genotype( ...
%       'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
%       "REM", ["Sigma","Beta","lGamma1","lGamma2"], []);

    if nargin < 4 || isempty(outPng)
        outPng = sprintf('EEG_PSD_AllMice_AmbTemp/%s_delta_baseline_vs_ambtemp.png', state);
    end

    % ---------- Load & prepare table ----------
    T = readtable(csvFile);
    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = string(T.Condition);
    T.State     = string(T.State);
    T.Band      = string(T.Band);

    % Only baseline + ambtemp, chosen state
    mask = ismember(lower(T.Condition), ["baseline","ambtemp"]) & T.State == state;
    T = T(mask, :);
    if isempty(T)
        error('No rows for state %s with baseline/ambtemp in %s', state, csvFile);
    end

    % Normalise condition labels
    condLower = lower(T.Condition);
    T.Condition(condLower=="baseline") = "baseline";
    T.Condition(condLower=="ambtemp")  = "ambtemp";

    nb = numel(bandsUse);
    figure('Color','w','Position',[100 100 320*nb 450]);

    for bi = 1:nb
        band = bandsUse(bi);

        subplot(1, nb, bi); hold on;
        title(sprintf('%s – %s (amb - base)', state, band));

        % ----- subset for this band -----
        Tb = T(T.Band == band, :);
        if isempty(Tb)
            text(0.5,0.5,'No data','HorizontalAlignment','center');
            axis off;
            continue;
        end

        % ----- compute Delta per mouse -----
        mIDs = unique(Tb.MouseID);
        Delta = [];
        Geno  = strings(0,1);

        for i = 1:numel(mIDs)
            mid  = mIDs(i);
            rows = Tb(Tb.MouseID==mid, :);

            baseVals = rows.Power_dB(rows.Condition=="baseline");
            ambVals  = rows.Power_dB(rows.Condition=="ambtemp");

            if isempty(baseVals) || isempty(ambVals)
                continue;   % need both
            end

            baseMean = mean(baseVals,'omitnan');
            ambMean  = mean(ambVals,'omitnan');

            if isnan(baseMean) || isnan(ambMean)
                continue;
            end

            Delta(end+1,1) = ambMean - baseMean; %#ok<AGROW>
            Geno(end+1,1)  = rows.Genotype(1);   %#ok<AGROW>
        end

        if isempty(Delta)
            text(0.5,0.5,'No paired baseline+ambtemp for this band', ...
                 'HorizontalAlignment','center');
            axis off;
            continue;
        end

        % ----- split by genotype -----
        valsWT  = Delta(Geno=="WT");
        valsAPP = Delta(Geno=="APP");

        [mWT,  eWT]  = mean_sem(valsWT);
        [mAPP, eAPP] = mean_sem(valsAPP);

        x       = [1 2];
        labels  = {'WT','APP'};
        COL_WT  = [0.6 0.6 0.6];
        COL_APP = [0.392 0.584 0.929];

        % mean markers + SEM
        if ~isnan(mWT)
            plot(1, mWT, 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', COL_WT, 'MarkerEdgeColor','k');
            line([1 1], [mWT-eWT, mWT+eWT], 'Color','k','LineWidth',1.5);
        end
        if ~isnan(mAPP)
            plot(2, mAPP, 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', COL_APP, 'MarkerEdgeColor','k');
            line([2 2], [mAPP-eAPP, mAPP+eAPP], 'Color','k','LineWidth',1.5);
        end

        % individual points
        jitter = 0.08;
        plot_individuals(valsWT,  1, COL_WT,  jitter);
        plot_individuals(valsAPP, 2, COL_APP, jitter);

        set(gca,'XTick',x,'XTickLabel',labels);
        ylabel('\Delta Power (dB)  (ambtemp - baseline)');
        grid on;

        % ----- stats -----
        pWT   = safe_ttest0(valsWT);        % WT Δ vs 0
        pAPP  = safe_ttest0(valsAPP);       % APP Δ vs 0
        pDiff = safe_ttest2(valsWT,valsAPP);% WT vs APP Δ

        allVals = [valsWT; valsAPP];
        yMin = min(allVals);
        yMax = max(allVals);
        if isempty(yMin) || isnan(yMin), yMin = -1; end
        if isempty(yMax) || isnan(yMax), yMax = 1;  end
        yRange = max(yMax - yMin, 1);

        % WT vs 0 above WT
        add_star_single(1, yMax + 0.05*yRange, pWT);
        % APP vs 0 above APP
        add_star_single(2, yMax + 0.05*yRange, pAPP);
        % WT vs APP between groups
        add_star_pair([1 2], yMax + 0.20*yRange, pDiff);

        ylim([yMin - 0.1*yRange, yMax + 0.35*yRange]);

        fprintf('\n%s – %s:\n', state, band);
        fprintf('  WT:   mean Δ = %.3f dB (n=%d), p(Δ≠0)=%.3g\n', ...
                mWT, numel(valsWT), pWT);
        fprintf('  APP:  mean Δ = %.3f dB (n=%d), p(Δ≠0)=%.3g\n', ...
                mAPP, numel(valsAPP), pAPP);
        fprintf('  WT vs APP change: p = %.3g\n', pDiff);
    end

    sgtitle(sprintf('%s bandpower change (ambtemp - baseline)', state));

    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    saveas(gcf, outPng);
    fprintf('Saved Δ bandpower figure: %s\n', outPng);
end

% ---------- helper functions ----------

function [m,e] = mean_sem(x)
    x = x(~isnan(x));
    if isempty(x)
        m = NaN; e = NaN;
    else
        m = mean(x);
        e = std(x) / sqrt(numel(x));
    end
end

function plot_individuals(vals, xpos, col, jitter)
    vals = vals(~isnan(vals));
    n = numel(vals);
    if n == 0, return; end
    x = xpos + (rand(n,1)-0.5)*jitter;
    scatter(x, vals, 25, col, 'filled', 'MarkerFaceAlpha',0.6);
end

function p = safe_ttest0(x)
    x = x(~isnan(x));
    if numel(x) >= 2
        [~,p] = ttest(x);   % test vs 0
    else
        p = NaN;
    end
end

function p = safe_ttest2(a,b)
    a = a(~isnan(a));
    b = b(~isnan(b));
    if numel(a) >= 2 && numel(b) >= 2
        [~,p] = ttest2(a,b);
    else
        p = NaN;
    end
end

function add_star_single(x, y, p)
    if isnan(p), return; end
    if     p < 0.001, stars = '***';
    elseif p < 0.01,  stars = '**';
    elseif p < 0.05,  stars = '*';
    else,             stars = 'n.s.';
    end
    text(x, y, stars, 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize',8);
end

function add_star_pair(xpair, y, p)
    if isnan(p), return; end
    if     p < 0.001, stars = '***';
    elseif p < 0.01,  stars = '**';
    elseif p < 0.05,  stars = '*';
    else,             stars = 'n.s.';
    end
    plot(xpair, [y y], 'k-', 'LineWidth',1);
    text(mean(xpair), y + 0.01*(abs(y)+1), stars, ...
         'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',8);
end
