function bandpower_ratio_baseline(csvFile, outPng)
% bandpower_ratio_baseline(csvFile, outPng)
%
% Computes and plots baseline band-power ratios:
%   - NREM: Sigma/Delta (as Sigma_dB - Delta_dB)
%   - REM : Theta/Delta (as Theta_dB - Delta_dB)
% for WT vs APP.
%
% INPUTS
%   csvFile : path to EEG_band_power_allmice.csv
%   outPng  : output PNG for the figure
%
% EXAMPLE
%   bandpower_ratio_baseline('EEG_PSD_AllMice/EEG_band_power_allmice.csv',[]);

    if nargin < 1 || isempty(csvFile)
        csvFile = 'EEG_PSD_AllMice/EEG_band_power_allmice.csv';
    end
    if nargin < 2 || isempty(outPng)
        outPng = 'EEG_PSD_AllMice/Bandpower_ratios_baseline.png';
    end

    % ----- Load table -----
    T = readtable(csvFile);
    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = string(T.Condition);
    T.State     = string(T.State);
    T.Band      = string(T.Band);

    % ----- Baseline only -----
    maskBase = lower(T.Condition) == "baseline";
    Tb = T(maskBase, :);
    if isempty(Tb)
        error('No baseline rows found in table.');
    end

    % Unique mouse × genotype combinations in baseline
    combo = unique(Tb(:, {'MouseID','Genotype'}));

    % Containers for NREM and REM ratio tables
    NREM_mouse   = [];
    NREM_geno    = [];
    NREM_ratio   = [];

    REM_mouse    = [];
    REM_geno     = [];
    REM_ratio    = [];

    % ----- Loop over mice -----
    for i = 1:height(combo)
        mid  = combo.MouseID(i);
        geno = combo.Genotype(i);

        maskMouse = Tb.MouseID == mid & Tb.Genotype == geno;

        % NREM Delta / Sigma
        maskNREM = maskMouse & Tb.State == "NREM";

        valsNREM_Delta = Tb.Power_dB(maskNREM & Tb.Band == "Delta");
        valsNREM_Sigma = Tb.Power_dB(maskNREM & Tb.Band == "Sigma");

        if ~isempty(valsNREM_Delta) && ~isempty(valsNREM_Sigma)
            dDelta = mean(valsNREM_Delta, 'omitnan');
            dSigma = mean(valsNREM_Sigma, 'omitnan');
            if ~isnan(dDelta) && ~isnan(dSigma)
                NREM_mouse   = [NREM_mouse; mid];
                NREM_geno    = [NREM_geno;  geno];
                NREM_ratio   = [NREM_ratio; dSigma - dDelta]; %#ok<AGROW>
            end
        end

        % REM Delta / Theta
        maskREM = maskMouse & Tb.State == "REM";

        valsREM_Delta = Tb.Power_dB(maskREM & Tb.Band == "Delta");
        valsREM_Theta = Tb.Power_dB(maskREM & Tb.Band == "Theta");

        if ~isempty(valsREM_Delta) && ~isempty(valsREM_Theta)
            rDelta = mean(valsREM_Delta, 'omitnan');
            rTheta = mean(valsREM_Theta, 'omitnan');
            if ~isnan(rDelta) && ~isnan(rTheta)
                REM_mouse    = [REM_mouse; mid];
                REM_geno     = [REM_geno;  geno];
                REM_ratio    = [REM_ratio; rTheta - rDelta]; %#ok<AGROW>
            end
        end
    end

    % Build ratio tables
    NREM = table(NREM_mouse, NREM_geno, NREM_ratio, ...
                 'VariableNames', {'MouseID','Genotype','SigmaMinusDelta_dB'});
    REM  = table(REM_mouse, REM_geno, REM_ratio, ...
                 'VariableNames', {'MouseID','Genotype','ThetaMinusDelta_dB'});

    % ----- Plot -----
    figure('Color','w','Position',[100 100 900 400]);

    % NREM panel
    subplot(1,2,1); hold on;
    title('NREM: Sigma/Delta (dB difference)');
    ylabel('Sigma_dB - Delta_dB');
    [pNREM, dNREM] = plotBoxWithStars(NREM, "SigmaMinusDelta_dB");

    % REM panel
    subplot(1,2,2); hold on;
    title('REM: Theta/Delta (dB difference)');
    ylabel('Theta_dB - Delta_dB');
    [pREM, dREM] = plotBoxWithStars(REM, "ThetaMinusDelta_dB");

    sgtitle('Baseline bandpower ratios (WT vs APP)');

    % Save
    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    saveas(gcf, outPng);

    fprintf('Saved bandpower ratio figure: %s\n', outPng);
    fprintf('NREM Sigma/Delta: p = %.4g, d = %.3f\n', pNREM, dNREM);
    fprintf('REM Theta/Delta : p = %.4g, d = %.3f\n', pREM, dREM);
end

% ---------- Helper: boxplot + individual points + t-test + Cohen''s d ----------
function [p, d] = plotBoxWithStars(Tsub, valueVar)

    genotypes = unique(Tsub.Genotype);
    if numel(genotypes) ~= 2
        error('Expected exactly 2 genotypes in the table.');
    end
    g1 = genotypes(1);
    g2 = genotypes(2);

    COL_1  = [0.6 0.6 0.6];          % WT-ish
    COL_2  = [0.392 0.584 0.929];    % APP-ish

    vals1 = Tsub.(valueVar)(Tsub.Genotype == g1);
    vals2 = Tsub.(valueVar)(Tsub.Genotype == g2);

    % Boxplot
    boxplot([vals1; vals2], ...
            [repmat({char(g1)},numel(vals1),1); ...
             repmat({char(g2)},numel(vals2),1)]);
    hold on;

    % Jittered scatter
    x1 = ones(size(vals1)) * 1;
    x2 = ones(size(vals2)) * 2;
    jitter = 0.05;

    scatter(x1 + (rand(size(x1))-0.5)*jitter, vals1, 40, COL_1, 'filled');
    scatter(x2 + (rand(size(x2))-0.5)*jitter, vals2, 40, COL_2, 'filled');

    set(gca,'XTick',[1 2],'XTickLabel',{char(g1), char(g2)});
    grid on;

    % t-test & effect size
    if numel(vals1) >= 2 && numel(vals2) >= 2
        [~, p] = ttest2(vals1, vals2);
    else
        p = NaN;
    end

    % Cohen's d (pooled SD)
    m1 = mean(vals1);
    m2 = mean(vals2);
    s1 = std(vals1);
    s2 = std(vals2);
    n1 = numel(vals1);
    n2 = numel(vals2);
    sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d = (m2 - m1) / sp;

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

    yMax   = max([vals1; vals2]);
    yMin   = min([vals1; vals2]);
    yRange = yMax - yMin;
    if yRange == 0
        yRange = 1;
    end
    yStar = yMax + 0.1*yRange;

    plot([1 2],[yStar yStar],'k-','LineWidth',1);
    text(1.5, yStar + 0.02*yRange, stars, ...
         'HorizontalAlignment','center', 'FontWeight','bold');
end
