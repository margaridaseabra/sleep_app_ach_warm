function plot_bandpower_geno_cond_anova(csvFile, state, bandsUse, outPng)
% plot_bandpower_geno_cond_anova
%
% For a given STATE (e.g. "NREM" or "REM") and list of BANDS, this:
%   - uses baseline + ambtemp conditions from the bandpower CSV
%   - fits: Power_dB ~ Genotype*Condition + (1|MouseID)
%   - does post-hocs:
%       WT vs APP @ baseline      (unpaired)
%       WT vs APP @ ambtemp       (unpaired)
%       baseline vs ambtemp in WT (paired)
%       baseline vs ambtemp in APP (paired)
%   - plots 4 groups per band with stars:
%       1: WT-baseline
%       2: WT-ambtemp
%       3: APP-baseline
%       4: APP-ambtemp
%
% EXAMPLE:
%   plot_bandpower_geno_cond_anova( ...
%       'EEG_PSD_AllMice_AmbTemp/EEG_band_power_allmice.csv', ...
%       "NREM", ["Sigma","Beta","lGamma1"], []);
%

    if nargin < 4 || isempty(outPng)
        outPng = sprintf('EEG_PSD_AllMice_AmbTemp/Bandpower_%s_baseline_vs_ambtemp.png', state);
    end

    % ---------- Load & prepare table ----------
    T = readtable(csvFile);

    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = string(T.Condition);
    T.State     = string(T.State);
    T.Band      = string(T.Band);

    % Only baseline + ambtemp for this STATE
    mask = ismember(lower(T.Condition), ["baseline","ambtemp"]) & T.State == state;
    T = T(mask, :);
    if isempty(T)
        error('No rows for state %s with baseline/ambtemp in %s', state, csvFile);
    end

    % Normalise condition labels
    condLower = lower(T.Condition);
    T.Condition(condLower == "baseline") = "baseline";
    T.Condition(condLower == "ambtemp")  = "ambtemp";

    % Categorical for model
    T.MouseID   = categorical(T.MouseID);
    T.Genotype  = categorical(T.Genotype);
    T.Condition = categorical(T.Condition);
    T.Band      = string(T.Band);  % keep bands as string

    nb = numel(bandsUse);
    figure('Color','w','Position',[100 100 300*nb 500]);

    for bi = 1:nb
        band = bandsUse(bi);

        subplot(1, nb, bi); hold on;
        title(sprintf('%s – %s', state, band));

        % ---------- Subset for this band ----------
        Tb = T(T.Band == band, :);
        if isempty(Tb)
            text(0.5,0.5,'No data','HorizontalAlignment','center');
            axis off;
            continue;
        end

        % ---------- Mixed model ANOVA ----------
        try
            lme = fitlme(Tb, 'Power_dB ~ Genotype*Condition + (1|MouseID)');
            anovTbl = anova(lme);
            fprintf('\n==============================================\n');
            fprintf('Mixed model for %s – %s\n', state, band);
            disp(anovTbl);
        catch ME
            warning('fitlme failed for %s-%s: %s', state, band, ME.message);
        end

        % ---------- Extract group values ----------
        genoLevels = categories(Tb.Genotype);
        if numel(genoLevels) ~= 2
            warning('Expected 2 genotypes, got %d', numel(genoLevels));
        end
        G_WT  = genoLevels{1};
        G_APP = genoLevels{2};

        vals_WT_base  = getVals(Tb, G_WT,  "baseline");
        vals_WT_amb   = getVals(Tb, G_WT,  "ambtemp");
        vals_APP_base = getVals(Tb, G_APP, "baseline");
        vals_APP_amb  = getVals(Tb, G_APP, "ambtemp");

        [m1,e1] = mean_sem(vals_WT_base);
        [m2,e2] = mean_sem(vals_WT_amb);
        [m3,e3] = mean_sem(vals_APP_base);
        [m4,e4] = mean_sem(vals_APP_amb);

        means = [m1 m2 m3 m4];
        sems  = [e1 e2 e3 e4];
        x     = 1:4;
        labels = { ...
            sprintf('%s-base', char(G_WT)), ...
            sprintf('%s-amb',  char(G_WT)), ...
            sprintf('%s-base', char(G_APP)), ...
            sprintf('%s-amb',  char(G_APP)) };

        COL_WT  = [0.6 0.6 0.6];
        COL_APP = [0.392 0.584 0.929];
        cols    = [COL_WT; COL_WT; COL_APP; COL_APP];

        % Mean ± SEM markers
        for i = 1:4
            if isnan(means(i)), continue; end
            plot(x(i), means(i), 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', cols(i,:), 'MarkerEdgeColor','k');
            line([x(i) x(i)], [means(i)-sems(i), means(i)+sems(i)], ...
                 'Color','k','LineWidth',1.5);
        end

        % Individual points
        jitter = 0.08;
        plot_individuals(vals_WT_base,  1, COL_WT,  jitter);
        plot_individuals(vals_WT_amb,   2, COL_WT,  jitter);
        plot_individuals(vals_APP_base, 3, COL_APP, jitter);
        plot_individuals(vals_APP_amb,  4, COL_APP, jitter);

        set(gca,'XTick',x,'XTickLabel',labels,'XTickLabelRotation',30);
        ylabel('Power (dB)');
        grid on;

        % ---------- Post-hoc tests + stars ----------
        % Genotype effects at each condition (unpaired)
        p_G_base = safe_ttest2(vals_WT_base, vals_APP_base);
        p_G_amb  = safe_ttest2(vals_WT_amb,  vals_APP_amb);

        % Condition effects within each genotype (paired by mouse)
        p_C_WT  = paired_by_mouse(Tb, G_WT);
        p_C_APP = paired_by_mouse(Tb, G_APP);

        allVals = [vals_WT_base; vals_WT_amb; vals_APP_base; vals_APP_amb];
        yMin    = min(allVals);
        yMax    = max(allVals);
        if isempty(yMin) || isnan(yMin), yMin = 0; end
        if isempty(yMax) || isnan(yMax), yMax = 1; end
        yRange  = max(yMax - yMin, 1);

        % WT vs APP @baseline between x=1 and 3
        add_star([1 3], yMax + 0.05*yRange, p_G_base);
        % WT vs APP @ambtemp between x=2 and 4
        add_star([2 4], yMax + 0.15*yRange, p_G_amb);
        % baseline vs ambtemp in WT between x=1 and 2
        add_star([1 2], yMax + 0.25*yRange, p_C_WT);
        % baseline vs ambtemp in APP between x=3 and 4
        add_star([3 4], yMax + 0.35*yRange, p_C_APP);

        ylim([yMin - 0.1*yRange, yMax + 0.5*yRange]);

        fprintf('Post-hocs for %s – %s:\n', state, band);
        fprintf('  WT vs APP @baseline: p = %.3g\n', p_G_base);
        fprintf('  WT vs APP @ambtemp : p = %.3g\n', p_G_amb);
        fprintf('  baseline vs ambtemp in WT : p = %.3g\n', p_C_WT);
        fprintf('  baseline vs ambtemp in APP: p = %.3g\n', p_C_APP);
    end

    sgtitle(sprintf('%s band power – baseline vs ambtemp (WT vs APP)', state));

    [outDir,~,~] = fileparts(outPng);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    saveas(gcf, outPng);
    fprintf('Saved figure: %s\n', outPng);
end

% ===== helpers (same file) =====

function vals = getVals(Tb, geno, cond)
    vals = Tb.Power_dB(Tb.Genotype == geno & Tb.Condition == cond);
end

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

function p = safe_ttest2(a,b)
    a = a(~isnan(a));
    b = b(~isnan(b));
    if numel(a) >= 2 && numel(b) >= 2
        [~,p] = ttest2(a,b);
    else
        p = NaN;
    end
end

function p = paired_by_mouse(Tb, geno)
    % paired baseline vs ambtemp comparison within a genotype
    Tg   = Tb(Tb.Genotype==geno, :);
    mice = categories(Tg.MouseID);
    diffVals = [];
    for i = 1:numel(mice)
        mID  = mice{i};
        Tb_m = Tg(Tg.MouseID==mID, :);
        v_base = Tb_m.Power_dB(Tb_m.Condition=="baseline");
        v_amb  = Tb_m.Power_dB(Tb_m.Condition=="ambtemp");
        if numel(v_base)==1 && numel(v_amb)==1
            diffVals(end+1,1) = v_amb - v_base; %#ok<AGROW>
        end
    end
    diffVals = diffVals(~isnan(diffVals));
    if numel(diffVals) >= 2
        [~,p] = ttest(diffVals);  % test mean change vs 0
    else
        p = NaN;
    end
end

function add_star(xpair, y, p)
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
