function STATS = plot_baseline_APP_vs_WT_REM_stats(ALL_SUMMARY, out_dir, varargin)
% plot_baseline_APP_vs_WT_REM_stats
% -------------------------------------------------------------------------
% Focused comparison of APP vs WT in *baseline* for REM-related metrics.
%
% - Filters ALL_SUMMARY to a single condition (default: "baseline")
% - Optionally removes outliers per group (WT / APP) using 1.5*IQR rule
% - Runs t-tests (or falls back to ranksum if needed)
% - Creates bar + scatter plots with significance stars
% - Saves one PNG per metric + a CSV with stats
%
% Usage:
%   STATS = plot_baseline_APP_vs_WT_REM_stats(ALL_SUMMARY, out_dir);
%
%   STATS = plot_baseline_APP_vs_WT_REM_stats(ALL_SUMMARY, out_dir, ...
%                 'baseline_cond', 'baseline', ...
%                 'remove_outliers', true);
%

    p = inputParser;
    addParameter(p,'baseline_cond',"baseline",@(x)ischar(x)||isstring(x));
    addParameter(p,'remove_outliers',true,@islogical);
    addParameter(p,'alpha',0.05,@(x)isscalar(x)&&x>0&&x<1);
    parse(p,varargin{:});
        % --- Color definitions (same as survival plot) ---
    WT_COLOR  = [0.5 0.5 0.5];           % grey
    APP_COLOR = [100 149 237] / 255;     % cornflower blue

    baseline_cond   = string(p.Results.baseline_cond);
    remove_outliers = p.Results.remove_outliers;
    alpha           = p.Results.alpha;

    if nargin < 2 || isempty(out_dir)
        out_dir = pwd;
    end

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    % Ensure string arrays
    cond = string(ALL_SUMMARY.cond);
    geno = string(ALL_SUMMARY.geno);

    mask_baseline = (cond == baseline_cond);
    if ~any(mask_baseline)
        error('No rows in ALL_SUMMARY with cond == "%s".', baseline_cond);
    end

    SUM = ALL_SUMMARY(mask_baseline, :);

    % Define metrics you care about for baseline comparison
    metrics = {
        'total_REM_dur_s',          'Total REM time (s)';
        'mean_REM_bout_len_s',      'Mean REM bout length (s)';
        'REM_frag_bouts_per_min_REM','REM fragmentation (bouts / min REM)';
        'propREM_long',             'Proportion of long REM bouts';
        'mean_n_MA_pre_normal',     'Mean # MAs before non-long REM';
        'mean_n_MA_pre_long',       'Mean # MAs before long REM';
        'prop_cluster_shortonly',   'Prop. REM clusters: short only';
        'prop_cluster_shortlong',   'Prop. REM clusters: short→long';
        'prop_cluster_longonly',    'Prop. REM clusters: long only';
        };

    rows = [];

    for iM = 1:size(metrics,1)
        varName = metrics{iM,1};
        yLabel  = metrics{iM,2};

        if ~ismember(varName, SUM.Properties.VariableNames)
            fprintf('Skipping %s (not present in ALL_SUMMARY).\n', varName);
            continue;
        end

        % --- Extract data for WT and APP ---
        y  = SUM.(varName);
        g  = geno(mask_baseline);

        yWT  = y(g == "WT");
        yAPP = y(g == "APP");

        % Drop NaNs
        yWT  = yWT(~isnan(yWT));
        yAPP = yAPP(~isnan(yAPP));

        if numel(yWT) == 0 || numel(yAPP) == 0
            fprintf('Skipping %s (no data for WT or APP in baseline).\n', varName);
            continue;
        end

        % --- Outlier removal (per group) ---
        nWT_orig  = numel(yWT);
        nAPP_orig = numel(yAPP);

        if remove_outliers
            [yWT, idxOutWT]   = remove_outliers_IQR(yWT);
            [yAPP, idxOutAPP] = remove_outliers_IQR(yAPP);

            if ~isempty(idxOutWT)
                fprintf('Removed %d WT outlier(s) for %s.\n', numel(idxOutWT), varName);
            end
            if ~isempty(idxOutAPP)
                fprintf('Removed %d APP outlier(s) for %s.\n', numel(idxOutAPP), varName);
            end
        else
            idxOutWT  = [];
            idxOutAPP = [];
        end

        nWT  = numel(yWT);
        nAPP = numel(yAPP);

        if nWT < 2 || nAPP < 2
            fprintf('Not enough data for robust stats in %s (WT n=%d, APP n=%d).\n', ...
                    varName, nWT, nAPP);
            p_val  = NaN;
            p_test = "none";
        else
            % Try t-test; if it fails, fall back to ranksum
            try
                [~, p_val] = ttest2(yWT, yAPP, 'Vartype','unequal');
                p_test = "ttest2";
            catch
                p_val  = ranksum(yWT, yAPP);
                p_test = "ranksum";
            end
        end

                meanWT  = mean(yWT,  'omitnan');
        meanAPP = mean(yAPP, 'omitnan');
        semWT   = std(yWT,  'omitnan') / sqrt(max(nWT,1));
        semAPP  = std(yAPP, 'omitnan') / sqrt(max(nAPP,1));

        % --- Plot ---
        fig = figure('Color','w','Units','normalized','Position',[0.25 0.25 0.4 0.5]);
        hold on;

        barData = [meanWT, meanAPP];
        bh = bar(1:2, barData);
        bh.FaceColor = 'flat';              % so we can color bars individually
        bh.CData(1,:) = WT_COLOR;           % bar 1 = WT
        bh.CData(2,:) = APP_COLOR;          % bar 2 = APP
        set(bh,'EdgeColor','k','LineWidth',1.2);

        % Error bars
        eb = errorbar(1:2, barData, [semWT, semAPP], '.');
        set(eb,'LineWidth',1.2);

        % Overlay individual points (jitter)
        jitter = 0.08;
        % WT points
        xWT = ones(size(yWT)) .* (1 - jitter + 2*jitter*rand(size(yWT)));
        scatter(xWT, yWT, 40, ...
            'MarkerFaceColor', WT_COLOR, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.8);

        % APP points
        xAPP = ones(size(yAPP)) .* (2 - jitter + 2*jitter*rand(size(yAPP)));
        scatter(xAPP, yAPP, 40, ...
            'MarkerFaceColor', APP_COLOR, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.8);

        % Axes labels & title
        xticks(1:2);
        xticklabels({'WT','APP'});
        ylabel(yLabel);
        title(sprintf('%s – baseline (%s)', varName, baseline_cond), 'Interpreter','none');

        set(gca,'FontSize',12,'LineWidth',1.2);
        grid on;

        % --- Determine a safe y-max from the actual data ---
        y_all = [yWT; yAPP];
        y_max = max(y_all, [], 'omitnan');

        if ~isfinite(y_max) || isempty(y_all)
            y_max = 1;          % fallback if everything is NaN/empty
        elseif y_max <= 0
            y_max = 1;          % avoid [0 0] limits
        end

        % --- Significance star ---
        if ~isnan(p_val)
            [starStr, starColor] = p_to_stars(p_val, alpha);
        else
            starStr   = 'n.s.';
            starColor = [0 0 0];
        end

        if ~isempty(starStr)
            line([1 2], [1 1]*y_max*1.10, 'Color','k','LineWidth',1.2);
            text(1.5, y_max*1.15, starStr, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'Color',starColor);
        end

        ylim([0, y_max*1.25]);


        % Save figure
        fname_fig = fullfile(out_dir, sprintf('baseline_%s_%s.png', ...
                               baseline_cond, varName));
        exportgraphics(fig, fname_fig, 'Resolution', 300);
        close(fig);

        % --- Collect stats row ---
        r = struct;
        r.metric            = string(varName);
        r.baseline_cond     = baseline_cond;
        r.n_WT              = nWT;
        r.n_APP             = nAPP;
        r.n_WT_orig         = nWT_orig;
        r.n_APP_orig        = nAPP_orig;
        r.mean_WT           = meanWT;
        r.mean_APP          = meanAPP;
        r.sem_WT            = semWT;
        r.sem_APP           = semAPP;
        r.p_value           = p_val;
        r.test_used         = p_test;
        r.alpha             = alpha;
        r.sig_star          = string(starStr);

        rows = [rows; r]; %#ok<AGROW>
    end

    if isempty(rows)
        STATS = table();
        warning('No metrics plotted / no stats computed.');
        return;
    end

    STATS = struct2table(rows);

    % Save stats CSV
    fname_csv = fullfile(out_dir, sprintf('baseline_APP_vs_WT_stats_%s.csv', baseline_cond));
    writetable(STATS, fname_csv);

    fprintf('Baseline APP vs WT stats written to %s\n', fname_csv);
end


% -------------------------------------------------------------------------
function [x_clean, idxOut] = remove_outliers_IQR(x)
% Remove outliers using 1.5*IQR rule.

    x = x(:);
    idxOut = [];

    if numel(x) < 4
        x_clean = x;
        return;
    end

    q1 = quantile(x, 0.25);
    q3 = quantile(x, 0.75);
    IQR = q3 - q1;
    lower = q1 - 1.5*IQR;
    upper = q3 + 1.5*IQR;

    idxOut = find(x < lower | x > upper);
    x_clean = x;
    x_clean(idxOut) = [];

end


% -------------------------------------------------------------------------
function [starStr, color] = p_to_stars(p, alpha)
% Map p-values to "*", "**", "***" and color.

    if isnan(p)
        starStr = '';
        color   = [0 0 0];
        return;
    end

    if p < 0.001
        starStr = '***';
        color   = [0.2 0.2 0.2];
    elseif p < 0.01
        starStr = '**';
        color   = [0.2 0.2 0.2];
    elseif p < alpha
        starStr = '*';
        color   = [0.2 0.2 0.2];
    else
        % non-significant – you can choose to label 'n.s.' or keep empty
        starStr = 'n.s.';   % or '' if you prefer no label
        color   = [0 0 0];
    end
end
