function OUT = make_baseline_total_bouts_per_hour_APPvsWT(rows_perhr, out_dir, states_to_plot)
% make_baseline_bouts_per_hour_by_state_APPvsWT_withIDs
% -------------------------------------------------------------------------
% For BASELINE recordings only:
%   For each requested STATE (WK, MA, NREM, REM), make a figure:
%
%       X: hour_idx (0, 1, 2, ...)
%       Y: bouts_per_h (for that state)
%
%   - Two curves: WT vs APP (mean ± SEM across mice)
%   - All individual mice shown as jittered dots with mouse ID labels.
%
% Statistics (per state):
%   - For each hour: WT vs APP bouts_per_h
%       * unpaired t-test (Welch)
%       * Mann-Whitney (ranksum)
%       * Cohen's d (APP - WT)
%   - BH-FDR correction across hours (based on t-test p-values)
%   - Stars plotted at hours where FDR-corrected p < 0.05.
%
% Input:
%   rows_perhr  : table from run_group_sleep_architecture (group_per_hour)
%   out_dir     : where to save figures
%   states_to_plot (optional) : string/cell array, default ["WK","MA","NREM","REM"]
%
% Output:
%   OUT : struct with info about states plotted, output files, and stats.
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

if nargin < 3 || isempty(states_to_plot)
    states_to_plot = ["WK","MA","NREM","REM"];
end

% normalize to string array
states_to_plot = string(states_to_plot(:)).';

P = rows_perhr;

% ---- 1) BASELINE only ----
cond_str    = lower(strtrim(P.condition));
is_baseline = cond_str == "baseline";
P = P(is_baseline, :);

if isempty(P)
    warning('No baseline rows found in rows_perhr (condition=="baseline"). Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline data');
    return;
end

OUT = struct();
OUT.success = true;
OUT.states_plotted = [];
OUT.files  = struct();
OUT.stats  = struct();   % will hold per-state stats tables

COL_WT  = [0.6 0.6 0.6];       % grey
COL_APP = [0.39 0.58 0.93];    % cornflower-ish blue

for s = 1:numel(states_to_plot)
    st = states_to_plot(s);

    % ---- 2) Filter to this state ----
    Pst = P(P.state == st, :);

    if isempty(Pst)
        warning('No baseline rows for state "%s". Skipping.', st);
        continue;
    end

    % ---- 3) We have per mouse/genotype/hour already: bouts_per_h ----
    G = Pst;

    hasWT  = any(G.genotype == "WT");
    hasAPP = any(G.genotype == "APP");

    if ~hasWT && ~hasAPP
        warning('No WT or APP data for state "%s". Skipping.', st);
        continue;
    end

    % All hours present for this state
    all_hours = unique(G.hour_idx);
    all_hours = sort(all_hours);

    nH = numel(all_hours);
    meanWT  = nan(1, nH); semWT  = nan(1, nH);
    meanAPP = nan(1, nH); semAPP = nan(1, nH);

    % stats arrays per hour
    p_t      = nan(1, nH);
    p_rs     = nan(1, nH);
    cohen_d  = nan(1, nH);
    nWT_vec  = nan(1, nH);
    nAPP_vec = nan(1, nH);
    max_y    = nan(1, nH);  % max value at that hour (for star position)

    for h = 1:nH
        hr = all_hours(h);

        valsWT  = [];
        valsAPP = [];

        if hasWT
            maskWT = (G.genotype == "WT") & (G.hour_idx == hr);
            valsWT = G.bouts_per_h(maskWT);
            valsWT = valsWT(~isnan(valsWT));
            if ~isempty(valsWT)
                meanWT(h) = mean(valsWT);
                semWT(h)  = std(valsWT) / sqrt(numel(valsWT));
            end
        end

        if hasAPP
            maskAPP = (G.genotype == "APP") & (G.hour_idx == hr);
            valsAPP = G.bouts_per_h(maskAPP);
            valsAPP = valsAPP(~isnan(valsAPP));
            if ~isempty(valsAPP)
                meanAPP(h) = mean(valsAPP);
                semAPP(h)  = std(valsAPP) / sqrt(numel(valsAPP));
            end
        end

        % --- Stats for this hour (if both groups present) ---
        if ~isempty(valsWT) && ~isempty(valsAPP)
            nWT_vec(h)  = numel(valsWT);
            nAPP_vec(h) = numel(valsAPP);

            % unpaired t-test (Welch)
            [~, p_t(h)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');

            % Mann–Whitney
            p_rs(h) = ranksum(valsWT, valsAPP);

            % Cohen's d (APP - WT)
            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
            cohen_d(h) = (m2 - m1) / sp;

            max_y(h) = max([valsWT; valsAPP]);
        elseif ~isempty(valsWT)
            max_y(h) = max(valsWT);
        elseif ~isempty(valsAPP)
            max_y(h) = max(valsAPP);
        end
    end

    % ---- 3b) Multiple-comparison correction over HOURS (BH-FDR on p_t) ----
    p_t_fdr = nan(size(p_t));
    valid   = ~isnan(p_t);
    pvals   = p_t(valid);
    if ~isempty(pvals)
        [sorted_p, sort_idx] = sort(pvals(:));
        m = numel(sorted_p);
        adj = sorted_p .* (m ./ (1:m)');        % BH adjustment
        % enforce monotonicity
        for i = m-1:-1:1
            adj(i) = min(adj(i), adj(i+1));
        end
        adj(adj>1) = 1;
        p_fdr_vals = nan(size(pvals));
        p_fdr_vals(sort_idx) = adj;
        p_t_fdr(valid) = p_fdr_vals;
    end

    % ---- 4) Plot for this state ----
    figure('Color','w'); hold on;
    use = (nWT_vec >= 3) & (nAPP_vec >= 3);
    all_hours = all_hours(use);
    meanWT    = meanWT(use);
    
    show_ids = false;  % <--- add at top as an option

   
    % Mean ± SEM lines
    if hasWT
        errorbar(all_hours, meanWT, semWT, '-o', ...
                 'Color', COL_WT, ...
                 'MarkerFaceColor', COL_WT, ...
                 'MarkerSize', 5, ...
                 'LineWidth', 1.2);
    end

    if hasAPP
        errorbar(all_hours, meanAPP, semAPP, '-o', ...
                 'Color', COL_APP, ...
                 'MarkerFaceColor', COL_APP, ...
                 'MarkerSize', 5, ...
                 'LineWidth', 1.2);
    end

    % ---- 5) Overlay per-mouse dots + IDs (jittered) ----
    jitterFrac = 0.25;
    y_offset   = 0.5;  % vertical offset for text above dot

    for h = 1:nH
        hr = all_hours(h);

        % WT
        if hasWT
            maskWT  = (G.genotype == "WT") & (G.hour_idx == hr);
            valsWT  = G.bouts_per_h(maskWT);
            mWT     = G.mouse(maskWT);

            if ~isempty(valsWT)
                xw = hr - 0.1 + (rand(size(valsWT)) - 0.5) * jitterFrac;
                plot(xw, valsWT, '.', 'Color', [0.3 0.3 0.3]);

                for j = 1:numel(valsWT)
                    thisID = char(mWT(j));
                    text(xw(j), valsWT(j) + y_offset, thisID, ...
                         'Rotation', 45, ...
                         'HorizontalAlignment','left', ...
                         'VerticalAlignment','bottom', ...
                         'FontSize',8, ...
                         'Color',[0.2 0.2 0.2]);


                end
            end
        end

        % APP
        if hasAPP
            maskAPP = (G.genotype == "APP") & (G.hour_idx == hr);
            valsAPP = G.bouts_per_h(maskAPP);
            mAPP    = G.mouse(maskAPP);

            if ~isempty(valsAPP)
                xa = hr + 0.1 + (rand(size(valsAPP)) - 0.5) * jitterFrac;
                plot(xa, valsAPP, '.', 'Color', [0.1 0.2 0.6]);

                for j = 1:numel(valsAPP)
                    thisID = char(mAPP(j));
                    text(xa(j), valsAPP(j) + y_offset, thisID, ...
                         'Rotation', 45, ...
                         'HorizontalAlignment','left', ...
                         'VerticalAlignment','bottom', ...
                         'FontSize',8, ...
                         'Color',[0.1 0.2 0.6]);
                end
            end
        end
    end

    xlabel('Hour from recording start');
    ylabel(sprintf('Bouts per hour (%s)', st));
    title(sprintf('Baseline: bouts per hour in %s (WT vs APP)', st));

    if hasWT && hasAPP
        legend({'WT mean±SEM','APP mean±SEM'}, 'Location','best');
    elseif hasWT
        legend({'WT mean±SEM'}, 'Location','best');
    elseif hasAPP
        legend({'APP mean±SEM'}, 'Location','best');
    end

    set(gca,'Box','off','FontSize',12);
    xlim([min(all_hours)-0.5, max(all_hours)+0.5]);

    % ---- 6) Add significance stars (per hour, using FDR p_t) ----
    if any(~isnan(max_y))
        global_max_y = max(max_y(~isnan(max_y)));
    else
        global_max_y = max(G.bouts_per_h, [], 'omitnan');
    end
    if isempty(global_max_y) || isnan(global_max_y)
        global_max_y = 1;
    end

    y_top  = global_max_y + 3;  % headroom
    ylim([0, y_top]);

    for h = 1:nH
        p_here = p_t_fdr(h);
        if isnan(p_here) || p_here >= 0.05
            continue;
        end

        % Decide number of stars
        if p_here < 0.001
            stars = '***';
        elseif p_here < 0.01
            stars = '**';
        else
            stars = '*';
        end

        y_star = max_y(h);
        if isnan(y_star)
            y_star = global_max_y * 0.8;
        end
        y_star = y_star + 1.0;  % a bit above points/means

        text(all_hours(h), y_star, stars, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',12, ...
             'FontWeight','bold');
    end

    % ---- 7) Save figure for this state ----
    fname = sprintf('baseline_bouts_per_hour_%s_APPvsWT_withIDs.png', lower(st));
    out_file = fullfile(out_dir, fname);
    saveas(gcf, out_file);

    OUT.states_plotted = [OUT.states_plotted, st];
    OUT.files.(matlab.lang.makeValidName(st)) = out_file;

    fprintf('✅ Baseline bouts/hour plot with IDs for %s saved to: %s\n', st, out_file);

    % ---- 8) Build and store stats table for this state ----
    stats_idx = ~isnan(p_t);  % hours with both groups present
    if any(stats_idx)
        stats_tbl = table( ...
            all_hours(stats_idx), ...
            nWT_vec(stats_idx)', ...
            nAPP_vec(stats_idx)', ...
            meanWT(stats_idx)', ...
            meanAPP(stats_idx)', ...
            p_t(stats_idx)', ...
            p_t_fdr(stats_idx)', ...
            p_rs(stats_idx)', ...
            cohen_d(stats_idx)', ...
            'VariableNames', {'Hour','nWT','nAPP','MeanWT','MeanAPP', ...
                              'p_ttest','p_ttest_FDR','p_ranksum','Cohen_d'});

        fprintf('\nBaseline bouts/hour (%s): WT vs APP stats per hour\n', st);
        disp(stats_tbl);

        OUT.stats.(matlab.lang.makeValidName(st)) = stats_tbl;
    else
        fprintf('\n[Stats %s] No hours with both WT and APP present. No tests performed.\n', st);
        OUT.stats.(matlab.lang.makeValidName(st)) = table();
    end
end
end
