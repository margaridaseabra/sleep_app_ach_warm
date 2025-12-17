function OUT = plot_REM_3h_baseline_vs_ambtemp_APPvsWT(ALL_SUMMARY_3h, out_dir, varargin)
% plot_REM_3h_baseline_vs_ambtemp_APPvsWT
% -------------------------------------------------------------------------
% For the 3 h CROPPED windows (baseline + ambtemp), compare REM metrics:
%
%   - total_REM_dur_s      (shown in minutes)
%   - mean_REM_bout_len_s
%   - REM_frag_bouts_per_min_REM
%   - propREM_long
%
% For each metric, makes one bar plot:
%   4 bars: WT-baseline, WT-ambtemp, APP-baseline, APP-ambtemp
%
% Stats:
%   - Within-genotype paired t-test (baseline vs ambtemp), WT & APP
%       • BH-FDR across metrics per genotype
%       • stars between the 2 bars (per genotype)
%   - Mixed-effects model per metric:
%       metric ~ cond * geno + (1|mouse)
%       => p(cond), p(geno), p(cond×geno) printed on the plot
%
% -------------------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});

minN  = p.Results.minNperGroupForStats;
useFDR = p.Results.useFDR;

T = ALL_SUMMARY_3h;
T.cond = lower(strtrim(string(T.cond)));
T.geno = string(T.geno);
T.mouse = string(T.mouse);

% keep only baseline + ambtemp
T = T(T.cond=="baseline" | T.cond=="ambtemp", :);
if isempty(T)
    warning('No baseline + ambtemp rows in ALL_SUMMARY_3h.');
    OUT = struct('success',false); return;
end

% restrict to mice that have both conditions (for clean paired tests)
G = groupsummary(T, {'mouse','geno'});
G.Properties.VariableNames(end) = {'nCond'};

haveBoth = G.nCond >= 2;
pair_keys = G.mouse(haveBoth);
T = T(ismember(T.mouse, pair_keys), :);

if height(T) < 2
    warning('Too few recordings with both baseline & ambtemp.');
    OUT = struct('success',false); return;
end

% Metrics to plot: {field, y_label, scale_factor, nice_name}
metrics = {
    'total_REM_dur_s',          'Total REM time (min, 3 h window)', 1/60,  'totalREMmin';
    'mean_REM_bout_len_s',      'Mean REM bout length (s)',         1,      'meanREMlen';
    'REM_frag_bouts_per_min_REM','REM fragmentation (bouts/min REM)',1,     'REMfrag';
    'propREM_long',             'Prop. of REM bouts that are long', 1,      'propLong';
    };

genotypes = ["WT","APP"];

% For FDR across metrics:
p_within = nan(numel(genotypes), size(metrics,1));

OUT = struct();
OUT.metric_stats = struct();

for m = 1:size(metrics,1)
    field    = metrics{m,1};
    y_label  = metrics{m,2};
    scale    = metrics{m,3};
    short_id = metrics{m,4};

    if ~ismember(field, T.Properties.VariableNames)
        warning('Metric %s missing in ALL_SUMMARY_3h. Skipping.', field);
        continue;
    end

    Ti = T;
    Ti.Y = double(Ti.(field)) * scale;
    Ti = Ti(~isnan(Ti.Y), :);
    if isempty(Ti), continue; end

    % wide per mouse × geno: baseline / ambtemp
    Twide = unstack(Ti(:,{'mouse','geno','cond','Y'}), 'Y','cond');
    if ~all(ismember({'baseline','ambtemp'}, Twide.Properties.VariableNames))
        warning('Metric %s: not all mice have both baseline & ambtemp.', field);
    end

    % containers for bar means / SEMs
    meanMat = nan(4,1);  % WT_b, WT_a, APP_b, APP_a
    semMat  = nan(4,1);
    paired_stats = {};   % for output table

    % within-genotype paired tests
    for gi = 1:numel(genotypes)
        gtype = genotypes(gi);
        Tg = Twide(Twide.geno==gtype, :);
        if isempty(Tg), continue; end

        if ~all(ismember({'baseline','ambtemp'}, Tg.Properties.VariableNames))
            continue;
        end

        base = Tg.baseline;
        amb  = Tg.ambtemp;
        ok   = ~isnan(base) & ~isnan(amb);
        base = base(ok); amb = amb(ok);
        nPairs = numel(base);
        if nPairs < minN, continue; end

        [~, p_pair] = ttest(base, amb);

        m_base = mean(base,'omitnan');
        m_amb  = mean(amb,'omitnan');
        se_base = std(base,'omitnan')/sqrt(nPairs);
        se_amb  = std(amb,'omitnan')/sqrt(nPairs);
        delta   = amb - base;

        p_within(gi,m) = p_pair;

        paired_stats = [paired_stats; ...
            {char(gtype), nPairs, m_base, m_amb, mean(delta,'omitnan'), p_pair}]; %#ok<AGROW>

        if gtype=="WT"
            meanMat(1) = m_base;  semMat(1) = se_base;
            meanMat(2) = m_amb;   semMat(2) = se_amb;
        else
            meanMat(3) = m_base;  semMat(3) = se_base;
            meanMat(4) = m_amb;   semMat(4) = se_amb;
        end
    end

    % Mixed-effects ANOVA: Y ~ cond * geno + (1|mouse)
    Ti2 = Ti;
    Ti2.cond  = categorical(Ti2.cond, {'baseline','ambtemp'});
    Ti2.geno  = categorical(Ti2.geno);
    Ti2.mouse = categorical(Ti2.mouse);

    Fg=NaN; pg=NaN; Fc=NaN; pc=NaN; Fi=NaN; pi=NaN;

    try
        lme = fitlme(Ti2, 'Y ~ cond * geno + (1|mouse)', ...
                     'DummyVarCoding','effects');
        a = anova(lme);


        terms = string(a.Term);

        rowCond = find(terms=="cond", 1);
        rowGen  = find(terms=="geno", 1);
        rowInt  = find(contains(terms,"cond") & contains(terms,"geno"), 1);


        if ~isempty(rowCond), Fc = a.FStat(rowCond); pc = a.pValue(rowCond); end
        if ~isempty(rowGen),  Fg = a.FStat(rowGen);  pg = a.pValue(rowGen);  end
        if ~isempty(rowInt),  Fi = a.FStat(rowInt);  pi = a.pValue(rowInt);  end
    catch ME
        warning('Mixed-effects for %s failed: %s', field, ME.message);
    end

    % ---------- FDR across metrics (filled later) ----------
    % (we just collected p_within; now correct)
    % We'll handle after loop; for now use raw p's for plotting.

    % ---------- build stats table for this metric ----------
    if isempty(paired_stats)
        within_tbl = table();
    else
        within_tbl = cell2table(paired_stats, ...
            'VariableNames', {'Genotype','n_pairs', ...
                              'Mean_baseline','Mean_ambtemp', ...
                              'Mean_delta','p_paired'});
    end

    anova_tbl = table(Fg,pg,Fc,pc,Fi,pi, ...
        'VariableNames', {'F_genotype','p_genotype', ...
                          'F_condition','p_condition', ...
                          'F_interaction','p_interaction'});

    OUT.metric_stats.(short_id).within = within_tbl;
    OUT.metric_stats.(short_id).anova  = anova_tbl;

    % ---------- plotting ----------
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.39 0.58 0.93];

    figure('Color','w','Position',[200 200 600 400]); hold on;
    x = 1:4;
    labels = {'WT base','WT amb','APP base','APP amb'};
    barWidth = 0.6;

    for k = 1:4
        if isnan(meanMat(k)), continue; end
        if k<=2
            col = COL_WT;
        else
            col = COL_APP;
        end
        if mod(k,2)==1
            % baseline: solid
            bar(x(k), meanMat(k), barWidth, 'FaceColor',col,'EdgeColor','none');
        else
            % ambtemp: outline
            bar(x(k), meanMat(k), barWidth, 'FaceColor','none', ...
                'EdgeColor',col,'LineWidth',1.3);
        end
        if ~isnan(semMat(k))
            errorbar(x(k), meanMat(k), semMat(k), 'k','LineStyle','none');
        end
    end

    all_means = meanMat(~isnan(meanMat));
    if isempty(all_means)
        ymax = 1;
    else
        ymax = max(all_means) + 0.3*max(all_means);
    end
    ylim([0 ymax]);

    % (temporarily use raw p's; we'll update FDR text after loop)
    % stars: we need raw p's now
    for gi = 1:numel(genotypes)
        gtype = genotypes(gi);
        if isempty(within_tbl), continue; end
        r = strcmp(within_tbl.Genotype, gtype);
        if ~any(r), continue; end
        p_raw = within_tbl.p_paired(r);
        if isnan(p_raw) || p_raw >= 0.05, continue; end

        if p_raw < 0.001, stars='***';
        elseif p_raw < 0.01, stars='**';
        else, stars='*';
        end

        if gtype=="WT"
            x1 = x(1); x2 = x(2);
            y_star = ymax*0.9;
        else
            x1 = x(3); x2 = x(4);
            y_star = ymax*0.78;
        end

        line([x1 x2],[y_star y_star],'Color','k','LineWidth',1);
        text(mean([x1 x2]), y_star+0.02*ymax, stars, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
    end

    txt = sprintf('Cond p=%.3f | Geno p=%.3f | C×G p=%.3f', pc, pg, pi);
    text(2.5, ymax*0.65, txt, ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontSize',9);

    set(gca,'XTick',x,'XTickLabel',labels,'FontSize',11);
    ylabel(y_label);
    title(sprintf('3 h REM metric: %s', field), 'Interpreter','none');
    set(gca,'Box','off');

    fig_file = fullfile(out_dir, sprintf('REM3h_%s_baseline_vs_ambtemp_APPvsWT.png', short_id));
    saveas(gcf, fig_file);

    OUT.metric_stats.(short_id).fig_file = fig_file;

    % save stats tables
    if ~isempty(within_tbl)
        writetable(within_tbl, ...
            fullfile(out_dir, sprintf('REM3h_%s_within_stats.csv', short_id)));
    end
    writetable(anova_tbl, ...
        fullfile(out_dir, sprintf('REM3h_%s_anova_stats.csv', short_id)));
end

% ---------- FDR across metrics for within-genotype tests ----------
if useFDR
    nM = size(metrics,1);
    for gi = 1:numel(genotypes)
        pv = p_within(gi,:);
        valid = ~isnan(pv);
        if ~any(valid), continue; end
        pvals = pv(valid);
        [sp,idx] = sort(pvals(:));
        m = numel(sp);
        adj = sp .* (m./(1:m))';
        for k = m-1:-1:1
            adj(k) = min(adj(k), adj(k+1));
        end
        adj(adj>1) = 1;
        pF = nan(size(pv));
        pF(valid) = adj;
        p_within(gi,:) = pF;
    end

    % store in OUT for later inspection
    OUT.p_within_FDR = p_within;
end

OUT.success = true;
fprintf('✅ REM 3 h bar plots + stats saved in %s\n', out_dir);
end
