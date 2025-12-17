function OUT = make_baseline_state_percent_APPvsWT_prism(rows_overall, out_dir, split_by_sex, varargin)
% make_baseline_state_percent_APPvsWT_prism
% ---------------------------------------------------------------
% Baseline-only plot of:
%   Y: % of total recording time spent in each state
%   X: states (WK, MA, NREM, REM, SLEEP)
%
% SLEEP is a composite state:
%   SLEEP_dur = NREM + REM   (per mouse)   [MA treated as wake]
%   %SLEEP    = 100 * SLEEP_dur / (WK+MA+NREM+REM)
%
% Two modes:
%   split_by_sex = false (default):
%       - Bars: WT vs APP (2 bars per state)
%       - Stats: unpaired t-test + ranksum, BH-FDR corrected p, stars on plot
%
%   split_by_sex = true:
%       - Bars: WT boys, WT girls, APP boys, APP girls (up to 4 bars per state)
%       - Boys = odd-numbered mouse IDs, Girls = even-numbered
%       - Currently NO stats (descriptive only).
%
% Optional name–value pairs:
%   'showIDs'             : true/false, show mouse IDs next to dots (default: true)
%   'useFDRforStars'      : true/false, use FDR p-values for stars (default: true)
%   'minNperGroupForStats': minimum n per genotype to draw a star (default: 3)
%
% OUTPUT
%   OUT.success             : logical
%   OUT.out_file            : PNG path
%   OUT.state_order         : states in x-axis order
%   OUT.split_by_sex        : logical
%   OUT.stats               : table with means, variability and p-values (genotype-only mode)
%
%   Prism-friendly outputs:
%   OUT.prism_long_table    : table with columns [MouseID, Genotype, Sex, State, PctTime]
%   OUT.prism_long_csv      : path to CSV with same content
%   OUT.prism_summary_table : per-state summary table (means/SEM etc.)
%   OUT.prism_summary_csv   : path to summary CSV
% ---------------------------------------------------------------

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end
if nargin < 3 || isempty(split_by_sex)
    split_by_sex = false;
end

% ---------- optional args ----------
p = inputParser;
addParameter(p,'showIDs',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'useFDRforStars',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minNperGroupForStats',3,@(x)isscalar(x)&&x>=1);
parse(p, varargin{:});

showIDs              = p.Results.showIDs;
useFDRforStars       = p.Results.useFDRforStars;
minNperGroupForStats = p.Results.minNperGroupForStats;

T = rows_overall;

% ---- 1) Keep only BASELINE recordings ----
cond_str    = lower(strtrim(T.condition));
is_baseline = cond_str == "baseline";
T = T(is_baseline, :);

if isempty(T)
    warning('No baseline rows found in rows_overall (condition=="baseline"). Nothing to plot.');
    OUT = struct('success', false, 'msg', 'no baseline data');
    return;
end

% ---- 2) "Base" states we use as primitives for SLEEP ----
base_desired = ["WK","MA","NREM","REM"];
avail_states = unique(T.state);
base_states  = base_desired(ismember(base_desired, avail_states));

if isempty(base_states)
    warning('No WK/MA/NREM/REM states found in baseline rows.');
    OUT = struct('success', false, 'msg', 'no states found');
    return;
end

% ---- 3) Per-mouse, per-state durations for the base states ----
Gbase = groupsummary(T, {'mouse','genotype','state'}, 'sum', 'total_dur_s');
% last column is sum_total_dur_s -> rename:
Gbase.Properties.VariableNames{end} = 'dur_s';
Gbase = Gbase(ismember(Gbase.state, base_states), :);

% ---- 4) Total recording time per mouse across all base states ----
Gtot = groupsummary(Gbase, {'mouse','genotype'}, 'sum', 'dur_s');
% last column is sum_dur_s -> rename:
Gtot.Properties.VariableNames{end} = 'total_dur_allstates';

% ---- 5) Build composite SLEEP state = NREM + REM per mouse ----
% MA is treated as wake (NOT included in SLEEP)
sleep_components = ["NREM","REM"];
Gsleep_src = Gbase(ismember(Gbase.state, sleep_components), :);

if ~isempty(Gsleep_src)
    Gsleep = groupsummary(Gsleep_src, {'mouse','genotype'}, 'sum', 'dur_s');
    Gsleep.Properties.VariableNames{end} = 'dur_s';   % sum_dur_s -> dur_s
    Gsleep.state = repmat("SLEEP", height(Gsleep), 1);
else
    Gsleep = Gbase([],:); % empty with same vars, if no NREM/REM
end

% ---- 6) Combine base states + SLEEP into one table ----
G = [Gbase; Gsleep];

% ---- 7) Attach total_dur_allstates and compute % time ----
G = innerjoin(G, Gtot(:, {'mouse','genotype','total_dur_allstates'}), ...
              'Keys', {'mouse','genotype'});

G.pct = 100 * (G.dur_s ./ G.total_dur_allstates);

% ---- 8) Decide final x-axis order including SLEEP ----
desired_order = ["WK","MA","NREM","REM","SLEEP"];
avail_states2 = unique(G.state);
state_order   = desired_order(ismember(desired_order, avail_states2));

if isempty(state_order)
    warning('No usable states after adding SLEEP.');
    OUT = struct('success', false, 'msg', 'no states after SLEEP build');
    return;
end

% ---- 9) Infer sex from mouse ID: boys=odd, girls=even ----
sex = strings(height(G),1);
for i = 1:height(G)
    m = char(G.mouse(i));
    d = regexp(m, '\d+', 'match', 'once');   % extract first number
    n = str2double(d);
    if isnan(n)
        sex(i) = "NA";
    elseif mod(n,2) == 1
        sex(i) = "boy";
    else
        sex(i) = "girl";
    end
end
G.sex = sex;

% =========================================================
%  9b) Prism export: long table (one row per mouse × state)
% =========================================================
prism_long_tbl = table;
prism_long_tbl.MouseID  = G.mouse;
prism_long_tbl.Genotype = G.genotype;
prism_long_tbl.Sex      = G.sex;
prism_long_tbl.State    = G.state;
prism_long_tbl.PctTime  = G.pct;

% sort for nicer appearance in Prism
prism_long_tbl = sortrows(prism_long_tbl, {'Genotype','Sex','State','MouseID'});

prism_long_csv = fullfile(out_dir, 'baseline_state_percent_APPvsWT_PrismLong.csv');
try
    writetable(prism_long_tbl, prism_long_csv);
catch ME
    warning('make_baseline_state_percent_APPvsWT_prism:PrismLongWriteFailed', ...
            'Could not write Prism long CSV (%s): %s', prism_long_csv, ME.message);
    prism_long_csv = '';
end

% will fill this only in the relevant branch below
prism_summary_tbl = table();
prism_summary_csv = '';

% ---- Colors ----
COL_WT       = [0.6 0.6 0.6];
COL_APP      = [0.39 0.58 0.93];
COL_WT_boy   = [0.25 0.25 0.25];
COL_WT_girl  = [0.75 0.75 0.75];
COL_APP_boy  = [0.05 0.25 0.65];
COL_APP_girl = [0.65 0.75 0.95];

% ---- 10) Plot depending on split_by_sex ----
figure('Color','w'); hold on;
x = 1:numel(state_order);
jitterFrac = 0.25;
y_offset   = 0.8;    % vertical offset for text labels

if ~split_by_sex
    % ==========================
    %  MODE 1: genotype only
    % ==========================
    hasWT  = any(G.genotype == "WT");
    hasAPP = any(G.genotype == "APP");

    if ~hasWT && ~hasAPP
        warning('No WT or APP genotypes found in baseline data.');
        OUT = struct('success', false, 'msg', 'no WT/APP found');
        return;
    end

    nStates = numel(state_order);
    meanWT  = nan(1, nStates); semWT  = nan(1, nStates);
    meanAPP = nan(1, nStates); semAPP = nan(1, nStates);

    % variability + stats arrays
    sdWT    = nan(1, nStates);
    sdAPP   = nan(1, nStates);
    cvWT    = nan(1, nStates);
    cvAPP   = nan(1, nStates);
    p_var   = nan(1, nStates);   % F-test p-value for variance equality

    p_t      = nan(1, nStates);
    p_rs     = nan(1, nStates);
    cohen_d  = nan(1, nStates);
    nWT_vec  = nan(1, nStates);
    nAPP_vec = nan(1, nStates);
    max_y    = nan(1, nStates);   % for placing stars

    for si = 1:nStates
        st = state_order(si);

        valsWT  = [];
        valsAPP = [];

        if hasWT
            maskWT = (G.genotype == "WT") & (G.state == st);
            valsWT = G.pct(maskWT);
            valsWT = valsWT(~isnan(valsWT));
            if ~isempty(valsWT)
                meanWT(si) = mean(valsWT);
                sdWT(si)   = std(valsWT);
                semWT(si)  = sdWT(si) / sqrt(numel(valsWT));
                cvWT(si)   = sdWT(si) / meanWT(si);
            end
        end

        if hasAPP
            maskAPP = (G.genotype == "APP") & (G.state == st);
            valsAPP = G.pct(maskAPP);
            valsAPP = valsAPP(~isnan(valsAPP));
            if ~isempty(valsAPP)
                meanAPP(si) = mean(valsAPP);
                sdAPP(si)   = std(valsAPP);
                semAPP(si)  = sdAPP(si) / sqrt(numel(valsAPP));
                cvAPP(si)   = sdAPP(si) / meanAPP(si);
            end
        end

        % --- stats for this state (if we have both groups) ---
        if ~isempty(valsWT) && ~isempty(valsAPP)
            nWT_vec(si)  = numel(valsWT);
            nAPP_vec(si) = numel(valsAPP);

            % unpaired t-test (Welch)
            [~, p_t(si)] = ttest2(valsWT, valsAPP, 'Vartype','unequal');

            % Mann–Whitney
            p_rs(si) = ranksum(valsWT, valsAPP);

            % Cohen's d (APP - WT)
            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            n1 = numel(valsWT); n2 = numel(valsAPP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
            cohen_d(si) = (m2 - m1) / sp;

            % F-test for equal variances (exploratory)
            if n1 >= 3 && n2 >= 3
                try
                    [~, pF] = vartest2(valsWT, valsAPP, 'Tail','both');
                    p_var(si) = pF;
                catch
                    p_var(si) = NaN;
                end
            end

            % for star placement
            max_y(si) = max([valsWT; valsAPP]);
        end
    end

    % ---- Multiple-comparison correction over STATES (BH-FDR on p_t) ----
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

    % ---- Bars ----
    barWidth = 0.4;
    if hasWT
        bar(x - barWidth/2, meanWT, barWidth, ...
            'FaceColor', COL_WT, 'EdgeColor', 'none');
    end
    if hasAPP
        bar(x + barWidth/2, meanAPP, barWidth, ...
            'FaceColor', COL_APP, 'EdgeColor', 'none');
    end

    % ---- Error bars ----
    if hasWT
        errorbar(x - barWidth/2, meanWT, semWT, 'k', ...
                 'LineStyle','none','LineWidth',1);
    end
    if hasAPP
        errorbar(x + barWidth/2, meanAPP, semAPP, 'k', ...
                 'LineStyle','none','LineWidth',1);
    end

    % ---- Dots + optional labels ----
    for si = 1:nStates
        st = state_order(si);

        if hasWT
            maskWT  = (G.genotype == "WT") & (G.state == st);
            valsWT  = G.pct(maskWT);
            mWT     = G.mouse(maskWT);
            xw      = (x(si) - barWidth/2) + (rand(size(valsWT)) - 0.5) * barWidth * jitterFrac;

            plot(xw, valsWT, '.', 'Color', [0.3 0.3 0.3]);
            if showIDs
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

        if hasAPP
            maskAPP = (G.genotype == "APP") & (G.state == st);
            valsAPP = G.pct(maskAPP);
            mAPP    = G.mouse(maskAPP);
            xa      = (x(si) + barWidth/2) + (rand(size(valsAPP)) - 0.5) * barWidth * jitterFrac;

            plot(xa, valsAPP, '.', 'Color', [0.1 0.2 0.6]);
            if showIDs
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

    % ---- Significance bars + stars ----
    if any(~isnan(p_t) | ~isnan(p_t_fdr))
        global_max_y = max(max_y(~isnan(max_y)));
        if isempty(global_max_y) || isnan(global_max_y)
            global_max_y = max(G.pct,[],'omitnan');
        end
        y_top = global_max_y + 10;  % leave some headroom for stars
        ylim([0, y_top]);

        for si = 1:nStates
            % choose which p to use for stars
            if useFDRforStars
                p_here = p_t_fdr(si);
            else
                p_here = p_t(si);
            end

            if isnan(p_here) || p_here >= 0.05
                continue;
            end
            if nWT_vec(si) < minNperGroupForStats || nAPP_vec(si) < minNperGroupForStats
                continue;
            end

            % decide number of stars
            if p_here < 0.001
                stars = '***';
            elseif p_here < 0.01
                stars = '**';
            else
                stars = '*';
            end

            % horizontal line between WT and APP bars
            y_star = max_y(si) + 4;   % a bit above highest dot for that state
            if isnan(y_star)
                y_star = global_max_y * 0.9;
            end
            line([x(si)-barWidth/2, x(si)+barWidth/2], [y_star, y_star], ...
                 'Color','k','LineWidth',1.0);
            text(x(si), y_star + 1.5, stars, ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',14, ...
                 'FontWeight','bold');
        end
    end

    % ---- Print & store stats table ----
    stats_idx = ~isnan(p_t);  % states where we had both groups
    if any(stats_idx)
        sd_ratio = sdAPP ./ sdWT;

        stats_tbl = table( ...
            state_order(stats_idx)', ...
            nWT_vec(stats_idx)', ...
            nAPP_vec(stats_idx)', ...
            meanWT(stats_idx)', ...
            meanAPP(stats_idx)', ...
            sdWT(stats_idx)', ...
            sdAPP(stats_idx)', ...
            cvWT(stats_idx)', ...
            cvAPP(stats_idx)', ...
            sd_ratio(stats_idx)', ...
            p_t(stats_idx)', ...
            p_t_fdr(stats_idx)', ...
            p_rs(stats_idx)', ...
            p_var(stats_idx)', ...
            cohen_d(stats_idx)', ...
            'VariableNames', {'State','nWT','nAPP', ...
                              'MeanWT','MeanAPP', ...
                              'SD_WT','SD_APP', ...
                              'CV_WT','CV_APP', ...
                              'SD_ratio_APP_vs_WT', ...
                              'p_ttest','p_ttest_FDR', ...
                              'p_ranksum','p_var_Ftest','Cohen_d'});
        fprintf('\nBaseline %%time in state: WT vs APP stats (per state)\n');
        disp(stats_tbl);

        OUT.stats = stats_tbl;
    else
        fprintf('\n[Stats] No states with both WT and APP present. No tests performed.\n');
    end

    % =====================================================
    %  Prism summary table (per-state means/SEM etc.)
    % =====================================================
    prism_summary_tbl = table( ...
        state_order(:), ...
        meanWT(:), semWT(:), nWT_vec(:), ...
        meanAPP(:), semAPP(:), nAPP_vec(:), ...
        'VariableNames', {'State', ...
                          'Mean_WT','SEM_WT','nWT', ...
                          'Mean_APP','SEM_APP','nAPP'});
    prism_summary_csv = fullfile(out_dir, 'baseline_state_percent_APPvsWT_PrismSummary.csv');
    try
        writetable(prism_summary_tbl, prism_summary_csv);
    catch ME
        warning('make_baseline_state_percent_APPvsWT_prism:PrismSummaryWriteFailed', ...
                'Could not write Prism summary CSV (%s): %s', prism_summary_csv, ME.message);
        prism_summary_csv = '';
    end

    ttl     = 'Baseline: % time in each state (WT vs APP, incl. SLEEP)';
    out_png = 'baseline_state_percent_APPvsWT.png';

else
    % ==========================
    %  MODE 2: genotype + sex
    % ==========================
    group_defs = struct( ...
        'name',  {'WT_boys','WT_girls','APP_boys','APP_girls'}, ...
        'gen',   {"WT","WT","APP","APP"}, ...
        'sex',   {"boy","girl","boy","girl"}, ...
        'label', {'WT boys','WT girls','APP boys','APP girls'}, ...
        'color', {COL_WT_boy, COL_WT_girl, COL_APP_boy, COL_APP_girl} );

    nGroups = numel(group_defs);
    nStates = numel(state_order);
    meanMat = nan(nGroups, nStates);
    semMat  = nan(nGroups, nStates);
    hasGroup = false(1, nGroups);

    % compute means/SEMs
    for g = 1:nGroups
        for si = 1:nStates
            st = state_order(si);
            mask = (G.genotype == group_defs(g).gen) & ...
                   (G.sex      == group_defs(g).sex) & ...
                   (G.state    == st);
            vals = G.pct(mask);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                meanMat(g,si) = mean(vals);
                semMat(g,si)  = std(vals) / sqrt(numel(vals));
                hasGroup(g)   = true;
            end
        end
    end

    barWidth = 0.18;
    offsets  = ((1:nGroups) - (nGroups+1)/2) * barWidth;
    hBar     = gobjects(1, nGroups);

    % Bars + error bars
    for g = 1:nGroups
        if ~hasGroup(g), continue; end
        this_x = x + offsets(g);
        hBar(g) = bar(this_x, meanMat(g,:), barWidth, ...
                      'FaceColor', group_defs(g).color, ...
                      'EdgeColor', 'none');
        errorbar(this_x, meanMat(g,:), semMat(g,:), 'k', ...
                 'LineStyle','none','LineWidth',1);
    end

    % Dots + optional labels
    max_y = nan(1, nStates);
    for si = 1:nStates
        st = state_order(si);

        for g = 1:nGroups
            if ~hasGroup(g), continue; end

            mask = (G.genotype == group_defs(g).gen) & ...
                   (G.sex      == group_defs(g).sex) & ...
                   (G.state    == st);
            vals = G.pct(mask);
            mID  = G.mouse(mask);

            if isempty(vals), continue; end

            base_x = x(si) + offsets(g);
            xd = base_x + (rand(size(vals)) - 0.5) * barWidth * jitterFrac;

            plot(xd, vals, '.', 'Color', group_defs(g).color);

            if showIDs
                for j = 1:numel(vals)
                    thisID = char(mID(j));
                    text(xd(j), vals(j) + y_offset, thisID, ...
                         'Rotation', 45, ...
                         'HorizontalAlignment','left', ...
                         'VerticalAlignment','bottom', ...
                         'FontSize',8, ...
                         'Color', group_defs(g).color);
                end
            end

            max_y(si) = max([max_y(si), max(vals)]);
        end
    end

    % Legend only for groups that exist
    legLabels = {group_defs(hasGroup).label};
    legHandles = hBar(hasGroup);
    if ~isempty(legLabels)
        legend(legHandles, legLabels, ...
               'Location','northoutside', ...
               'Orientation','horizontal');
    end

    % increase y-limit a bit for readability
    if any(~isnan(max_y))
        ylim([0, max(max_y,[],'omitnan') + 10]);
    end

    % =====================================================
    %  Prism summary table for genotype+sex mode
    %    (rows: State × Group, with Mean/SEM)
    % =====================================================
    rows_state  = strings(0,1);
    rows_group  = strings(0,1);
    rows_gen    = strings(0,1);
    rows_sex    = strings(0,1);
    rows_mean   = [];
    rows_sem    = [];

    for g = 1:nGroups
        if ~hasGroup(g), continue; end
        for si = 1:nStates
            rows_state(end+1,1) = state_order(si); %#ok<AGROW>
            rows_group(end+1,1) = group_defs(g).label; %#ok<AGROW>
            rows_gen(end+1,1)   = group_defs(g).gen;   %#ok<AGROW>
            rows_sex(end+1,1)   = group_defs(g).sex;   %#ok<AGROW>
            rows_mean(end+1,1)  = meanMat(g,si);       %#ok<AGROW>
            rows_sem(end+1,1)   = semMat(g,si);        %#ok<AGROW>
        end
    end

    if ~isempty(rows_state)
        prism_summary_tbl = table( ...
            rows_state, rows_group, rows_gen, rows_sex, ...
            rows_mean, rows_sem, ...
            'VariableNames', {'State','GroupLabel','Genotype','Sex', ...
                              'Mean_Pct','SEM_Pct'});
        prism_summary_csv = fullfile(out_dir, 'baseline_state_percent_APPvsWT_bySex_PrismSummary.csv');
        try
            writetable(prism_summary_tbl, prism_summary_csv);
        catch ME
            warning('make_baseline_state_percent_APPvsWT_prism:PrismSummaryBySexWriteFailed', ...
                    'Could not write Prism by-sex summary CSV (%s): %s', ...
                    prism_summary_csv, ME.message);
            prism_summary_csv = '';
        end
    end

    ttl     = 'Baseline: % time in each state (WT vs APP, by sex, incl. SLEEP)';
    out_png = 'baseline_state_percent_APPvsWT_bySex.png';
end

% ---- 11) Common axes formatting, labels, save ----
xticks(x);
xticklabels(state_order);
xlabel('State');
ylabel('% time in state');
title(ttl);

set(gca,'Box','off','FontSize',12);

out_file = fullfile(out_dir, out_png);
saveas(gcf, out_file);

% finalize OUT
if ~exist('OUT','var') || ~isstruct(OUT)
    OUT = struct();
end
OUT.success           = true;
OUT.out_file          = out_file;
OUT.state_order       = state_order;
OUT.split_by_sex      = logical(split_by_sex);

% Prism outputs
OUT.prism_long_table    = prism_long_tbl;
OUT.prism_long_csv      = prism_long_csv;
OUT.prism_summary_table = prism_summary_tbl;
OUT.prism_summary_csv   = prism_summary_csv;

fprintf('✅ Baseline %%time state plot saved to: %s\n', out_file);
if ~isempty(prism_long_csv)
    fprintf('📄 Prism LONG CSV (mouse × state) saved to: %s\n', prism_long_csv);
end
if ~isempty(prism_summary_csv)
    fprintf('📄 Prism SUMMARY CSV saved to: %s\n', prism_summary_csv);
end
end
