function STATS = run_REM_group_stats(ALL_SUMMARY, out_dir)
% run_REM_group_stats
% -------------------------------------------------------------------------
% Simple APP vs WT stats per condition for key REM metrics.
%
% Outputs:
%   STATS : table with columns:
%           cond, metric, n_WT, n_APP, mean_WT, mean_APP, p_ttest2

if nargin < 2 || isempty(out_dir)
    out_dir = pwd;
end

if isempty(ALL_SUMMARY)
    warning('ALL_SUMMARY is empty, no stats computed.');
    STATS = table();
    return;
end

conds = unique(ALL_SUMMARY.cond(:));
genos = categorical(ALL_SUMMARY.geno);
metrics = {
    'n_REM_total'
    'propREM_long'
    'total_REM_dur_s'
    'mean_REM_bout_len_s'
    'REM_frag_bouts_per_min_REM'
    'prop_cluster_shortonly'
    'prop_cluster_shortlong'
    'prop_cluster_longonly'
    };

rows = [];

for iC = 1:numel(conds)
    cond_i = conds(iC);

    mask_cond = (ALL_SUMMARY.cond == cond_i);

    for iM = 1:numel(metrics)
        mName = metrics{iM};

        if ~ismember(mName, ALL_SUMMARY.Properties.VariableNames)
            continue; % metric not available
        end

        % Extract values per genotype
        mask_WT  = mask_cond & (genos == 'WT');
        mask_APP = mask_cond & (genos == 'APP');

        xWT  = ALL_SUMMARY.(mName)(mask_WT);
        xAPP = ALL_SUMMARY.(mName)(mask_APP);

        xWT  = xWT(~isnan(xWT));
        xAPP = xAPP(~isnan(xAPP));

        if numel(xWT) < 2 || numel(xAPP) < 2
            p = NaN;
        else
            [~, p] = ttest2(xWT, xAPP, 'Vartype','unequal');
        end

        r.cond       = cond_i;
        r.metric     = string(mName);
        r.n_WT       = numel(xWT);
        r.n_APP      = numel(xAPP);
        r.mean_WT    = mean(xWT, 'omitnan');
        r.mean_APP   = mean(xAPP, 'omitnan');
        r.p_ttest2   = p;

        rows = [rows; r]; %#ok<AGROW>
    end
end

if isempty(rows)
    STATS = table();
    warning('No stats rows created (check genotypes/conditions).');
    return;
end

STATS = struct2table(rows);

% Save to CSV
fname = fullfile(out_dir, 'REM_group_stats_APP_vs_WT.csv');
writetable(STATS, fname);

fprintf('REM stats written to %s\n', fname);
end
