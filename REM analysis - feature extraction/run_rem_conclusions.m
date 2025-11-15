function run_rem_conclusions(ROOT)
if nargin<1, ROOT = pwd; end

% 1) Gather all bout-level rows
ALL = collect_rem_features(ROOT);

% 2) Per-mouse medians (robust)
M = aggregate_per_mouse(ALL);

% ---- helper to pick the first existing column for a band, regardless of suffix
pickBand = @(tbl, base) pick_band_var(tbl, base);

% Build shortlist of endpoints using real column names that exist
want = {'dur_s','ibi_prev_s','theta_peak_Hz','theta_peak_sd','theta_bw_Hz', ...
        'emg_burst_pct','pac_theta_hgamma_MI'};

theta_col  = pickBand(M,'theta');
hgamma_col = pickBand(M,'hgamma');

if ~isempty(theta_col),  want = [want, {theta_col}];  end
if ~isempty(hgamma_col), want = [want, {hgamma_col}]; end

% Filter to columns that actually exist in M
vars = want( ismember(want, M.Properties.VariableNames) );

% 3) Mixed-effects (or nonparametric fallback)
OUT = mixed_effects_summary(M, vars);

% 4) Survival / hazard
KM = km_survival_plot_rem(ALL);

% 5) Save tidy CSVs for figures/tables
outd = fullfile(ROOT,'_rem_conclusions'); if ~exist(outd,'dir'), mkdir(outd); end
writetable(ALL, fullfile(outd,'ALL_bout_level.csv'));
writetable(M,   fullfile(outd,'PER_MOUSE_median.csv'));
save(fullfile(outd,'mixed_models.mat'),'OUT');
save(fullfile(outd,'km_survival.mat'),'KM');

fprintf('[OK] REM conclusions pack saved under %s\n', outd);
end

% -------- helper (suffix-agnostic band column picker) ----------
function name = pick_band_var(tbl, base)
cands = { base, [base '_db'], [base '_z'], [base '_zrem'], [base '_none'], [base '_raw'] };
name = '';
for i = 1:numel(cands)
    if any(strcmpi(tbl.Properties.VariableNames, cands{i}))
        name = cands{i};
        return
    end
end
% as a final fallback, take the first column that starts with "base_"
vn = string(tbl.Properties.VariableNames);
hit = startsWith(lower(vn), lower(base+"_"));
if any(hit), name = char(vn(find(hit,1,'first'))); end
end
