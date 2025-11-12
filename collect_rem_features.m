function ALL = collect_rem_features(root)
% Recursively gather *_REM_TFR_features.csv into one table with mouse/cond parsed.
% Robust to band column suffixes: *_db, *_z, *_zrem, *_none (and plain names).

if nargin<1, root = pwd; end
files = dir(fullfile(root,'**','*_REM_TFR_features.csv'));
assert(~isempty(files),'No *_REM_TFR_features.csv found under %s',root);

rows = {};
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    T = readtable(f);

    % --- parse mouse & condition from filename ---
    [mouse, cond] = parse_mouse_cond_matlab(files(k).name);
    T.source_file = repmat(string(files(k).name), height(T),1);
    T.mouse       = repmat(string(mouse),           height(T),1);
    T.condition   = repmat(string(cond),            height(T),1);

    % --- robust getters for band columns (handle *_db, *_z, *_zrem, *_none, or plain) ---
    theta = get_band(T,'theta');
    delta = get_band(T,'delta');
    beta  = get_band(T,'beta');
    hgam  = get_band(T,'hgamma');  % high gamma
    lgam  = get_band(T,'lgamma');  % low gamma

    % Derived features (guarded: produce NaN if inputs missing)
    T.theta_to_delta = safe_ratio(theta, delta);
    T.gamma_to_beta  = safe_ratio(hgam,  beta);   % use high-gamma vs beta
    T.lgamma_to_beta = safe_ratio(lgam,  beta);

    % Simple proxy for REM density based on gap to previous bout
    if ismember('ibi_prev_s', T.Properties.VariableNames)
        T.rem_density_bout = 1 ./ max(T.ibi_prev_s, eps);
    else
        T.rem_density_bout = nan(height(T),1);
    end

    rows{end+1} = T; %#ok<AGROW>
end

ALL = vertcat(rows{:});   % long table

end

% ---------- helpers ----------
function v = get_band(T, base)
% Return numeric column for a band regardless of suffix. If none, NaN.
cands = band_name_candidates(base);
v = nan(height(T),1);
for i = 1:numel(cands)
    name = cands{i};
    hit = strcmpi(T.Properties.VariableNames, name);
    if any(hit)
        v = T{:, find(hit,1,'first')};
        return
    end
end
% Fallback: try prefix match like 'theta_' or exact 'theta'
vn = string(T.Properties.VariableNames);
hit = startsWith(lower(vn), lower(base+"_"));
if any(hit)
    v = T{:, find(hit,1,'first')};
end
end

function C = band_name_candidates(base)
% Order matters; include common spellings
C = {
    char(base), ...
    [char(base) '_db'], ...
    [char(base) '_z'], ...
    [char(base) '_zrem'], ...
    [char(base) '_none'], ...
    [char(base) '_raw'] ...
};
end

function r = safe_ratio(a,b)
try
    r = a ./ max(b, eps);
    r(~isfinite(r)) = NaN;
catch
    n = max([numel(a), numel(b), 1]);
    r = nan(n,1);
end
end
