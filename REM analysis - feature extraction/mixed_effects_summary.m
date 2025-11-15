function OUT = mixed_effects_summary(M, varnames, varargin)
% If fitlme is installed AND licensed, use it; otherwise do nonparametrics.
% Usage: OUT = mixed_effects_summary(M, vars)  % quiet by default

p = inputParser;
addParameter(p,'verbose',false,@islogical);
parse(p,varargin{:});
verbose = p.Results.verbose;

OUT = struct();

% Check availability AND license for Statistics Toolbox
hasFun = exist('fitlme','file')==2;
hasLic = false;
try
    hasLic = license('test','Statistics_Toolbox');
catch, hasLic = false;
end
useLME = hasFun && hasLic;

conds = categorical(M.condition);
groups = categories(conds);

for v = varnames(:)'
    vname = v{1};
    if ~ismember(vname, M.Properties.VariableNames), continue; end
    y = M.(vname);

    % Group summaries (median/mean/IQR)
    S = struct();
    for i = 1:numel(groups)
        mask = conds == groups{i};
        S.(groups{i}).n      = sum(mask);
        S.(groups{i}).median = median(y(mask),'omitnan');
        S.(groups{i}).mean   = mean(y(mask),'omitnan');
        S.(groups{i}).iqr    = iqr(y(mask));
    end

    if useLME
        try
            T = table(y, categorical(M.condition), categorical(M.mouse), ...
                      'VariableNames',{'y','cond','mouse'});
            lme = fitlme(T,'y ~ cond + (1|mouse)');
            CI  = coefCI(lme);
            tbl = dataset2table(lme.Coefficients); %#ok<DSET2TBL>
            OUT.(vname) = struct('mode','lme','summaries',S,'lme',lme,'coeff',tbl,'CI',CI);
            continue
        catch
            if verbose
                warning('fitlme failed for %s, falling back to nonparametric.', vname);
            end
        end
    end

    % --- Nonparametric fallback (quiet) ---
    if numel(groups)==2
        a = y(conds==groups{1}); b = y(conds==groups{2});
        P = safe_ranksum(a,b);
        eff = hedges_g(a,b);
        OUT.(vname) = struct('mode','nonparametric','summaries',S, ...
                             'test','Wilcoxon rank-sum','p',P,'effect',eff);
    else
        P = safe_kruskal(y, conds);
        OUT.(vname) = struct('mode','nonparametric','summaries',S, ...
                             'test','Kruskal–Wallis','p',P,'effect',NaN);
    end
end
end

% ----- helpers -----
function P = safe_ranksum(a,b)
P = NaN; a=a(isfinite(a)); b=b(isfinite(b));
try, P = ranksum(a,b); end
end

function P = safe_kruskal(y,g)
P = NaN; y=y(isfinite(y)); g=g(isfinite(y));
try, P = kruskalwallis(y,g,'off'); end
end

function g = hedges_g(a,b)
a=a(isfinite(a)); b=b(isfinite(b));
na=numel(a); nb=numel(b);
if na<2 || nb<2, g = NaN; return; end
sa = var(a,1); sb = var(b,1);
sp = sqrt(((na-1)*sa + (nb-1)*sb) / max(na+nb-2,1));
if sp==0 || ~isfinite(sp), g = NaN; return; end
d = (mean(a)-mean(b)) / sp;
J = 1 - 3/(4*(na+nb)-9);
g = J*d;
end
