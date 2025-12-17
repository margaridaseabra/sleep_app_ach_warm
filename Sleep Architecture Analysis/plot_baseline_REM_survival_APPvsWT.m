function SURV_BASE = plot_baseline_REM_survival_APPvsWT(ALL_REM, out_dir, varargin)
% plot_baseline_REM_survival_APPvsWT
% -------------------------------------------------------------------------
% Baseline-only survival analysis of REM bout durations:
%   - filters ALL_REM to baseline condition (or user-specified)
%   - plots WT vs APP survival curves (Kaplan–Meier style)
%   - performs a log-rank test (WT vs APP)
%
% Usage:
%   SURV_BASE = plot_baseline_REM_survival_APPvsWT(ALL_REM, out_dir);
%
% Options (name-value):
%   'baseline_cond' : which condition label to use (default "baseline")
%   'min_bouts'     : minimum # of bouts per group to run stats (default 2)

    p = inputParser;
    addParameter(p,'baseline_cond',"baseline",@(x)ischar(x)||isstring(x));
    addParameter(p,'min_bouts',2,@(x)isscalar(x)&&x>=1);
    parse(p,varargin{:});

    % --- Color definitions (match other function) ---
    WT_COLOR       = [0 0 0];                    % black
    APP_COLOR      = [100 149 237] / 255;       % cornflower blue

    baseline_cond = string(p.Results.baseline_cond);
    min_bouts     = p.Results.min_bouts;

    if nargin < 2 || isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    if isempty(ALL_REM)
        warning('ALL_REM is empty, skipping baseline survival analysis.');
        SURV_BASE = table();
        return;
    end

    cond   = string(ALL_REM.cond);
    geno   = string(ALL_REM.geno);
    dur_s  = ALL_REM.dur_s;

    % ---------- 1) Filter to baseline only ----------
    mask_base = (cond == baseline_cond);
    if ~any(mask_base)
        warning('No REM bouts found for baseline condition "%s".', baseline_cond);
        SURV_BASE = table();
        return;
    end

    dur_base  = dur_s(mask_base);
    geno_base = geno(mask_base);

    % Split by genotype
    dWT  = dur_base(geno_base == "WT");
    dAPP = dur_base(geno_base == "APP");

    dWT  = dWT(~isnan(dWT) & dWT > 0);
    dAPP = dAPP(~isnan(dAPP) & dAPP > 0);

    nWT  = numel(dWT);
    nAPP = numel(dAPP);

    if nWT < min_bouts || nAPP < min_bouts
        warning('Not enough REM bouts in baseline (WT n=%d, APP n=%d).', nWT, nAPP);
        p_logrank    = NaN;
        chi2_logrank = NaN;
    else
        % ---------- 2) Log-rank test (WT vs APP) ----------
        time = [dWT; dAPP];
        grp  = [zeros(nWT,1); ones(nAPP,1)];  % 0 = WT, 1 = APP

        [p_logrank, chi2_logrank] = logrank_two_groups(time, grp);
    end

    % ---------- 3) Plot survival curves (baseline only) ----------
    figure('Color','w','Position',[200 200 700 450]); % match other function
    hold on;

    if nWT >= 1
        [tWT, SWT] = empirical_survival(dWT);
        % Use stairs + same style as other function (solid, 1.5)
        stairs(tWT, SWT, ...
            'LineWidth', 1.5, ...
            'LineStyle','-', ...
            'Color', WT_COLOR, ...
            'DisplayName','WT');
    end

    if nAPP >= 1
        [tAPP, SAPP] = empirical_survival(dAPP);
        stairs(tAPP, SAPP, ...
            'LineWidth', 1.5, ...
            'LineStyle','-', ...
            'Color', APP_COLOR, ...
            'DisplayName','APP');
    end

    % Axis labels + title style to match
    xlabel('REM bout duration (s)');
    ylabel('Survival S(t) = 1 - F(t)');
    title(sprintf('REM bout survival – baseline (%s)', baseline_cond), ...
          'Interpreter','none');

    if nWT >= 1 || nAPP >= 1
        legend('Location','southwest');
    end

    % Match axes style: no grid, box off, font size 11
    set(gca,'Box','off','FontSize',11);

    % Limits: mimic style of second function (0 to 95th percentile, 0–1)
    all_d = [dWT; dAPP];
    all_d = all_d(~isnan(all_d));
    if isempty(all_d)
        xmax = 1;
    else
        xmax = prctile(all_d, 95);
        if ~isfinite(xmax) || xmax <= 0
            xmax = 1;
        end
    end
    xlim([0, xmax]);
    ylim([0, 1]);

    % Annotate log-rank p-value (font size aligned with axes)
    if ~isnan(p_logrank)
        txt = sprintf('log-rank p = %.3g', p_logrank);
    else
        txt = 'log-rank p = n/a';
    end

    xl = xlim; yl = ylim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.90*diff(yl), ...
        txt, 'FontSize',11, 'FontWeight','bold');

    % Save using same style as the other function
    fname = fullfile(out_dir, sprintf('REM_survival_baseline_%s.png', baseline_cond));
    saveas(gcf, fname);

    % ---------- 4) Collect stats ----------
    SURV_BASE = table;
    SURV_BASE.baseline_cond = baseline_cond;
    SURV_BASE.n_WT_bouts    = nWT;
    SURV_BASE.n_APP_bouts   = nAPP;
    SURV_BASE.p_logrank     = p_logrank;
    SURV_BASE.chi2_logrank  = chi2_logrank;

    % Save to CSV for the thesis
    csv_name = fullfile(out_dir, sprintf('REM_survival_baseline_logrank_%s.csv', baseline_cond));
    writetable(SURV_BASE, csv_name);
    fprintf('Baseline survival log-rank stats written to %s\n', csv_name);
end


% =====================================================================
function [t, S] = empirical_survival(d)
% empirical_survival
% Compute Kaplan–Meier-like survival curve S(t) from durations d (>0).

d = sort(d(:));
t = unique(d);
N = numel(d);
S = zeros(size(t));

for k = 1:numel(t)
    S(k) = sum(d >= t(k)) / N;
end
end


% =====================================================================
function [p, chi2_stat] = logrank_two_groups(time, grp)
% logrank_two_groups
% ---------------------------------------------------------------------
% Log-rank test comparing two survival curves:
%   time : [N x 1] durations (positive)
%   grp  : [N x 1] group labels (0 or 1)
%
% Assumes all events are observed (no censoring), which matches your
% REM bout durations (all bouts end).

time = time(:);
grp  = grp(:);

if numel(time) ~= numel(grp)
    error('time and grp must have same length.');
end

% Sort unique event times
t_unique = unique(time);
O1 = 0;  % observed events in group 1
E1 = 0;  % expected events in group 1
V1 = 0;  % variance

for j = 1:numel(t_unique)
    tj = t_unique(j);

    % Risk set: subjects with time >= tj
    at_risk = (time >= tj);
    n1j = sum(grp(at_risk) == 1);
    n0j = sum(grp(at_risk) == 0);
    nj  = n1j + n0j;

    % Events at time tj
    events = (time == tj);
    d1j = sum(grp(events) == 1);
    d0j = sum(grp(events) == 0);
    dj  = d1j + d0j;

    if nj <= 1 || dj == 0
        continue;
    end

    % Expected events in group 1
    E1j = dj * (n1j / nj);

    % Variance contribution
    V1j = (n1j * n0j * dj * (nj - dj)) / (nj^2 * (nj - 1));

    O1 = O1 + d1j;
    E1 = E1 + E1j;
    V1 = V1 + V1j;
end

if V1 <= 0
    chi2_stat = NaN;
    p = NaN;
    return;
end

chi2_stat = (O1 - E1)^2 / V1;

if exist('chi2cdf','file')
    p = 1 - chi2cdf(chi2_stat, 1);
else
    p = NaN;
end
end
