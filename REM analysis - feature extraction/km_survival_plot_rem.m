function KM = km_survival_plot_rem(ALL, outdir)
% Pure-MATLAB Kaplan–Meier survival for REM bout durations by condition.
% No Statistics Toolbox required.

if nargin<2, outdir = fullfile(pwd,'_rem_survival'); end
if ~exist(outdir,'dir'), mkdir(outdir); end

conds = unique(string(ALL.condition));
KM = struct();

figure('Color','w'); hold on;
for i = 1:numel(conds)
    c = conds(i);
    d = double(ALL.dur_s(string(ALL.condition)==c));
    d = d(isfinite(d) & d>0);
    if isempty(d), continue; end

    % --- KM estimate (no censoring) ---
    d = sort(d(:));
    [t,~,ic] = unique(d);           % event times
    ni = accumarray(ic, 1);         % events at time t_i
    % risk set just before t_i:
    % r_i = number still at risk = total - number with time < t_i
    cum_events_before = [0; cumsum(ni(1:end-1))];
    r = numel(d) - cum_events_before;
    S = cumprod(1 - ni./r);         % KM survivor

    KM.(c).t = t; KM.(c).S = S;

    % Plot as stairs
    stairs([0; t],[1; S],'LineWidth',1.5,'DisplayName',char(c));
end
xlabel('REM bout duration (s)'); ylabel('Survival S(t)');
title('REM bout survival by condition'); grid on; legend('Location','best');
saveas(gcf, fullfile(outdir,'km_survival.png'));
end
