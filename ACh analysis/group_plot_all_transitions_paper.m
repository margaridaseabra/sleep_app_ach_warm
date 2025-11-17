function group_plot_all_transitions_paper(GROUP)
% Make paper-style transition figures (Wake, NREM, REM onsets)
% Rows = conditions, columns: traces (left) + peak bar (right).

trans_list  = {'Wake_onset','NREM_onset','REM_onset'};
title_list  = {'WAKE ONSETS','NREM ONSETS','REM ONSETS'};

for i = 1:numel(trans_list)
    plot_one_transition(GROUP, trans_list{i}, title_list{i});
end
end

% ======================================================================
function plot_one_transition(GROUP, trans_name, big_title)

sessions = GROUP.sessions;
conds    = unique(string({sessions.cond}),'stable');   % baseline / ambtemp / drugs
nCond    = numel(conds);

% figure out genotype labels and order (WT first if present)
geno_vals = unique(string({sessions.geno}),'stable');
if any(geno_vals=="WT")
    genoWT  = "WT";
    genoMut = geno_vals(geno_vals~="WT");
else
    genoWT  = geno_vals(1);
    genoMut = geno_vals(2:end);
end
if isempty(genoMut)
    error('Need at least two genotypes (WT + APP).');
end
genoMut = genoMut(1);   % assume one mutant group (APP/PS1)

% colours (WT grey, APP cornflower blue)
COL_WT  = [0.6 0.6 0.6];           % grey
COL_APP = [0.392 0.584 0.929];     % cornflower blue


figure('Name',sprintf('%s – %s',big_title,trans_name),'Color','w');

for c = 1:nCond
    cond = conds(c);

    WT_traces   = [];
    APP_traces  = [];
    WT_peaks    = [];
    APP_peaks   = [];
    t_ref       = [];

    % ---------- collect per-session means for this condition ----------
    for k = 1:numel(sessions)
        if string(sessions(k).cond) ~= cond
            continue;
        end

        OUT = GROUP.out{k};

        % find this transition in OUT.transitions
        idxT = [];
        if isfield(OUT,'transitions')
            idxT = find(strcmp({OUT.transitions.name}, trans_name));
        end
        if isempty(idxT); continue; end

        Ttr = OUT.transitions(idxT);
        t_this = Ttr.t_rel(:);
        mean_trace = Ttr.mean(:);
        peak_mean  = mean(Ttr.peaks,'omitnan');

        % ---- define common time base & resample if needed ----
        if isempty(t_ref)
            t_ref = t_this;      % first session defines reference
        elseif numel(t_this) ~= numel(t_ref) || any(abs(t_this - t_ref) > 1e-6)
            % resample onto reference grid
            mean_trace = interp1(t_this, mean_trace, t_ref, 'linear', 'extrap');
        end

        g = string(sessions(k).geno);
        if g == genoWT
            WT_traces = [WT_traces mean_trace]; %#ok<AGROW>
            WT_peaks  = [WT_peaks  peak_mean];  %#ok<AGROW>
        elseif g == genoMut
            APP_traces = [APP_traces mean_trace]; %#ok<AGROW>
            APP_peaks  = [APP_peaks  peak_mean];  %#ok<AGROW>
        end
    end

    % ---------------- LEFT: traces -----------------------
    ax1 = subplot(nCond,2,2*c-1); hold(ax1,'on');
    hWTline  = [];
    hAPPline = [];

    if ~isempty(WT_traces)
        m  = mean(WT_traces,2);
        se = std(WT_traces,0,2)/sqrt(size(WT_traces,2));
        fill_between(t_ref, m-se, m+se, COL_WT, 0.2);
        hWTline = plot(t_ref, m, 'Color',COL_WT,'LineWidth',1.8);
    end

    if ~isempty(APP_traces)
        m  = mean(APP_traces,2);
        se = std(APP_traces,0,2)/sqrt(size(APP_traces,2));
        fill_between(t_ref, m-se, m+se, COL_APP, 0.2);
        hAPPline = plot(t_ref, m, 'Color',COL_APP,'LineWidth',1.8);
    end

    yline(0,'k:'); xline(0,'k--');
    xlabel('Time from onset (s)');
    ylabel('dF/F (baseline norm.)');
    title(sprintf('%s – %s', cond, strrep(trans_name,'_',' ')), ...
          'Interpreter','none');
    box off;

    if c == 1
        % legend only from the line handles (ignore fills)
        legHandles = [hWTline hAPPline];
        legLabels  = {char(genoWT), char(genoMut)};
        % remove any empties in case one group is missing for this condition
        keep = ~cellfun(@isempty, num2cell(legHandles));
        legHandles = legHandles(keep);
        legLabels  = legLabels(keep);
        if ~isempty(legHandles)
            legend(legHandles, legLabels, ...
                   'Location','northwest','Box','off');
        end
    end

    % ---------------- RIGHT: peak bar plot ----------------
    ax2 = subplot(nCond,2,2*c); hold(ax2,'on');

    xs       = [1 2];
    mu_peaks = [mean(WT_peaks,'omitnan')  mean(APP_peaks,'omitnan')];
    se_peaks = [std(WT_peaks,'omitnan')/max(1,sqrt(sum(~isnan(WT_peaks)))) ...
                std(APP_peaks,'omitnan')/max(1,sqrt(sum(~isnan(APP_peaks))))];

    if ~all(isnan(mu_peaks))
        bar(1, mu_peaks(1), 0.6, 'FaceColor',COL_WT);
        bar(2, mu_peaks(2), 0.6, 'FaceColor',COL_APP);
        errorbar(xs, mu_peaks, se_peaks, 'k','LineStyle','none','CapSize',8);
    end

    % scatter individual sessions
    jitter = 0.06;
    if ~isempty(WT_peaks)
        scatter(1 + (rand(size(WT_peaks))-0.5)*2*jitter, WT_peaks, ...
                25, COL_WT,'filled');
    end
    if ~isempty(APP_peaks)
        scatter(2 + (rand(size(APP_peaks))-0.5)*2*jitter, APP_peaks, ...
                25, COL_APP,'filled');
    end

    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{char(genoWT), char(genoMut)});
    ylabel('Peak dF/F');
    title(sprintf('%s – peak dF/F', cond), 'Interpreter','none');
    box off;
end

sgtitle(big_title);

end

% ----------------------------------------------------------------------
function fill_between(x,y1,y2,col,alphaVal)
% small helper for shaded SEM
x   = x(:);
y1  = y1(:);
y2  = y2(:);
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
     'EdgeColor','none','FaceAlpha',alphaVal);
end
