function group_plot_wake_onsets(GROUP)

conds = unique(string({GROUP.sessions.cond}),'stable');
COL_WT  = [0.2 0.2 0.2];
COL_APP = [1.0 0.5 0];

figure('Name','Group Wake onsets','Color','w');

for c = 1:numel(conds)
    cond = conds(c);

    WT_traces = [];
    APP_traces = [];
    t_ref = [];
    WT_peaks  = []; APP_peaks  = [];
    WT_slopes = []; APP_slopes = [];

    for k = 1:numel(GROUP.sessions)
        if string(GROUP.sessions(k).cond) ~= cond; continue; end
        OUT = GROUP.out{k};
        idxW = find(strcmp({OUT.transitions.name},'Wake_onset'));
        if isempty(idxW); continue; end
        Tw = OUT.transitions(idxW);

        if isempty(t_ref); t_ref = Tw.t_rel(:); end

        % use per-session mean trace so each session contributes equally
        if strcmpi(GROUP.sessions(k).geno,'WT')
            WT_traces = [WT_traces Tw.mean(:)];
            WT_peaks  = [WT_peaks  mean(Tw.peaks,'omitnan')];
            WT_slopes = [WT_slopes mean(Tw.slopes,'omitnan')];
        else
            APP_traces = [APP_traces Tw.mean(:)];
            APP_peaks  = [APP_peaks  mean(Tw.peaks,'omitnan')];
            APP_slopes = [APP_slopes mean(Tw.slopes,'omitnan')];
        end
    end

    r = c;

    % (A) traces
    subplot(numel(conds),2,2*r-1); hold on;
    if ~isempty(WT_traces)
        m = mean(WT_traces,2); se = std(WT_traces,0,2)/sqrt(size(WT_traces,2));
        fill_between(t_ref, m-se, m+se, COL_WT, 0.2); 
        plot(t_ref, m, 'Color',COL_WT,'LineWidth',1.8);
    end
    if ~isempty(APP_traces)
        m = mean(APP_traces,2); se = std(APP_traces,0,2)/sqrt(size(APP_traces,2));
        fill_between(t_ref, m-se, m+se, COL_APP, 0.2); 
        plot(t_ref, m, 'Color',COL_APP,'LineWidth',1.8);
    end
    yline(0,'k:'); xline(0,'k--');
    xlabel('Time from wake onset (s)');
    ylabel('dF/F (baseline norm.)');
    title(sprintf('Wake onsets – %s',cond));
    legend({'WT','APP/PS1'},'Location','northwest');
    box off;

    % (B) bar plot for peak OR slope – here I show peak, slope stacked
    subplot(numel(conds),2,2*r); hold on;

    % jittered scatter of per-session peak means
    scatter(0.8*ones(size(WT_peaks)), WT_peaks, 30, COL_WT,'filled');
    scatter(1.2*ones(size(APP_peaks)), APP_peaks, 30, COL_APP,'filled');

    % group means + error
    muP = [mean(WT_peaks,'omitnan') mean(APP_peaks,'omitnan')];
    seP = [std(WT_peaks,'omitnan')/sqrt(max(1,numel(WT_peaks))) ...
           std(APP_peaks,'omitnan')/sqrt(max(1,numel(APP_peaks)))];
    errorbar([1 2], muP, seP,'k','LineStyle','none','CapSize',8);

    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','APP'});
    ylabel('Peak dF/F');
    title(sprintf('Wake onset peaks – %s',cond));
    box off;

    % if you'd prefer slope bars instead, just change to WT_slopes/APP_slopes
end

sgtitle('ACh at wake onsets – WT vs APP/PS1 across conditions');

end
