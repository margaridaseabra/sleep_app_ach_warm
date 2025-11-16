function group_plot_psd_NREM(GROUP)

conds = unique(string({GROUP.sessions.cond}),'stable');   % ['baseline','ambtemp','drug']
genos = {'WT','APPPS1'};

COL_WT   = [0.2 0.2 0.2];
COL_APP  = [1.0 0.5 0];   % orange-ish

figure('Name','Group NREM ACh PSD','Color','w');

nC = numel(conds);
for c = 1:nC
    cond = conds(c);

    % collect PSD curves per genotype
    PSD_WT  = [];
    PSD_APP = [];
    f_ref   = [];

    pow_WT = []; pow_APP = [];
    pf_WT  = []; pf_APP  = [];
    pa_WT  = []; pa_APP  = [];

    for k = 1:numel(GROUP.sessions)
        if string(GROUP.sessions(k).cond) ~= cond
            continue;
        end
        s = GROUP.sessions(k);
        OUT = GROUP.out{k};
        if isempty(OUT.psd.NREM.f); continue; end

        f  = OUT.psd.NREM.f;
        ps = OUT.psd.NREM.psd;
        ps = ps / max(ps);   % normalize like your single-session plot

        if isempty(f_ref); f_ref = f; end
        % (assume all same length / frequencies)

        if strcmpi(s.geno,'WT')
            PSD_WT  = [PSD_WT  ps(:)];
            pow_WT  = [pow_WT  OUT.psd.NREM.band_power];
            pf_WT   = [pf_WT   OUT.psd.NREM.peak_freq];
            pa_WT   = [pa_WT   OUT.psd.NREM.peak_amp];
        else
            PSD_APP = [PSD_APP ps(:)];
            pow_APP = [pow_APP OUT.psd.NREM.band_power];
            pf_APP  = [pf_APP  OUT.psd.NREM.peak_freq];
            pa_APP  = [pa_APP  OUT.psd.NREM.peak_amp];
        end
    end

    % --- row index for this condition ---
    r = c;

    % (1) PSD curves
    subplot(nC,3,3*(r-1)+1);
    hold on;
    if ~isempty(PSD_WT)
        m = mean(PSD_WT,2); se = std(PSD_WT,0,2)/sqrt(size(PSD_WT,2));
        fill_between(f_ref, m-se, m+se, COL_WT, 0.2);
        plot(f_ref, m,'Color',COL_WT,'LineWidth',1.8);
    end
    if ~isempty(PSD_APP)
        m = mean(PSD_APP,2); se = std(PSD_APP,0,2)/sqrt(size(PSD_APP,2));
        fill_between(f_ref, m-se, m+se, COL_APP, 0.2);
        plot(f_ref, m,'Color',COL_APP,'LineWidth',1.8);
    end
    xlabel('Frequency (Hz)');
    ylabel('ACh power (A.U.)');
    title(sprintf('NREM PSD – %s',cond));
    legend({'WT','APP/PS1'},'Location','northeast');
    box off;

    % (2) band power bar
    subplot(nC,3,3*(r-1)+2);
    hold on;
    bar(1, mean(pow_WT,'omitnan'),'FaceColor',COL_WT);
    bar(2, mean(pow_APP,'omitnan'),'FaceColor',COL_APP);
    errorbar([1 2], ...
        [mean(pow_WT,'omitnan') mean(pow_APP,'omitnan')], ...
        [std(pow_WT,'omitnan')/sqrt(max(1,numel(pow_WT))) ...
         std(pow_APP,'omitnan')/sqrt(max(1,numel(pow_APP)))], ...
        'k','LineStyle','none','CapSize',8);
    xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'WT','APP'});
    ylabel('NREM power');
    title('ACh PSD power');
    box off;

    % (3) peak freq & amp combined as dots
    subplot(nC,3,3*(r-1)+3); hold on;
    scatter(ones(size(pf_WT))*0.9, pf_WT, 30, COL_WT,'filled');
    scatter(ones(size(pf_APP))*1.1, pf_APP, 30, COL_APP,'filled');
    scatter(ones(size(pa_WT))*1.9, pa_WT, 30, COL_WT,'o');
    scatter(ones(size(pa_APP))*2.1, pa_APP, 30, COL_APP,'o');
    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Peak Hz','Peak amp'});
    ylabel('Value');
    title('Peak metrics (dots = mice)');
    box off;
end

sgtitle('NREM ACh PSD – WT vs APP/PS1 across conditions');
end

% simple shaded area helper
function fill_between(x,y1,y2,col,alphaVal)
fill([x; flipud(x)], [y1; flipud(y2)], col, ...
     'EdgeColor','none','FaceAlpha',alphaVal);
end
