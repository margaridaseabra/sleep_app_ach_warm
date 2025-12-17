function OUT = EEG_bandpower_baseline_APPvsWT(T, ...
    mouseVar, genoVar, condVar, stateVar, bandVars, out_dir)

subdir = fullfile(out_dir, 'baseline_APPvsWT');
if ~isfolder(subdir), mkdir(subdir); end

TT = T;

cond = lower(strtrim(string(TT.(condVar))));
geno = string(TT.(genoVar));
state= string(TT.(stateVar));
mouse= string(TT.(mouseVar));

TT.(condVar)  = cond;
TT.(genoVar)  = geno;
TT.(stateVar) = state;
TT.(mouseVar) = mouse;

% Baseline only
TT = TT(cond=="baseline", :);
if isempty(TT)
    warning('No baseline rows in EEG bandpower table.');
    OUT = struct('success',false); return;
end

% Group per mouse × geno × state (mean bandpower per mouse, per state)
groupVars = {mouseVar, genoVar, stateVar};
statVars  = bandVars;
G = groupsummary(TT, groupVars, 'mean', statVars);

% Clean variable names: mean_band -> band
for k = 1:numel(statVars)
    old = ['mean_' statVars{k}];
    if ismember(old, G.Properties.VariableNames)
        G.Properties.VariableNames{old} = statVars{k};
    end
end

states = unique(G.(stateVar));
OUT = struct();
OUT.per_mouse = G;
OUT.stats = struct();

for s = 1:numel(states)
    st = states(s);
    Gi = G(G.(stateVar)==st, :);
    if isempty(Gi), continue; end
    
    for b = 1:numel(bandVars)
        var = bandVars{b};
        if ~ismember(var, Gi.Properties.VariableNames)
            continue;
        end
        
        Y = Gi.(var);
        ok = ~isnan(Y);
        Gi2 = Gi(ok,:);
        Y   = Y(ok);
        if height(Gi2) < 3, continue; end
        
        % LME: band ~ geno + (1|mouse)
        Gi2.(mouseVar) = categorical(Gi2.(mouseVar));
        Gi2.(genoVar)  = categorical(Gi2.(genoVar));
        
        try
            formula = sprintf('%s ~ %s + (1|%s)', var, genoVar, mouseVar);
            lme = fitlme(Gi2, formula, 'DummyVarCoding','effects');
            a   = anova(lme);
            rowG= strcmp(a.Term, genoVar);
            pG  = a.pValue(rowG);
        catch ME
            warning('LME failed (%s, %s): %s', st, var, ME.message);
            pG = NaN;
        end
        
        % Unpaired WT vs APP + Cohen d
        WT  = Gi2.(var)(Gi2.(genoVar)=="WT");
        APP = Gi2.(var)(Gi2.(genoVar)=="APP");
        WT  = WT(~isnan(WT)); APP = APP(~isnan(APP));
        p_t = NaN; d = NaN;
        if numel(WT)>=3 && numel(APP)>=3
            [~, p_t] = ttest2(WT, APP, 'Vartype','unequal');
            m1 = mean(WT); m2 = mean(APP);
            s1 = std(WT);  s2 = std(APP);
            n1 = numel(WT); n2 = numel(APP);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2)/max(1,(n1+n2-2)));
            d  = (m2 - m1)/sp;
        end
        
        % Plot
        figure('Color','w','Position',[200 200 450 350]); hold on;
        COL_WT  = [0.6 0.6 0.6];
        COL_APP = [0.39 0.58 0.93];
        
        means = [mean(WT,'omitnan'), mean(APP,'omitnan')];
        sems  = [std(WT,'omitnan')/sqrt(max(1,numel(WT))), ...
                 std(APP,'omitnan')/sqrt(max(1,numel(APP)))];
        
        bar(1, means(1), 0.6, 'FaceColor',COL_WT,'EdgeColor','none');
        bar(2, means(2), 0.6, 'FaceColor',COL_APP,'EdgeColor','none');
        errorbar(1:2, means, sems, 'k','LineStyle','none');
        
        xlim([0.5 2.5]);
        set(gca,'XTick',1:2,'XTickLabel',{'WT','APP'});
        ylabel(sprintf('%s (baseline, %s)', var, st), 'Interpreter','none');
        
        yTop = max(means + sems) * 1.3;
        if isempty(yTop) || isnan(yTop), yTop = 1; end
        ylim([0 yTop]);
        
        txt = sprintf('LME geno p=%.3f | t-test p=%.3f | d=%.2f', pG, p_t, d);
        text(1.5, yTop*0.9, txt, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top','FontSize',9);
        
        title(sprintf('Baseline %s – %s', st, var), 'Interpreter','none');
        set(gca,'Box','off','FontSize',11);
        
        fname = fullfile(subdir, sprintf('EEGbase_%s_%s_APPvsWT.png', ...
                        char(st), var));
        saveas(gcf, fname);
        
        key = sprintf('%s_%s', st, var);
        OUT.stats.(matlab.lang.makeValidName(key)) = struct( ...
            'state',st, 'band',var, ...
            'p_geno_LME',pG, 'p_ttest',p_t, 'd',d, ...
            'nWT',numel(WT), 'nAPP',numel(APP), ...
            'fig_file',fname);
    end
end

OUT.success = true;
save(fullfile(subdir,'baseline_APPvsWT_stats.mat'),'OUT');
end
