function OUT = run_EEG_bandpower_group_stats_manual(csvFile, outDir)
% run_EEG_bandpower_group_stats_manual
% -------------------------------------------------------------------------
% Group-level EEG band power analysis (baseline or ambtemp)
% using manually specified variable names.
%
% INPUTS
%   csvFile : EEG_band_power_allmice.csv (from run_all_mice_eeg_psd_auto)
%   outDir  : output folder for plots and stats
%
% OUTPUT
%   OUT.tables  : summary tables per state and band
%   OUT.figures : paths to plots
% -------------------------------------------------------------------------

if nargin < 2 || isempty(outDir)
    outDir = fullfile(pwd, 'EEG_PSD_Group');
end
if ~isfolder(outDir), mkdir(outDir); end

% ---- load table ----
T = readtable(csvFile);
disp('Loaded EEG bandpower table:');
disp(T.Properties.VariableNames');

% ---- manually assign variable names ----
mouseVar = 'MouseID';        % if your column is called 'Mouse'
genoVar  = 'Genotype';     % if your column is 'Genotype'
condVar  = 'Condition';    % <-- THIS is the one causing the error
stateVar = 'State';        % if your column is 'State'


% ---- band power variables (edit if needed) ----
bandVars = {'delta_rel','theta_rel','sigma_rel','beta_rel','gamma_rel'}; 
% check with T.Properties.VariableNames to adjust

% ---- sanitize ----
T.(genoVar)  = string(T.(genoVar));
T.(condVar)  = lower(string(T.(condVar)));
T.(stateVar) = upper(string(T.(stateVar)));
T.(mouseVar) = string(T.(mouseVar));

% ---- summarize ----
conds = unique(T.(condVar));
states = unique(T.(stateVar));
bands = intersect(bandVars, T.Properties.VariableNames);

OUT = struct();
OUT.tables = struct();
OUT.figures = struct();

for s = 1:numel(states)
    st = states(s);
    Ts = T(T.(stateVar)==st,:);
    if isempty(Ts), continue; end

    fprintf('\n=== %s ===\n', st);

    for b = 1:numel(bands)
        band = bands{b};
        Y = Ts.(band);
        if all(isnan(Y)), continue; end

        % Group means
        G = groupsummary(Ts, {genoVar,condVar}, 'mean', band);

        % ---- plotting ----
        f = figure('Color','w','Position',[200 200 500 400]); hold on;
        cats = categorical(G.(genoVar) + "_" + G.(condVar));
        bar(cats, G.("mean_"+band));
        ylabel(sprintf('%s power (rel)', band));
        title(sprintf('%s – %s', st, band));
        set(gca,'Box','off','FontSize',11);
        xtickangle(30);

        % save
        figFile = fullfile(outDir, sprintf('EEG_%s_%s_bar.png', st, band));
        saveas(f, figFile);
        OUT.figures.(sprintf('%s_%s',st,band)) = figFile;
        close(f);

        % ---- stats: 2×2 ANOVA (genotype × condition) ----
        try
            L = Ts(:,{mouseVar,genoVar,condVar});
            L.Y = Ts.(band);
            L.(genoVar) = categorical(L.(genoVar));
            L.(condVar) = categorical(L.(condVar));
            L.(mouseVar) = categorical(L.(mouseVar));

            lme = fitlme(L,'Y ~ condition * genotype + (1|mouse)');
            a = anova(lme);

            OUT.tables.(sprintf('%s_%s',st,band)) = a;
            writetable(a, fullfile(outDir, sprintf('EEG_%s_%s_anova.csv', st, band)));
        catch ME
            warning('ANOVA failed for %s %s: %s', st, band, ME.message);
        end
    end
end

fprintf('\n✅ Group EEG band-power stats saved in %s\n', outDir);
end
