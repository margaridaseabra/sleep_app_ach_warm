function OUT = run_EEG_bandpower_APPvsWT_simple(bandCsv, out_dir, varargin)
% run_EEG_bandpower_APPvsWT_simple
% -------------------------------------------------------------------------
% Simple APP vs WT EEG band-power stats + plots for ONE condition.
%
% INPUT
%   bandCsv  : path to EEG_band_power_allmice.csv
%              (from run_all_mice_eeg_psd_auto)
%
%   out_dir  : folder for stats + plots
%
% NAME–VALUE
%   'condition' : which Condition to analyze (default 'baseline')
%                 e.g. 'baseline' or 'ambtemp'
%   'useFDR'    : apply BH-FDR across bands (per state) for stars (default true)
%   'minN'      : minimum # of mice per genotype to run stats (default 3)
%
% CSV columns expected:
%   MouseID, Genotype, Condition, State, Band, Power_dB
%
% OUTPUT
%   OUT.success
%   OUT.state_order
%   OUT.band_order
%   OUT.perState.(StateName).stats   : table with per-band stats
%   OUT.perState.(StateName).fig_file: PNG path
% -------------------------------------------------------------------------

% ---------- args ----------
if nargin < 2 || isempty(out_dir)
    [p,~,~] = fileparts(bandCsv);
    out_dir = fullfile(p, 'EEG_group_APPvsWT');
end
if ~isfolder(out_dir), mkdir(out_dir); end

p = inputParser;
addParameter(p,'condition','baseline',@(s)ischar(s)||isstring(s));
addParameter(p,'useFDR',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'minN',3,@(x)isscalar(x)&&x>=1);
parse(p,varargin{:});

condSelect = lower(string(p.Results.condition));
useFDR     = p.Results.useFDR;
minN       = p.Results.minN;

% ---------- load & normalize table ----------
T = readtable(bandCsv);

% Expected columns:
%   MouseID, Genotype, Condition, State, Band, Power_dB
mustHave = {'MouseID','Genotype','Condition','State','Band','Power_dB'};
missing  = setdiff(mustHave, T.Properties.VariableNames);
if ~isempty(missing)
    error('EEG bandpower CSV is missing columns: %s', strjoin(missing, ', '));
end

MouseID   = string(T.MouseID);
Genotype  = string(T.Genotype);
Condition = lower(strtrim(string(T.Condition)));
State     = string(T.State);
Band      = string(T.Band);
Power_dB  = double(T.Power_dB);

% keep only requested condition
mask_cond = Condition == condSelect;
T = T(mask_cond,:);
MouseID   = MouseID(mask_cond);
Genotype  = Genotype(mask_cond);
State     = State(mask_cond);
Band      = Band(mask_cond);
Power_dB  = Power_dB(mask_cond);

if isempty(T)
    warning('No rows for condition "%s" in %s', condSelect, bandCsv);
    OUT = struct('success',false); 
    return;
end

% keep only WT + APP
mask_geno = (Genotype=="WT" | Genotype=="APP");
MouseID   = MouseID(mask_geno);
Genotype  = Genotype(mask_geno);
State     = State(mask_geno);
Band      = Band(mask_geno);
Power_dB  = Power_dB(mask_geno);

if isempty(MouseID)
    warning('No WT or APP rows for condition "%s".', condSelect);
    OUT = struct('success',false);
    return;
end

% Build a per-mouse table (average across any duplicate rows)
Tmouse = table(MouseID, Genotype, State, Band, Power_dB, ...
               'VariableNames', {'Mouse','Genotype','State','Band','Power_dB'});

G = groupsummary(Tmouse, {'Mouse','Genotype','State','Band'}, 'mean','Power_dB');
G.Properties.VariableNames{end} = 'Power_dB';  % rename mean_Power_dB -> Power_dB

% ---------- determine state & band order ----------
allStates  = unique(G.State);
allStates  = allStates(:);                      % column
prefStates = ["NREM","REM","WAKE"].';           % preferred order (column)

state_order = prefStates(ismember(prefStates, allStates));
state_order = [state_order; allStates(~ismember(allStates, state_order))];

allBands  = unique(G.Band);
allBands  = allBands(:);                        % column
prefBands = ["delta","theta","sigma","beta","gamma"].'; % preferred order

band_order = prefBands(ismember(prefBands, allBands));
band_order = [band_order; allBands(~ismember(allBands, band_order))];

OUT = struct();
OUT.success     = true;
OUT.state_order = state_order;
OUT.band_order  = band_order;
OUT.perState    = struct();

COL_WT  = [0.6 0.6 0.6];
COL_APP = [0.39 0.58 0.93];

% ---------- loop states ----------
for si = 1:numel(state_order)
    st = state_order(si);
    Gst = G(G.State==st, :);
    if isempty(Gst)
        continue;
    end

    fprintf('\n=== %s (condition: %s) ===\n', st, condSelect);

    nBands = numel(band_order);

    meanWT   = nan(1,nBands);
    semWT    = nan(1,nBands);
    meanAPP  = nan(1,nBands);
    semAPP   = nan(1,nBands);
    nWT_vec  = nan(1,nBands);
    nAPP_vec = nan(1,nBands);
    p_t      = nan(1,nBands);
    d_cohen  = nan(1,nBands);

    stats_rows = {};

    for bi = 1:nBands
        bd = band_order(bi);
        Gb = Gst(Gst.Band==bd, :);
        if isempty(Gb), continue; end

        valsWT  = Gb.Power_dB(Gb.Genotype=="WT");
        valsAPP = Gb.Power_dB(Gb.Genotype=="APP");

        valsWT  = valsWT(~isnan(valsWT));
        valsAPP = valsAPP(~isnan(valsAPP));

        nWT = numel(valsWT);
        nAPP = numel(valsAPP);

        if nWT>0
            meanWT(bi) = mean(valsWT);
            semWT(bi)  = std(valsWT)/sqrt(nWT);
        end
        if nAPP>0
            meanAPP(bi) = mean(valsAPP);
            semAPP(bi)  = std(valsAPP)/sqrt(nAPP);
        end

        nWT_vec(bi)  = nWT;
        nAPP_vec(bi) = nAPP;

        p_val = NaN;
        d_val = NaN;

        if nWT>=minN && nAPP>=minN
            [~, p_val] = ttest2(valsWT, valsAPP, 'Vartype','unequal');

            m1 = mean(valsWT);  m2 = mean(valsAPP);
            s1 = std(valsWT);   s2 = std(valsAPP);
            sp = sqrt(((nWT-1)*s1^2 + (nAPP-1)*s2^2) / max(1,(nWT+nAPP-2)));
            if sp>0
                d_val = (m2 - m1)/sp;   % APP - WT
            end
        end

        p_t(bi)     = p_val;
        d_cohen(bi) = d_val;

        stats_rows = [stats_rows; ...
            {char(bd), nWT, nAPP, ...
             meanWT(bi), meanAPP(bi), ...
             semWT(bi), semAPP(bi), ...
             p_val, d_val}]; %#ok<AGROW>
    end

    % ---------- FDR across bands (per state) ----------
    p_t_fdr = nan(size(p_t));
    if useFDR
        valid = ~isnan(p_t);
        pv = p_t(valid);
        if ~isempty(pv)
            [sp, idx] = sort(pv(:));
            m = numel(sp);
            adj = sp .* (m./(1:m))';
            for k = m-1:-1:1
                adj(k) = min(adj(k), adj(k+1));
            end
            adj(adj>1) = 1;
            adj_full = nan(size(pv));
            adj_full(idx) = adj;
            p_t_fdr(valid) = adj_full;
        end
    end

    % ---------- build stats table ----------
    if isempty(stats_rows)
        stats_tbl = table();
    else
        stats_tbl = cell2table(stats_rows, ...
            'VariableNames', {'Band','nWT','nAPP', ...
                              'MeanWT','MeanAPP', ...
                              'SEM_WT','SEM_APP', ...
                              'p_ttest','Cohen_d'});
        % add FDR column
        pF = nan(height(stats_tbl),1);
        for r = 1:height(stats_tbl)
            bd = string(stats_tbl.Band(r));
            bi = find(band_order==bd,1);
            if ~isempty(bi)
                pF(r) = p_t_fdr(bi);
            end
        end
        stats_tbl.p_ttest_FDR = pF;
    end

    OUT.perState.(matlab.lang.makeValidName(st)).stats = stats_tbl;

    % ---------- plotting ----------
    fig = figure('Color','w','Position',[200 200 700 420]); hold on;

    x = 1:nBands;
    barWidth = 0.35;

    % WT
    bar(x - barWidth/2, meanWT, barWidth, ...
        'FaceColor',COL_WT, 'EdgeColor','none');
    errorbar(x - barWidth/2, meanWT, semWT, 'k','LineStyle','none');

    % APP
    bar(x + barWidth/2, meanAPP, barWidth, ...
        'FaceColor',COL_APP, 'EdgeColor','none');
    errorbar(x + barWidth/2, meanAPP, semAPP, 'k','LineStyle','none');

    all_means = [meanWT(~isnan(meanWT)), meanAPP(~isnan(meanAPP))];
    if isempty(all_means)
        ymax = 1;
    else
        ymax = max(all_means) + 0.25*max(all_means);
    end
    ylim([min([all_means,0]) - 1, ymax]);


    % stars
    for bi = 1:nBands
        if isnan(p_t(bi)), continue; end
        if useFDR
            p_use = p_t_fdr(bi);
        else
            p_use = p_t(bi);
        end
        if isnan(p_use) || p_use >= 0.05
            continue;
        end

        if p_use < 0.001, stars = '***';
        elseif p_use < 0.01, stars = '**';
        else, stars = '*';
        end

        y_star = max(meanWT(bi)+semWT(bi), meanAPP(bi)+semAPP(bi));
        if isnan(y_star), y_star = ymax*0.8; end
        y_star = y_star + 0.05*ymax;

        text(x(bi), y_star, stars, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
    end

    xticks(x);
    xticklabels(band_order);
    xlabel('Band');
    ylabel('Power (dB)');
    title(sprintf('%s – %s (APP vs WT)', st, condSelect));

    legend({'WT','APP'}, 'Location','northoutside', ...
           'Orientation','horizontal');
    set(gca,'Box','off','FontSize',11);

    fig_file = fullfile(out_dir, ...
        sprintf('EEG_bandpower_%s_%s_APPvsWT.png', lower(st), condSelect));
    saveas(fig, fig_file);

    OUT.perState.(matlab.lang.makeValidName(st)).fig_file = fig_file;

    % also save CSV
    if ~isempty(stats_tbl)
        csv_out = fullfile(out_dir, ...
            sprintf('EEG_bandpower_%s_%s_APPvsWT_stats.csv', lower(st), condSelect));
        writetable(stats_tbl, csv_out);
        fprintf('  → stats saved: %s\n', csv_out);
    end
end

fprintf('\n✅ EEG band-power APP vs WT stats done for condition "%s".\n', condSelect);
end
