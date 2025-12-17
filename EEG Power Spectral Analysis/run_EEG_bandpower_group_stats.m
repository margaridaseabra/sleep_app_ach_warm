function OUT = run_EEG_bandpower_group_stats(bandpowerCsv, out_dir)
% run_EEG_bandpower_group_stats
% -------------------------------------------------------------------------
% Designed for your long-format EEG_band_power_allmice.csv:
%
% Columns (required):
%   MouseID   : mouse identifier (string/char/numeric)
%   Genotype  : 'WT' / 'APP'
%   Condition : 'baseline','ambtemp','drugs', ...
%   State     : 'Wake','NREM','REM', ...
%   Band      : 'Delta','Theta','Sigma','Beta','lGamma1','lGamma2','hGamma',...
%   Power_dB  : band power in dB
%
% This function currently does:
%   A) Baseline-only APP vs WT stats per State x Band
%      - mean Power_dB per mouse
%      - t-test APP vs WT
%      - ranksum
%      - Cohen's d
%      - BH-FDR across bands within each state
%      - bar plots (WT vs APP) per state with stars on sig. bands
%
% OUTPUT:
%   OUT.baseline.stats   : table with all state x band stats
%   OUT.baseline.figures : struct with fig file paths per state
% -------------------------------------------------------------------------

    if nargin < 2 || isempty(out_dir)
        [p,~,~] = fileparts(bandpowerCsv);
        out_dir = fullfile(p, 'EEG_group_stats');
    end
    if ~isfolder(out_dir), mkdir(out_dir); end

    T = readtable(bandpowerCsv);

    % ---------- normalize types ----------
    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = lower(strtrim(string(T.Condition)));
    T.State     = string(T.State);
    T.Band      = string(T.Band);
    T.Power_dB  = double(T.Power_dB);

    OUT = struct();
    [OUT.baseline.stats, OUT.baseline.figures] = local_baseline_APPvsWT(T, out_dir);

    save(fullfile(out_dir,'EEG_group_stats.mat'), 'OUT');
    fprintf('✅ EEG band-power group stats saved in %s\n', out_dir);
end

% =====================================================================
function [stats_tbl, fig_struct] = local_baseline_APPvsWT(T, out_dir)
% Baseline-only APP vs WT, per State x Band

    Tb = T(T.Condition=="baseline", :);
    if isempty(Tb)
        warning('No baseline rows in EEG band-power table.');
        stats_tbl  = table();
        fig_struct = struct();
        return;
    end

    states = unique(Tb.State);
    bands_all = unique(Tb.Band);

    % nice band order if present
    band_order_pref = ["Delta","Theta","Sigma","Beta","lGamma1","lGamma2","hGamma"];
    band_order = band_order_pref(ismember(band_order_pref, bands_all));
    band_order = [band_order, setdiff(bands_all, band_order_pref, 'stable')];

    rows = {};
    fig_struct = struct();

    for si = 1:numel(states)
        st = states(si);
        Ts = Tb(Tb.State==st, :);
        if isempty(Ts), continue; end

        p_vec = nan(numel(band_order),1);  % for FDR
        have_row = false(numel(band_order),1);

        for bi = 1:numel(band_order)
            bd = band_order(bi);
            Tsb = Ts(Ts.Band==bd, :);
            if isempty(Tsb), continue; end

            % average per mouse x genotype
            G = groupsummary(Tsb, {'MouseID','Genotype'}, 'mean', 'Power_dB');
            if ~ismember("Genotype", G.Properties.VariableNames)
                % older MATLAB sometimes renames; guard
                if ismember("Group_Genotype", G.Properties.VariableNames)
                    G.Genotype = G.Group_Genotype;
                else
                    error('Could not find Genotype column after groupsummary.');
                end
            end
            meanPow = G.("mean_Power_dB");

            valsWT  = meanPow(G.Genotype=="WT");
            valsAPP = meanPow(G.Genotype=="APP");

            valsWT  = valsWT(~isnan(valsWT));
            valsAPP = valsAPP(~isnan(valsAPP));

            nWT  = numel(valsWT);
            nAPP = numel(valsAPP);

            if nWT < 2 || nAPP < 2
                continue;
            end

            % t-test & ranksum
            [~, p_t] = ttest2(valsWT, valsAPP, 'Vartype','unequal');
            p_rs     = ranksum(valsWT, valsAPP);

            mWT   = mean(valsWT,'omitnan');
            mAPP  = mean(valsAPP,'omitnan');
            sWT   = std(valsWT,'omitnan');
            sAPP  = std(valsAPP,'omitnan');
            sp    = sqrt(((nWT-1)*sWT^2 + (nAPP-1)*sAPP^2) / max(1,(nWT+nAPP-2)));
            d_eff = (mAPP - mWT) / sp;

            rows = [rows; ...
                {char(st), char(bd), nWT, nAPP, mWT, mAPP, p_t, p_rs, d_eff}]; %#ok<AGROW>

            p_vec(bi) = p_t;
            have_row(bi) = true;
        end

        % ------------------ FDR per state ------------------
        idx_valid = find(have_row & ~isnan(p_vec));
        pF = nan(size(p_vec));
        if ~isempty(idx_valid)
            pF(idx_valid) = local_fdr(p_vec(idx_valid));
        end

        % store temporarily in fig_struct so we can use for plotting
        fig_struct.(matlab.lang.makeValidName(st)).p_raw = p_vec;
        fig_struct.(matlab.lang.makeValidName(st)).p_FDR = pF;
        fig_struct.(matlab.lang.makeValidName(st)).bands = band_order(:);
    end

    if isempty(rows)
        warning('No valid APP vs WT baseline comparisons (too few mice per group?).');
        stats_tbl  = table();
        fig_struct = struct();
        return;
    end

    stats_tbl = cell2table(rows, ...
        'VariableNames', {'State','Band','nWT','nAPP', ...
                          'MeanWT_dB','MeanAPP_dB', ...
                          'p_ttest','p_ranksum','Cohen_d'});
    % Convert State and Band columns from cellstr to string for easy comparison
    stats_tbl.State = string(stats_tbl.State);
    stats_tbl.Band  = string(stats_tbl.Band);


    % ---------- attach FDR to stats_tbl ----------
    pF_all = nan(height(stats_tbl),1);
    for r = 1:height(stats_tbl)
        st = string(stats_tbl.State(r));
        bd = string(stats_tbl.Band(r));
        f  = fig_struct.(matlab.lang.makeValidName(st));
        bi = find(f.bands == bd, 1);
        if ~isempty(bi)
            pF_all(r) = f.p_FDR(bi);
        end
    end
    stats_tbl.p_ttest_FDR = pF_all;

    % save CSV
    csv_out = fullfile(out_dir, 'EEG_bandpower_baseline_APPvsWT_stats.csv');
    writetable(stats_tbl, csv_out);
    fprintf('Baseline APP vs WT stats saved to: %s\n', csv_out);

    % ---------- plots per state ----------
    COL_WT  = [0.6 0.6 0.6];
    COL_APP = [0.39 0.58 0.93];

    uniqueStates = unique(stats_tbl.State);
    for si = 1:numel(uniqueStates)
        st = uniqueStates(si);
        sub = stats_tbl(stats_tbl.State==st, :);
        f  = fig_struct.(matlab.lang.makeValidName(st));

        % keep band order & indices
        [~, idx] = ismember(sub.Band, f.bands);
        [~, order] = sort(idx);
        sub = sub(order,:);

        bands = string(sub.Band);
        nb    = numel(bands);

        x = 1:nb;
        mWT   = sub.MeanWT_dB;
        mAPP  = sub.MeanAPP_dB;
        p_raw = sub.p_ttest;
        p_fdr = sub.p_ttest_FDR;

        fig = figure('Color','w','Position',[150 150 800 400]); hold on;

        barWidth = 0.35;
        xWT  = x - barWidth/2;
        xAPP = x + barWidth/2;

        % WT bars
        bar(xWT, mWT, barWidth, 'FaceColor',COL_WT, 'EdgeColor','none');
        % APP bars
        bar(xAPP, mAPP, barWidth, 'FaceColor',COL_APP, 'EdgeColor','none');

        % no SEMs here (we could add them later if you like)

        ymax = max([mWT(:); mAPP(:)],[],'omitnan');
        if isempty(ymax) || isnan(ymax), ymax = 1; end
        ylim([ymax-20, ymax+5]); % dB, zoom into top (optional tweak)

        % stars above band where FDR < 0.05
        for bi = 1:nb
            p_use = p_fdr(bi);
            if isnan(p_use) || p_use >= 0.05
                continue;
            end
            stars = local_p2stars(p_use);
            y_star = ymax + 0.5;
            text(x(bi), y_star, stars, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',11,'FontWeight','bold');
        end

        set(gca,'XTick',x,'XTickLabel',bands, ...
                'XTickLabelRotation',45,'FontSize',11);
        ylabel('Band power (dB)');
        title(sprintf('Baseline band power: APP vs WT (%s)', st));

        legend({'WT','APP'}, 'Location','northoutside','Orientation','horizontal');
        set(gca,'Box','off');

        fpath = fullfile(out_dir, sprintf('EEG_baseline_APPvsWT_%s.png', lower(st)));
        saveas(fig, fpath);
        fig_struct.(matlab.lang.makeValidName(st)).fig_file = fpath;
    end
end

% =====================================================================
function p_adj = local_fdr(p)
% BH-FDR for a vector p (no NaNs)
    p = p(:);
    m = numel(p);
    [sp, idx] = sort(p);
    adj = sp .* (m ./ (1:m)');
    for k = m-1:-1:1
        adj(k) = min(adj(k), adj(k+1));
    end
    adj(adj>1) = 1;
    p_adj = nan(size(p));
    p_adj(idx) = adj;
end

% =====================================================================
function stars = local_p2stars(p)
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = '';
    end
end
