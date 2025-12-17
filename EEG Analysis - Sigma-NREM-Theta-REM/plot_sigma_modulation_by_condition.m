function plot_sigma_modulation_by_condition(outRoot)
% plot_sigma_modulation_by_condition(outRoot)
%
% Reads:
%   - outRoot/SigmaTheta_modulation_allmice.csv
%   - outRoot/<Genotype>/<Condition>/<MouseID>/SigmaTheta_OUT.mat
%
% For EACH Condition, makes a figure with:
%   (1) mean sigma modulation PSD (WT vs APP, 0–0.15 Hz, A.U.)
%   (2) sigma_mod_power bar plot WT vs APP
%   (3) sigma_mod_peak_f  bar plot WT vs APP
%   (4) sigma_mod_peak_amp bar plot WT vs APP
%
% One PNG per condition is saved in outRoot.

    if nargin < 1 || isempty(outRoot)
        outRoot = 'SigmaTheta_ModAnalysis';
    end

    % ---------- Load group metrics ----------
    csvFile = fullfile(outRoot, 'SigmaTheta_modulation_allmice.csv');
    tbl = readtable(csvFile);

    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);

    % Genotypes and colours
    genotypes = ["WT","APP"];
    colMap.WT  = [0.6 0.6 0.6];
    colMap.APP = [0.392 0.584 0.929];

    % Conditions (keep original order)
    [~, idxFirst] = unique(tbl.Condition, 'stable');
    conditions    = tbl.Condition(sort(idxFirst))';

    % Common modulation frequency axis for group PSD
    F_common = (0.02:0.001:0.15)';   % 0–0.15 Hz, good resolution

    for ci = 1:numel(conditions)
        cond = conditions(ci);

        fig = figure('Color','w', 'Name', cond, ...
                     'Position',[100 100 1500 400]);

        % --------- 1) Group sigma modulation PSD ----------
        subplot(1,4,1); hold on;
        hLines = gobjects(1, numel(genotypes));

        for gi = 1:numel(genotypes)
            g = genotypes(gi);

            % Find all mice of this genotype & condition
            mask = tbl.Condition == cond & tbl.Genotype == g;
            sub  = tbl(mask, :);
            if isempty(sub)
                continue;
            end

            allPSD = [];

            for r = 1:height(sub)
                mouseID = sub.MouseID(r);
                sesDir  = fullfile(outRoot, char(g), char(cond), char(mouseID));
                matFile = fullfile(sesDir, 'SigmaTheta_OUT.mat');
                if ~exist(matFile,'file')
                    warning('Missing SigmaTheta_OUT for %s %s %s', ...
                            mouseID, g, cond);
                    continue;
                end
                S = load(matFile);
                OUT = S.OUT;

                f  = OUT.sigma.mod_f(:);
                ps = OUT.sigma.mod_psd_plot(:);          % raw PSD of envelope
                
                % Interpolate onto common freq axis (start at mod_f_min to avoid extrap):
                F_common = (0.01:0.001:0.15)';

                ps_common = interp1(f, ps, F_common, 'linear', NaN);

                allPSD(:, end+1) = ps_common;
            end

            if isempty(allPSD)
                continue;
            end

            mPSD = mean(allPSD, 2, 'omitnan');
            nSub = sum(~isnan(allPSD(1,:)));
            sPSD = std(allPSD, 0, 2, 'omitnan') ./ max(nSub,1).^0.5;

            c = colMap.(g);

            upper = mPSD + sPSD;
            lower = mPSD - sPSD;

            fill([F_common; flipud(F_common)], ...
                 [upper; flipud(lower)], ...
                 c, 'FaceAlpha',0.2, 'EdgeColor','none');
            hLines(gi) = plot(F_common, mPSD, 'Color', c*0.9, 'LineWidth',2);
        end

        xlim([0 0.15]);
        xlabel('Modulation frequency (Hz)');
        ylabel('Sigma power (A.U.)');
        title(sprintf('%s – sigma modulation', cond));
        grid on;
        legend(hLines(isgraphics(hLines)), genotypes, 'Location','northeast');

        % --------- 2–4) Bar plots of metrics ----------
        metrics = {'sigma_mod_power','sigma_mod_peak_f','sigma_mod_peak_amp'};
        ylabels = {'Modulation power (A.U.)', ...
                   'Peak modulation frequency (Hz)', ...
                   'Peak modulation amplitude (A.U.)'};
        titles  = {'Sigma modulation power', ...
                   'Sigma peak modulation frequency', ...
                   'Sigma peak modulation amplitude'};

        for m = 1:3
            subplot(1,4,1+m); hold on;

            M = nan(1, numel(genotypes));
            E = nan(1, numel(genotypes));

            % Keep individual points for scatter
            xPts = []; yPts = []; gIdx = [];

            for gi = 1:numel(genotypes)
                g = genotypes(gi);
                mask = tbl.Condition == cond & tbl.Genotype == g;
                vals = tbl.(metrics{m})(mask);

                if ~isempty(vals)
                    M(gi) = mean(vals, 'omitnan');
                    E(gi) = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));

                    xPts = [xPts; gi*ones(numel(vals),1)]; %#ok<AGROW>
                    yPts = [yPts; vals(:)];              %#ok<AGROW>
                    gIdx = [gIdx; gi*ones(numel(vals),1)]; %#ok<AGROW>
                end
            end

            % Bars (WT = grey, APP = cornflower blue)
            bh = bar(1:numel(genotypes), M);
            set(gca,'XTick',1:numel(genotypes), ...
                    'XTickLabel', genotypes);

            bh.FaceColor = 'flat';
            bh.CData = [
                colMap.WT;   % bar at x=1 (WT)
                colMap.APP;  % bar at x=2 (APP)
            ];

            % Error bars
            errorbar(1:numel(genotypes), M, E, 'k', ...
                'LineStyle','none','LineWidth',1);

            % Scatter individual mice
            for gi = 1:numel(genotypes)
                jitter = (rand(sum(gIdx==gi),1)-0.5)*0.15;
                plot(gi + jitter, yPts(gIdx==gi), 'k.', 'MarkerSize',10);
            end

            ylabel(ylabels{m});
            title(titles{m});
            box off; grid on;
        end
        sgtitle(sprintf('Condition: %s', cond), 'FontWeight','bold');

        % --------- Save figure ----------
        safeCond = regexprep(cond, '[^a-zA-Z0-9]', '');
        outPng   = fullfile(outRoot, sprintf('SigmaMod_%s.png', safeCond));
        saveas(fig, outPng);
        fprintf('Saved sigma modulation figure for %s: %s\n', cond, outPng);
        

    end
end
