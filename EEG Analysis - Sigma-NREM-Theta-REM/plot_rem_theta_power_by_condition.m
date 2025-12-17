function plot_rem_theta_power_by_condition(outRoot)
% plot_rem_theta_power_by_condition(outRoot)
%
% Uses:
%   outRoot/SigmaTheta_modulation_allmice.csv
%
% For EACH condition, makes a figure with:
%   (1) REM theta absolute power (WT vs APP)
%   (2) REM theta relative power (WT vs APP)
%   (3) REM theta peak frequency (WT vs APP)
%
% One PNG per condition is saved in outRoot.

    if nargin < 1 || isempty(outRoot)
        outRoot = 'SigmaTheta_ModAnalysis';
    end

    % ---------- Load group metrics ----------
    csvFile = fullfile(outRoot, 'SigmaTheta_modulation_allmice.csv');
    tbl = readtable(csvFile);

    % Expect these columns to exist:
    %   MouseID, Genotype, Condition,
    %   theta_abs_power, theta_rel_power, theta_peak_freq
    tbl.MouseID   = string(tbl.MouseID);
    tbl.Genotype  = string(tbl.Genotype);
    tbl.Condition = string(tbl.Condition);

    genotypes = ["WT","APP"];
    colMap.WT  = [0.6 0.6 0.6];
    colMap.APP = [0.392 0.584 0.929];

    % Conditions in original order
    [~, idxFirst] = unique(tbl.Condition, 'stable');
    conditions    = tbl.Condition(sort(idxFirst))';

    metrics = {'theta_abs_power','theta_rel_power','theta_peak_freq'};
    ylabels = {'REM theta power (a.u.)', ...
               'REM theta power (relative)', ...
               'REM theta peak freq (Hz)'};
    titles  = {'REM theta absolute power', ...
               'REM theta relative power', ...
               'REM theta peak frequency'};

    for ci = 1:numel(conditions)
        cond = conditions(ci);

        fig = figure('Color','w', 'Name', cond, ...
                     'Position',[100 100 1200 400]);

        for m = 1:3
            subplot(1,3,m); hold on;

            M = nan(1, numel(genotypes));
            E = nan(1, numel(genotypes));

            xPts = []; yPts = []; gIdx = [];

            for gi = 1:numel(genotypes)
                g = genotypes(gi);

                mask = tbl.Condition == cond & tbl.Genotype == g;
                vals = tbl.(metrics{m})(mask);

                if ~isempty(vals)
                    M(gi) = mean(vals, 'omitnan');
                    E(gi) = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));

                    xPts = [xPts; gi*ones(numel(vals),1)]; %#ok<AGROW>
                    yPts = [yPts; vals(:)];               %#ok<AGROW>
                    gIdx = [gIdx; gi*ones(numel(vals),1)]; %#ok<AGROW>
                end
            end

            % Bars WT vs APP
            bh = bar(1:numel(genotypes), M);
            set(gca,'XTick',1:numel(genotypes), ...
                    'XTickLabel', genotypes);

            bh.FaceColor = 'flat';
            bh.CData = [
                colMap.WT;
                colMap.APP;
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

        sgtitle(sprintf('Condition: %s – REM theta metrics', cond), ...
                'FontWeight','bold');

        % Save figure
        safeCond = regexprep(cond, '[^a-zA-Z0-9]', '');
        outPng   = fullfile(outRoot, sprintf('REMtheta_%s.png', safeCond));
        saveas(fig, outPng);
        fprintf('Saved REM theta metrics figure for %s: %s\n', cond, outPng);
    end
end
