function plot_summary_nrem_pac(ALL_RESULTS, animals, conditions)
% PLOT_SUMMARY_NREM_PAC
% Plot NREM PAC (spindle-SO coupling) across animals

animal_ids = {animals.id};
nAnimals = numel(animals);
nCond = numel(conditions);

% Extract ROI MI values (SO-spindle: 0.5-1.5 Hz × 11-16 Hz)
phase_roi = [0.5 1.5];
amp_roi = [11 16];

mi_roi = nan(nAnimals, nCond);

for a = 1:nAnimals
    animal_id = animal_ids{a};
    for c = 1:nCond
        cond = conditions{c};
        
        if isfield(ALL_RESULTS.(animal_id).(cond), 'NREM_PAC')
            pac_data = ALL_RESULTS.(animal_id).(cond).NREM_PAC;
            
            if ~isempty(pac_data.MI)
                pf = pac_data.phase_freqs;
                af = pac_data.amp_freqs;
                
                ip = pf >= phase_roi(1) & pf <= phase_roi(2);
                ia = af >= amp_roi(1) & af <= amp_roi(2);
                
                roi_vals = pac_data.MI(ip, ia);
                mi_roi(a, c) = mean(roi_vals(:), 'omitnan');
            end
        end
    end
end

% Plot
figure('Color','w','Position',[140 140 1200 400]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Plot 1: MI by condition
nexttile; hold on;
mean_mi = nanmean(mi_roi, 1);
sem_mi = nanstd(mi_roi, 0, 1) ./ sqrt(sum(~isnan(mi_roi), 1));

bar(1:nCond, mean_mi, 'FaceColor', [0.6 0.8 0.9], 'EdgeColor', 'k');
errorbar(1:nCond, mean_mi, sem_mi, 'k.', 'LineWidth', 1.5, 'CapSize', 10);

% Individual animals
cols = lines(nAnimals);
for a = 1:nAnimals
    plot(1:nCond, mi_roi(a,:), 'o-', 'Color', cols(a,:), ...
         'MarkerFaceColor', cols(a,:), 'LineWidth', 1.5);
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('SO-Spindle PAC (MI)');
title('NREM PAC Across Conditions');
legend([{'Group Mean'}, animal_ids], 'Location', 'best', 'Box', 'off');
grid on; box on;

% Plot 2: Heatmap of all animals × conditions
nexttile;
imagesc(mi_roi');
colormap(turbo);
colorbar;
set(gca, 'XTick', 1:nAnimals, 'XTickLabel', animal_ids, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nCond, 'YTickLabel', conditions);
xlabel('Animal');
ylabel('Condition');
title('PAC Strength (MI)');

sgtitle('NREM PAC Summary (SO-Spindle)', 'FontSize', 14, 'FontWeight', 'bold');
end

