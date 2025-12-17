function plot_summary_rem_pac(ALL_RESULTS, animals, conditions)
% PLOT_SUMMARY_REM_PAC  
% Plot REM PAC (theta-gamma coupling) across animals

% Similar structure to NREM PAC
animal_ids = {animals.id};
nAnimals = numel(animals);
nCond = numel(conditions);

% Extract ROI MI values (theta-gamma: 7-9 Hz × 40-70 Hz)
phase_roi = [7 9];
amp_roi = [40 70];

mi_roi = nan(nAnimals, nCond);

for a = 1:nAnimals
    animal_id = animal_ids{a};
    for c = 1:nCond
        cond = conditions{c};
        
        if isfield(ALL_RESULTS.(animal_id).(cond), 'REM_PAC')
            pac_data = ALL_RESULTS.(animal_id).(cond).REM_PAC;
            
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
figure('Color','w','Position',[160 160 800 400]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile; hold on;
mean_mi = nanmean(mi_roi, 1);
sem_mi = nanstd(mi_roi, 0, 1) ./ sqrt(sum(~isnan(mi_roi), 1));

bar(1:nCond, mean_mi, 'FaceColor', [0.9 0.7 0.5], 'EdgeColor', 'k');
errorbar(1:nCond, mean_mi, sem_mi, 'k.', 'LineWidth', 1.5);

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('Theta-Gamma PAC (MI)');
title('REM PAC Across Conditions');
grid on; box on;

nexttile;
imagesc(mi_roi');
colormap(turbo);
colorbar;
set(gca, 'XTick', 1:nAnimals, 'XTickLabel', animal_ids, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nCond, 'YTickLabel', conditions);
title('PAC Strength');

sgtitle('REM PAC Summary (Theta-Gamma)', 'FontSize', 14, 'FontWeight', 'bold');
end
