function plot_summary_rem_plv(ALL_RESULTS, animals, conditions)
% PLOT_SUMMARY_REM_PLV
% Plot REM PLV results across all animals and conditions

animal_ids = {animals.id};
nAnimals = numel(animals);
nCond = numel(conditions);

% Organize data (same structure as NREM)
data_by_condition = cell(1, nCond);

for c = 1:nCond
    cond = conditions{c};
    all_plv = [];
    all_dur = [];
    animal_idx = [];
    
    for a = 1:nAnimals
        animal_id = animal_ids{a};
        
        if isfield(ALL_RESULTS.(animal_id).(cond), 'REM_PLV')
            plv_data = ALL_RESULTS.(animal_id).(cond).REM_PLV;
            plv_vals = plv_data.plv_all;
            dur_vals = plv_data.dur_all;
            
            valid = ~isnan(plv_vals);
            all_plv = [all_plv; plv_vals(valid)];
            all_dur = [all_dur; dur_vals(valid)];
            animal_idx = [animal_idx; a*ones(sum(valid),1)];
        end
    end
    
    data_by_condition{c} = struct('plv', all_plv, 'dur', all_dur, 'animal', animal_idx);
end

% Create figure (similar structure to NREM)
figure('Color','w','Position',[120 120 1400 500]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Plot 1: PLV by condition
nexttile; hold on;
mean_plv = nan(1, nCond);
sem_plv = nan(1, nCond);
for c = 1:nCond
    if ~isempty(data_by_condition{c}.plv)
        plv_vals = data_by_condition{c}.plv;
        mean_plv(c) = mean(plv_vals);
        sem_plv(c) = std(plv_vals) / sqrt(length(plv_vals));
    end
end

bar(1:nCond, mean_plv, 'FaceColor', [0.9 0.6 0.6], 'EdgeColor', 'k');
errorbar(1:nCond, mean_plv, sem_plv, 'k.', 'LineWidth', 1.5, 'CapSize', 10);

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('REM PLV (theta EEG-ACh)');
title('REM PLV Across Conditions');
grid on; box on;

% Plot 2: Compare NREM vs REM PLV
nexttile; hold on;
% Would need NREM data passed in - simplified version here
bar(1:nCond, mean_plv, 'FaceColor', [0.9 0.6 0.6]);
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('PLV');
title('REM PLV (Theta-ACh Coupling)');
grid on; box on;

% Plot 3: Duration distribution
nexttile; hold on;
for c = 1:nCond
    if ~isempty(data_by_condition{c}.dur)
        histogram(data_by_condition{c}.dur, 'BinWidth', 10, ...
                 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
end
xlabel('REM Bout Duration (s)');
ylabel('Count');
title('REM Bout Duration Distribution');
legend(conditions, 'Location', 'best');
grid on; box on;

sgtitle('REM PLV Summary - All Animals', 'FontSize', 14, 'FontWeight', 'bold');
end
