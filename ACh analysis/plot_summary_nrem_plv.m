function plot_summary_nrem_plv(ALL_RESULTS, animals, conditions)
% PLOT_SUMMARY_NREM_PLV
% Plot NREM PLV results across all animals and conditions

animal_ids = {animals.id};
genotypes = {animals.genotype};
nAnimals = numel(animals);
nCond = numel(conditions);

% Organize data
data_by_condition = cell(1, nCond);
data_by_genotype_cond = struct();

for c = 1:nCond
    cond = conditions{c};
    all_plv = [];
    all_dur = [];
    animal_idx = [];
    genotype_idx = [];
    
    for a = 1:nAnimals
        animal_id = animal_ids{a};
        geno = genotypes{a};
        
        if isfield(ALL_RESULTS.(animal_id).(cond), 'NREM_PLV')
            plv_data = ALL_RESULTS.(animal_id).(cond).NREM_PLV;
            plv_vals = plv_data.plv_all;
            dur_vals = plv_data.dur_all;
            
            valid = ~isnan(plv_vals);
            all_plv = [all_plv; plv_vals(valid)];
            all_dur = [all_dur; dur_vals(valid)];
            animal_idx = [animal_idx; a*ones(sum(valid),1)];
            
            % Store by genotype
            if ~isfield(data_by_genotype_cond, geno)
                data_by_genotype_cond.(geno).(cond) = [];
            end
            data_by_genotype_cond.(geno).(cond) = [data_by_genotype_cond.(geno).(cond); plv_vals(valid)];
        end
    end
    
    data_by_condition{c} = struct('plv', all_plv, 'dur', all_dur, 'animal', animal_idx);
end

% Create figure
figure('Color','w','Position',[100 100 1400 500]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Plot 1: PLV by condition (all animals)
nexttile; hold on;
mean_plv = nan(1, nCond);
sem_plv = nan(1, nCond);
for c = 1:nCond
    plv_vals = data_by_condition{c}.plv;
    mean_plv(c) = mean(plv_vals);
    sem_plv(c) = std(plv_vals) / sqrt(length(plv_vals));
end

bar(1:nCond, mean_plv, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
errorbar(1:nCond, mean_plv, sem_plv, 'k.', 'LineWidth', 1.5, 'CapSize', 10);

% Add individual animal means
cols = lines(nAnimals);
for a = 1:nAnimals
    animal_means = nan(1, nCond);
    for c = 1:nCond
        data = data_by_condition{c};
        animal_data = data.plv(data.animal == a);
        if ~isempty(animal_data)
            animal_means(c) = mean(animal_data);
        end
    end
    plot(1:nCond, animal_means, 'o-', 'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), 'LineWidth', 1.5);
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('NREM PLV (delta EEG-ACh)');
title('NREM PLV Across Conditions');
legend([{'Group Mean'}, animal_ids], 'Location', 'best', 'Box', 'off');
grid on; box on;

% Plot 2: PLV by genotype
nexttile; hold on;
genotype_names = fieldnames(data_by_genotype_cond);
nGeno = numel(genotype_names);
x_offset = [-0.15, 0.15];
cols_geno = [0.8 0.3 0.3; 0.3 0.3 0.8];

for g = 1:nGeno
    geno = genotype_names{g};
    means = nan(1, nCond);
    sems = nan(1, nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        if isfield(data_by_genotype_cond.(geno), cond)
            vals = data_by_genotype_cond.(geno).(cond);
            means(c) = mean(vals);
            sems(c) = std(vals) / sqrt(length(vals));
        end
    end
    
    x_pos = (1:nCond) + x_offset(g);
    bar(x_pos, means, 0.3, 'FaceColor', cols_geno(g,:), 'EdgeColor', 'k');
    errorbar(x_pos, means, sems, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('NREM PLV');
title('NREM PLV by Genotype');
legend(genotype_names, 'Location', 'best');
grid on; box on;

% Plot 3: PLV vs Duration correlation
nexttile; hold on;
for c = 1:nCond
    scatter(data_by_condition{c}.dur, data_by_condition{c}.plv, 50, 'filled', ...
            'MarkerFaceAlpha', 0.5);
end

xlabel('Bout Duration (s)');
ylabel('NREM PLV');
title('PLV vs Bout Duration');
legend(conditions, 'Location', 'best');
grid on; box on;

% Overall correlation
all_plv = vertcat(data_by_condition{:}.plv);
all_dur = vertcat(data_by_condition{:}.dur);
[r, p] = corr(all_dur, all_plv, 'rows', 'complete');
text(0.05, 0.95, sprintf('r = %.3f, p = %.3f', r, p), ...
     'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');

sgtitle('NREM PLV Summary - All Animals', 'FontSize', 14, 'FontWeight', 'bold');
end
