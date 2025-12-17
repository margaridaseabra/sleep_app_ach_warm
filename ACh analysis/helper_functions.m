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


function plot_genotype_comparison(ALL_RESULTS, animals, conditions)
% PLOT_GENOTYPE_COMPARISON
% Compare APP vs WT across all metrics

genotypes = {animals.genotype};
unique_genos = unique(genotypes);
nGeno = numel(unique_genos);
nCond = numel(conditions);

% Initialize storage
geno_data = struct();
for g = 1:nGeno
    geno_data.(unique_genos{g}) = struct();
end

% Collect data by genotype
animal_ids = {animals.id};
for a = 1:numel(animals)
    geno = genotypes{a};
    animal_id = animal_ids{a};
    
    for c = 1:nCond
        cond = conditions{c};
        
        % NREM PLV
        if isfield(ALL_RESULTS.(animal_id).(cond), 'NREM_PLV')
            plv_data = ALL_RESULTS.(animal_id).(cond).NREM_PLV.plv_all;
            if ~isfield(geno_data.(geno), 'NREM_PLV')
                geno_data.(geno).NREM_PLV = cell(1, nCond);
            end
            geno_data.(geno).NREM_PLV{c} = [geno_data.(geno).NREM_PLV{c}; plv_data(~isnan(plv_data))];
        end
        
        % Add REM PLV, NREM PAC, REM PAC similarly...
    end
end

% Plot comparison
figure('Color','w','Position',[180 180 1200 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% NREM PLV comparison
nexttile; hold on;
x_offset = [-0.15, 0.15];
cols_geno = [0.8 0.3 0.3; 0.3 0.3 0.8];

for g = 1:nGeno
    geno = unique_genos{g};
    means = nan(1, nCond);
    sems = nan(1, nCond);
    
    for c = 1:nCond
        if isfield(geno_data.(geno), 'NREM_PLV') && ~isempty(geno_data.(geno).NREM_PLV{c})
            vals = geno_data.(geno).NREM_PLV{c};
            means(c) = mean(vals);
            sems(c) = std(vals) / sqrt(length(vals));
        end
    end
    
    x_pos = (1:nCond) + x_offset(g);
    bar(x_pos, means, 0.3, 'FaceColor', cols_geno(g,:));
    errorbar(x_pos, means, sems, 'k.', 'LineWidth', 1.5);
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions);
ylabel('NREM PLV');
title('NREM PLV: Genotype Comparison');
legend(unique_genos, 'Location', 'best');
grid on; box on;

% Add other comparisons...
nexttile;
text(0.5, 0.5, 'Add REM PLV comparison here', ...
     'HorizontalAlignment', 'center', 'FontSize', 12);

nexttile;
text(0.5, 0.5, 'Add NREM PAC comparison here', ...
     'HorizontalAlignment', 'center', 'FontSize', 12);

nexttile;
text(0.5, 0.5, 'Add REM PAC comparison here', ...
     'HorizontalAlignment', 'center', 'FontSize', 12);

sgtitle('Genotype × Condition Interaction', 'FontSize', 16, 'FontWeight', 'bold');
end