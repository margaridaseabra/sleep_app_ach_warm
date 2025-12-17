%% SIMPLE PLOTTING SCRIPT - Works with your data
% Run this after loading: load('all_animals_plv_pac_results.mat')

clearvars -except ALL_RESULTS animals conditions CODES
close all;

% Load if not already in workspace
if ~exist('ALL_RESULTS', 'var')
    load('all_animals_plv_pac_results.mat');
end

animal_ids = {animals.id};
genotypes = {animals.genotype};
nAnimals = numel(animals);
nCond = numel(conditions);

fprintf('\n=== Creating Summary Plots ===\n');
fprintf('Animals: %d\n', nAnimals);
fprintf('Conditions: %d\n', nCond);
for i = 1:nAnimals
    fprintf('  %d. %s (%s)\n', i, animal_ids{i}, genotypes{i});
end
fprintf('\n');

%% ============================================================
%  FIGURE 1: NREM & REM PLV SUMMARY
% ============================================================

figure('Color','w','Position',[50 50 1600 900]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

% --- Panel 1: NREM PLV by condition ---
nexttile; hold on;
cols = lines(nAnimals);
for a = 1:nAnimals
    plv_means = nan(1,nCond);
    plv_sems = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
           isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
            plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
            plv_data = plv_data(~isnan(plv_data));
            if ~isempty(plv_data)
                plv_means(c) = mean(plv_data);
                plv_sems(c) = std(plv_data) / sqrt(length(plv_data));
            end
        end
    end
    
    errorbar(1:nCond, plv_means, plv_sems, 'o-', 'LineWidth', 2, ...
             'MarkerSize', 10, 'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), ...
             'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Mean NREM PLV');
title('NREM PLV (Delta EEG-ACh)');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 2: REM PLV by condition ---
nexttile; hold on;
for a = 1:nAnimals
    plv_means = nan(1,nCond);
    plv_sems = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
           isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'REM_PLV')
            plv_data = ALL_RESULTS.(animal_ids{a}).(cond).REM_PLV.plv_all;
            plv_data = plv_data(~isnan(plv_data));
            if ~isempty(plv_data)
                plv_means(c) = mean(plv_data);
                plv_sems(c) = std(plv_data) / sqrt(length(plv_data));
            end
        end
    end
    
    errorbar(1:nCond, plv_means, plv_sems, 'o-', 'LineWidth', 2, ...
             'MarkerSize', 10, 'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), ...
             'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Mean REM PLV');
title('REM PLV (Theta EEG-ACh)');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 3: NREM PAC by condition ---
nexttile; hold on;
phase_roi = [2.0 3.0];  % Delta (avoid 1.5 Hz that has NaN)
amp_roi = [11 16];      % Spindle

for a = 1:nAnimals
    mi_means = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
           isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PAC')
            pac_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PAC;
            
            if ~isempty(pac_data.MI) && ~isempty(pac_data.phase_freqs)
                pf = pac_data.phase_freqs;
                af = pac_data.amp_freqs;
                ip = pf >= phase_roi(1) & pf <= phase_roi(2);
                ia = af >= amp_roi(1) & af <= amp_roi(2);
                
                if any(ip) && any(ia)
                    mi_means(c) = mean(pac_data.MI(ip, ia), 'all', 'omitnan');
                end
            end
        end
    end
    
    plot(1:nCond, mi_means, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), ...
         'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Mean PAC (MI)');
title('NREM PAC (Delta-Spindle)');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 4: NREM bout count ---
nexttile; hold on;
for a = 1:nAnimals
    n_bouts = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
           isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
            n_bouts(c) = numel(ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all);
        end
    end
    
    plot(1:nCond, n_bouts, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'Color', cols(a,:), 'MarkerFaceColor', cols(a,:), ...
         'DisplayName', animal_ids{a});
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Number of Bouts');
title('NREM Bout Count');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 5: Genotype comparison - NREM PLV ---
nexttile; hold on;
unique_genos = unique(genotypes);
nGeno = numel(unique_genos);
x_offset = linspace(-0.15, 0.15, nGeno);
geno_colors = [0.8 0.3 0.3; 0.3 0.3 0.8];

for g = 1:nGeno
    geno = unique_genos{g};
    means = nan(1,nCond);
    sems = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        all_plv = [];
        
        % Collect data from all animals of this genotype
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
                    plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
                    plv_data = plv_data(~isnan(plv_data));
                    all_plv = [all_plv; plv_data];
                end
            end
        end
        
        if ~isempty(all_plv)
            means(c) = mean(all_plv);
            sems(c) = std(all_plv) / sqrt(length(all_plv));
        end
    end
    
    x_pos = (1:nCond) + x_offset(g);
    bar(x_pos, means, 0.3, 'FaceColor', geno_colors(g,:), 'EdgeColor', 'k', ...
        'DisplayName', geno);
    errorbar(x_pos, means, sems, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('NREM PLV');
title('Genotype Comparison - NREM');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 6: PLV distribution ---
nexttile; hold on;
for g = 1:nGeno
    geno = unique_genos{g};
    all_plv = [];
    
    for c = 1:nCond
        cond = conditions{c};
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
                    plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
                    all_plv = [all_plv; plv_data(~isnan(plv_data))];
                end
            end
        end
    end
    
    histogram(all_plv, 'BinWidth', 0.02, 'FaceColor', geno_colors(g,:), ...
              'FaceAlpha', 0.5, 'DisplayName', geno);
end
xlabel('NREM PLV');
ylabel('Count');
title('PLV Distribution by Genotype');
legend('Location', 'best');
grid on; box on;

% --- Panel 7: PAC comparison by genotype ---
nexttile; hold on;
for g = 1:nGeno
    geno = unique_genos{g};
    means = nan(1,nCond);
    
    for c = 1:nCond
        cond = conditions{c};
        all_mi = [];
        
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PAC')
                    pac_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PAC;
                    
                    if ~isempty(pac_data.MI)
                        pf = pac_data.phase_freqs;
                        af = pac_data.amp_freqs;
                        ip = pf >= phase_roi(1) & pf <= phase_roi(2);
                        ia = af >= amp_roi(1) & af <= amp_roi(2);
                        
                        if any(ip) && any(ia)
                            mi_val = mean(pac_data.MI(ip, ia), 'all', 'omitnan');
                            if ~isnan(mi_val)
                                all_mi = [all_mi; mi_val];
                            end
                        end
                    end
                end
            end
        end
        
        if ~isempty(all_mi)
            means(c) = mean(all_mi);
        end
    end
    
    x_pos = (1:nCond) + x_offset(g);
    bar(x_pos, means, 0.3, 'FaceColor', geno_colors(g,:), 'EdgeColor', 'k', ...
        'DisplayName', geno);
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('PAC (MI)');
title('Genotype Comparison - PAC');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

% --- Panel 8: Summary stats table ---
nexttile;
axis off;

% Create summary text
summary_text = {};
summary_text{end+1} = 'Summary Statistics:';
summary_text{end+1} = '';

for g = 1:nGeno
    geno = unique_genos{g};
    summary_text{end+1} = sprintf('--- %s ---', geno);
    
    for c = 1:nCond
        cond = conditions{c};
        all_plv = [];
        
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
                    plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
                    all_plv = [all_plv; plv_data(~isnan(plv_data))];
                end
            end
        end
        
        if ~isempty(all_plv)
            summary_text{end+1} = sprintf('%s: n=%d, PLV=%.4f±%.4f', ...
                                          cond, length(all_plv), ...
                                          mean(all_plv), std(all_plv));
        end
    end
    summary_text{end+1} = '';
end

text(0.1, 0.9, summary_text, 'VerticalAlignment', 'top', ...
     'FontName', 'Courier', 'FontSize', 9);

sgtitle('Sleep Analysis Summary: Multi-Animal Comparison', ...
        'FontSize', 16, 'FontWeight', 'bold');

fprintf('Figure 1 created successfully!\n');

%% ============================================================
%  PRINT DETAILED STATISTICS
% ============================================================

fprintf('\n=== DETAILED STATISTICS ===\n\n');

for g = 1:nGeno
    geno = unique_genos{g};
    fprintf('GENOTYPE: %s\n', geno);
    fprintf('%-10s | N Bouts | Mean PLV | SEM PLV  | Mean PAC | Valid %%\n', 'Condition');
    fprintf('-----------|---------|----------|----------|----------|---------\n');
    
    for c = 1:nCond
        cond = conditions{c};
        
        % Collect PLV data
        all_plv = [];
        n_total = 0;
        n_valid = 0;
        
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PLV')
                    plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
                    n_total = n_total + numel(plv_data);
                    n_valid = n_valid + sum(~isnan(plv_data));
                    all_plv = [all_plv; plv_data(~isnan(plv_data))];
                end
            end
        end
        
        % Collect PAC data
        all_mi = [];
        for a = 1:nAnimals
            if strcmp(genotypes{a}, geno)
                if isfield(ALL_RESULTS.(animal_ids{a}), cond) && ...
                   isfield(ALL_RESULTS.(animal_ids{a}).(cond), 'NREM_PAC')
                    pac_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PAC;
                    if ~isempty(pac_data.MI)
                        pf = pac_data.phase_freqs;
                        af = pac_data.amp_freqs;
                        ip = pf >= 2.0 & pf <= 3.0;
                        ia = af >= 11 & af <= 16;
                        if any(ip) && any(ia)
                            mi_val = mean(pac_data.MI(ip, ia), 'all', 'omitnan');
                            if ~isnan(mi_val)
                                all_mi = [all_mi; mi_val];
                            end
                        end
                    end
                end
            end
        end
        
        mean_plv = mean(all_plv);
        sem_plv = std(all_plv) / sqrt(length(all_plv));
        mean_pac = mean(all_mi);
        valid_pct = 100 * n_valid / n_total;
        
        fprintf('%-10s | %7d | %8.4f | %8.4f | %8.4f | %6.1f%%\n', ...
                cond, n_total, mean_plv, sem_plv, mean_pac, valid_pct);
    end
    fprintf('\n');
end

fprintf('=== ANALYSIS COMPLETE ===\n');
fprintf('Generated 1 figure with 8 panels\n\n');