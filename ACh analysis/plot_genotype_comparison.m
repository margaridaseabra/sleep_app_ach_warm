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