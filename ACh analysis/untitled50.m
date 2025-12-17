%% QUICK PLOTS - Run this now!
load('all_animals_plv_pac_results.mat')

animal_ids = {animals.id};
genotypes = {animals.genotype};
nAnimals = numel(animals);
conditions = {'baseline', 'ambtemp', 'drugs'};

%% Figure 1: NREM PLV Summary
figure('Color','w','Position',[100 100 1400 600]);
tiledlayout(2,3,'TileSpacing','compact');

% Panel A: Mean PLV by condition
nexttile; hold on;
for a = 1:nAnimals
    plv_means = nan(1,3);
    for c = 1:3
        cond = conditions{c};
        plv_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
        plv_means(c) = mean(plv_data, 'omitnan');
    end
    
    plot(1:3, plv_means, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:3, 'XTickLabel', conditions);
ylabel('Mean NREM PLV');
title('NREM PLV by Condition');
legend('Location', 'best');
grid on; box on;

% Panel B: REM PLV by condition
nexttile; hold on;
for a = 1:nAnimals
    plv_means = nan(1,3);
    for c = 1:3
        cond = conditions{c};
        plv_data = ALL_RESULTS.(animal_ids{a}).(cond).REM_PLV.plv_all;
        plv_means(c) = mean(plv_data, 'omitnan');
    end
    
    plot(1:3, plv_means, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:3, 'XTickLabel', conditions);
ylabel('Mean REM PLV');
title('REM PLV by Condition');
legend('Location', 'best');
grid on; box on;

% Panel C: NREM PAC by condition
nexttile; hold on;
phase_roi = [2.0 3.0];  % Delta (avoid 1.5 Hz that failed)
amp_roi = [11 16];      % Spindle

for a = 1:nAnimals
    mi_means = nan(1,3);
    for c = 1:3
        cond = conditions{c};
        pac_data = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PAC;
        
        if ~isempty(pac_data.MI)
            pf = pac_data.phase_freqs;
            af = pac_data.amp_freqs;
            ip = pf >= phase_roi(1) & pf <= phase_roi(2);
            ia = af >= amp_roi(1) & af <= amp_roi(2);
            mi_means(c) = mean(pac_data.MI(ip, ia), 'all', 'omitnan');
        end
    end
    
    plot(1:3, mi_means, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'DisplayName', sprintf('%s (%s)', animal_ids{a}, genotypes{a}));
end
set(gca, 'XTick', 1:3, 'XTickLabel', conditions);
ylabel('Mean Delta-Spindle PAC');
title('NREM PAC (Delta-Spindle)');
legend('Location', 'best');
grid on; box on;

% Panel D: Genotype comparison - NREM PLV
nexttile; hold on;
app_data = [];
wt_data = [];
cond_labels = [];

for c = 1:3
    cond = conditions{c};
    for a = 1:nAnimals
        plv_vals = ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all;
        plv_vals = plv_vals(~isnan(plv_vals));
        
        if strcmp(genotypes{a}, 'APP')
            app_data = [app_data; plv_vals];
            cond_labels = [cond_labels; c*ones(size(plv_vals))];
        else
            wt_data = [wt_data; plv_vals];
        end
    end
end

% Bar plot by genotype
x_pos = [1 1.3 2 2.3 3 3.3];
app_means = [mean(app_data(cond_labels==1)), mean(app_data(cond_labels==2)), mean(app_data(cond_labels==3))];
wt_means = [mean(wt_data(cond_labels==1)), mean(wt_data(cond_labels==2)), mean(wt_data(cond_labels==3))];

bar([1 2 3]-0.15, app_means, 0.3, 'FaceColor', [0.8 0.3 0.3]);
bar([1 2 3]+0.15, wt_means, 0.3, 'FaceColor', [0.3 0.3 0.8]);

set(gca, 'XTick', 1:3, 'XTickLabel', conditions);
ylabel('NREM PLV');
title('Genotype Comparison');
legend({'APP', 'WT'}, 'Location', 'best');
grid on; box on;

% Panel E: Distribution of PLV values
nexttile; hold on;
all_app = [];
all_wt = [];
for c = 1:3
    for a = 1:nAnimals
        plv_vals = ALL_RESULTS.(animal_ids{a}).(conditions{c}).NREM_PLV.plv_all;
        plv_vals = plv_vals(~isnan(plv_vals));
        if strcmp(genotypes{a}, 'APP')
            all_app = [all_app; plv_vals];
        else
            all_wt = [all_wt; plv_vals];
        end
    end
end

histogram(all_app, 'BinWidth', 0.02, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5);
histogram(all_wt, 'BinWidth', 0.02, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.5);
xlabel('NREM PLV');
ylabel('Count');
title('PLV Distribution by Genotype');
legend({'APP', 'WT'});
grid on; box on;

% Panel F: Number of bouts
nexttile; hold on;
for a = 1:nAnimals
    n_bouts = nan(1,3);
    for c = 1:3
        cond = conditions{c};
        n_bouts(c) = numel(ALL_RESULTS.(animal_ids{a}).(cond).NREM_PLV.plv_all);
    end
    plot(1:3, n_bouts, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'DisplayName', animal_ids{a});
end
set(gca, 'XTick', 1:3, 'XTickLabel', conditions);
ylabel('Number of NREM Bouts');
title('Bout Count by Condition');
legend('Location', 'best');
grid on; box on;

sgtitle('Sleep Analysis Summary: APP vs WT', 'FontSize', 16, 'FontWeight', 'bold');
```

---

## 📊 **Key Results from Your Output**

### **🧬 Genotype × Condition Interaction:**

| Metric | Genotype | Baseline | Ambtemp | Drugs | Pattern |
|--------|----------|----------|---------|-------|---------|
| **NREM PLV** | APP | 0.0585 | 0.0589 | 0.0643 | ↗ Increase |
| | WT | 0.0710 | 0.0662 | 0.0428 | ↘ **Decrease** |
| **REM PLV** | APP | 0.0041 | 0.0034 | 0.0031 | ↘ Slight decrease |
| | WT | 0.0031 | 0.0045 | 0.0019 | ↕ Variable |

### **🎯 Main Finding:**

**OPPOSITE patterns between genotypes!**

**APP Mouse:**
- NREM PLV increases from baseline (0.058) → ambtemp (0.059) → drugs (0.064)
- Stable or slightly improving coupling

**WT Mouse:**
- NREM PLV decreases from baseline (0.071) → ambtemp (0.066) → drugs (0.043)
- **40% decrease in drugs condition!**

---

## ⚠️ **Data Quality Check**

**WT Drugs condition has issues:**
```
NREM: 22/24 valid (91.7%) - 2 bouts failed
REM: 3/4 valid (75.0%) - 1 bout failed