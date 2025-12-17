#!/usr/bin/env python3
"""
Aggregate PLV/PAC results across all recordings and add metadata for analysis.

Part 1: Baseline comparison (WT vs APP) with t-tests
Part 2: Multi-condition comparison (baseline, ambtemp, drugs) with RM-ANOVA and time blocks
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ===========================
# Configuration
# ===========================

RESULTS_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/plv_pac_results"
OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/plv_pac_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Color scheme
COLORS_BASELINE = {
    'WT': '#808080',    # Grey
    'APP': '#6495ED'    # Blue
}

COLORS_CONDITIONS = {
    'WT': {
        'baseline': '#A9A9A9',    # Light grey
        'ambtemp': '#808080',     # Medium grey
        'drugs': '#505050'        # Dark grey
    },
    'APP': {
        'baseline': '#87CEEB',    # Light blue
        'ambtemp': '#6495ED',     # Medium blue
        'drugs': '#4169E1'        # Dark blue (royal blue)
    }
}

# ===========================
# Metadata Extraction
# ===========================

def parse_filename(filename):
    """Extract metadata from filename."""
    parts = filename.split('_')
    
    metadata = {
        'date': parts[0] if len(parts) > 0 else 'unknown',
        'condition': parts[1] if len(parts) > 1 else 'unknown',
        'mouse_id': parts[2] if len(parts) > 2 else 'unknown',
        'genotype': parts[3] if len(parts) > 3 else 'unknown',
    }
    
    return metadata


def assign_time_blocks(df, block_hours=[3, 6]):
    """Assign each bout to a time block based on its start time."""
    df['start_h'] = df['start_s'] / 3600.0
    
    conditions = [
        (df['start_h'] < block_hours[0]),
        (df['start_h'] >= block_hours[0]) & (df['start_h'] < block_hours[1]),
        (df['start_h'] >= block_hours[1])
    ]
    
    choices = [
        f'0-{block_hours[0]}h',
        f'{block_hours[0]}-{block_hours[1]}h',
        f'{block_hours[1]}h+'
    ]
    
    df['time_block'] = np.select(conditions, choices, default='unknown')
    
    return df


def load_and_aggregate_results(results_dir, state='NREM', condition_filter=None, use_time_blocks=False):
    """Load all CSV files for a given state and aggregate."""
    pattern = os.path.join(results_dir, f"*_{state}_all_bouts_plv_pac.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        print(f"⚠️  No result files found matching: {pattern}")
        return None
    
    print(f"Found {len(csv_files)} result files for {state}")
    
    dfs = []
    for csv_path in csv_files:
        df = pd.read_csv(csv_path)
        
        base = os.path.basename(csv_path).replace(f"_{state}_all_bouts_plv_pac.csv", "")
        metadata = parse_filename(base)
        
        for key, value in metadata.items():
            df[key] = value
        
        # Normalize condition names
        df['condition'] = df['condition'].str.lower()
        df['condition'] = df['condition'].replace({
            'amb': 'ambtemp',
            'ambtemp': 'ambtemp',
            'baseline': 'baseline',
            'drug': 'drugs',
            'drugs': 'drugs'
        })
        
        dfs.append(df)
    
    df_all = pd.concat(dfs, ignore_index=True)
    
    # Filter by condition if specified
    if condition_filter:
        df_all = df_all[df_all['condition'].isin(condition_filter)].copy()
    
    print(f"Total bouts: {len(df_all)}")
    print(f"Genotypes: {df_all['genotype'].unique()}")
    print(f"Conditions: {df_all['condition'].unique()}")
    
    if use_time_blocks:
        df_all = assign_time_blocks(df_all)
        print(f"Time blocks: {df_all['time_block'].unique()}")
    
    return df_all


def compute_per_animal_averages(df, group_by_time_block=False):
    """Compute average PLV/PAC values per animal."""
    group_cols = ['mouse_id', 'genotype', 'condition', 'state']
    if group_by_time_block and 'time_block' in df.columns:
        group_cols.append('time_block')
    
    animal_avg = df.groupby(group_cols).agg({
        'plv_eeg_ach': 'mean',
        'pac_eeg_theta_gamma': 'mean',
        'pac_eeg_phase_ach_amp': 'mean',
        'duration_s': 'mean',
        'bout_index': 'count'
    }).reset_index()
    
    animal_avg.rename(columns={'bout_index': 'n_bouts'}, inplace=True)
    
    return animal_avg


# ===========================
# PART 1: BASELINE ANALYSIS
# ===========================

def analyze_baseline(df_baseline, state='NREM'):
    """
    Part 1: Compare WT vs APP in baseline condition only.
    Uses t-tests and grey/blue color scheme.
    """
    print("\n" + "="*70)
    print("PART 1: BASELINE COMPARISON (WT vs APP)")
    print("="*70)
    
    # Compute per-animal averages
    df_animal = compute_per_animal_averages(df_baseline, group_by_time_block=False)
    
    print(f"\nBaseline animals:")
    print(f"  WT: {len(df_animal[df_animal['genotype']=='WT'])}")
    print(f"  APP: {len(df_animal[df_animal['genotype']=='APP'])}")
    
    # Statistical tests
    measures = ['plv_eeg_ach', 'pac_eeg_theta_gamma', 'pac_eeg_phase_ach_amp']
    measure_names = ['PLV (EEG-ACh)', 'PAC (EEG θ→γ)', 'PAC (EEG→ACh)']
    
    print("\n--- T-TEST RESULTS ---")
    results = {}
    for meas, name in zip(measures, measure_names):
        wt_data = df_animal[df_animal['genotype']=='WT'][meas].dropna()
        app_data = df_animal[df_animal['genotype']=='APP'][meas].dropna()
        
        if len(wt_data) >= 2 and len(app_data) >= 2:
            t_stat, p_val = stats.ttest_ind(wt_data, app_data)
            
            results[meas] = {
                'test': 't-test',
                't_statistic': t_stat,
                'p_value': p_val,
                'wt_mean': wt_data.mean(),
                'wt_std': wt_data.std(),
                'wt_n': len(wt_data),
                'app_mean': app_data.mean(),
                'app_std': app_data.std(),
                'app_n': len(app_data)
            }
            
            print(f"\n{name}:")
            print(f"  WT:  {results[meas]['wt_mean']:.4f} ± {results[meas]['wt_std']:.4f} (n={results[meas]['wt_n']})")
            print(f"  APP: {results[meas]['app_mean']:.4f} ± {results[meas]['app_std']:.4f} (n={results[meas]['app_n']})")
            print(f"  t={t_stat:.3f}, p={p_val:.4f} {'***' if p_val<0.001 else '**' if p_val<0.01 else '*' if p_val<0.05 else 'ns'}")
    
    # Save statistics
    stats_df = pd.DataFrame(results).T
    stats_path = os.path.join(OUTPUT_DIR, f'{state}_baseline_statistics.csv')
    stats_df.to_csv(stats_path)
    print(f"\n✓ Saved statistics: {stats_path}")
    
    # Create plot
    plot_baseline_comparison(df_animal, results, state)
    
    return df_animal, results


def plot_baseline_comparison(df_animal, stats_results, state='NREM'):
    """Create baseline comparison plot (WT vs APP)."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    measures = [
        ('plv_eeg_ach', 'Phase-Locking Value\n(PLV)'),
        ('pac_eeg_theta_gamma', 'Phase-Amplitude Coupling\n(EEG θ→γ)'),
        ('pac_eeg_phase_ach_amp', 'Cross-Signal PAC\n(EEG θ→ACh)')
    ]
    
    for ax, (meas, label) in zip(axes, measures):
        # Bar plot
        genotypes = ['WT', 'APP']
        means = [stats_results[meas]['wt_mean'], stats_results[meas]['app_mean']]
        stds = [stats_results[meas]['wt_std'], stats_results[meas]['app_std']]
        colors = [COLORS_BASELINE['WT'], COLORS_BASELINE['APP']]
        
        bars = ax.bar(genotypes, means, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax.errorbar(genotypes, means, yerr=stds, fmt='none', ecolor='black', capsize=5, linewidth=1.5)
        
        # Add individual points
        for i, geno in enumerate(genotypes):
            data = df_animal[df_animal['genotype']==geno][meas].values
            x = np.random.normal(i, 0.04, len(data))
            ax.scatter(x, data, color='black', alpha=0.6, s=30, zorder=3)
        
        # Statistics annotation
        p_val = stats_results[meas]['p_value']
        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
        
        y_max = max(means[0] + stds[0], means[1] + stds[1])
        y_pos = y_max * 1.15
        
        ax.plot([0, 1], [y_pos, y_pos], 'k-', linewidth=1.5)
        ax.text(0.5, y_pos * 1.02, sig, ha='center', va='bottom', fontsize=12, fontweight='bold')
        ax.text(0.5, y_pos * 1.08, f'p={p_val:.4f}', ha='center', va='bottom', fontsize=9)
        
        ax.set_ylabel(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Genotype', fontsize=11, fontweight='bold')
        ax.set_ylim(bottom=0, top=y_pos * 1.15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    fig.suptitle(f'{state} State - Baseline Comparison: WT vs APP\n'
                f'Independent t-tests', 
                fontsize=14, fontweight='bold')
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_baseline_comparison.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


# ===========================
# PART 2: MULTI-CONDITION ANALYSIS
# ===========================

def analyze_conditions_with_timeblocks(df_multi, state='NREM'):
    """
    Part 2: Compare baseline, ambtemp, drugs across genotypes and time blocks.
    Uses RM-ANOVA and progressive color darkening.
    """
    print("\n" + "="*70)
    print("PART 2: MULTI-CONDITION COMPARISON WITH TIME BLOCKS")
    print("="*70)
    
    # Compute per-animal averages per time block
    df_animal = compute_per_animal_averages(df_multi, group_by_time_block=True)
    
    print(f"\nData summary:")
    print(df_animal.groupby(['genotype', 'condition', 'time_block']).size())
    
    # Perform RM-ANOVA
    try:
        from statsmodels.formula.api import ols
        from statsmodels.stats.anova import anova_lm
        
        measures = ['plv_eeg_ach', 'pac_eeg_theta_gamma', 'pac_eeg_phase_ach_amp']
        measure_names = ['PLV (EEG-ACh)', 'PAC (EEG θ→γ)', 'PAC (EEG→ACh)']
        
        anova_results = {}
        
        print("\n--- RM-ANOVA RESULTS ---")
        print("Design: Genotype (between) × Condition × Time Block (within)\n")
        
        for meas, name in zip(measures, measure_names):
            print(f"\n{name}:")
            
            # Three-way ANOVA
            formula = f"{meas} ~ C(genotype) * C(condition) * C(time_block)"
            model = ols(formula, data=df_animal).fit()
            anova_table = anova_lm(model, typ=2)
            
            print(anova_table)
            
            anova_results[meas] = {
                'table': anova_table,
                'genotype_p': anova_table.loc['C(genotype)', 'PR(>F)'],
                'condition_p': anova_table.loc['C(condition)', 'PR(>F)'],
                'time_block_p': anova_table.loc['C(time_block)', 'PR(>F)'],
                'geno_x_cond_p': anova_table.loc['C(genotype):C(condition)', 'PR(>F)'],
                'geno_x_time_p': anova_table.loc['C(genotype):C(time_block)', 'PR(>F)'],
                'cond_x_time_p': anova_table.loc['C(condition):C(time_block)', 'PR(>F)'],
                'three_way_p': anova_table.loc['C(genotype):C(condition):C(time_block)', 'PR(>F)']
            }
            
            print(f"\nMain effects:")
            print(f"  Genotype:    p={anova_results[meas]['genotype_p']:.6f} {'***' if anova_results[meas]['genotype_p']<0.001 else '**' if anova_results[meas]['genotype_p']<0.01 else '*' if anova_results[meas]['genotype_p']<0.05 else 'ns'}")
            print(f"  Condition:   p={anova_results[meas]['condition_p']:.6f} {'***' if anova_results[meas]['condition_p']<0.001 else '**' if anova_results[meas]['condition_p']<0.01 else '*' if anova_results[meas]['condition_p']<0.05 else 'ns'}")
            print(f"  Time Block:  p={anova_results[meas]['time_block_p']:.6f} {'***' if anova_results[meas]['time_block_p']<0.001 else '**' if anova_results[meas]['time_block_p']<0.01 else '*' if anova_results[meas]['time_block_p']<0.05 else 'ns'}")
            print(f"\nInteractions:")
            print(f"  G × C: p={anova_results[meas]['geno_x_cond_p']:.6f} {'***' if anova_results[meas]['geno_x_cond_p']<0.001 else '**' if anova_results[meas]['geno_x_cond_p']<0.01 else '*' if anova_results[meas]['geno_x_cond_p']<0.05 else 'ns'}")
            print(f"  G × T: p={anova_results[meas]['geno_x_time_p']:.6f} {'***' if anova_results[meas]['geno_x_time_p']<0.001 else '**' if anova_results[meas]['geno_x_time_p']<0.01 else '*' if anova_results[meas]['geno_x_time_p']<0.05 else 'ns'}")
            print(f"  C × T: p={anova_results[meas]['cond_x_time_p']:.6f} {'***' if anova_results[meas]['cond_x_time_p']<0.001 else '**' if anova_results[meas]['cond_x_time_p']<0.01 else '*' if anova_results[meas]['cond_x_time_p']<0.05 else 'ns'}")
            print(f"  G × C × T: p={anova_results[meas]['three_way_p']:.6f} {'***' if anova_results[meas]['three_way_p']<0.001 else '**' if anova_results[meas]['three_way_p']<0.01 else '*' if anova_results[meas]['three_way_p']<0.05 else 'ns'}")
        
        # Save ANOVA results
        anova_summary_path = os.path.join(OUTPUT_DIR, f'{state}_multicondition_ANOVA.txt')
        with open(anova_summary_path, 'w') as f:
            f.write(f"RM-ANOVA Results - {state} State\n")
            f.write(f"Design: Genotype × Condition × Time Block\n")
            f.write(f"{'='*70}\n\n")
            
            for meas, name in zip(measures, measure_names):
                f.write(f"\n{name}:\n")
                f.write(f"{anova_results[meas]['table']}\n\n")
        
        print(f"\n✓ Saved ANOVA summary: {anova_summary_path}")
        
    except ImportError:
        print("⚠️  statsmodels not installed. Run: pip install statsmodels")
        anova_results = None
    except Exception as e:
        print(f"⚠️  ANOVA failed: {e}")
        anova_results = None
    
    # Create plots
    plot_conditions_by_timeblock(df_animal, anova_results, state)
    plot_interaction_heatmaps(df_animal, state)
    
    return df_animal, anova_results


def plot_conditions_by_timeblock(df_animal, anova_results, state='NREM'):
    """Create bar plots for each time block showing all conditions."""
    measures = [
        ('plv_eeg_ach', 'Phase-Locking Value (PLV)'),
        ('pac_eeg_theta_gamma', 'PAC (EEG θ→γ)'),
        ('pac_eeg_phase_ach_amp', 'PAC (EEG θ→ACh)')
    ]
    
    time_blocks = sorted(df_animal['time_block'].unique())
    conditions = ['baseline', 'ambtemp', 'drugs']
    
    for meas, label in measures:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        for ax, time_block in zip(axes, time_blocks):
            df_block = df_animal[df_animal['time_block'] == time_block]
            
            # Plot bars
            x_pos = 0
            x_labels = []
            x_positions = []
            
            for condition in conditions:
                for genotype in ['WT', 'APP']:
                    df_group = df_block[
                        (df_block['genotype'] == genotype) &
                        (df_block['condition'] == condition)
                    ]
                    
                    if len(df_group) > 0:
                        mean_val = df_group[meas].mean()
                        sem_val = stats.sem(df_group[meas])
                        color = COLORS_CONDITIONS[genotype][condition]
                        
                        ax.bar(x_pos, mean_val, 0.35, yerr=sem_val,
                              color=color, alpha=0.8, capsize=4,
                              edgecolor='black', linewidth=1)
                        
                        # Add points
                        x_jitter = np.random.normal(x_pos, 0.02, len(df_group))
                        ax.scatter(x_jitter, df_group[meas].values,
                                 color='black', alpha=0.5, s=20, zorder=3)
                        
                        x_pos += 1
                
                x_labels.append(condition.capitalize())
                x_positions.append(x_pos - 1.5)
                x_pos += 0.3
            
            ax.set_xticks(x_positions)
            ax.set_xticklabels(x_labels)
            ax.set_ylabel(label if time_block == time_blocks[0] else '', fontsize=11, fontweight='bold')
            ax.set_xlabel('Condition', fontsize=11, fontweight='bold')
            ax.set_title(f'{time_block}', fontsize=12, fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylim(bottom=0)
            ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#A9A9A9', edgecolor='black', label='WT'),
            Patch(facecolor='#87CEEB', edgecolor='black', label='APP')
        ]
        axes[0].legend(handles=legend_elements, loc='upper left')
        
        # Add ANOVA results
        if anova_results and meas in anova_results:
            anova_text = f"RM-ANOVA: "
            anova_text += f"G: {'***' if anova_results[meas]['genotype_p']<0.001 else '**' if anova_results[meas]['genotype_p']<0.01 else '*' if anova_results[meas]['genotype_p']<0.05 else 'ns'}, "
            anova_text += f"C: {'***' if anova_results[meas]['condition_p']<0.001 else '**' if anova_results[meas]['condition_p']<0.01 else '*' if anova_results[meas]['condition_p']<0.05 else 'ns'}, "
            anova_text += f"T: {'***' if anova_results[meas]['time_block_p']<0.001 else '**' if anova_results[meas]['time_block_p']<0.01 else '*' if anova_results[meas]['time_block_p']<0.05 else 'ns'}, "
            anova_text += f"G×C: {'***' if anova_results[meas]['geno_x_cond_p']<0.001 else '**' if anova_results[meas]['geno_x_cond_p']<0.01 else '*' if anova_results[meas]['geno_x_cond_p']<0.05 else 'ns'}"
            
            fig.text(0.5, 0.02, anova_text, ha='center', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        fig.suptitle(f'{state} State - {label} Across Conditions and Time Blocks',
                    fontsize=14, fontweight='bold')
        fig.tight_layout(rect=[0, 0.05, 1, 0.98])
        
        out_path = os.path.join(OUTPUT_DIR, f'{state}_multicondition_{meas}.png')
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {out_path}")
        plt.close()


def plot_interaction_heatmaps(df_animal, state='NREM'):
    """Create heatmaps showing mean values across conditions and time blocks."""
    measures = [
        ('plv_eeg_ach', 'PLV (EEG-ACh)'),
        ('pac_eeg_theta_gamma', 'PAC (EEG θ→γ)'),
        ('pac_eeg_phase_ach_amp', 'PAC (EEG θ→ACh)')
    ]
    
    for meas, label in measures:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        for ax, genotype in zip(axes, ['WT', 'APP']):
            df_geno = df_animal[df_animal['genotype'] == genotype]
            
            # Pivot table for heatmap
            pivot = df_geno.pivot_table(
                values=meas,
                index='condition',
                columns='time_block',
                aggfunc='mean'
            )
            
            # Reorder
            pivot = pivot.reindex(['baseline', 'ambtemp', 'drugs'])
            pivot = pivot[sorted(pivot.columns)]
            
            # Plot heatmap
            sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlBu_r',
                       ax=ax, cbar_kws={'label': label})
            ax.set_title(f'{genotype}', fontsize=13, fontweight='bold')
            ax.set_ylabel('Condition', fontsize=11, fontweight='bold')
            ax.set_xlabel('Time Block', fontsize=11, fontweight='bold')
        
        fig.suptitle(f'{state} State - {label} Heatmap',
                    fontsize=14, fontweight='bold')
        fig.tight_layout()
        
        out_path = os.path.join(OUTPUT_DIR, f'{state}_heatmap_{meas}.png')
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {out_path}")
        plt.close()


# ===========================
# Main Analysis
# ===========================

def main():
    print("="*70)
    print("PLV/PAC TWO-PART ANALYSIS")
    print("="*70)
    
    # Choose state
    print("\nAvailable states:")
    print("  NREM (or 1)")
    print("  REM  (or 2)")
    print("  Wake (or 0)")
    state_input = input("Enter state to analyze (default NREM): ").strip() or "NREM"
    
    state_map = {
        '0': 'Wake',
        '1': 'NREM',
        '2': 'REM',
        'wake': 'Wake',
        'nrem': 'NREM',
        'rem': 'REM'
    }
    
    state = state_map.get(state_input.lower(), state_input)
    
    # ===== PART 1: BASELINE ANALYSIS =====
    print("\n" + "="*70)
    print("Loading baseline data...")
    print("="*70)
    
    df_baseline = load_and_aggregate_results(
        RESULTS_DIR,
        state=state,
        condition_filter=['baseline'],
        use_time_blocks=False
    )
    
    if df_baseline is not None and len(df_baseline) > 0:
        df_baseline_animal, baseline_stats = analyze_baseline(df_baseline, state)
        
        # Save baseline data
        baseline_path = os.path.join(OUTPUT_DIR, f'{state}_baseline_per_animal.csv')
        df_baseline_animal.to_csv(baseline_path, index=False)
        print(f"✓ Saved baseline data: {baseline_path}")
    else:
        print("⚠️  No baseline data found!")
    
    # ===== PART 2: MULTI-CONDITION ANALYSIS =====
    print("\n" + "="*70)
    print("Loading multi-condition data...")
    print("="*70)
    
    df_multi = load_and_aggregate_results(
        RESULTS_DIR,
        state=state,
        condition_filter=['baseline', 'ambtemp', 'drugs'],
        use_time_blocks=True
    )
    
    if df_multi is not None and len(df_multi) > 0:
        df_multi_animal, anova_stats = analyze_conditions_with_timeblocks(df_multi, state)
        
        # Save multi-condition data
        multi_path = os.path.join(OUTPUT_DIR, f'{state}_multicondition_per_animal.csv')
        df_multi_animal.to_csv(multi_path, index=False)
        print(f"✓ Saved multi-condition data: {multi_path}")
    else:
        print("⚠️  No multi-condition data found!")
    
    print("\n" + "="*70)
    print(f"✓ ANALYSIS COMPLETE!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()