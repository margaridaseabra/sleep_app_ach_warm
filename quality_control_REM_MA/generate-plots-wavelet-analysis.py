#!/usr/bin/env python3
"""
Generate presentation-ready plots for WT vs APP comparison.

Focus:
1. Baseline comparison (WT vs APP, all time points)
2. Condition comparison (WT vs APP, matched time blocks only)
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import combinations

# ===========================
# Configuration
# ===========================

RESULTS_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/plv_pac_results"
OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/presentation_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

# Color scheme: WT (grays), APP (blues) with increasing darkness
COLORS_CONDITION = {
    'WT': {
        'baseline': '#B0B0B0',    # Light gray
        'ambtemp': '#707070',      # Medium gray
        'drugs': '#303030'         # Dark gray
    },
    'APP': {
        'baseline': '#A8C5E6',    # Light cornflower blue
        'ambtemp': '#6495ED',      # Cornflower blue (standard)
        'drugs': '#4169B3'         # Dark cornflower blue
    }
}

# For backward compatibility and baseline-only plots
COLORS_GENOTYPE = {
    'WT': COLORS_CONDITION['WT']['baseline'],
    'APP': COLORS_CONDITION['APP']['baseline'],
    'wt': COLORS_CONDITION['WT']['baseline'],
    'app': COLORS_CONDITION['APP']['baseline']
}

# ===========================
# Data Loading
# ===========================

def parse_filename(filename):
    """Extract metadata from filename."""
    parts = filename.split('_')
    return {
        'date': parts[0] if len(parts) > 0 else 'unknown',
        'condition': parts[1] if len(parts) > 1 else 'unknown',
        'mouse_id': parts[2] if len(parts) > 2 else 'unknown',
        'genotype': parts[3] if len(parts) > 3 else 'unknown',
    }


def assign_time_blocks(df, block_hours=[3, 6]):
    """Assign time blocks to bouts."""
    df['start_h'] = df['start_s'] / 3600.0
    
    conditions = [
        (df['start_h'] < block_hours[0]),
        (df['start_h'] >= block_hours[0]) & (df['start_h'] < block_hours[1]),
        (df['start_h'] >= block_hours[1])
    ]
    
    choices = [
        f'0-{block_hours[0]}h',
        f'{block_hours[0]}-{block_hours[1]}h',
        f'Washout'
    ]
    
    df['time_block'] = np.select(conditions, choices, default='unknown')
    return df


def load_data(state='NREM'):
    """Load all data for a given state."""
    pattern = os.path.join(RESULTS_DIR, f"*_{state}_all_bouts_plv_pac.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        print(f"⚠️  No files found for {state}")
        return None
    
    print(f"Loading {len(csv_files)} files for {state}...")
    
    dfs = []
    for csv_path in csv_files:
        df = pd.read_csv(csv_path)
        base = os.path.basename(csv_path).replace(f"_{state}_all_bouts_plv_pac.csv", "")
        metadata = parse_filename(base)
        
        for key, value in metadata.items():
            df[key] = value
        
        dfs.append(df)
    
    df_all = pd.concat(dfs, ignore_index=True)
    df_all = assign_time_blocks(df_all)
    
    # Normalize genotype names
    df_all['genotype'] = df_all['genotype'].str.upper()
    df_all['condition'] = df_all['condition'].str.lower()
    
    print(f"Total bouts: {len(df_all)}")
    print(f"Genotypes: {df_all['genotype'].unique()}")
    print(f"Conditions: {df_all['condition'].unique()}")
    print(f"Time blocks: {df_all['time_block'].unique()}")
    
    return df_all


def compute_per_animal_averages(df, group_by_time_block=False):
    """Compute per-animal averages."""
    group_cols = ['mouse_id', 'genotype', 'condition']  # Changed 'file' to 'mouse_id'
    if group_by_time_block:
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
# Statistical Functions
# ===========================

def compare_two_groups_paired(df, measure, group_col, group1, group2, subject_col='mouse_id'):
    """
    Perform paired t-test for repeated measures design.
    
    Parameters
    ----------
    df : DataFrame
        Data with measurements
    measure : str
        Column name of dependent variable
    group_col : str
        Column name defining groups (e.g., 'genotype', 'time_block')
    group1, group2 : str
        Names of the two groups to compare
    subject_col : str
        Column identifying individual subjects (animals)
    
    Returns
    -------
    dict : Test results
    """
    # Get data for each group
    df1 = df[df[group_col] == group1].set_index(subject_col)[measure]
    df2 = df[df[group_col] == group2].set_index(subject_col)[measure]
    
    # Find common subjects (animals present in both groups)
    common_subjects = df1.index.intersection(df2.index)
    
    if len(common_subjects) < 3:
        print(f"⚠️  Only {len(common_subjects)} paired subjects, cannot perform paired t-test")
        return None
    
    # Extract paired data
    data1 = df1.loc[common_subjects]
    data2 = df2.loc[common_subjects]
    
    # Paired t-test
    stat, p_val = stats.ttest_rel(data1, data2)
    
    # Effect size (Cohen's d for paired samples)
    diff = data1 - data2
    cohens_d = diff.mean() / diff.std() if diff.std() > 0 else 0
    
    return {
        'test': 'Paired t-test',
        'statistic': stat,
        'p_value': p_val,
        'mean1': data1.mean(),
        'mean2': data2.mean(),
        'std1': data1.std(),
        'std2': data2.std(),
        'mean_diff': diff.mean(),
        'std_diff': diff.std(),
        'n_pairs': len(common_subjects),
        'cohens_d': cohens_d
    }


def compare_two_groups_unpaired(data1, data2):
    """Unpaired t-test (for baseline where animals may differ)."""
    stat, p_val = stats.ttest_ind(data1, data2)
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt((data1.std()**2 + data2.std()**2) / 2)
    cohens_d = (data1.mean() - data2.mean()) / pooled_std if pooled_std > 0 else 0
    
    return {
        'test': 'Independent t-test',
        'statistic': stat,
        'p_value': p_val,
        'mean1': data1.mean(),
        'mean2': data2.mean(),
        'std1': data1.std(),
        'std2': data2.std(),
        'n1': len(data1),
        'n2': len(data2),
        'cohens_d': cohens_d
    }


def three_way_repeated_measures_anova(df, measure='plv_eeg_ach', subject_col='mouse_id'):
    """
    Perform three-way repeated measures ANOVA: genotype (between) × condition × time_block (within).
    
    Note: This is a mixed design:
    - Between-subjects factor: genotype (WT vs APP)
    - Within-subjects factors: condition, time_block
    """
    try:
        from statsmodels.stats.anova import AnovaRM
    except ImportError:
        print("⚠️  statsmodels not installed. Run: pip install statsmodels")
        return None
    
    # Check if we have all factors
    n_genotypes = len(df['genotype'].unique())
    n_conditions = len(df['condition'].unique())
    n_blocks = len(df['time_block'].unique())
    
    print(f"  Found: {n_genotypes} genotypes, {n_conditions} conditions, {n_blocks} time blocks")
    
    if n_genotypes < 2 or n_conditions < 2 or n_blocks < 2:
        print(f"⚠️  Need at least 2 levels per factor")
        return None
    
    # Check for repeated measures
    n_obs_per_subject = df.groupby(subject_col).size()
    expected_obs = n_conditions * n_blocks
    complete_subjects = n_obs_per_subject[n_obs_per_subject == expected_obs].index
    
    if len(complete_subjects) < 3:
        print(f"⚠️  Only {len(complete_subjects)} subjects with complete data")
        return None
    
    print(f"  Using {len(complete_subjects)} subjects with complete data across all conditions/blocks")
    
    # Filter to complete subjects only
    df_complete = df[df[subject_col].isin(complete_subjects)].copy()
    
    # Create combined within-subjects factor (condition × time_block)
    df_complete['cond_time'] = df_complete['condition'] + '_' + df_complete['time_block']
    
    try:
        # Perform repeated measures ANOVA with one within-subject factor
        # We'll do separate analyses for main effects and interactions
        
        # Main effect of condition (collapsed across time blocks)
        df_cond = df_complete.groupby([subject_col, 'genotype', 'condition'])[measure].mean().reset_index()
        model_cond = AnovaRM(df_cond, measure, subject_col, within=['condition'], 
                            aggregate_func='mean')
        result_cond = model_cond.fit()
        
        # Main effect of time block (collapsed across conditions)
        df_time = df_complete.groupby([subject_col, 'genotype', 'time_block'])[measure].mean().reset_index()
        model_time = AnovaRM(df_time, measure, subject_col, within=['time_block'],
                            aggregate_func='mean')
        result_time = model_time.fit()
        
        # Interaction: condition × time_block (full model)
        model_full = AnovaRM(df_complete, measure, subject_col, within=['cond_time'],
                            aggregate_func='mean')
        result_full = model_full.fit()
        
        # Extract results
        result = {
            'test': 'Repeated Measures ANOVA (within-subjects)',
            'n_subjects': len(complete_subjects),
            'condition_table': result_cond.anova_table,
            'time_block_table': result_time.anova_table,
            'full_table': result_full.anova_table,
            'condition_p': result_cond.anova_table.loc['condition', 'Pr > F'],
            'condition_F': result_cond.anova_table.loc['condition', 'F Value'],
            'time_block_p': result_time.anova_table.loc['time_block', 'Pr > F'],
            'time_block_F': result_time.anova_table.loc['time_block', 'F Value'],
        }
        
        # Between-subjects effect (genotype) - use simple comparison
        wt_means = df_complete[df_complete['genotype']=='WT'].groupby(subject_col)[measure].mean()
        app_means = df_complete[df_complete['genotype']=='APP'].groupby(subject_col)[measure].mean()
        
        if len(wt_means) > 0 and len(app_means) > 0:
            t_stat, p_val = stats.ttest_ind(wt_means, app_means)
            result['genotype_p'] = p_val
            result['genotype_t'] = t_stat
        
        return result
        
    except Exception as e:
        print(f"⚠️  Repeated measures ANOVA failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def format_rm_anova_results(result, measure_name):
    """Format repeated measures ANOVA results for display."""
    if result is None:
        return "RM-ANOVA could not be computed"
    
    def sig_marker(p):
        return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    
    text = f"Repeated Measures ANOVA: {measure_name}\n"
    text += "=" * 50 + "\n"
    text += f"n = {result['n_subjects']} subjects with complete data\n\n"
    
    if 'genotype_p' in result:
        text += "Between-Subjects Effect:\n"
        text += f"  Genotype: t={result['genotype_t']:.2f}, p={result['genotype_p']:.4f} {sig_marker(result['genotype_p'])}\n\n"
    
    text += "Within-Subjects Effects:\n"
    text += f"  Condition:  F={result['condition_F']:.2f}, p={result['condition_p']:.4f} {sig_marker(result['condition_p'])}\n"
    text += f"  Time Block: F={result['time_block_F']:.2f}, p={result['time_block_p']:.4f} {sig_marker(result['time_block_p'])}\n"
    
    return text


def add_significance_bar(ax, x1, x2, y, p_val, height_offset=0.05):
    """Add significance bar to plot."""
    sig_marker = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
    
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    bar_height = y + height_offset * y_range
    
    ax.plot([x1, x1, x2, x2], [bar_height, bar_height + 0.02*y_range, 
            bar_height + 0.02*y_range, bar_height], 'k-', linewidth=1.5)
    
    ax.text((x1 + x2) / 2, bar_height + 0.03*y_range, sig_marker,
           ha='center', va='bottom', fontsize=14, fontweight='bold')


# ===========================
# PLOT 1: Baseline Comparison (WT vs APP) - T-TEST
# ===========================

def plot_baseline_comparison(df_all, state='NREM'):
    """
    Compare WT vs APP in BASELINE condition.
    Uses baseline colors (light gray for WT, light blue for APP).
    """
    # Filter baseline only
    df_baseline = df_all[df_all['condition'] == 'baseline'].copy()
    
    if len(df_baseline) == 0:
        print("⚠️  No baseline data found!")
        return
    
    # Compute per-animal averages (no time blocking)
    df_animal = compute_per_animal_averages(df_baseline, group_by_time_block=False)
    
    # Check genotypes
    genotypes = sorted(df_animal['genotype'].unique())
    if len(genotypes) != 2 or set(genotypes) != {'WT', 'APP'}:
        print(f"⚠️  Expected WT and APP, found: {genotypes}")
        return
    
    print(f"\n=== BASELINE COMPARISON (WT vs APP - T-TEST) ===")
    
    # Check if we have paired data
    wt_animals = set(df_animal[df_animal['genotype']=='WT']['mouse_id'])
    app_animals = set(df_animal[df_animal['genotype']=='APP']['mouse_id'])
    common_animals = wt_animals.intersection(app_animals)
    
    if len(common_animals) >= 3:
        print(f"Using paired design: {len(common_animals)} animals in both WT and APP groups")
        use_paired = True
    else:
        print(f"⚠️  Only {len(common_animals)} animals in both groups")
        print(f"WT animals: {len(wt_animals)}, APP animals: {len(app_animals)}")
        print("Using unpaired t-test")
        use_paired = False
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    measures = [
        ('plv_eeg_ach', 'Phase-Locking Value (PLV)', 'EEG-ACh Coupling'),
        ('pac_eeg_theta_gamma', 'Phase-Amplitude Coupling (PAC)', 'EEG Theta → Gamma'),
        ('pac_eeg_phase_ach_amp', 'Cross-Signal PAC', 'EEG Theta → ACh Amplitude')
    ]
    
    results = []
    
    for ax, (meas, title, subtitle) in zip(axes, measures):
        if use_paired:
            # Paired t-test
            result = compare_two_groups_paired(df_animal, meas, 'genotype', 'WT', 'APP')
            if result is None:
                continue
            
            # Get paired data for plotting
            df_wt = df_animal[df_animal['genotype']=='WT'].set_index('mouse_id')[meas]
            df_app = df_animal[df_animal['genotype']=='APP'].set_index('mouse_id')[meas]
            common = df_wt.index.intersection(df_app.index)
            wt_data = df_wt.loc[common]
            app_data = df_app.loc[common]
            
        else:
            # Unpaired t-test
            wt_data = df_animal[df_animal['genotype'] == 'WT'][meas].dropna()
            app_data = df_animal[df_animal['genotype'] == 'APP'][meas].dropna()
            result = compare_two_groups_unpaired(wt_data, app_data)
        
        results.append((meas, result))
        
        # Bar plot with individual points
        x_positions = [0, 1]
        means = [result['mean1'], result['mean2']]
        
        if use_paired:
            sems = [result['std1']/np.sqrt(result['n_pairs']), 
                    result['std2']/np.sqrt(result['n_pairs'])]
        else:
            sems = [result['std1']/np.sqrt(result['n1']), 
                    result['std2']/np.sqrt(result['n2'])]
        
        # Bars with baseline colors
        colors = [COLORS_CONDITION['WT']['baseline'], COLORS_CONDITION['APP']['baseline']]
        ax.bar(x_positions, means, yerr=sems, 
              color=colors, alpha=0.85, capsize=5, width=0.6, 
              edgecolor='black', linewidth=1.5)
        
        # Individual points
        if use_paired:
            # Connect paired points with lines
            for i in range(len(wt_data)):
                ax.plot([0, 1], [wt_data.iloc[i], app_data.iloc[i]], 
                       'gray', alpha=0.3, linewidth=1)
            ax.scatter([0]*len(wt_data), wt_data, alpha=0.7, s=50, 
                      color='black', edgecolors='white', linewidths=0.5, zorder=3)
            ax.scatter([1]*len(app_data), app_data, alpha=0.7, s=50, 
                      color='black', edgecolors='white', linewidths=0.5, zorder=3)
        else:
            x_wt = np.random.normal(0, 0.04, size=len(wt_data))
            x_app = np.random.normal(1, 0.04, size=len(app_data))
            ax.scatter(x_wt, wt_data, alpha=0.7, s=50, 
                      color='black', edgecolors='white', linewidths=0.5, zorder=3)
            ax.scatter(x_app, app_data, alpha=0.7, s=50, 
                      color='black', edgecolors='white', linewidths=0.5, zorder=3)
        
        # Significance bar
        max_val = max(wt_data.max(), app_data.max())
        add_significance_bar(ax, 0, 1, max_val, result['p_value'])
        
        # Labels
        ax.set_xticks(x_positions)
        ax.set_xticklabels(['WT', 'APP'])
        ax.set_ylabel(title, fontweight='bold')
        ax.set_title(f"{title}\n{subtitle}", fontsize=12, fontweight='bold')
        
        # Stats text
        if use_paired:
            stat_text = (f"Paired t-test\n"
                        f"t({result['n_pairs']-1}) = {result['statistic']:.3f}\n"
                        f"p = {result['p_value']:.4f}\n"
                        f"Cohen's d = {result['cohens_d']:.2f}\n"
                        f"Mean diff = {result['mean_diff']:.3f}±{result['std_diff']:.3f}\n"
                        f"WT: {result['mean1']:.3f}±{result['std1']:.3f}\n"
                        f"APP: {result['mean2']:.3f}±{result['std2']:.3f}\n"
                        f"n = {result['n_pairs']} pairs")
        else:
            stat_text = (f"Independent t-test\n"
                        f"t({result['n1']+result['n2']-2}) = {result['statistic']:.3f}\n"
                        f"p = {result['p_value']:.4f}\n"
                        f"Cohen's d = {result['cohens_d']:.2f}\n"
                        f"WT: {result['mean1']:.3f}±{result['std1']:.3f} (n={result['n1']})\n"
                        f"APP: {result['mean2']:.3f}±{result['std2']:.3f} (n={result['n2']})")
        
        ax.text(0.98, 0.02, stat_text, transform=ax.transAxes,
               fontsize=8, va='bottom', ha='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    test_type = "Paired t-test" if use_paired else "Independent t-test"
    fig.suptitle(f'{state} State - BASELINE: WT vs APP Comparison ({test_type})\n'
                f'All recording time points combined',
                fontsize=16, fontweight='bold')
    fig.tight_layout()
    
    # Save
    suffix = 'paired' if use_paired else 'unpaired'
    out_path = os.path.join(OUTPUT_DIR, f'{state}_BASELINE_WT_vs_APP_{suffix}_ttest.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()
    
    # Print statistics summary
    print(f"\nStatistical Results ({test_type}):")
    for meas, result in results:
        if result:
            print(f"\n{meas}:")
            if use_paired:
                print(f"  t({result['n_pairs']-1}) = {result['statistic']:.4f}, p = {result['p_value']:.6f}")
                print(f"  Mean difference: {result['mean_diff']:.4f} ± {result['std_diff']:.4f}")
            else:
                print(f"  t({result['n1']+result['n2']-2}) = {result['statistic']:.4f}, p = {result['p_value']:.6f}")
            print(f"  Effect size (Cohen's d): {result['cohens_d']:.3f}")
            print(f"  WT:  {result['mean1']:.4f} ± {result['std1']:.4f}")
            print(f"  APP: {result['mean2']:.4f} ± {result['std2']:.4f}")


# ===========================
# PLOT 2: ALL CONDITIONS BY TIME BLOCK
# ===========================

def plot_matched_time_blocks(df_all, state='NREM'):
    """
    Compare WT vs APP across ALL conditions (baseline, ambtemp, drugs) organized by time blocks.
    Each time block shows all conditions side-by-side.
    
    Color scheme:
    - WT: Shades of gray (lighter for baseline, darker for ambtemp, darkest for drugs)
    - APP: Shades of cornflower blue (lighter for baseline, darker for ambtemp, darkest for drugs)
    """
    # Compute per-animal averages with time blocking
    df_animal = compute_per_animal_averages(df_all, group_by_time_block=True)
    
    # Define conditions and colors
    conditions = ['baseline', 'ambtemp', 'drugs']
    time_blocks = ['0-3h', '3-6h', 'Washout']
    
    # Color scheme: WT (grays), APP (blues)
    COLORS_CONDITION = {
        'WT': {
            'baseline': '#B0B0B0',    # Light gray
            'ambtemp': '#707070',      # Medium gray
            'drugs': '#303030'         # Dark gray
        },
        'APP': {
            'baseline': '#A8C5E6',    # Light cornflower blue
            'ambtemp': '#6495ED',      # Cornflower blue (standard)
            'drugs': '#4169B3'         # Dark cornflower blue
        }
    }
    
    print(f"\n=== ALL CONDITIONS BY TIME BLOCK ===")
    print(f"Conditions: {conditions}")
    print(f"Time blocks: {time_blocks}")
    
    measures = [
        ('plv_eeg_ach', 'PLV (EEG-ACh Coupling)'),
        ('pac_eeg_theta_gamma', 'PAC (EEG Theta → Gamma)'),
        ('pac_eeg_phase_ach_amp', 'Cross-PAC (EEG → ACh)')
    ]
    
    # Perform RM-ANOVA for overall effects (excluding baseline for within-subjects)
    print("\n=== Overall Repeated Measures ANOVA (ambtemp & drugs only) ===")
    df_treatment = df_all[df_all['condition'] != 'baseline'].copy()
    df_animal_treatment = compute_per_animal_averages(df_treatment, group_by_time_block=True)
    
    rm_anova_results = {}
    for meas, meas_label in measures:
        print(f"\nComputing RM-ANOVA for {meas_label}...")
        result = three_way_repeated_measures_anova(df_animal_treatment, measure=meas)
        rm_anova_results[meas] = result
        
        if result:
            print(format_rm_anova_results(result, meas_label))
    
    # Create plots
    for meas, meas_label in measures:
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        
        all_pairwise_results = []
        
        for ax_idx, time_block in enumerate(time_blocks):
            ax = axes[ax_idx]
            
            # Filter data for this time block
            df_block = df_animal[df_animal['time_block'] == time_block].copy()
            
            if len(df_block) == 0:
                ax.text(0.5, 0.5, f'No data for {time_block}',
                       ha='center', va='center', transform=ax.transAxes)
                continue
            
            # Prepare x-axis positions
            # Each condition gets 2 bars (WT, APP) with some spacing between conditions
            n_conditions = len(conditions)
            bar_width = 0.35
            group_spacing = 0.2  # Space between condition groups
            
            x_positions = []
            x_labels = []
            x_tick_positions = []
            
            current_x = 0
            for cond_idx, condition in enumerate(conditions):
                # Check if condition exists in data
                if condition not in df_block['condition'].values:
                    continue
                
                # WT bar position
                x_positions.append((condition, 'WT', current_x))
                # APP bar position
                x_positions.append((condition, 'APP', current_x + bar_width))
                
                # Label position (center of WT and APP bars)
                x_tick_positions.append(current_x + bar_width/2)
                x_labels.append(condition.capitalize())
                
                # Move to next condition group
                current_x += 2 * bar_width + group_spacing
            
            # Plot bars
            for condition, genotype, x_pos in x_positions:
                df_subset = df_block[
                    (df_block['condition'] == condition) & 
                    (df_block['genotype'] == genotype)
                ]
                
                if len(df_subset) == 0:
                    continue
                
                vals = df_subset[meas].values
                
                if len(vals) > 0:
                    # Plot bar
                    color = COLORS_CONDITION[genotype][condition]
                    ax.bar(x_pos, vals.mean(), bar_width,
                          yerr=stats.sem(vals), color=color,
                          alpha=0.85, capsize=4, edgecolor='black', linewidth=1.5,
                          label=f'{genotype} {condition}' if ax_idx == 0 else '')
                    
                    # Individual points
                    x_jitter = np.random.normal(x_pos, 0.02, size=len(vals))
                    ax.scatter(x_jitter, vals, alpha=0.7, s=50, 
                             color='black', edgecolors='white', linewidths=0.5, zorder=3)
            
            # Add significance bars for WT vs APP within each condition
            current_x = 0
            for condition in conditions:
                if condition not in df_block['condition'].values:
                    continue
                
                # Get data
                df_cond = df_block[df_block['condition'] == condition]
                wt_vals = df_cond[df_cond['genotype'] == 'WT'][meas].values
                app_vals = df_cond[df_cond['genotype'] == 'APP'][meas].values
                
                # Unpaired t-test
                if len(wt_vals) >= 2 and len(app_vals) >= 2:
                    result = compare_two_groups_unpaired(pd.Series(wt_vals), pd.Series(app_vals))
                    
                    all_pairwise_results.append({
                        'condition': condition,
                        'time_block': time_block,
                        'measure': meas,
                        'p_value': result['p_value'],
                        't_stat': result['statistic'],
                        'cohens_d': result['cohens_d'],
                        'n_wt': result['n1'],
                        'n_app': result['n2'],
                        'mean_wt': result['mean1'],
                        'mean_app': result['mean2'],
                        'std_wt': result['std1'],
                        'std_app': result['std2']
                    })
                    
                    # Add significance marker (NO CORRECTION)
                    if result['p_value'] < 0.05:
                        max_val = max(wt_vals.max() if len(wt_vals)>0 else 0, 
                                    app_vals.max() if len(app_vals)>0 else 0)
                        y_pos = max_val + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
                        sig = '***' if result['p_value'] < 0.001 else '**' if result['p_value'] < 0.01 else '*'
                        
                        # Draw bar between WT and APP for this condition
                        x1 = current_x
                        x2 = current_x + bar_width
                        ax.plot([x1, x1, x2, x2], 
                               [y_pos, y_pos + 0.015 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 
                                y_pos + 0.015 * (ax.get_ylim()[1] - ax.get_ylim()[0]), y_pos], 
                               'k-', linewidth=1.5)
                        ax.text((x1 + x2) / 2, y_pos + 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 
                               sig, ha='center', va='bottom', fontsize=12, fontweight='bold')
                
                current_x += 2 * bar_width + group_spacing
            
            # Format axes
            ax.set_xticks(x_tick_positions)
            ax.set_xticklabels(x_labels, fontsize=11, fontweight='bold')
            ax.set_ylabel(meas_label, fontweight='bold', fontsize=12)
            ax.set_title(f'{time_block}', fontsize=14, fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(axis='both', which='major', labelsize=10)
            
            # Add legend only to first subplot
            if ax_idx == 0:
                from matplotlib.patches import Patch
                legend_elements = []
                for genotype in ['WT', 'APP']:
                    for condition in conditions:
                        color = COLORS_CONDITION[genotype][condition]
                        legend_elements.append(
                            Patch(facecolor=color, label=f'{genotype} {condition.capitalize()}', 
                                 alpha=0.85, edgecolor='black')
                        )
                ax.legend(handles=legend_elements, loc='upper left', frameon=True, 
                         fontsize=9, ncol=1)
        
        # Add RM-ANOVA results as text
        rm_result = rm_anova_results.get(meas)
        if rm_result:
            def sig(p):
                return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            
            anova_text = f"RM-ANOVA (ambtemp & drugs, n={rm_result['n_subjects']} subjects): "
            if 'genotype_p' in rm_result:
                anova_text += f"Genotype (between): {sig(rm_result['genotype_p'])} | "
            anova_text += f"Condition (within): {sig(rm_result['condition_p'])} | "
            anova_text += f"Time (within): {sig(rm_result['time_block_p'])}\n"
            anova_text += f"Pairwise comparisons: Independent t-tests (uncorrected)"
            
            fig.text(0.5, 0.02, anova_text, ha='center', va='bottom',
                    fontsize=10, bbox=dict(boxstyle='round', facecolor='lightyellow', 
                                         alpha=0.8, edgecolor='gray'))
        
        fig.suptitle(f'{state} State - {meas_label}\n'
                    f'WT vs APP across all conditions (organized by time block)',
                    fontsize=16, fontweight='bold')
        fig.tight_layout(rect=[0, 0.06, 1, 1])
        
        # Save
        meas_short = meas.replace('plv_eeg_ach', 'PLV').replace('pac_eeg_theta_gamma', 'PAC_EEG').replace('pac_eeg_phase_ach_amp', 'PAC_Cross')
        out_path = os.path.join(OUTPUT_DIR, f'{state}_ALL_CONDITIONS_{meas_short}.png')
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {out_path}")
        plt.close()
        
        # Save pairwise results to CSV
        if all_pairwise_results:
            csv_path = os.path.join(OUTPUT_DIR, f'{state}_{meas}_pairwise_ttests_all_conditions.csv')
            df_results = pd.DataFrame(all_pairwise_results)
            df_results = df_results.sort_values(['time_block', 'condition'])
            df_results.to_csv(csv_path, index=False)
            print(f"✓ Saved pairwise results: {csv_path}")
            
            # Print significant results
            sig_results = df_results[df_results['p_value'] < 0.05]
            if len(sig_results) > 0:
                print(f"\n  Significant differences (p < 0.05):")
                for _, row in sig_results.iterrows():
                    print(f"    {row['condition']} - {row['time_block']}: "
                          f"t({row['n_wt']+row['n_app']-2})={row['t_stat']:.2f}, "
                          f"p={row['p_value']:.4f}, d={row['cohens_d']:.2f}")
    
    # Save RM-ANOVA tables
    for meas, result in rm_anova_results.items():
        if result and 'condition_table' in result:
            csv_path = os.path.join(OUTPUT_DIR, f'{state}_{meas}_RM_ANOVA_tables.csv')
            with open(csv_path, 'w') as f:
                f.write("=== Repeated Measures ANOVA Results ===\n")
                f.write("Note: Baseline excluded (between-subjects comparison)\n")
                f.write(f"n = {result['n_subjects']} subjects with complete data\n\n")
                
                if 'genotype_p' in result:
                    f.write(f"Between-Subjects Factor (Genotype):\n")
                    f.write(f"  t = {result['genotype_t']:.4f}\n")
                    f.write(f"  p = {result['genotype_p']:.6f}\n\n")
                
                f.write("Within-Subjects Factors:\n\n")
                f.write("Condition Effect (ambtemp vs drugs):\n")
                result['condition_table'].to_csv(f)
                f.write("\n\nTime Block Effect:\n")
                result['time_block_table'].to_csv(f)
            print(f"✓ Saved RM-ANOVA tables: {csv_path}")


# ===========================
# PLOT 3: Time Course Summary (All Genotypes/Conditions)
# ===========================

def plot_time_course_summary(df_all, state='NREM'):
    """
    Show time course of coupling measures across all conditions.
    Uses consistent color scheme: grays for WT, blues for APP, darker for drugs.
    """
    df_animal = compute_per_animal_averages(df_all, group_by_time_block=True)
    
    measures = [
        ('plv_eeg_ach', 'PLV (EEG-ACh)'),
        ('pac_eeg_theta_gamma', 'PAC (EEG θ→γ)'),
        ('pac_eeg_phase_ach_amp', 'PAC (EEG→ACh)')
    ]
    
    conditions = sorted(df_animal['condition'].unique())
    time_blocks = ['0-3h', '3-6h', 'Washout']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for ax, (meas, meas_label) in zip(axes, measures):
        for condition in conditions:
            for genotype in ['WT', 'APP']:
                df_subset = df_animal[
                    (df_animal['condition'] == condition) &
                    (df_animal['genotype'] == genotype)
                ]
                
                if len(df_subset) == 0:
                    continue
                
                # Compute means per time block
                means = []
                sems = []
                for tb in time_blocks:
                    data = df_subset[df_subset['time_block'] == tb][meas]
                    if len(data) > 0:
                        means.append(data.mean())
                        sems.append(data.sem())
                    else:
                        means.append(np.nan)
                        sems.append(np.nan)
                
                # Get color for this genotype/condition combination
                color = COLORS_CONDITION[genotype][condition]
                
                # Line style varies by condition for clarity
                linestyle = '-' if condition == 'baseline' else '--' if condition == 'ambtemp' else ':'
                label = f'{genotype} {condition.capitalize()}'
                
                x = np.arange(len(time_blocks))
                ax.errorbar(x, means, yerr=sems, marker='o', label=label,
                           color=color, linestyle=linestyle, linewidth=2.5,
                           markersize=8, capsize=4, alpha=0.9)
        
        ax.set_xticks(range(len(time_blocks)))
        ax.set_xticklabels(time_blocks)
        ax.set_xlabel('Time Block', fontweight='bold')
        ax.set_ylabel(meas_label, fontweight='bold')
        ax.set_title(meas_label, fontsize=13, fontweight='bold')
        ax.legend(fontsize=8, loc='best', frameon=True, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    fig.suptitle(f'{state} State - Time Course Summary\n'
                f'Changes in coupling across recording duration',
                fontsize=16, fontweight='bold')
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_TIME_COURSE_SUMMARY.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


def plot_animal_trajectories(df_all, state='NREM', measure='plv_eeg_ach'):
    """
    Plot individual animal trajectories across time blocks.
    Uses condition-specific colors with appropriate shading.
    """
    df_animal = compute_per_animal_averages(df_all, group_by_time_block=True)
    
    conditions = sorted(df_animal['condition'].unique())
    time_blocks = ['0-3h', '3-6h', 'Washout']
    genotypes = sorted(df_animal['genotype'].unique())
    
    measure_labels = {
        'plv_eeg_ach': 'PLV (EEG-ACh Coupling)',
        'pac_eeg_theta_gamma': 'PAC (EEG Theta → Gamma)',
        'pac_eeg_phase_ach_amp': 'Cross-PAC (EEG → ACh)'
    }
    
    fig, axes = plt.subplots(len(genotypes), len(conditions), 
                            figsize=(6*len(conditions), 5*len(genotypes)))
    
    if len(genotypes) == 1:
        axes = axes.reshape(1, -1)
    if len(conditions) == 1:
        axes = axes.reshape(-1, 1)
    
    for i, genotype in enumerate(genotypes):
        for j, condition in enumerate(conditions):
            ax = axes[i, j]
            
            df_subset = df_animal[
                (df_animal['genotype'] == genotype) & 
                (df_animal['condition'] == condition)
            ]
            
            if len(df_subset) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                       transform=ax.transAxes)
                ax.set_title(f'{genotype} - {condition.capitalize()}')
                continue
            
            # Get base color for this genotype/condition
            base_color = COLORS_CONDITION[genotype][condition]
            
            # Plot each animal's trajectory
            animals = sorted(df_subset['mouse_id'].unique())
            
            # Create color variations around base color for different animals
            from matplotlib.colors import to_rgb, to_hex
            import colorsys
            
            base_rgb = to_rgb(base_color)
            base_hsv = colorsys.rgb_to_hsv(*base_rgb)
            
            colors = []
            for idx in range(len(animals)):
                # Vary the value (brightness) slightly for each animal
                offset = (idx - len(animals)/2) * 0.1 / len(animals)
                new_v = max(0.2, min(1.0, base_hsv[2] + offset))
                new_rgb = colorsys.hsv_to_rgb(base_hsv[0], base_hsv[1], new_v)
                colors.append(new_rgb)
            
            for animal, color in zip(animals, colors):
                df_animal_data = df_subset[df_subset['mouse_id'] == animal].copy()
                
                # Sort by time block
                df_animal_data['tb_idx'] = df_animal_data['time_block'].apply(
                    lambda x: time_blocks.index(x) if x in time_blocks else -1
                )
                df_animal_data = df_animal_data.sort_values('tb_idx')
                
                # Plot line
                x = df_animal_data['tb_idx'].values
                y = df_animal_data[measure].values
                
                ax.plot(x, y, marker='o', label=animal, linewidth=2.5, 
                       markersize=9, color=color, alpha=0.85,
                       markeredgecolor='white', markeredgewidth=0.5)
            
            # Formatting
            ax.set_xticks(range(len(time_blocks)))
            ax.set_xticklabels(time_blocks, rotation=15)
            ax.set_xlabel('Time Block', fontweight='bold')
            ax.set_ylabel(measure_labels.get(measure, measure), fontweight='bold')
            ax.set_title(f'{genotype} - {condition.capitalize()}', 
                        fontsize=12, fontweight='bold',
                        color=base_color)
            ax.legend(fontsize=8, loc='best', frameon=True, title='Mouse ID')
            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
    
    fig.suptitle(f'{state} State - Individual Animal Trajectories\n'
                f'{measure_labels.get(measure, measure)}',
                fontsize=16, fontweight='bold')
    fig.tight_layout()
    
    meas_short = measure.replace('plv_eeg_ach', 'PLV').replace('pac_eeg_theta_gamma', 'PAC_EEG').replace('pac_eeg_phase_ach_amp', 'PAC_Cross')
    out_path = os.path.join(OUTPUT_DIR, f'{state}_animal_trajectories_{meas_short}.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


# ===========================
# DATA AVAILABILITY PLOTS
# ===========================

def plot_data_availability(df_all, state='NREM'):
    """
    Visualize which animals have data for each condition and time block.
    Shows a heatmap of data availability.
    """
    print(f"\n=== DATA AVAILABILITY SUMMARY ===")
    
    # Get per-animal data with time blocks
    df_animal = compute_per_animal_averages(df_all, group_by_time_block=True)
    
    # Create summary
    animals = sorted(df_animal['mouse_id'].unique())
    conditions = sorted(df_animal['condition'].unique())
    time_blocks = ['0-3h', '3-6h', 'Washout']
    genotypes = sorted(df_animal['genotype'].unique())
    
    print(f"\nTotal animals: {len(animals)}")
    for geno in genotypes:
        geno_animals = df_animal[df_animal['genotype'] == geno]['mouse_id'].unique()
        print(f"  {geno}: {len(geno_animals)} animals - {list(geno_animals)}")
    
    print(f"\nConditions: {conditions}")
    print(f"Time blocks: {time_blocks}")
    
    # Create figure with subplots for each genotype
    fig, axes = plt.subplots(1, len(genotypes), figsize=(8*len(genotypes), 6))
    if len(genotypes) == 1:
        axes = [axes]
    
    for ax, genotype in zip(axes, genotypes):
        df_geno = df_animal[df_animal['genotype'] == genotype].copy()
        geno_animals = sorted(df_geno['mouse_id'].unique())
        
        # Create matrix: rows=animals, columns=condition×time_block combinations
        cond_time_combos = []
        for cond in conditions:
            for tb in time_blocks:
                cond_time_combos.append(f"{cond}\n{tb}")
        
        # Build availability matrix (1=has data, 0=no data)
        matrix = np.zeros((len(geno_animals), len(cond_time_combos)))
        
        for i, animal in enumerate(geno_animals):
            for j, (cond, tb) in enumerate([(c, t) for c in conditions for t in time_blocks]):
                has_data = len(df_geno[
                    (df_geno['mouse_id'] == animal) & 
                    (df_geno['condition'] == cond) & 
                    (df_geno['time_block'] == tb)
                ]) > 0
                matrix[i, j] = 1 if has_data else 0
        
        # Plot heatmap
        im = ax.imshow(matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
        
        # Set ticks
        ax.set_xticks(range(len(cond_time_combos)))
        ax.set_xticklabels(cond_time_combos, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(geno_animals)))
        ax.set_yticklabels(geno_animals, fontsize=9)
        
        # Add grid
        ax.set_xticks(np.arange(len(cond_time_combos)) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(geno_animals)) - 0.5, minor=True)
        ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
        
        # Add text annotations
        for i in range(len(geno_animals)):
            for j in range(len(cond_time_combos)):
                if matrix[i, j] == 1:
                    # Count number of bouts for this animal/condition/block
                    cond, tb = [(c, t) for c in conditions for t in time_blocks][j]
                    n_bouts = df_geno[
                        (df_geno['mouse_id'] == geno_animals[i]) & 
                        (df_geno['condition'] == cond) & 
                        (df_geno['time_block'] == tb)
                    ]['n_bouts'].sum()
                    ax.text(j, i, f'{int(n_bouts)}', ha='center', va='center',
                           color='black', fontsize=7, fontweight='bold')
        
        ax.set_title(f'{genotype} Animals\nData Availability (numbers = bout count)',
                    fontsize=13, fontweight='bold')
        ax.set_xlabel('Condition × Time Block', fontweight='bold')
        ax.set_ylabel('Mouse ID', fontweight='bold')
    
    # Add colorbar
    cbar = fig.colorbar(im, ax=axes, orientation='horizontal', pad=0.1, aspect=30)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['No Data', 'Has Data'])
    
    fig.suptitle(f'{state} State - Data Availability Matrix\n'
                f'Green = data available, Red = missing data, Numbers = bout count',
                fontsize=16, fontweight='bold')
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_data_availability.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()
    
    # Print detailed summary
    print("\n=== DETAILED ANIMAL SUMMARY ===")
    for genotype in genotypes:
        print(f"\n{genotype} Animals:")
        df_geno = df_animal[df_animal['genotype'] == genotype]
        for animal in sorted(df_geno['mouse_id'].unique()):
            df_animal_data = df_geno[df_geno['mouse_id'] == animal]
            conds_present = sorted(df_animal_data['condition'].unique())
            blocks_present = sorted(df_animal_data['time_block'].unique())
            total_bouts = df_animal_data['n_bouts'].sum()
            print(f"  {animal}:")
            print(f"    Conditions: {conds_present}")
            print(f"    Time blocks: {blocks_present}")
            print(f"    Total bouts: {int(total_bouts)}")
    
    # Check which animals have complete data
    print("\n=== COMPLETE DATA CHECK (for RM-ANOVA) ===")
    n_conditions = len(df_animal['condition'].unique())
    n_blocks = len(time_blocks)
    expected_combos = n_conditions * n_blocks
    
    for genotype in genotypes:
        df_geno = df_animal[df_animal['genotype'] == genotype]
        print(f"\n{genotype}:")
        for animal in sorted(df_geno['mouse_id'].unique()):
            df_animal_data = df_geno[df_geno['mouse_id'] == animal]
            n_combos = len(df_animal_data)
            complete = "✓ COMPLETE" if n_combos == expected_combos else f"✗ INCOMPLETE ({n_combos}/{expected_combos})"
            print(f"  {animal}: {complete}")
    
    # Animals with complete data across all conditions/blocks
    animal_combo_counts = df_animal.groupby('mouse_id').size()
    complete_animals = animal_combo_counts[animal_combo_counts == expected_combos].index.tolist()
    print(f"\nAnimals with COMPLETE data (all {expected_combos} condition×block combos):")
    for animal in complete_animals:
        geno = df_animal[df_animal['mouse_id'] == animal]['genotype'].iloc[0]
        print(f"  {animal} ({geno})")
    print(f"Total: {len(complete_animals)} animals")


# ===========================
# Main
# ===========================

def main():
    print("="*70)
    print("PRESENTATION PLOTS GENERATOR")
    print("Focus: WT vs APP Comparison")
    print("="*70)
    
    # Choose state
    print("\nAvailable states: NREM, REM, Wake")
    state_input = input("Enter state (default NREM): ").strip() or "NREM"
    
    state_map = {'0': 'Wake', '1': 'NREM', '2': 'REM',
                 'wake': 'Wake', 'nrem': 'NREM', 'rem': 'REM'}
    state = state_map.get(state_input.lower(), state_input)
    
    # Load data
    print(f"\nLoading {state} data...")
    df_all = load_data(state=state)
    
    if df_all is None or len(df_all) == 0:
        print("No data found!")
        return
    
    # Generate plots
    print(f"\n{'='*70}")
    print("GENERATING PLOTS")
    print(f"{'='*70}\n")
    
    print("0. Data availability check...")
    plot_data_availability(df_all, state=state)
    
    print("\n0b. Individual animal trajectories...")
    for measure in ['plv_eeg_ach', 'pac_eeg_theta_gamma', 'pac_eeg_phase_ach_amp']:
        plot_animal_trajectories(df_all, state=state, measure=measure)
    
    print("\n1. All conditions organized by time block...")
    plot_matched_time_blocks(df_all, state=state)
    
    print("\n2. Time course summary...")
    plot_time_course_summary(df_all, state=state)
    
    print(f"\n{'='*70}")
    print(f"✓ All plots saved to: {OUTPUT_DIR}")
    print(f"{'='*70}")
    
    # Summary
    print("\nGenerated files:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if f.endswith('.png'):
            print(f"  • {f}")


if __name__ == "__main__":
    main()