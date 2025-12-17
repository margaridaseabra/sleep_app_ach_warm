#!/usr/bin/env python3
"""
Simple bar plots for Mean Power and Variability across conditions and time windows
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Import from main script
import sys
sys.path.append('/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA')

# ===========================
# Configuration
# ===========================

OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/sigma_theta_plots"

# Colors - progressive darkness for conditions
COLORS = {
    'WT': {
        'baseline': '#A9A9A9',    # Light grey
        'ambtemp': '#808080',     # Medium grey
        'drugs': '#505050'        # Dark grey
    },
    'APP': {
        'baseline': '#87CEEB',    # Light blue
        'ambtemp': '#6495ED',     # Medium blue (cornflower)
        'drugs': '#4169E1'        # Dark blue (royal blue)
    }
}

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9


def plot_simple_barplots(df_envelope, state='NREM', band_name='Sigma'):
    """
    Create simple bar plots:
    - Top row: Mean Power (3 panels, one per time window)
    - Bottom row: Variability (3 panels, one per time window)
    """
    df_state = df_envelope[df_envelope['state'] == state].copy()
    
    if len(df_state) == 0:
        print(f"⚠️  No {state} data found")
        return
    
    print(f"\n{'='*60}")
    print(f"BAR PLOTS: {band_name} during {state}")
    print(f"{'='*60}")
    
    # Create figure with 2 rows x 3 columns
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    condition_order = ['baseline', 'ambtemp', 'drugs']
    time_windows = ['0-3h', '3-6h', 'washout']
    
    # ===== TOP ROW: MEAN POWER =====
    for col_idx, time_window in enumerate(time_windows):
        ax = axes[0, col_idx]
        
        df_window = df_state[df_state['time_window'] == time_window]
        
        if len(df_window) == 0:
            ax.text(0.5, 0.5, f'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{time_window}', fontweight='bold')
            continue
        
        # Plot bars for each condition
        bar_positions = []
        bar_idx = 0
        bar_width = 0.35
        
        for condition in condition_order:
            for genotype in ['WT', 'APP']:
                df_group = df_window[
                    (df_window['genotype'] == genotype) &
                    (df_window['condition'] == condition)
                ]
                
                if len(df_group) > 0:
                    mean_val = np.mean(df_group['mean_power'])
                    sem_val = stats.sem(df_group['mean_power'])
                    color = COLORS[genotype][condition]
                    
                    ax.bar(bar_idx, mean_val, bar_width,
                          yerr=sem_val, color=color, alpha=0.8,
                          capsize=4, edgecolor='black', linewidth=1)
                    
                    bar_idx += 1
            
            bar_positions.append(bar_idx - 1)
            bar_idx += 0.3  # spacing between conditions
        
        # Style
        ax.set_title(f'{time_window}', fontweight='bold', fontsize=13)
        ax.set_xticks([0.5, 2.15, 3.8])
        ax.set_xticklabels(['Baseline', 'Ambtemp', 'Drugs'], rotation=0, ha='center')
        
        if col_idx == 0:
            ax.set_ylabel(f'Mean {band_name} Power (A.U.)', fontweight='bold')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add legend to first panel
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#A9A9A9', edgecolor='black', label='WT'),
        Patch(facecolor='#87CEEB', edgecolor='black', label='APP')
    ]
    axes[0, 0].legend(handles=legend_elements, loc='upper left')
    
    # ===== BOTTOM ROW: VARIABILITY =====
    for col_idx, time_window in enumerate(time_windows):
        ax = axes[1, col_idx]
        
        df_window = df_state[df_state['time_window'] == time_window]
        
        if len(df_window) == 0:
            ax.text(0.5, 0.5, f'No data', ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Plot bars for each condition
        bar_idx = 0
        bar_width = 0.35
        
        for condition in condition_order:
            for genotype in ['WT', 'APP']:
                df_group = df_window[
                    (df_window['genotype'] == genotype) &
                    (df_window['condition'] == condition)
                ]
                
                if len(df_group) > 0:
                    mean_val = np.mean(df_group['std_power'])
                    sem_val = stats.sem(df_group['std_power'])
                    color = COLORS[genotype][condition]
                    
                    ax.bar(bar_idx, mean_val, bar_width,
                          yerr=sem_val, color=color, alpha=0.8,
                          capsize=4, edgecolor='black', linewidth=1)
                    
                    bar_idx += 1
            
            bar_idx += 0.3
        
        # Style
        ax.set_xticks([0.5, 2.15, 3.8])
        ax.set_xticklabels(['Baseline', 'Ambtemp', 'Drugs'], rotation=0, ha='center')
        
        if col_idx == 0:
            ax.set_ylabel(f'{band_name} Power Variability (SD)', fontweight='bold')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add legend to first panel
    axes[1, 0].legend(handles=legend_elements, loc='upper left')
    
    # Overall title
    fig.suptitle(f'{state} - {band_name} Power Across Conditions and Time Windows',
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{band_name}_simple_barplots.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()


# ===========================
# Main - Run from original script
# ===========================

if __name__ == "__main__":
    # This should be called after running the main script
    print("\n" + "="*70)
    print("Run this after the main script to generate simple bar plots")
    print("="*70)