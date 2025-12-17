#!/usr/bin/env python3
"""
Power Spectral Density (PSD) Analysis of Sigma/Theta Envelope Oscillations

Analyzes slow oscillations in sigma/theta power during natural sleep:
1. Extract sigma power over time during NREM (envelope)
2. Compute PSD of sigma power time series to find cyclic patterns
3. Compare across conditions (baseline, ambtemp, drugs) and time windows
"""

import os
import glob
import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy import stats, signal
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# ===========================
# Configuration
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"
OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/sigma_theta_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Bout selection
MIN_BOUT_LEN_SEC = 5.0  # minimum bout length in seconds

# Time windows (in seconds)
TIME_WINDOWS = {
    '0-3h': (0, 3 * 3600),
    '3-6h': (3 * 3600, 6 * 3600),
    'washout': (6 * 3600, float('inf'))
}

# Frequency bands for power extraction
SIGMA_RANGE = (11, 16)  # Hz - NREM sigma
THETA_RANGE = (4, 10)   # Hz - REM theta

# Envelope PSD parameters
WINDOW_SIZE_SEC = 4.0  # Window size for computing band power
WINDOW_STEP_SEC = 0.5  # Step size for sliding window
ENVELOPE_PSD_RANGE = (0, 0.15)  # Frequency range for envelope PSD

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

# ===========================
# Helper Functions
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


def load_mat_signal(mat_path):
    """Load EEG and sampling rate from .mat file."""
    try:
        data = sio.loadmat(mat_path, squeeze_me=True)
        eeg = np.asarray(data["eeg"], dtype=float)
        fs = float(np.squeeze(data["eeg_frequency"]))
        return eeg, fs
    except Exception as e:
        print(f"⚠️  Error loading {mat_path}: {e}")
        return None, None


def load_scores_1hz(csv_path):
    """Load per-second scoring from CSV."""
    try:
        df = pd.read_csv(csv_path)
        num_cols = list(df.select_dtypes(include="number").columns)
        
        if not num_cols:
            raise ValueError(f"No numeric columns found in {csv_path}")
        
        # Find score column
        score_col = None
        for c in num_cols:
            lc = c.lower()
            if "score" in lc or "state" in lc:
                score_col = c
                break
        
        if score_col is None:
            score_col = min(num_cols, key=lambda c: df[c].nunique())
        
        states = df[score_col].to_numpy().astype(int)
        
        # Find time column
        time_col = None
        for c in num_cols:
            lc = c.lower()
            if c == score_col:
                continue
            if "time" in lc or "sec" in lc or "second" in lc:
                time_col = c
                break
        
        if time_col is not None:
            time_s = df[time_col].to_numpy().astype(float)
        else:
            time_s = np.arange(len(states), dtype=float)
        
        if np.any(np.diff(time_s) <= 0):
            time_s = np.arange(len(states), dtype=float)
        
        return time_s, states
        
    except Exception as e:
        print(f"⚠️  Error loading {csv_path}: {e}")
        return None, None


def find_state_bouts(states, min_len_sec=MIN_BOUT_LEN_SEC):
    """Find contiguous bouts for each state."""
    bouts_by_state = {s: [] for s in np.unique(states)}
    in_bout = False
    start = None
    current_state = None
    
    for i, s in enumerate(states):
        if not in_bout:
            in_bout = True
            start = i
            current_state = s
        else:
            if s != current_state:
                end = i
                length = end - start
                if length >= min_len_sec:
                    bouts_by_state[current_state].append((start, end))
                start = i
                current_state = s
    
    if in_bout:
        end = len(states)
        length = end - start
        if length >= min_len_sec:
            bouts_by_state[current_state].append((start, end))
    
    return bouts_by_state


def get_time_window_label(start_s):
    """Determine which time window a bout belongs to."""
    for label, (window_start, window_end) in TIME_WINDOWS.items():
        if window_start <= start_s < window_end:
            return label
    return None


# ===========================
# Band Power Envelope Extraction
# ===========================

def extract_band_power_envelope(eeg_signal, srate, band_range, window_size_sec=4.0, step_sec=0.5):
    """Extract band power over time using sliding window."""
    window_samples = int(window_size_sec * srate)
    step_samples = int(step_sec * srate)
    
    n_windows = int((len(eeg_signal) - window_samples) / step_samples) + 1
    
    time_points = []
    power_envelope = []
    
    for i in range(n_windows):
        start_idx = i * step_samples
        end_idx = start_idx + window_samples
        
        if end_idx > len(eeg_signal):
            break
        
        window_data = eeg_signal[start_idx:end_idx]
        
        # Compute PSD for this window
        freqs, psd = signal.welch(
            window_data,
            fs=srate,
            nperseg=min(window_samples, len(window_data)),
            scaling='density'
        )
        
        # Extract power in band
        band_mask = (freqs >= band_range[0]) & (freqs <= band_range[1])
        band_power = np.trapz(psd[band_mask], freqs[band_mask])
        
        time_points.append(start_idx / srate + window_size_sec / 2)
        power_envelope.append(band_power)
    
    return np.array(time_points), np.array(power_envelope)


def compute_envelope_psd(power_envelope, envelope_srate, freq_range=(0, 0.15)):
    """Compute PSD of the power envelope to find slow oscillations."""
    # Detrend to remove DC component
    power_envelope_detrended = signal.detrend(power_envelope, type='linear')
    
    # Use longer windows for better low-frequency resolution
    nperseg = min(len(power_envelope_detrended), int(60 * envelope_srate))
    noverlap = nperseg // 2
    
    freqs, psd = signal.welch(
        power_envelope_detrended,
        fs=envelope_srate,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling='density',
        window='hann'
    )
    
    # Filter to specified range
    mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    return freqs[mask], psd[mask]


# ===========================
# Load All Bouts
# ===========================

def load_all_bouts():
    """Load all REM and NREM bouts from all conditions."""
    signal_files = glob.glob(os.path.join(SIGNAL_DIR, "*.mat"))
    
    if not signal_files:
        print(f"⚠️  No .mat files found in {SIGNAL_DIR}")
        return None
    
    print(f"\nFound {len(signal_files)} total .mat files")
    
    # Exclude group/combined files
    valid_files = [f for f in signal_files if not any(
        exclude in os.path.basename(f).lower() 
        for exclude in ['group', 'all', 'concat', 'combined']
    )]
    
    print(f"Processing {len(valid_files)} individual recording files...")
    
    all_bouts = []
    
    for mat_path in valid_files:
        base_name = os.path.basename(mat_path).replace('.mat', '')
        metadata = parse_filename(base_name)
        
        # Normalize condition names
        condition = metadata['condition'].lower()
        if 'baseline' in condition:
            condition = 'baseline'
        elif 'ambtemp' in condition or 'amb' in condition:
            condition = 'ambtemp'
        elif 'drug' in condition:
            condition = 'drugs'
        else:
            continue  # Skip unknown conditions
        
        print(f"\n{'='*60}")
        print(f"Processing: {base_name}")
        print(f"  Condition: {condition}")
        print(f"  Genotype: {metadata['genotype']}")
        
        # Find score file
        score_path = os.path.join(SCORES_DIR, base_name + '_scored_scores_1Hz.csv')
        
        if not os.path.exists(score_path):
            print(f"  ⚠️  No score file found")
            continue
        
        # Load signal
        eeg, fs = load_mat_signal(mat_path)
        if eeg is None:
            continue
        
        # Load scores
        time_s, states = load_scores_1hz(score_path)
        if time_s is None:
            continue
        
        # Trim states
        total_sec = len(eeg) / fs
        valid_mask = time_s < total_sec
        time_s = time_s[valid_mask]
        states = states[valid_mask]
        
        print(f"  Signal: {len(eeg)/fs:.1f}s @ {fs}Hz")
        
        # Find bouts
        bouts_by_state = find_state_bouts(states, min_len_sec=MIN_BOUT_LEN_SEC)
        
        # Extract NREM and REM bouts
        for state_num, state_name in [(1, 'NREM'), (2, 'REM')]:
            bouts_idx = bouts_by_state.get(state_num, [])
            print(f"  Found {len(bouts_idx)} {state_name} bouts")
            
            for bout_num, (start_idx, end_idx) in enumerate(bouts_idx):
                start_s = time_s[start_idx]
                end_s = time_s[end_idx - 1] + 1.0
                duration_s = end_s - start_s
                
                # Determine time window
                time_window = get_time_window_label(start_s)
                if time_window is None:
                    continue
                
                start_sample = int(start_s * fs)
                end_sample = int(end_s * fs)
                
                eeg_bout = eeg[start_sample:end_sample]
                
                if len(eeg_bout) < fs:
                    continue
                
                bout_info = {
                    'mouse_id': metadata['mouse_id'],
                    'genotype': metadata['genotype'].upper(),
                    'condition': condition,
                    'time_window': time_window,
                    'state': state_name,
                    'bout_number': bout_num,
                    'start_s': start_s,
                    'duration_s': duration_s,
                    'n_samples': len(eeg_bout),
                    'srate': fs,
                    'eeg_signal': eeg_bout
                }
                all_bouts.append(bout_info)
    
    if not all_bouts:
        return None
    
    df_bouts = pd.DataFrame(all_bouts)
    
    print(f"\n{'='*60}")
    print("DATA SUMMARY")
    print(f"{'='*60}")
    print(f"Total bouts: {len(df_bouts)}")
    print(f"\nBy condition × genotype × state:")
    print(df_bouts.groupby(['condition', 'genotype', 'state']).size())
    print(f"\nBy time_window × genotype × state:")
    print(df_bouts.groupby(['time_window', 'genotype', 'state']).size())
    
    return df_bouts


# ===========================
# Compute Envelope Features
# ===========================

def compute_envelope_features(df_bouts):
    """Compute band power envelope and its PSD for each bout."""
    print("\nComputing band power envelopes and envelope PSDs...")
    
    envelope_srate = 1.0 / WINDOW_STEP_SEC
    
    results = []
    
    for idx, row in df_bouts.iterrows():
        if idx % 100 == 0:
            print(f"  Processing bout {idx+1}/{len(df_bouts)}...")
        
        eeg = row['eeg_signal']
        srate = row['srate']
        state = row['state']
        
        # Choose band based on state
        if state == 'NREM':
            band_range = SIGMA_RANGE
        else:
            band_range = THETA_RANGE
        
        # Extract power envelope
        time_points, power_env = extract_band_power_envelope(
            eeg, srate, band_range,
            window_size_sec=WINDOW_SIZE_SEC,
            step_sec=WINDOW_STEP_SEC
        )
        
        if len(power_env) < 10:
            continue
        
        # Compute PSD of envelope
        env_freqs, env_psd = compute_envelope_psd(
            power_env,
            envelope_srate,
            freq_range=ENVELOPE_PSD_RANGE
        )
        
        result = {
            **row.to_dict(),
            'power_envelope': power_env,
            'envelope_time': time_points,
            'envelope_psd_freqs': env_freqs,
            'envelope_psd_values': env_psd,
            'mean_power': np.mean(power_env),
            'std_power': np.std(power_env)
        }
        
        results.append(result)
    
    df_envelope = pd.DataFrame(results)
    print(f"✓ Computed envelope features for {len(df_envelope)} bouts")
    
    return df_envelope


# ===========================
# Statistical Analysis Functions
# ===========================

def perform_mixed_anova(df_data, dependent_var='mean_power'):
    """
    Perform mixed-way ANOVA.
    
    Parameters:
      df_data: DataFrame with columns for dependent variable and grouping factors
      dependent_var: name of dependent variable column
      
    Returns:
      anova_table
    """
    # Create formula for 3-way ANOVA
    formula = f'{dependent_var} ~ C(genotype) * C(condition) * C(time_window)'
    
    # Fit model
    model = ols(formula, data=df_data).fit()
    anova_table = anova_lm(model, typ=2)
    
    return anova_table


def add_significance_bar(ax, x1, x2, y, p_value, height_offset=0.05):
    """Add significance bar with stars to plot."""
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    bar_height = y + height_offset * y_range
    
    ax.plot([x1, x2], [bar_height, bar_height], 'k-', linewidth=1.5)
    
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'
    
    ax.text((x1 + x2) / 2, bar_height + 0.02 * y_range, sig_text,
            ha='center', va='bottom', fontsize=10, fontweight='bold')


# ===========================
# Plotting Functions
# ===========================

def plot_condition_comparison_with_stats(df_envelope, state='NREM', band_name='Sigma'):
    """
    Create comprehensive plot with:
    - Top row: 3 panels of PSD comparisons (one per time window)
    - Bottom row: Bar plots for mean power and variability with statistics
    """
    df_state = df_envelope[df_envelope['state'] == state].copy()
    
    if len(df_state) == 0:
        print(f"⚠️  No {state} data found")
        return
    
    print(f"\n{'='*60}")
    print(f"CONDITION COMPARISON: {band_name} during {state}")
    print(f"{'='*60}")
    
    # Create figure with 2 rows: top for PSDs, bottom for bar plots
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.5, 1], wspace=0.3, hspace=0.4)
    
    # Top row: 3 PSD panels
    axes_psd = [fig.add_subplot(gs[0, i:i+2]) for i in [0, 2, 3]]
    
    # Bottom row: 2 bar plot panels
    ax_mean_power = fig.add_subplot(gs[1, 0:2])
    ax_variability = fig.add_subplot(gs[1, 3:5])
    
    # Common frequency grid
    common_freqs = np.linspace(0, 0.15, 100)
    
    # Order of conditions
    condition_order = ['baseline', 'ambtemp', 'drugs']
    
    # ===== TOP ROW: PSD COMPARISONS =====
    for ax_idx, (time_window, ax) in enumerate(zip(['0-3h', '3-6h', 'washout'], axes_psd)):
        df_window = df_state[df_state['time_window'] == time_window]
        
        if len(df_window) == 0:
            ax.text(0.5, 0.5, f'No data for {time_window}', 
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        print(f"\n{time_window}:")
        
        # Plot each genotype-condition combination
        for genotype in ['WT', 'APP']:
            for condition in condition_order:
                df_group = df_window[
                    (df_window['genotype'] == genotype) & 
                    (df_window['condition'] == condition)
                ]
                
                if len(df_group) == 0:
                    continue
                
                # Interpolate PSDs
                all_psds_interp = []
                for _, row in df_group.iterrows():
                    freqs = row['envelope_psd_freqs']
                    psd_vals = row['envelope_psd_values']
                    
                    if len(freqs) < 2 or len(psd_vals) < 2:
                        continue
                    
                    psd_interp = np.interp(common_freqs, freqs, psd_vals, 
                                          left=psd_vals[0], right=psd_vals[-1])
                    all_psds_interp.append(psd_interp)
                
                if len(all_psds_interp) == 0:
                    continue
                
                psd_array = np.array(all_psds_interp)
                mean_psd = np.mean(psd_array, axis=0)
                sem_psd = stats.sem(psd_array, axis=0)
                
                # Get color
                color = COLORS[genotype][condition]
                
                # Create label
                label = f'{genotype}-{condition} (n={len(all_psds_interp)})'
                
                # Plot
                ax.plot(common_freqs, mean_psd, color=color, linewidth=2.5,
                       label=label, alpha=0.9)
                ax.fill_between(common_freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                               color=color, alpha=0.2)
                
                print(f"  {genotype}-{condition}: n={len(all_psds_interp)}")
        
        # Style panel
        ax.set_xlabel('Frequency (Hz)', fontweight='bold')
        if ax_idx == 0:
            ax.set_ylabel(f'{band_name} power (A.U.)', fontweight='bold')
        ax.set_title(f'{time_window}', fontweight='bold', fontsize=13)
        ax.legend(loc='upper right', frameon=True, fancybox=False, 
                 edgecolor='black', fontsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(0, 0.15)
        ax.set_ylim(bottom=0)
        ax.set_xticks([0, 0.05, 0.10, 0.15])
    
    # ===== BOTTOM ROW: BAR PLOTS WITH STATISTICS =====
    
    print("\n" + "="*60)
    print("STATISTICAL ANALYSIS")
    print("="*60)
    
    # Storage for p-values to export
    pvalue_records = []
    
    # Mean Power Analysis
    print("\n--- MEAN POWER ANALYSIS ---")
    anova_mean = perform_mixed_anova(df_state, dependent_var='mean_power')
    print("\nMixed-way ANOVA for Mean Power:")
    print(anova_mean)
    
    # Store ANOVA results
    for effect in anova_mean.index:
        pvalue_records.append({
            'state': state,
            'band': band_name,
            'measure': 'mean_power',
            'test_type': 'ANOVA',
            'comparison': effect,
            'F_statistic': anova_mean.loc[effect, 'F'] if 'F' in anova_mean.columns else np.nan,
            'p_value': anova_mean.loc[effect, 'PR(>F)'] if 'PR(>F)' in anova_mean.columns else np.nan,
            'n_WT': np.nan,
            'n_APP': np.nan
        })
    
    # Bar plot for Mean Power
    x_labels = []
    x_positions = []
    bar_positions_by_group = {}
    bar_idx = 0
    bar_width = 0.35
    
    group_idx = 0
    for time_window in ['0-3h', '3-6h', 'washout']:
        for condition in condition_order:
            wt_pos = None
            app_pos = None
            
            for geno_idx, genotype in enumerate(['WT', 'APP']):
                df_group = df_state[
                    (df_state['genotype'] == genotype) & 
                    (df_state['condition'] == condition) &
                    (df_state['time_window'] == time_window)
                ]
                
                if len(df_group) > 0:
                    mean_val = np.mean(df_group['mean_power'])
                    sem_val = stats.sem(df_group['mean_power'])
                    color = COLORS[genotype][condition]
                    
                    ax_mean_power.bar(bar_idx, mean_val, bar_width, 
                                     yerr=sem_val, color=color, alpha=0.8,
                                     capsize=4, edgecolor='black', linewidth=1)
                    
                    if genotype == 'WT':
                        wt_pos = bar_idx
                    else:
                        app_pos = bar_idx
                    
                    if geno_idx == 0:
                        x_positions.append(bar_idx + bar_width/2)
                        x_labels.append(f"{time_window}\n{condition}")
                    
                    bar_idx += 1
            
            bar_positions_by_group[group_idx] = {
                'wt_pos': wt_pos,
                'app_pos': app_pos,
                'time_window': time_window,
                'condition': condition
            }
            group_idx += 1
            
            bar_idx += 0.3
    
    ax_mean_power.set_xlabel('Time Window - Condition', fontweight='bold')
    ax_mean_power.set_ylabel(f'Mean {band_name} Power (A.U.)', fontweight='bold')
    ax_mean_power.set_title(f'Mean {band_name} Power', fontweight='bold', fontsize=13)
    ax_mean_power.set_xticks(x_positions)
    ax_mean_power.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    ax_mean_power.spines['top'].set_visible(False)
    ax_mean_power.spines['right'].set_visible(False)
    ax_mean_power.set_ylim(bottom=0)
    ax_mean_power.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#808080', edgecolor='black', label='WT'),
        Patch(facecolor='#6495ED', edgecolor='black', label='APP')
    ]
    ax_mean_power.legend(handles=legend_elements, loc='upper left')
    
    # Add significance bars and collect p-values (Mean Power)
    print("\n--- Pairwise t-tests (WT vs APP) for Mean Power ---")
    for idx, info in bar_positions_by_group.items():
        if info['wt_pos'] is not None and info['app_pos'] is not None:
            df_wt = df_state[
                (df_state['genotype'] == 'WT') &
                (df_state['condition'] == info['condition']) &
                (df_state['time_window'] == info['time_window'])
            ]
            df_app = df_state[
                (df_state['genotype'] == 'APP') &
                (df_state['condition'] == info['condition']) &
                (df_state['time_window'] == info['time_window'])
            ]
            
            if len(df_wt) >= 2 and len(df_app) >= 2:
                wt_vals = df_wt['mean_power'].values
                app_vals = df_app['mean_power'].values
                
                stat, p_val = stats.ttest_ind(wt_vals, app_vals)
                
                print(f"{info['time_window']} - {info['condition']}: t={stat:.3f}, p={p_val:.4f}")
                
                # Store p-value
                pvalue_records.append({
                    'state': state,
                    'band': band_name,
                    'measure': 'mean_power',
                    'test_type': 't-test',
                    'comparison': f"WT_vs_APP_{info['condition']}_{info['time_window']}",
                    'F_statistic': stat,
                    'p_value': p_val,
                    'n_WT': len(df_wt),
                    'n_APP': len(df_app)
                })
                
                # Get y-position for significance bar
                wt_mean = np.mean(wt_vals)
                wt_sem = stats.sem(wt_vals)
                app_mean = np.mean(app_vals)
                app_sem = stats.sem(app_vals)
                
                y_max = max(wt_mean + wt_sem, app_mean + app_sem)
                
                add_significance_bar(ax_mean_power, info['wt_pos'], info['app_pos'], 
                                   y_max, p_val, height_offset=0.03)
    
    # Variability Analysis
    print("\n--- VARIABILITY ANALYSIS ---")
    anova_var = perform_mixed_anova(df_state, dependent_var='std_power')
    print("\nMixed-way ANOVA for Variability:")
    print(anova_var)
    
    # Store ANOVA results
    for effect in anova_var.index:
        pvalue_records.append({
            'state': state,
            'band': band_name,
            'measure': 'std_power',
            'test_type': 'ANOVA',
            'comparison': effect,
            'F_statistic': anova_var.loc[effect, 'F'] if 'F' in anova_var.columns else np.nan,
            'p_value': anova_var.loc[effect, 'PR(>F)'] if 'PR(>F)' in anova_var.columns else np.nan,
            'n_WT': np.nan,
            'n_APP': np.nan
        })
    
    # Bar plot for Variability
    bar_positions_by_group_var = {}
    bar_idx = 0
    group_idx = 0
    
    for time_window in ['0-3h', '3-6h', 'washout']:
        for condition in condition_order:
            wt_pos = None
            app_pos = None
            
            for genotype in ['WT', 'APP']:
                df_group = df_state[
                    (df_state['genotype'] == genotype) & 
                    (df_state['condition'] == condition) &
                    (df_state['time_window'] == time_window)
                ]
                
                if len(df_group) > 0:
                    mean_val = np.mean(df_group['std_power'])
                    sem_val = stats.sem(df_group['std_power'])
                    color = COLORS[genotype][condition]
                    
                    ax_variability.bar(bar_idx, mean_val, bar_width, 
                                      yerr=sem_val, color=color, alpha=0.8,
                                      capsize=4, edgecolor='black', linewidth=1)
                    
                    if genotype == 'WT':
                        wt_pos = bar_idx
                    else:
                        app_pos = bar_idx
                    
                    bar_idx += 1
            
            bar_positions_by_group_var[group_idx] = {
                'wt_pos': wt_pos,
                'app_pos': app_pos,
                'time_window': time_window,
                'condition': condition
            }
            group_idx += 1
            
            bar_idx += 0.3
    
    ax_variability.set_xlabel('Time Window - Condition', fontweight='bold')
    ax_variability.set_ylabel(f'{band_name} Power Variability (SD)', fontweight='bold')
    ax_variability.set_title(f'{band_name} Power Variability', fontweight='bold', fontsize=13)
    ax_variability.set_xticks(x_positions)
    ax_variability.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    ax_variability.legend(handles=legend_elements, loc='upper left')
    ax_variability.spines['top'].set_visible(False)
    ax_variability.spines['right'].set_visible(False)
    ax_variability.set_ylim(bottom=0)
    ax_variability.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add significance bars and collect p-values (Variability)
    print("\n--- Pairwise t-tests (WT vs APP) for Variability ---")
    for idx, info in bar_positions_by_group_var.items():
        if info['wt_pos'] is not None and info['app_pos'] is not None:
            df_wt = df_state[
                (df_state['genotype'] == 'WT') &
                (df_state['condition'] == info['condition']) &
                (df_state['time_window'] == info['time_window'])
            ]
            df_app = df_state[
                (df_state['genotype'] == 'APP') &
                (df_state['condition'] == info['condition']) &
                (df_state['time_window'] == info['time_window'])
            ]
            
            if len(df_wt) >= 2 and len(df_app) >= 2:
                wt_vals = df_wt['std_power'].values
                app_vals = df_app['std_power'].values
                
                stat, p_val = stats.ttest_ind(wt_vals, app_vals)
                
                print(f"{info['time_window']} - {info['condition']}: t={stat:.3f}, p={p_val:.4f}")
                
                # Store p-value
                pvalue_records.append({
                    'state': state,
                    'band': band_name,
                    'measure': 'std_power',
                    'test_type': 't-test',
                    'comparison': f"WT_vs_APP_{info['condition']}_{info['time_window']}",
                    'F_statistic': stat,
                    'p_value': p_val,
                    'n_WT': len(df_wt),
                    'n_APP': len(df_app)
                })
                
                # Get y-position for significance bar
                wt_mean = np.mean(wt_vals)
                wt_sem = stats.sem(wt_vals)
                app_mean = np.mean(app_vals)
                app_sem = stats.sem(app_vals)
                
                y_max = max(wt_mean + wt_sem, app_mean + app_sem)
                
                add_significance_bar(ax_variability, info['wt_pos'], info['app_pos'], 
                                   y_max, p_val, height_offset=0.03)
    
    # Overall title
    fig.suptitle(f'{state} - {band_name} Power Analysis Across Conditions and Time', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{band_name}_full_analysis.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()
    
    # Export p-values to CSV
    df_pvalues = pd.DataFrame(pvalue_records)
    pvalue_path = os.path.join(OUTPUT_DIR, f'{state}_{band_name}_statistics.csv')
    df_pvalues.to_csv(pvalue_path, index=False)
    print(f"✓ Saved statistics: {pvalue_path}")
    
    return df_pvalues


def plot_condition_comparison(df_envelope, state='NREM', band_name='Sigma'):
    """Wrapper to call the full analysis with statistics."""
    plot_condition_comparison_with_stats(df_envelope, state, band_name)


def plot_simple_barplots(df_envelope, state='NREM', band_name='Sigma'):
    """
    Create simple bar plots without statistics:
    - Top row: Mean Power (3 panels, one per time window)
    - Bottom row: Variability (3 panels, one per time window)
    """
    df_state = df_envelope[df_envelope['state'] == state].copy()
    
    if len(df_state) == 0:
        print(f"⚠️  No {state} data found")
        return
    
    print(f"\n{'='*60}")
    print(f"SIMPLE BAR PLOTS: {band_name} during {state}")
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
        
        # Plot bars
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
            
            bar_idx += 0.3
        
        ax.set_title(f'{time_window}', fontweight='bold', fontsize=13)
        ax.set_xticks([0.5, 2.15, 3.8])
        ax.set_xticklabels(['Baseline', 'Ambtemp', 'Drugs'], rotation=0)
        
        if col_idx == 0:
            ax.set_ylabel(f'Mean {band_name} Power (A.U.)', fontweight='bold')
        
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
    axes[0, 0].legend(handles=legend_elements, loc='upper left')
    
    # ===== BOTTOM ROW: VARIABILITY =====
    for col_idx, time_window in enumerate(time_windows):
        ax = axes[1, col_idx]
        
        df_window = df_state[df_state['time_window'] == time_window]
        
        if len(df_window) == 0:
            continue
        
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
        
        ax.set_xticks([0.5, 2.15, 3.8])
        ax.set_xticklabels(['Baseline', 'Ambtemp', 'Drugs'], rotation=0)
        
        if col_idx == 0:
            ax.set_ylabel(f'{band_name} Power Variability (SD)', fontweight='bold')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    axes[1, 0].legend(handles=legend_elements, loc='upper left')
    
    fig.suptitle(f'{state} - {band_name} Power Across Conditions and Time Windows',
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{band_name}_simple_barplots.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()


# ===========================
# Main Function
# ===========================

def main():
    print("="*70)
    print("ENVELOPE PSD ANALYSIS - CONDITION COMPARISON")
    print("Comparing baseline, ambtemp, drugs across time windows")
    print("="*70)
    
    # Step 1: Load bouts
    df_bouts = load_all_bouts()
    
    if df_bouts is None:
        print("\n⚠️  No bouts found. Exiting.")
        return
    
    # Step 2: Compute envelope features
    df_envelope = compute_envelope_features(df_bouts)
    
    # Step 3: Create comparison plots
    print("\n" + "="*70)
    print("Creating condition comparison plots...")
    print("="*70)
    
    # NREM - Sigma
    plot_condition_comparison(df_envelope, state='NREM', band_name='Sigma')
    
    # REM - Theta
    plot_condition_comparison(df_envelope, state='REM', band_name='Theta')
    
    # Add simple bar plots without statistics
    print("\n" + "="*70)
    print("Creating simple bar plots...")
    print("="*70)
    
    plot_simple_barplots(df_envelope, state='NREM', band_name='Sigma')
    plot_simple_barplots(df_envelope, state='REM', band_name='Theta')
    
    print("\n" + "="*70)
    print("✓ ANALYSIS COMPLETE")
    print(f"Plots saved to: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()