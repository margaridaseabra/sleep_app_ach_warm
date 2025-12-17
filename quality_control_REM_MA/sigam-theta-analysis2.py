#!/usr/bin/env python3
"""
Power Spectral Density (PSD) Analysis of Sigma/Theta Envelope Oscillations

Analyzes slow oscillations in sigma/theta power during natural sleep:
1. Extract sigma power over time during NREM (envelope)
2. Compute PSD of sigma power time series to find cyclic patterns
3. Similarly for theta power during REM

This reveals how sigma/theta power oscillates at slow frequencies.
"""

import os
import glob
import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy import stats, signal

# ===========================
# Configuration
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"
OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/sigma_theta_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Bout selection
MIN_BOUT_LEN_SEC = 5.0  # minimum bout length in seconds

# Frequency bands for power extraction
SIGMA_RANGE = (11, 16)  # Hz - NREM sigma
THETA_RANGE = (4, 10)   # Hz - REM theta

# Envelope PSD parameters
WINDOW_SIZE_SEC = 4.0  # Window size for computing band power
WINDOW_STEP_SEC = 0.5  # Step size for sliding window
ENVELOPE_PSD_RANGE = (0, 0.2)  # Frequency range for envelope PSD (slow oscillations)

# Colors
COLOR_WT = '#808080'  # Grey
COLOR_APP = '#6495ED'  # Cornflower blue

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 10

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


# ===========================
# Band Power Envelope Extraction
# ===========================

def extract_band_power_envelope(eeg_signal, srate, band_range, window_size_sec=4.0, step_sec=0.5):
    """
    Extract band power over time using sliding window.
    
    Returns:
      time_points : array of time points (seconds)
      power_envelope : array of band power at each time point
    """
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


def compute_envelope_psd(power_envelope, envelope_srate, freq_range=(0, 0.2)):
    """
    Compute PSD of the power envelope to find slow oscillations.
    
    Parameters:
      power_envelope : time series of band power
      envelope_srate : sampling rate of the envelope (1/step_sec)
      freq_range : frequency range for output
      
    Returns:
      freqs, psd
    """
    # Detrend to remove DC component
    power_envelope_detrended = signal.detrend(power_envelope, type='linear')
    
    # Compute PSD with appropriate window length
    # Use longer windows for better frequency resolution at low frequencies
    nperseg = min(len(power_envelope_detrended), int(30 * envelope_srate))  # 30 second windows
    noverlap = nperseg // 2
    
    freqs, psd = signal.welch(
        power_envelope_detrended,
        fs=envelope_srate,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling='density',
        window='hann'
    )
    
    # Filter to positive frequencies only and specified range
    mask = (freqs > 0) & (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    return freqs[mask], psd[mask]


# ===========================
# Load Baseline Bouts
# ===========================

def load_baseline_bouts():
    """Load all REM and NREM bouts from baseline recordings."""
    signal_files = glob.glob(os.path.join(SIGNAL_DIR, "*.mat"))
    
    if not signal_files:
        print(f"⚠️  No .mat files found in {SIGNAL_DIR}")
        return None
    
    # Filter to baseline only
    baseline_files = [f for f in signal_files if 'baseline' in os.path.basename(f).lower()]
    
    # Exclude group files
    baseline_files = [f for f in baseline_files if not any(
        exclude in os.path.basename(f).lower() 
        for exclude in ['group', 'all', 'concat', 'combined']
    )]
    
    print(f"Found {len(baseline_files)} baseline files...")
    
    all_bouts = []
    
    for mat_path in baseline_files:
        base_name = os.path.basename(mat_path).replace('.mat', '')
        metadata = parse_filename(base_name)
        
        if metadata['condition'].lower() != 'baseline':
            continue
        
        print(f"\nProcessing: {base_name}")
        
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
                
                start_sample = int(start_s * fs)
                end_sample = int(end_s * fs)
                
                eeg_bout = eeg[start_sample:end_sample]
                
                if len(eeg_bout) < fs:
                    continue
                
                bout_info = {
                    'mouse_id': metadata['mouse_id'],
                    'genotype': metadata['genotype'].upper(),
                    'condition': metadata['condition'].lower(),
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
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total bouts: {len(df_bouts)}")
    print(f"\nBy genotype × state:")
    print(df_bouts.groupby(['genotype', 'state']).size())
    
    return df_bouts


# ===========================
# Compute Envelope Features
# ===========================

def compute_envelope_features(df_bouts):
    """
    Compute band power envelope and its PSD for each bout.
    """
    print("\nComputing band power envelopes and envelope PSDs...")
    
    envelope_srate = 1.0 / WINDOW_STEP_SEC  # Sampling rate of envelope
    
    results = []
    
    for idx, row in df_bouts.iterrows():
        if idx % 50 == 0:
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
        
        if len(power_env) < 10:  # Skip if too short
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
# Plotting Functions
# ===========================

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
            ha='center', va='bottom', fontsize=12, fontweight='bold')


def plot_envelope_psd_analysis(df_envelope, state='NREM', band_name='Sigma', band_range=SIGMA_RANGE):
    """
    Create publication-style plot for envelope PSD analysis.
    
    - Left: PSD of band power envelope (slow oscillations)
    - Middle: Mean band power
    - Right: Variability of band power (std)
    """
    df_state = df_envelope[df_envelope['state'] == state].copy()
    
    if len(df_state) == 0:
        print(f"⚠️  No {state} data found")
        return
    
    print(f"\n{'='*60}")
    print(f"ENVELOPE PSD ANALYSIS: {band_name} Oscillations during Natural {state} Sleep")
    print(f"{'='*60}")
    
    sleep_type = "Natural NREM Sleep" if state == 'NREM' else "Natural REM Sleep"
    
    # Create figure
    fig = plt.figure(figsize=(15, 4))
    gs = fig.add_gridspec(1, 3, width_ratios=[2, 1, 1], wspace=0.4)
    
    # ===== LEFT: Envelope PSD =====
    ax_spectrum = fig.add_subplot(gs[0])
    
    for genotype, color in [('WT', COLOR_WT), ('APP', COLOR_APP)]:
        df_geno = df_state[df_state['genotype'] == genotype]
        
        if len(df_geno) == 0:
            continue
        
        # Average envelope PSD across all bouts
        all_psds = []
        common_freqs = None
        
        for _, row in df_geno.iterrows():
            freqs = row['envelope_psd_freqs']
            psd_vals = row['envelope_psd_values']
            
            # Skip if invalid
            if len(freqs) == 0 or len(psd_vals) == 0:
                continue
            
            if common_freqs is None:
                common_freqs = freqs
            
            # Only add if frequencies match
            if len(freqs) == len(common_freqs) and np.allclose(freqs, common_freqs):
                all_psds.append(psd_vals)
        
        if len(all_psds) == 0:
            continue
        
        psd_array = np.array(all_psds)
        mean_psd = np.mean(psd_array, axis=0)
        sem_psd = stats.sem(psd_array, axis=0)
        
        # Plot with cleaner styling
        ax_spectrum.plot(common_freqs, mean_psd, color=color, linewidth=2.5,
                        label=f'{genotype} (n={len(all_psds)} bouts)', alpha=0.9)
        ax_spectrum.fill_between(common_freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                                color=color, alpha=0.2)
        
        print(f"\n{genotype} envelope PSD:")
        print(f"  Frequency range: {common_freqs[0]:.4f} - {common_freqs[-1]:.4f} Hz")
        print(f"  PSD range: {np.min(mean_psd):.2e} - {np.max(mean_psd):.2e}")
    
    ax_spectrum.set_xlabel('Frequency (Hz)', fontweight='bold')
    ax_spectrum.set_ylabel(f'{band_name} Power (A.U.)', fontweight='bold')
    ax_spectrum.set_title(f'{sleep_type} - {band_name} Power Oscillations', fontweight='bold', fontsize=13)
    ax_spectrum.legend(loc='upper right', frameon=True)
    ax_spectrum.spines['top'].set_visible(False)
    ax_spectrum.spines['right'].set_visible(False)
    ax_spectrum.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Set x-axis to start from 0
    ax_spectrum.set_xlim(left=0, right=ENVELOPE_PSD_RANGE[1])
    
    # Use scientific notation if values are very small
    ax_spectrum.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    
    # ===== MIDDLE: Mean Band Power =====
    ax_power = fig.add_subplot(gs[1])
    
    wt_power = df_state[df_state['genotype'] == 'WT']['mean_power'].values
    app_power = df_state[df_state['genotype'] == 'APP']['mean_power'].values
    
    x_pos = [0, 1]
    means = [np.mean(wt_power), np.mean(app_power)]
    sems = [stats.sem(wt_power), stats.sem(app_power)]
    colors = [COLOR_WT, COLOR_APP]
    
    bars = ax_power.bar(x_pos, means, 0.6, yerr=sems, color=colors,
                       alpha=0.8, capsize=5, edgecolor='black', linewidth=1.5)
    
    if len(wt_power) >= 2 and len(app_power) >= 2:
        stat, p_val = stats.ttest_ind(wt_power, app_power)
        print(f"\nMean {band_name} Power:")
        print(f"  WT: {np.mean(wt_power):.4e} ± {stats.sem(wt_power):.4e}")
        print(f"  APP: {np.mean(app_power):.4e} ± {stats.sem(app_power):.4e}")
        print(f"  t-test: t={stat:.3f}, p={p_val:.4f}")
        
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        add_significance_bar(ax_power, 0, 1, y_max, p_val)
    
    ax_power.set_xticks(x_pos)
    ax_power.set_xticklabels(['WT', 'APP'], fontweight='bold')
    ax_power.set_ylabel(f'Mean {band_name} Power (A.U.)', fontweight='bold')
    ax_power.set_title(f'Mean {band_name} Power', fontweight='bold', fontsize=13)
    ax_power.spines['top'].set_visible(False)
    ax_power.spines['right'].set_visible(False)
    ax_power.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    ax_power.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    
    # ===== RIGHT: Power Variability =====
    ax_std = fig.add_subplot(gs[2])
    
    wt_std = df_state[df_state['genotype'] == 'WT']['std_power'].values
    app_std = df_state[df_state['genotype'] == 'APP']['std_power'].values
    
    means = [np.mean(wt_std), np.mean(app_std)]
    sems = [stats.sem(wt_std), stats.sem(app_std)]
    
    bars = ax_std.bar(x_pos, means, 0.6, yerr=sems, color=colors,
                     alpha=0.8, capsize=5, edgecolor='black', linewidth=1.5)
    
    if len(wt_std) >= 2 and len(app_std) >= 2:
        stat, p_val = stats.ttest_ind(wt_std, app_std)
        print(f"\n{band_name} Power Variability (SD):")
        print(f"  WT: {np.mean(wt_std):.4e} ± {stats.sem(wt_std):.4e}")
        print(f"  APP: {np.mean(app_std):.4e} ± {stats.sem(app_std):.4e}")
        print(f"  t-test: t={stat:.3f}, p={p_val:.4f}")
        
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        add_significance_bar(ax_std, 0, 1, y_max, p_val)
    
    ax_std.set_xticks(x_pos)
    ax_std.set_xticklabels(['WT', 'APP'], fontweight='bold')
    ax_std.set_ylabel(f'{band_name} Power Variability (SD)', fontweight='bold')
    ax_std.set_title(f'{band_name} Power Variability', fontweight='bold', fontsize=13)
    ax_std.spines['top'].set_visible(False)
    ax_std.spines['right'].set_visible(False)
    ax_std.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    ax_std.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{band_name}_envelope_PSD_analysis.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()


# ===========================
# Main Function
# ===========================

def main():
    print("="*70)
    print("ENVELOPE PSD ANALYSIS")
    print("Analyzing slow oscillations in Sigma/Theta power")
    print("="*70)
    
    # Step 1: Load bouts
    print("\nStep 1: Loading natural sleep bouts from baseline recordings...")
    df_bouts = load_baseline_bouts()
    
    if df_bouts is None:
        print("\n⚠️  No bouts found. Exiting.")
        return
    
    # Step 2: Compute envelope features
    print("\nStep 2: Computing band power envelopes and envelope PSDs...")
    df_envelope = compute_envelope_features(df_bouts)
    
    # Step 3: Create plots
    print("\n" + "="*70)
    print("Step 3: Creating envelope PSD analysis plots...")
    print("="*70)
    
    # NREM - Sigma envelope oscillations
    plot_envelope_psd_analysis(df_envelope, state='NREM', band_name='Sigma', band_range=SIGMA_RANGE)
    
    # REM - Theta envelope oscillations
    plot_envelope_psd_analysis(df_envelope, state='REM', band_name='Theta', band_range=THETA_RANGE)
    
    print("\n" + "="*70)
    print("✓ ENVELOPE PSD ANALYSIS COMPLETE")
    print(f"Plots saved to: {OUTPUT_DIR}")
    print("="*70)
    
    print("\nAnalysis Summary:")
    print(f"  • Analyzed natural sleep bouts ≥{MIN_BOUT_LEN_SEC}s")
    print(f"  • NREM Sigma band: {SIGMA_RANGE[0]}-{SIGMA_RANGE[1]} Hz")
    print(f"  • REM Theta band: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz")
    print(f"  • Window size: {WINDOW_SIZE_SEC}s, step: {WINDOW_STEP_SEC}s")
    print(f"  • Envelope PSD range: {ENVELOPE_PSD_RANGE[0]}-{ENVELOPE_PSD_RANGE[1]} Hz")
    print(f"  • Shows slow oscillations in band power (cyclic patterns)")


if __name__ == "__main__":
    main()