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
ENVELOPE_PSD_RANGE = (0, 0.15)  # Frequency range for envelope PSD (0-0.15 Hz like paper)

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


def compute_envelope_psd(power_envelope, envelope_srate, freq_range=(0, 0.15)):
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
    
    # Use longer windows for better low-frequency resolution
    # This ensures we get frequencies closer to 0 Hz
    nperseg = min(len(power_envelope_detrended), int(60 * envelope_srate))  # 60 second windows
    noverlap = nperseg // 2
    
    freqs, psd = signal.welch(
        power_envelope_detrended,
        fs=envelope_srate,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling='density',
        window='hann'
    )
    
    # Filter to specified range (include 0 Hz now)
    mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    return freqs[mask], psd[mask]


# ===========================
# Load Baseline Bouts
# ===========================

def load_baseline_bouts():
    """Load all REM and NREM bouts from baseline recordings ONLY."""
    signal_files = glob.glob(os.path.join(SIGNAL_DIR, "*.mat"))
    
    if not signal_files:
        print(f"⚠️  No .mat files found in {SIGNAL_DIR}")
        return None
    
    print(f"\nFound {len(signal_files)} total .mat files")
    
    # Filter to baseline only - CHECK MULTIPLE PATTERNS
    baseline_files = []
    for f in signal_files:
        fname_lower = os.path.basename(f).lower()
        # Check if baseline is in the filename
        if 'baseline' in fname_lower:
            # Exclude group/combined files
            if not any(exclude in fname_lower for exclude in ['group', 'all', 'concat', 'combined']):
                baseline_files.append(f)
                print(f"  ✓ Baseline file: {os.path.basename(f)}")
    
    print(f"\n{'='*60}")
    print(f"FILTERED TO {len(baseline_files)} BASELINE FILES")
    print(f"{'='*60}")
    
    if len(baseline_files) == 0:
        print("⚠️  No baseline files found! Check your filename patterns.")
        return None
    
    all_bouts = []
    
    for mat_path in baseline_files:
        base_name = os.path.basename(mat_path).replace('.mat', '')
        metadata = parse_filename(base_name)
        
        # Double-check condition
        if 'baseline' not in metadata['condition'].lower():
            print(f"  ⚠️  Skipping {base_name} - not baseline in metadata")
            continue
        
        print(f"\n{'='*60}")
        print(f"Processing: {base_name}")
        print(f"  Condition: {metadata['condition']}")
        print(f"  Genotype: {metadata['genotype']}")
        print(f"  Mouse ID: {metadata['mouse_id']}")
        
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
                    'condition': 'baseline',  # Force to baseline
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
    
    # Verify all are baseline
    print(f"\n{'='*60}")
    print("BASELINE VERIFICATION")
    print(f"{'='*60}")
    print(f"Conditions in data: {df_bouts['condition'].unique()}")
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
    Create publication-style plot matching the paper format.
    """
    df_state = df_envelope[df_envelope['state'] == state].copy()
    
    if len(df_state) == 0:
        print(f"⚠️  No {state} data found")
        return
    
    print(f"\n{'='*60}")
    print(f"ENVELOPE PSD ANALYSIS: {band_name} Oscillations during Natural {state} Sleep")
    print(f"{'='*60}")
    
    # Create figure
    fig = plt.figure(figsize=(15, 4))
    gs = fig.add_gridspec(1, 3, width_ratios=[2, 1, 1], wspace=0.4)
    
    # ===== LEFT: Envelope PSD =====
    ax_spectrum = fig.add_subplot(gs[0])
    
    # Create common frequency grid for interpolation
    common_freqs = np.linspace(0, 0.15, 100)
    
    for genotype, color in [('WT', COLOR_WT), ('APP', COLOR_APP)]:
        df_geno = df_state[df_state['genotype'] == genotype]
        
        if len(df_geno) == 0:
            continue
        
        # Interpolate all PSDs to common frequency grid
        all_psds_interp = []
        
        for _, row in df_geno.iterrows():
            freqs = row['envelope_psd_freqs']
            psd_vals = row['envelope_psd_values']
            
            if len(freqs) < 2 or len(psd_vals) < 2:
                continue
            
            # Interpolate to common grid
            psd_interp = np.interp(common_freqs, freqs, psd_vals, left=psd_vals[0], right=psd_vals[-1])
            all_psds_interp.append(psd_interp)
        
        if len(all_psds_interp) == 0:
            print(f"  ⚠️  No valid {genotype} bouts after interpolation")
            continue
        
        psd_array = np.array(all_psds_interp)
        mean_psd = np.mean(psd_array, axis=0)
        sem_psd = stats.sem(psd_array, axis=0)
        
        # Plot
        ax_spectrum.plot(common_freqs, mean_psd, color=color, linewidth=2.5,
                        label=f'{genotype} (n={len(all_psds_interp)})', alpha=0.9)
        ax_spectrum.fill_between(common_freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                                color=color, alpha=0.2)
        
        print(f"\n{genotype}:")
        print(f"  Total bouts in state: {len(df_geno)}")
        print(f"  Valid bouts after interpolation: {len(all_psds_interp)}")
        print(f"  PSD max: {np.max(mean_psd):.3e}")
    
    # Match paper style
    ax_spectrum.set_xlabel('Frequency (Hz)', fontweight='bold')
    ax_spectrum.set_ylabel(f'{band_name} power (A.U.)', fontweight='bold')
    ax_spectrum.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black')
    ax_spectrum.spines['top'].set_visible(False)
    ax_spectrum.spines['right'].set_visible(False)
    ax_spectrum.set_xlim(0, 0.15)
    ax_spectrum.set_ylim(bottom=0)
    ax_spectrum.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
    
    # ===== MIDDLE: Mean Band Power =====
    ax_power = fig.add_subplot(gs[1])
    
    wt_power = df_state[df_state['genotype'] == 'WT']['mean_power'].values
    app_power = df_state[df_state['genotype'] == 'APP']['mean_power'].values
    
    print(f"\nMean power sample sizes:")
    print(f"  WT: {len(wt_power)} bouts")
    print(f"  APP: {len(app_power)} bouts")
    
    x_pos = [0, 1]
    means = [np.mean(wt_power), np.mean(app_power)]
    sems = [stats.sem(wt_power), stats.sem(app_power)]
    colors = [COLOR_WT, COLOR_APP]
    
    bars = ax_power.bar(x_pos, means, 0.6, yerr=sems, color=colors,
                       alpha=0.8, capsize=5, edgecolor='black', linewidth=1.5)
    
    if len(wt_power) >= 2 and len(app_power) >= 2:
        stat, p_val = stats.ttest_ind(wt_power, app_power)
        print(f"\nMean {band_name} Power:")
        print(f"  WT: {np.mean(wt_power):.3e} ± {stats.sem(wt_power):.3e}")
        print(f"  APP: {np.mean(app_power):.3e} ± {stats.sem(app_power):.3e}")
        print(f"  t-test: t={stat:.3f}, p={p_val:.4f}")
        
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        add_significance_bar(ax_power, 0, 1, y_max, p_val)
    
    ax_power.set_xticks(x_pos)
    ax_power.set_xticklabels(['WT', 'APP/PS1'], fontweight='bold')
    ax_power.set_ylabel(f'Mean {band_name} Power (A.U.)', fontweight='bold')
    ax_power.spines['top'].set_visible(False)
    ax_power.spines['right'].set_visible(False)
    ax_power.set_ylim(bottom=0)
    
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
        print(f"  WT: {np.mean(wt_std):.3e} ± {stats.sem(wt_std):.3e}")
        print(f"  APP: {np.mean(app_std):.3e} ± {stats.sem(app_std):.3e}")
        print(f"  t-test: t={stat:.3f}, p={p_val:.4f}")
        
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        add_significance_bar(ax_std, 0, 1, y_max, p_val)
    
    ax_std.set_xticks(x_pos)
    ax_std.set_xticklabels(['WT', 'APP/PS1'], fontweight='bold')
    ax_std.set_ylabel(f'{band_name} Power Variability (SD)', fontweight='bold')
    ax_std.spines['top'].set_visible(False)
    ax_std.spines['right'].set_visible(False)
    ax_std.set_ylim(bottom=0)
    
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'baseline_{state}_{band_name}_envelope_PSD.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {out_path}")
    plt.close()


# ===========================
# Main Function
# ===========================

def main():
    print("="*70)
    print("ENVELOPE PSD ANALYSIS - BASELINE ONLY")
    print("Analyzing slow oscillations in Sigma/Theta power")
    print("="*70)
    
    # Step 1: Load bouts
    df_bouts = load_baseline_bouts()
    
    if df_bouts is None:
        print("\n⚠️  No bouts found. Exiting.")
        return
    
    # Step 2: Compute envelope features
    df_envelope = compute_envelope_features(df_bouts)
    
    # Step 3: Create plots
    print("\n" + "="*70)
    print("Creating plots...")
    print("="*70)
    
    # NREM - Sigma
    plot_envelope_psd_analysis(df_envelope, state='NREM', band_name='Sigma', band_range=SIGMA_RANGE)
    
    # REM - Theta
    plot_envelope_psd_analysis(df_envelope, state='REM', band_name='Theta', band_range=THETA_RANGE)
    
    print("\n" + "="*70)
    print("✓ ANALYSIS COMPLETE")
    print(f"Plots saved to: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()