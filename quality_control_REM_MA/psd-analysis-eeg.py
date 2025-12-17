#!/usr/bin/env python3
"""
Generate presentation-ready PSD and band power plots for WT vs APP comparison.

Analyzes:
1. Power Spectral Density (PSD) across frequencies
2. Band power (delta, theta, alpha, beta, gamma)
3. Comparisons: baseline (WT vs APP) and conditions × time blocks
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, signal
from itertools import combinations
import mne
import scipy.io as sio

# ===========================
# Configuration
# ===========================

# DATA PATHS - UPDATE THESE

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"
OUTPUT_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/psd_plots"
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

# Frequency bands (Hz)
FREQ_BANDS = {
    'delta': (0.5, 4),
    'theta': (4, 8),
    'alpha': (8, 12),
    'sigma': (9, 15),      # ADDED SIGMA BAND
    'beta': (12, 30),
    'gamma': (30, 80)
}

# Sampling rate (adjust based on your data)
SRATE = 256  # Hz - ADJUST THIS TO YOUR ACTUAL SAMPLING RATE

# Notch filter settings
NOTCH_FREQ = 50  # Hz - frequency to notch out (powerline interference)
NOTCH_Q = 30  # Quality factor (higher = narrower notch)

# Recordings to EXCLUDE from notch filtering (will NOT be notched)
EXCLUDE_FROM_NOTCH = [
    'baseline_mouse1',
    'baseline_mouse2', 
    'ambtemp_mouse1',
    'baseline_mouse8',
    'baseline_mouse4'
]

# ===========================
# Signal Preprocessing
# ===========================

def apply_notch_filter(signal_data, srate, notch_freq=50, quality_factor=30):
    """
    Apply notch filter to remove powerline interference.
    
    Parameters
    ----------
    signal_data : array
        Time series data
    srate : float
        Sampling rate in Hz
    notch_freq : float
        Frequency to notch out (e.g., 50 Hz or 60 Hz)
    quality_factor : float
        Q factor - higher values = narrower notch
        
    Returns
    -------
    filtered_signal : array
        Notch-filtered signal
    """
    from scipy.signal import iirnotch, filtfilt
    
    # Ensure signal is 1D
    signal_data = np.array(signal_data).flatten()
    
    # Check if signal is long enough
    if len(signal_data) < 3 * srate:
        print(f"    ⚠️  Signal too short for notch filter ({len(signal_data)} samples), skipping...")
        return signal_data
    
    try:
        # Design notch filter
        # Normalize frequency to Nyquist frequency
        nyquist_freq = srate / 2.0
        w0 = notch_freq / nyquist_freq
        
        # Check if notch frequency is valid
        if w0 <= 0 or w0 >= 1:
            print(f"    ⚠️  Invalid notch frequency {notch_freq} Hz for sampling rate {srate} Hz")
            return signal_data
        
        b, a = iirnotch(w0, quality_factor)
        
        # Apply filter (zero-phase filtering)
        filtered_signal = filtfilt(b, a, signal_data)
        
        print(f"    ✓ Notch filter applied successfully at {notch_freq} Hz")
        
        return filtered_signal
        
    except Exception as e:
        print(f"    ⚠️  Error applying notch filter: {e}")
        return signal_data


def should_apply_notch(metadata):
    """
    Determine if notch filter should be applied based on metadata.
    
    Parameters
    ----------
    metadata : dict
        Contains 'condition' and 'mouse_id' keys
        
    Returns
    -------
    apply_notch : bool
        True if notch should be applied, False otherwise
    """
    condition = metadata.get('condition', '').lower()
    mouse_id = metadata.get('mouse_id', '').lower()
    
    # Create identifier string
    identifier = f"{condition}_{mouse_id}"
    
    # Check if this recording should be excluded from notching
    for exclude_pattern in EXCLUDE_FROM_NOTCH:
        if exclude_pattern.lower() in identifier:
            return False
    
    return True


# ===========================
# Data Loading and Segmentation
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


def load_mat_file(mat_path):
    """
    Load signal data from .mat file.
    
    Returns
    -------
    signals : dict
        Dictionary with keys like 'EEG', 'ACh', 'srate', etc.
    """
    try:
        mat_data = sio.loadmat(mat_path)
        
        # Check what variables are in the .mat file
        mat_keys = [k for k in mat_data.keys() if not k.startswith('__')]
        print(f"    Available variables in .mat: {mat_keys}")
        
        # Try different possible variable names
        eeg = None
        ach = None
        srate = SRATE  # Default to global setting
        timestamps = None
        
        # EEG signal
        for key in ['EEG', 'eeg', 'EEG_signal', 'eeg_signal', 'EEG1', 'signal']:
            if key in mat_data:
                eeg = mat_data[key].flatten()
                print(f"    Found EEG as '{key}', length: {len(eeg)}")
                break
        
        # ACh signal
        for key in ['ACh', 'ach', 'ACh_signal', 'ach_signal', 'Ach', 'ACH']:
            if key in mat_data:
                ach = mat_data[key].flatten()
                print(f"    Found ACh as '{key}', length: {len(ach)}")
                break
        
        # Sampling rate - ADD eeg_frequency to the list
        for key in ['eeg_frequency', 'srate', 'Fs', 'fs', 'sampling_rate', 'SamplingRate', 'EEG_frequency']:
            if key in mat_data:
                srate_val = mat_data[key]
                if hasattr(srate_val, 'flatten'):
                    srate = float(srate_val.flatten()[0])
                else:
                    srate = float(srate_val)
                print(f"    Found sampling rate from '{key}': {srate} Hz")
                break
        
        if srate == SRATE:
            print(f"    Using default sampling rate: {srate} Hz")
        
        # Timestamps
        for key in ['timestamps', 'time', 't', 'Time']:
            if key in mat_data:
                timestamps = mat_data[key].flatten()
                print(f"    Found timestamps, length: {len(timestamps)}")
                break
        
        if eeg is None:
            print(f"    ⚠️  Could not find EEG signal in {mat_keys}")
            return None
        
        signals = {
            'eeg': eeg,
            'ach': ach,
            'srate': srate,
            'timestamps': timestamps,
        }
        
        return signals
        
    except Exception as e:
        print(f"⚠️  Error loading {mat_path}: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_score_file(csv_path):
    """
    Load sleep scoring from .csv file.
    
    Returns
    -------
    scores : DataFrame
        With columns: state (and optionally epoch, time)
    """
    try:
        scores = pd.read_csv(csv_path)
        
        # Check what columns are available
        print(f"    CSV columns: {list(scores.columns)}")
        
        # Try to find state column
        state_col = None
        for col in ['state', 'State', 'STATE', 'score', 'Score', 'label', 'Label']:
            if col in scores.columns:
                state_col = col
                break
        
        if state_col is None:
            # Try first column
            state_col = scores.columns[0]
            print(f"    ⚠️  No 'state' column found, using first column: '{state_col}'")
        
        # Rename to standard 'state' column
        if state_col != 'state':
            scores = scores.rename(columns={state_col: 'state'})
        
        print(f"    Found {len(scores)} epochs")
        print(f"    State values: {scores['state'].unique()}")
        
        return scores
        
    except Exception as e:
        print(f"⚠️  Error loading {csv_path}: {e}")
        import traceback
        traceback.print_exc()
        return None


def segment_signal_by_state(signals, scores, target_state='NREM', epoch_duration=4):
    """
    Extract bouts of specific sleep state from continuous signal.
    """
    # Map state names to numbers if needed
    state_map = {'Wake': 0, 'NREM': 1, 'REM': 2, 
                 'wake': 0, 'nrem': 1, 'rem': 2,
                 'W': 0, 'N': 1, 'R': 2,
                 0: 0, 1: 1, 2: 2}
    
    if isinstance(target_state, str):
        target_state_num = state_map.get(target_state, target_state)
    else:
        target_state_num = target_state
    
    # Find consecutive epochs of target state
    state_array = scores['state'].values
    
    # Convert states to numbers if they're strings
    try:
        state_array = np.array([state_map.get(s, s) for s in state_array])
    except:
        pass
    
    srate = signals['srate']
    samples_per_epoch = int(epoch_duration * srate)
    
    bouts = []
    in_bout = False
    bout_start_epoch = 0
    
    for i, state in enumerate(state_array):
        if state == target_state_num:
            if not in_bout:
                in_bout = True
                bout_start_epoch = i
        else:
            if in_bout:
                # End current bout
                bout_end_epoch = i
                
                # Extract signal segment
                start_sample = bout_start_epoch * samples_per_epoch
                end_sample = bout_end_epoch * samples_per_epoch
                
                eeg_bout = signals['eeg'][start_sample:end_sample] if signals['eeg'] is not None else None
                ach_bout = signals['ach'][start_sample:end_sample] if signals['ach'] is not None else None
                
                bout = {
                    'eeg_signal': eeg_bout,
                    'ach_signal': ach_bout,
                    'start_s': bout_start_epoch * epoch_duration,
                    'duration_s': (bout_end_epoch - bout_start_epoch) * epoch_duration,
                    'n_epochs': bout_end_epoch - bout_start_epoch
                }
                
                bouts.append(bout)
                in_bout = False
    
    # Handle case where recording ends in target state
    if in_bout:
        bout_end_epoch = len(state_array)
        start_sample = bout_start_epoch * samples_per_epoch
        end_sample = min(bout_end_epoch * samples_per_epoch, len(signals['eeg']))
        
        eeg_bout = signals['eeg'][start_sample:end_sample] if signals['eeg'] is not None else None
        ach_bout = signals['ach'][start_sample:end_sample] if signals['ach'] is not None else None
        
        bout = {
            'eeg_signal': eeg_bout,
            'ach_signal': ach_bout,
            'start_s': bout_start_epoch * epoch_duration,
            'duration_s': (bout_end_epoch - bout_start_epoch) * epoch_duration,
            'n_epochs': bout_end_epoch - bout_start_epoch
        }
        bouts.append(bout)
    
    return bouts


def load_signal_data(state='NREM'):
    """
    Load and segment signal data from .mat files using .csv scores.
    NOW: Applies notch filter to FULL recording BEFORE segmentation.
    THEN: Segments into bouts based on user-specified state (NREM, REM, or Wake).
    
    Parameters
    ----------
    state : str or int
        Sleep state to extract: 'NREM', 'REM', 'Wake' (or 1, 2, 0)
    
    Returns DataFrame with columns: mouse_id, genotype, condition, time_block,
                                    bout_index, eeg_signal, ach_signal, start_s
    """
    # Find all .mat signal files
    signal_files = glob.glob(os.path.join(SIGNAL_DIR, "*.mat"))
    
    if not signal_files:
        print(f"⚠️  No .mat files found in {SIGNAL_DIR}")
        return None
    
    print(f"Found {len(signal_files)} signal files...")
    print(f"Target state: {state}")
    
    # DIAGNOSTIC: Print exclusion list
    print(f"\n=== NOTCH FILTER EXCLUSION LIST ===")
    print(f"Files matching these patterns will NOT be notched:")
    for pattern in EXCLUDE_FROM_NOTCH:
        print(f"  - {pattern}")
    print(f"{'='*40}\n")
    
    # Find all score files
    score_files = glob.glob(os.path.join(SCORES_DIR, "*.csv"))
    print(f"Found {len(score_files)} score files...")
    
    # Create a mapping of (mouse_id, condition, genotype) -> score_path
    score_map = {}
    for score_path in score_files:
        score_name = os.path.basename(score_path).lower()
        
        # Extract mouse_id
        mouse_id = None
        for part in score_name.replace('.csv', '').split('_'):
            if 'mouse' in part:
                mouse_id = part
                break
        
        # Extract genotype
        genotype = None
        if 'wt' in score_name:
            genotype = 'WT'
        elif 'app' in score_name:
            genotype = 'APP'
        
        # Extract condition
        condition = None
        if 'baseline' in score_name:
            condition = 'baseline'
        elif 'ambtemp' in score_name or 'amb' in score_name:
            condition = 'ambtemp'
        elif 'drug' in score_name:
            condition = 'drugs'
        
        if mouse_id and genotype and condition:
            key = (mouse_id, condition, genotype)
            score_map[key] = score_path
            print(f"  Mapped: {key} -> {os.path.basename(score_path)}")
    
    print(f"\nCreated mapping for {len(score_map)} score files")
    
    all_data = []
    notch_summary = []  # Track notch decisions
    
    for mat_path in signal_files:
        # Extract metadata from filename
        base_name = os.path.basename(mat_path).replace('.mat', '')
        metadata = parse_filename(base_name)
        
        print(f"\n{'='*70}")
        print(f"Processing: {base_name}")
        print(f"{'='*70}")
        print(f"  Metadata: {metadata}")
        
        # Normalize metadata
        mouse_id = metadata['mouse_id'].lower()
        genotype = metadata['genotype'].upper()
        condition = metadata['condition'].lower()
        
        if 'baseline' in condition:
            condition = 'baseline'
        elif 'ambtemp' in condition or 'amb' in condition:
            condition = 'ambtemp'
        elif 'drug' in condition:
            condition = 'drugs'
        
        # CREATE IDENTIFIER STRING FOR DEBUGGING
        identifier = f"{condition}_{mouse_id}"
        print(f"  Identifier string: '{identifier}'")
        
        # Try to find matching score file
        score_path = None
        
        # Try exact match first
        key = (mouse_id, condition, genotype)
        if key in score_map:
            score_path = score_map[key]
            print(f"  ✓ Found exact match: {os.path.basename(score_path)}")
        else:
            # Try flexible matching
            print(f"  Looking for: mouse_id={mouse_id}, condition={condition}, genotype={genotype}")
            for (map_mouse, map_cond, map_geno), path in score_map.items():
                if (map_mouse == mouse_id and 
                    map_cond == condition and 
                    map_geno == genotype):
                    score_path = path
                    print(f"  ✓ Found match: {os.path.basename(score_path)}")
                    break
        
        if score_path is None:
            # Try even more flexible matching (just mouse_id and genotype)
            print(f"  ⚠️  No exact match, trying flexible matching...")
            for (map_mouse, map_cond, map_geno), path in score_map.items():
                if map_mouse == mouse_id and map_geno == genotype:
                    score_path = path
                    print(f"  ⚠️  Using: {os.path.basename(score_path)} (condition mismatch: {map_cond} vs {condition})")
                    break
        
        if score_path is None:
            print(f"  ⚠️  No matching score file found, SKIPPING")
            continue
        
        # Check if notch filter should be applied
        apply_notch = should_apply_notch(metadata)
        
        # DETAILED NOTCH DECISION LOGGING
        print(f"\n  → NOTCH FILTER DECISION:")
        print(f"     Identifier: '{identifier}'")
        print(f"     Checking against exclusion list...")
        
        matched_exclusion = False
        for exclude_pattern in EXCLUDE_FROM_NOTCH:
            if exclude_pattern.lower() in identifier:
                print(f"     ✗ MATCHED exclusion pattern: '{exclude_pattern}'")
                matched_exclusion = True
                break
        
        if not matched_exclusion:
            print(f"     ✓ NOT in exclusion list")
        
        print(f"     Final decision: apply_notch = {apply_notch}")
        
        if apply_notch:
            print(f"  ✓✓✓ WILL APPLY {NOTCH_FREQ} Hz notch filter TO FULL RECORDING ✓✓✓")
        else:
            print(f"  ⚠️⚠️⚠️ SKIPPING notch filter (excluded recording) ⚠️⚠️⚠️")
        
        # Track decision
        notch_summary.append({
            'filename': base_name,
            'identifier': identifier,
            'genotype': genotype,
            'condition': condition,
            'mouse_id': mouse_id,
            'apply_notch': apply_notch,
            'reason': 'excluded' if not apply_notch else 'will_notch'
        })
        
        # Load signals and scores
        signals = load_mat_file(mat_path)
        scores = load_score_file(score_path)
        
        if signals is None or scores is None:
            continue
        
        # IMPORTANT: Store the actual sampling rate from the file
        actual_srate = signals['srate']
        print(f"\n  → Sampling rate: {actual_srate} Hz")
        print(f"  → Full EEG length: {len(signals['eeg'])} samples ({len(signals['eeg'])/actual_srate:.1f} seconds)")
        
        # ========================================
        # APPLY NOTCH TO FULL RECORDING (BEFORE SEGMENTATION)
        # ========================================
        if apply_notch and signals['eeg'] is not None:
            print(f"\n  ⚡⚡⚡ APPLYING NOTCH TO FULL EEG RECORDING ⚡⚡⚡")
            print(f"      Frequency: {NOTCH_FREQ} Hz")
            print(f"      Sampling rate: {actual_srate} Hz")
            print(f"      Signal length: {len(signals['eeg'])} samples")
            
            signals['eeg'] = apply_notch_filter(
                signals['eeg'], 
                actual_srate,
                notch_freq=NOTCH_FREQ,
                quality_factor=NOTCH_Q
            )
            print(f"  ✓ FULL EEG recording notch-filtered")
        else:
            if not apply_notch:
                print(f"  ⊘ NOTCH SKIPPED for full EEG (file in exclusion list)")
        
        # Apply notch filter to full ACh recording if present
        if apply_notch and signals['ach'] is not None:
            print(f"\n  ⚡⚡⚡ APPLYING NOTCH TO FULL ACh RECORDING ⚡⚡⚡")
            print(f"      Frequency: {NOTCH_FREQ} Hz")
            print(f"      Signal length: {len(signals['ach'])} samples")
            
            signals['ach'] = apply_notch_filter(
                signals['ach'], 
                actual_srate,
                notch_freq=NOTCH_FREQ,
                quality_factor=NOTCH_Q
            )
            print(f"  ✓ FULL ACh recording notch-filtered")
        
        # ========================================
        # NOW SEGMENT INTO BOUTS (BASED ON USER-SPECIFIED STATE)
        # ========================================
        print(f"\n  → Segmenting into {state} bouts...")
        bouts = segment_signal_by_state(signals, scores, target_state=state)
        print(f"  ✓ Found {len(bouts)} {state} bouts")
        
        # Add to dataset (notch already applied to full recording)
        for bout_idx, bout in enumerate(bouts):
            entry = {
                'mouse_id': metadata['mouse_id'],
                'genotype': genotype,
                'condition': condition,
                'bout_index': bout_idx,
                'eeg_signal': bout['eeg_signal'],
                'ach_signal': bout['ach_signal'],
                'start_s': bout['start_s'],
                'duration_s': bout['duration_s'],
                'srate': actual_srate,
                'notch_applied': apply_notch,
                'state': state  # ADD STATE TO TRACK WHICH STATE WAS ANALYZED
            }
            all_data.append(entry)
    
    # PRINT NOTCH SUMMARY TABLE
    print("\n" + "="*80)
    print("NOTCH FILTER SUMMARY")
    print("="*80)
    df_notch = pd.DataFrame(notch_summary)
    print("\nFiles that WERE notched (full recording):")
    notched_files = df_notch[df_notch['apply_notch'] == True]
    if len(notched_files) > 0:
        for _, row in notched_files.iterrows():
            print(f"  ✓ {row['filename']} ({row['genotype']}, {row['condition']})")
    else:
        print("  (None)")
    
    print("\nFiles that were NOT notched (excluded):")
    excluded_files = df_notch[df_notch['apply_notch'] == False]
    if len(excluded_files) > 0:
        for _, row in excluded_files.iterrows():
            print(f"  ✗ {row['filename']} ({row['genotype']}, {row['condition']}) - EXCLUDED")
    else:
        print("  (None)")
    print("="*80)
    
    if not all_data:
        print("\n⚠️  No valid data loaded!")
        return None
    
    df = pd.DataFrame(all_data)
    
    # Assign time blocks
    df = assign_time_blocks(df)
    
    # Report sampling rates found
    unique_srates = df['srate'].unique()
    print(f"\n{'='*60}")
    print(f"✓ SUCCESS! Total {state} bouts loaded: {len(df)}")
    print(f"  • Sampling rates found: {unique_srates} Hz")
    
    # Report notch filter statistics
    notched = df['notch_applied'].sum()
    total = len(df)
    print(f"  • Notch filter applied: {notched} bouts ({100*notched/total:.1f}%)")
    print(f"  • No notch (excluded): {total-notched} bouts ({100*(total-notched)/total:.1f}%)")
    print(f"\nGenotypes: {df['genotype'].unique()}")
    print(f"Conditions: {df['condition'].unique()}")
    print(f"Animals: {sorted(df['mouse_id'].unique())}")
    print(f"Time blocks: {df['time_block'].unique()}")
    
    # BASELINE-ONLY SUMMARY (FOR ANY STATE)
    df_baseline = df[df['condition'] == 'baseline']
    if len(df_baseline) > 0:
        print(f"\n{'='*60}")
        print(f"BASELINE {state.upper()} RECORDINGS ONLY:")
        print(f"  • Total baseline {state} bouts: {len(df_baseline)}")
        print(f"  • WT animals: {sorted(df_baseline[df_baseline['genotype']=='WT']['mouse_id'].unique())}")
        print(f"  • APP animals: {sorted(df_baseline[df_baseline['genotype']=='APP']['mouse_id'].unique())}")
        print(f"  • WT n={df_baseline[df_baseline['genotype']=='WT']['mouse_id'].nunique()}")
        print(f"  • APP n={df_baseline[df_baseline['genotype']=='APP']['mouse_id'].nunique()}")
        
        # Check notch status for baseline
        baseline_notched = df_baseline['notch_applied'].sum()
        baseline_total = len(df_baseline)
        print(f"  • Baseline {state} bouts notched: {baseline_notched}/{baseline_total} ({100*baseline_notched/baseline_total:.1f}%)")
    
    print(f"{'='*60}")
    
    return df


# ===========================
# PSD Computation
# ===========================

def compute_psd_welch(signal_data, srate=256, nperseg=None, noverlap=None, freq_range=(0.5, 100), convert_to_db=False):  # CHANGED: default False
    """
    Compute PSD using Welch's method.
    
    Parameters
    ----------
    signal_data : array
        Time series data
    srate : int
        Sampling rate in Hz
    nperseg : int, optional
        Length of each segment (default: 4 seconds)
    noverlap : int, optional
        Overlap between segments (default: 50%)
    freq_range : tuple
        (min_freq, max_freq) to return
    convert_to_db : bool
        If True, convert PSD to dB scale (10 * log10(PSD)) - DEFAULT FALSE
        
    Returns
    -------
    freqs : array
        Frequency values
    psd : array
        Power spectral density (in LINEAR scale by default)
    """
    if nperseg is None:
        nperseg = int(4 * srate)  # 4 second windows
    if noverlap is None:
        noverlap = nperseg // 2
    
    # Use scipy.signal.welch
    from scipy import signal as scipy_signal
    freqs, psd = scipy_signal.welch(signal_data, fs=srate, nperseg=nperseg, 
                                     noverlap=noverlap, scaling='density')
    
    # Convert to dB if requested (but default is NO)
    if convert_to_db:
        psd = 10 * np.log10(psd + 1e-20)
    
    # Filter to frequency range
    freq_mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    freqs = freqs[freq_mask]
    psd = psd[freq_mask]
    
    return freqs, psd


def compute_band_power(freqs, psd, band_range):
    """
    Compute total power in a frequency band.
    
    Parameters
    ----------
    freqs : array
        Frequency values
    psd : array
        Power spectral density in LINEAR scale
    band_range : tuple
        (min_freq, max_freq)
        
    Returns
    -------
    power : float
        Total power in band (integrated PSD in LINEAR scale)
    """
    mask = (freqs >= band_range[0]) & (freqs < band_range[1])
    if not np.any(mask):
        return 0.0
    
    # Integrate PSD in linear scale
    power = np.trapz(psd[mask], freqs[mask])
    
    return power  # Returns LINEAR scale power


def convert_power_to_db(power):
    """
    Convert power values to dB scale.
    
    Parameters
    ----------
    power : float or array
        Power value(s) in linear scale
        
    Returns
    -------
    power_db : float or array
        Power in dB scale (10 * log10(power))
    """
    return 10 * np.log10(power + 1e-20)  # Add small epsilon to avoid log(0)


def convert_to_db_for_plot(values):
    """
    Convert power values to dB scale for plotting only.
    Handles both single values and arrays.
    
    Parameters
    ----------
    values : float, array, or Series
        Power value(s) in linear scale
        
    Returns
    -------
    values_db : float or array
        Power in dB scale (10 * log10(values))
    """
    values_array = np.array(values)
    return 10 * np.log10(values_array + 1e-20)


# ===========================
# PSD Processing Functions
# ===========================

def process_all_bouts_psd(df, signal_col='eeg_signal'):
    """
    Compute PSD and band powers for all bouts.
    Stores everything in LINEAR scale.
    
    Returns
    -------
    df_psd : DataFrame
        Contains freqs, psd (LINEAR), and band powers (LINEAR) for each bout
    """
    results = []
    
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"  Processing bout {idx+1}/{len(df)}...")
        
        signal_data = row[signal_col]
        
        if signal_data is None or len(signal_data) < 256:
            continue
        
        # Use actual sampling rate from the bout
        srate = row.get('srate', SRATE)
        
        # Compute PSD in LINEAR scale
        freqs, psd = compute_psd_welch(signal_data, srate=srate, convert_to_db=False)
        
        # Compute band powers (in linear scale)
        band_powers = {}
        for band_name, band_range in FREQ_BANDS.items():
            band_powers[f'{band_name}_power'] = compute_band_power(freqs, psd, band_range)
        
        # Total power
        band_powers['total_power'] = compute_band_power(freqs, psd, (0.5, 100))
        
        # Relative powers (normalized by total)
        for band_name in FREQ_BANDS.keys():
            if band_powers['total_power'] > 0:
                band_powers[f'{band_name}_rel_power'] = (
                    band_powers[f'{band_name}_power'] / band_powers['total_power']
                )
            else:
                band_powers[f'{band_name}_rel_power'] = 0
        
        # Store results (ALL IN LINEAR SCALE)
        result = {
            'mouse_id': row['mouse_id'],
            'genotype': row['genotype'],
            'condition': row['condition'],
            'time_block': row['time_block'],
            'bout_index': row['bout_index'],
            'notch_applied': row.get('notch_applied', False),
            'srate': srate,
            'freqs': freqs,
            'psd': psd,  # LINEAR scale
            **band_powers  # LINEAR scale
        }
        
        results.append(result)
    
    return pd.DataFrame(results)


def compute_per_animal_psd(df_psd):
    """
    Average PSDs per animal (across bouts and time blocks).
    Interpolates all PSDs to a common frequency grid.
    
    Parameters
    ----------
    df_psd : DataFrame
        Output from process_all_bouts_psd with columns:
        mouse_id, genotype, condition, time_block, freqs, psd, band_powers
    
    Returns
    -------
    df_animal : DataFrame
        One row per (animal, condition, time_block) with mean PSD and band powers
        All PSDs use common frequency grid for consistent averaging
    """
    results = []
    
    # Define common frequency grid (0.5 to 100 Hz)
    common_freqs = np.linspace(0.5, 100, 200)  # 200 points for smooth interpolation
    
    print("\n  Averaging PSDs per animal...")
    
    # Group by animal, condition, time block
    for (mouse_id, genotype, condition, time_block), group in df_psd.groupby(
        ['mouse_id', 'genotype', 'condition', 'time_block']
    ):
        # Interpolate all PSDs to common frequency grid
        interpolated_psds = []
        
        for idx, row in group.iterrows():
            freqs = row['freqs']
            psd = row['psd']
            
            if len(freqs) < 2 or len(psd) < 2:
                continue
            
            # Interpolate to common grid
            psd_interp = np.interp(common_freqs, freqs, psd, 
                                  left=psd[0], right=psd[-1])
            interpolated_psds.append(psd_interp)
        
        if len(interpolated_psds) == 0:
            continue
        
        # Convert to array and average
        psds_array = np.array(interpolated_psds)
        mean_psd = np.mean(psds_array, axis=0)
        
        # Average band powers (already in linear scale)
        band_powers = {}
        for band_name in FREQ_BANDS.keys():
            band_powers[f'{band_name}_power'] = group[f'{band_name}_power'].mean()
            band_powers[f'{band_name}_rel_power'] = group[f'{band_name}_rel_power'].mean()
        
        band_powers['total_power'] = group['total_power'].mean()
        
        result = {
            'mouse_id': mouse_id,
            'genotype': genotype,
            'condition': condition,
            'time_block': time_block,
            'freqs': common_freqs,  # Use common frequency grid
            'mean_psd': mean_psd,    # Averaged PSD (linear scale)
            'n_bouts': len(group),
            **band_powers
        }
        
        results.append(result)
    
    df_animal = pd.DataFrame(results)
    
    print(f"\n✓ Computed per-animal PSDs:")
    print(f"  Total animals: {len(df_animal)}")
    print(f"  Genotypes: {df_animal['genotype'].unique()}")
    print(f"  Conditions: {df_animal['condition'].unique()}")
    print(f"  Time blocks: {df_animal['time_block'].unique()}")
    
    return df_animal


# ===========================
# Statistical Functions
# ===========================

def compare_two_groups_unpaired(group1, group2):
    """
    Perform unpaired t-test between two groups.
    
    Parameters
    ----------
    group1, group2 : Series or array
        Data for each group
        
    Returns
    -------
    results : dict
        Contains t_statistic, p_value, means, stds, cohens_d, n values
    """
    from scipy import stats
    
    # Remove NaNs
    g1 = group1.dropna() if hasattr(group1, 'dropna') else group1[~np.isnan(group1)]
    g2 = group2.dropna() if hasattr(group2, 'dropna') else group2[~np.isnan(group2)]
    
    if len(g1) < 2 or len(g2) < 2:
        return {
            'statistic': np.nan,
            'p_value': 1.0,
            'mean1': np.nan,
            'mean2': np.nan,
            'std1': np.nan,
            'std2': np.nan,
            'n1': len(g1),
            'n2': len(g2),
            'cohens_d': np.nan
        }
    
    # Perform t-test
    t_stat, p_val = stats.ttest_ind(g1, g2)
    
    # Calculate Cohen's d
    pooled_std = np.sqrt(((len(g1) - 1) * np.var(g1, ddof=1) + 
                          (len(g2) - 1) * np.var(g2, ddof=1)) / 
                         (len(g1) + len(g2) - 2))
    
    cohens_d = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0
    
    return {
        'statistic': t_stat,
        'p_value': p_val,
        'mean1': np.mean(g1),
        'mean2': np.mean(g2),
        'std1': np.std(g1, ddof=1),
        'std2': np.std(g2, ddof=1),
        'n1': len(g1),
        'n2': len(g2),
        'cohens_d': cohens_d
    }


def perform_rm_anova(df, value_col, within_factors, between_factor, subject_col='mouse_id'):
    """
    Perform repeated measures ANOVA using statsmodels.
    
    Parameters
    ----------
    df : DataFrame
        Data with columns for subject, factors, and values
    value_col : str
        Column name for dependent variable
    within_factors : list
        List of within-subject factor column names
    between_factor : str
        Between-subject factor column name
    subject_col : str
        Column name for subject identifier
        
    Returns
    -------
    results : dict
        ANOVA results including F-statistics, p-values, effect sizes
    """
    try:
        from statsmodels.formula.api import ols
        from statsmodels.stats.anova import anova_lm
        
        # Build formula
        factor_str = ' * '.join([f'C({f})' for f in within_factors])
        formula = f"{value_col} ~ C({between_factor}) * ({factor_str})"
        
        # Fit model
        model = ols(formula, data=df).fit()
        
        # ANOVA table
        anova_table = anova_lm(model, typ=2)
        
        # Extract p-values
        pvalues = {}
        for idx in anova_table.index:
            if 'Residual' not in idx:
                pvalues[idx] = anova_table.loc[idx, 'PR(>F)']
        
        return {
            'model': formula,
            'table': anova_table,
            'pvalues': pvalues,
            'summary': str(anova_table)
        }
        
    except Exception as e:
        print(f"    ⚠️  RM-ANOVA failed: {e}")
        return {'error': str(e)}


# ===========================
# PLOT 1: Baseline PSD Comparison
# ===========================

def plot_baseline_psd_comparison(df_animal, state='NREM', signal_type='EEG'):
    """
    Plot PSD curves and band powers for baseline WT vs APP.
    Data stored in linear scale, converted to dB for display only.
    ONLY INCLUDES BASELINE CONDITION RECORDINGS.
    """
    # FILTER FOR BASELINE ONLY
    df_baseline = df_animal[df_animal['condition'] == 'baseline'].copy()
    
    if len(df_baseline) == 0:
        print("⚠️  No baseline data found!")
        return
    
    # VERIFY: Check that we only have baseline
    unique_conditions = df_baseline['condition'].unique()
    if len(unique_conditions) > 1 or 'baseline' not in unique_conditions:
        print(f"⚠️  WARNING: Non-baseline data found: {unique_conditions}")
        print("   Re-filtering to baseline only...")
        df_baseline = df_baseline[df_baseline['condition'] == 'baseline'].copy()
    
    print(f"\n=== BASELINE PSD COMPARISON (WT vs APP) ===")
    print(f"  Conditions in data: {df_baseline['condition'].unique()}")
    
    # Count UNIQUE animals per genotype (BASELINE ONLY)
    n_wt = df_baseline[df_baseline['genotype'] == 'WT']['mouse_id'].nunique()
    n_app = df_baseline[df_baseline['genotype'] == 'APP']['mouse_id'].nunique()
    
    print(f"  Unique BASELINE animals: WT n={n_wt}, APP n={n_app}")
    
    # List the actual animals
    wt_animals = sorted(df_baseline[df_baseline['genotype'] == 'WT']['mouse_id'].unique())
    app_animals = sorted(df_baseline[df_baseline['genotype'] == 'APP']['mouse_id'].unique())
    print(f"  WT animals: {wt_animals}")
    print(f"  APP animals: {app_animals}")
    
    # Create figure
    fig = plt.figure(figsize=(18, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=0.3)
    
    ax_psd = fig.add_subplot(gs[0])
    ax_bands = fig.add_subplot(gs[1])
    
    # --- PSD Curves (convert to dB for plotting) ---
    for genotype in ['WT', 'APP']:
        df_geno = df_baseline[df_baseline['genotype'] == genotype]
        
        if len(df_geno) == 0:
            print(f"  ⚠️  No {genotype} baseline data!")
            continue
        
        # VERIFY: All should be baseline condition
        assert all(df_geno['condition'] == 'baseline'), f"Non-baseline data found in {genotype} group!"
        
        # Get all PSDs (in linear scale)
        psds = np.array([p for p in df_geno['mean_psd'].values])
        freqs = df_geno['freqs'].iloc[0]
        
        # Mean and SEM in linear scale
        mean_psd = np.mean(psds, axis=0)
        sem_psd = stats.sem(psds, axis=0)
        
        # ❌ THIS WAS WRONG - sem of linear values, then convert to dB doesn't work
        # ✅ CORRECT: Convert each PSD to dB first, THEN compute mean and SEM
        psds_db = convert_to_db_for_plot(psds)
        mean_psd_db = np.mean(psds_db, axis=0)
        sem_psd_db = stats.sem(psds_db, axis=0)
        
        color = COLORS_CONDITION[genotype]['baseline']
        n_animals = df_geno['mouse_id'].nunique()
        
        ax_psd.plot(freqs, mean_psd_db, label=f'{genotype} (n={n_animals})',
                   color=color, linewidth=2.5, alpha=0.9)
        ax_psd.fill_between(freqs, mean_psd_db - sem_psd_db, mean_psd_db + sem_psd_db,
                           color=color, alpha=0.2)
    
    ax_psd.set_xlabel('Frequency (Hz)', fontweight='bold')
    ax_psd.set_ylabel('Power (dB)', fontweight='bold')
    ax_psd.set_title(f'{signal_type} Power Spectral Density', fontsize=13, fontweight='bold')
    ax_psd.legend(loc='upper right', frameon=True)
    ax_psd.set_xlim([0, 100])
    ax_psd.grid(True, alpha=0.3)
    ax_psd.spines['top'].set_visible(False)
    ax_psd.spines['right'].set_visible(False)
    
    # --- Band Powers (convert to dB for display) ---
    band_names = list(FREQ_BANDS.keys())
    x_positions = np.arange(len(band_names))
    width = 0.35
    
    wt_data = df_baseline[df_baseline['genotype'] == 'WT']
    app_data = df_baseline[df_baseline['genotype'] == 'APP']
    
    # Get linear values and convert to dB for plotting
    wt_means_db = []
    wt_sems_db = []
    for band in band_names:
        vals_linear = wt_data[f'{band}_power'].values
        if len(vals_linear) > 0:
            vals_db = convert_to_db_for_plot(vals_linear)
            wt_means_db.append(np.mean(vals_db))
            wt_sems_db.append(stats.sem(vals_db))
        else:
            wt_means_db.append(0)
            wt_sems_db.append(0)
    
    app_means_db = []
    app_sems_db = []
    for band in band_names:
        vals_linear = app_data[f'{band}_power'].values
        if len(vals_linear) > 0:
            vals_db = convert_to_db_for_plot(vals_linear)
            app_means_db.append(np.mean(vals_db))
            app_sems_db.append(stats.sem(vals_db))
        else:
            app_means_db.append(0)
            app_sems_db.append(0)
    
    ax_bands.bar(x_positions - width/2, wt_means_db, width, yerr=wt_sems_db,
                color=COLORS_CONDITION['WT']['baseline'], alpha=0.85,
                capsize=4, edgecolor='black', linewidth=1.5, label='WT')
    ax_bands.bar(x_positions + width/2, app_means_db, width, yerr=app_sems_db,
                color=COLORS_CONDITION['APP']['baseline'], alpha=0.85,
                capsize=4, edgecolor='black', linewidth=1.5, label='APP')
    
    # Add significance markers (statistical test on LINEAR values, display in dB)
    print("\n  Baseline band power comparisons:")
    for i, band in enumerate(band_names):
        wt_vals_linear = wt_data[f'{band}_power'].values
        app_vals_linear = app_data[f'{band}_power'].values
        
        if len(wt_vals_linear) >= 2 and len(app_vals_linear) >= 2:
            # Perform t-test on LINEAR values
            result = compare_two_groups_unpaired(pd.Series(wt_vals_linear), pd.Series(app_vals_linear))
            
            # But report in dB for readability
            wt_mean_db = np.mean(convert_to_db_for_plot(wt_vals_linear))
            wt_std_db = np.std(convert_to_db_for_plot(wt_vals_linear))
            app_mean_db = np.mean(convert_to_db_for_plot(app_vals_linear))
            app_std_db = np.std(convert_to_db_for_plot(app_vals_linear))
            
            print(f"    {band.capitalize()}: WT={wt_mean_db:.2f}±{wt_std_db:.2f} dB, "
                  f"APP={app_mean_db:.2f}±{app_std_db:.2f} dB, "
                  f"p={result['p_value']:.4f}, d={result['cohens_d']:.3f}")
            
            if result['p_value'] < 0.05:
                max_val = max(wt_means_db[i] + wt_sems_db[i], app_means_db[i] + app_sems_db[i])
                y_pos = max_val * 1.1
                sig = '***' if result['p_value'] < 0.001 else '**' if result['p_value'] < 0.01 else '*'
                
                ax_bands.plot([i - width/2, i + width/2], [y_pos, y_pos], 'k-', linewidth=1.5)
                ax_bands.text(i, y_pos, sig, ha='center', va='bottom',
                            fontsize=12, fontweight='bold')
    
    ax_bands.set_xticks(x_positions)
    ax_bands.set_xticklabels([b.capitalize() for b in band_names], rotation=15)
    ax_bands.set_ylabel('Band Power (dB)', fontweight='bold')
    ax_bands.set_title('Band Power', fontsize=13, fontweight='bold')
    ax_bands.legend(loc='upper right', frameon=True)
    ax_bands.spines['top'].set_visible(False)
    ax_bands.spines['right'].set_visible(False)
    ax_bands.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add sample size info to title
    fig.suptitle(f'{state} State - BASELINE: WT vs APP Comparison\n'
                f'{signal_type} Power Spectrum Analysis (WT n={n_wt}, APP n={n_app})',
                fontsize=16, fontweight='bold')
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_BASELINE_{signal_type}_PSD_comparison.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


# ===========================
# PLOT 2: All Conditions by Time Block (Band Powers)
# ===========================

def plot_band_powers_by_condition(df_animal, state='NREM', signal_type='EEG'):
    """
    Plot band powers across all conditions organized by time blocks.
    NOW SHOWS POWER IN dB SCALE AND INCLUDES RM-ANOVA.
    """
    conditions = ['baseline', 'ambtemp', 'drugs']
    time_blocks = ['0-3h', '3-6h', 'Washout']
    band_names = list(FREQ_BANDS.keys())
    
    print(f"\n=== BAND POWER ANALYSIS BY CONDITION ===")
    
    # Create one figure per frequency band
    for band in band_names:
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        
        power_col = f'{band}_power'
        all_pairwise_results = []
        
        # Prepare data for RM-ANOVA (across all time blocks)
        anova_data = []
        for idx, row in df_animal.iterrows():
            if power_col in row and not pd.isna(row[power_col]):
                power_db = convert_power_to_db(row[power_col])
                anova_data.append({
                    'mouse_id': row['mouse_id'],
                    'genotype': row['genotype'],
                    'condition': row['condition'],
                    'time_block': row['time_block'],
                    'power_db': power_db
                })
        
        df_anova = pd.DataFrame(anova_data)
        
        # Perform RM-ANOVA
        if len(df_anova) > 0:
            print(f"\n  Performing RM-ANOVA for {band} band...")
            anova_results = perform_rm_anova(
                df_anova, 
                value_col='power_db',
                within_factors=['condition', 'time_block'],
                between_factor='genotype',
                subject_col='mouse_id'
            )
            
            if 'error' not in anova_results:
                print(f"    ✓ RM-ANOVA completed")
                print(f"    Model: {anova_results['model']}")
                if 'pvalues' in anova_results:
                    print(f"    Main effects and interactions:")
                    for effect, pval in anova_results['pvalues'].items():
                        if pval < 0.05:
                            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                            print(f"      {effect}: p={pval:.4f} {sig}")
            else:
                print(f"    ⚠️  {anova_results['error']}")
        
        for ax_idx, time_block in enumerate(time_blocks):
            ax = axes[ax_idx]
            
            df_block = df_animal[df_animal['time_block'] == time_block].copy()
            
            if len(df_block) == 0:
                ax.text(0.5, 0.5, f'No data for {time_block}',
                       ha='center', va='center', transform=ax.transAxes)
                continue
            
            # X-axis setup
            bar_width = 0.35
            group_spacing = 0.2
            
            x_positions = []
            x_labels = []
            x_tick_positions = []
            
            current_x = 0
            for condition in conditions:
                if condition not in df_block['condition'].values:
                    continue
                
                x_positions.append((condition, 'WT', current_x))
                x_positions.append((condition, 'APP', current_x + bar_width))
                
                x_tick_positions.append(current_x + bar_width/2)
                x_labels.append(condition.capitalize())
                
                current_x += 2 * bar_width + group_spacing
            
            # Plot bars (CONVERT TO dB)
            for condition, genotype, x_pos in x_positions:
                df_subset = df_block[
                    (df_block['condition'] == condition) &
                    (df_block['genotype'] == genotype)
                ]
                
                if len(df_subset) == 0:
                    continue
                
                vals = df_subset[power_col].values
                
                if len(vals) > 0:
                    # CONVERT TO dB
                    vals_db = convert_power_to_db(vals)
                    
                    color = COLORS_CONDITION[genotype][condition]
                    ax.bar(x_pos, vals_db.mean(), bar_width,
                          yerr=stats.sem(vals_db), color=color,
                          alpha=0.85, capsize=4, edgecolor='black', linewidth=1.5)
                    
                    x_jitter = np.random.normal(x_pos, 0.02, size=len(vals_db))
                    ax.scatter(x_jitter, vals_db, alpha=0.7, s=50,
                             color='black', edgecolors='white', linewidths=0.5, zorder=3)
            
            # Significance testing (using dB values)
            current_x = 0
            for condition in conditions:
                if condition not in df_block['condition'].values:
                    continue
                
                df_cond = df_block[df_block['condition'] == condition]
                wt_vals = df_cond[df_cond['genotype'] == 'WT'][power_col].values
                app_vals = df_cond[df_cond['genotype'] == 'APP'][power_col].values
                
                if len(wt_vals) >= 2 and len(app_vals) >= 2:
                    # Convert to dB for comparison
                    wt_vals_db = convert_power_to_db(wt_vals)
                    app_vals_db = convert_power_to_db(app_vals)
                    
                    result = compare_two_groups_unpaired(pd.Series(wt_vals_db), pd.Series(app_vals_db))
                    
                    all_pairwise_results.append({
                        'band': band,
                        'condition': condition,
                        'time_block': time_block,
                        'p_value': result['p_value'],
                        't_stat': result['statistic'],
                        'cohens_d': result['cohens_d'],
                        'n_wt': result['n1'],
                        'n_app': result['n2'],
                        'mean_wt_db': result['mean1'],
                        'mean_app_db': result['mean2']
                    })
                    
                    if result['p_value'] < 0.05:
                        max_val = max(wt_vals_db.max(), app_vals_db.max())
                        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
                        y_pos = max_val + 0.05 * y_range
                        sig = '***' if result['p_value'] < 0.001 else '**' if result['p_value'] < 0.01 else '*'
                        
                        x1 = current_x
                        x2 = current_x + bar_width
                        ax.plot([x1, x1, x2, x2],
                               [y_pos, y_pos + 0.015 * y_range,
                                y_pos + 0.015 * y_range, y_pos],
                               'k-', linewidth=1.5)
                        ax.text((x1 + x2) / 2, y_pos + 0.02 * y_range,
                               sig, ha='center', va='bottom', fontsize=12, fontweight='bold')
                
                current_x += 2 * bar_width + group_spacing
            
            # Format (CHANGED Y-LABEL TO dB)
            ax.set_xticks(x_tick_positions)
            ax.set_xticklabels(x_labels, fontsize=11, fontweight='bold')
            ax.set_ylabel(f'{band.capitalize()} Power (dB)', fontweight='bold', fontsize=12)
            ax.set_title(f'{time_block}', fontsize=14, fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
            
            # Legend
            if ax_idx == 0:
                from matplotlib.patches import Patch
                legend_elements = []
                for genotype in ['WT', 'APP']:
                    for cond in conditions:
                        color = COLORS_CONDITION[genotype][cond]
                        legend_elements.append(
                            Patch(facecolor=color, label=f'{genotype} {cond.capitalize()}',
                                 alpha=0.85, edgecolor='black')
                        )
                ax.legend(handles=legend_elements, loc='upper left', frameon=True,
                         fontsize=9, ncol=1)
        
        band_range = FREQ_BANDS[band]
        fig.suptitle(f'{state} State - {band.capitalize()} Band Power ({band_range[0]}-{band_range[1]} Hz)\n'
                    f'{signal_type}: WT vs APP across conditions (Power in dB)',
                    fontsize=16, fontweight='bold')
        fig.tight_layout(rect=[0, 0.02, 1, 1])
        
        out_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_{band}_power_by_condition.png')
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {out_path}")
        plt.close()
        
        # Save statistics
        if all_pairwise_results:
            csv_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_{band}_power_statistics.csv')
            pd.DataFrame(all_pairwise_results).to_csv(csv_path, index=False)
            print(f"✓ Saved statistics: {csv_path}")
        
        # Save RM-ANOVA results
        if len(df_anova) > 0 and 'error' not in anova_results:
            anova_csv_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_{band}_power_RMANOVA.txt')
            with open(anova_csv_path, 'w') as f:
                f.write(f"Repeated Measures ANOVA Results\n")
                f.write(f"Band: {band} ({band_range[0]}-{band_range[1]} Hz)\n")
                f.write(f"State: {state}\n")
                f.write(f"Signal: {signal_type}\n")
                f.write(f"\n{anova_results['summary']}\n")
            print(f"✓ Saved RM-ANOVA: {anova_csv_path}")


# ===========================
# PLOT 3: PSD Curves by Condition
# ===========================

def plot_psd_curves_by_condition(df_animal, state='NREM', signal_type='EEG'):
    """
    Plot PSD curves for each condition×genotype combination.
    """
    conditions = sorted(df_animal['condition'].unique())
    
    fig, axes = plt.subplots(1, len(conditions), figsize=(6*len(conditions), 5))
    if len(conditions) == 1:
        axes = [axes]
    
    for ax, condition in zip(axes, conditions):
        df_cond = df_animal[df_animal['condition'] == condition]
        
        for genotype in ['WT', 'APP']:
            df_geno = df_cond[df_cond['genotype'] == genotype]
            
            if len(df_geno) == 0:
                continue
            
            psds = np.array([p for p in df_geno['mean_psd'].values])
            freqs = df_geno['freqs'].iloc[0]
            
            # ✅ CONVERT TO dB FIRST, THEN compute mean/sem
            psds_db = convert_to_db_for_plot(psds)
            mean_psd = np.mean(psds_db, axis=0)
            sem_psd = stats.sem(psds_db, axis=0)
            
            color = COLORS_CONDITION[genotype][condition]
            
            ax.plot(freqs, mean_psd, label=f'{genotype} (n={len(df_geno)})',
                   color=color, linewidth=2.5, alpha=0.9)
            ax.fill_between(freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                           color=color, alpha=0.2)
        
        ax.set_xlabel('Frequency (Hz)', fontweight='bold')
        ax.set_ylabel('Power (dB)', fontweight='bold')
        ax.set_title(f'{condition.capitalize()}', fontsize=13, fontweight='bold')
        ax.legend(loc='upper right', frameon=True)
        ax.set_xlim([0, 100])
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    fig.suptitle(f'{state} State - {signal_type} Power Spectral Density by Condition',
                fontsize=16, fontweight='bold')
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_PSD_by_condition.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


# ===========================
# PLOT 4: Comprehensive PSD Comparison by Time Block
# ===========================

def plot_comprehensive_psd_by_timeblock(df_animal, state='NREM', signal_type='EEG'):
    """
    Create a 3-panel figure showing PSD curves for all genotype×condition combinations,
    one panel per time block (0-3h, 3-6h, Washout).
    
    Each panel shows 6 lines:
    - WT baseline, WT ambtemp, WT drugs
    - APP baseline, APP ambtemp, APP drugs
    """
    print(f"\n=== COMPREHENSIVE PSD BY TIME BLOCK ===")
    
    time_blocks = ['0-3h', '3-6h', 'Washout']
    conditions = ['baseline', 'ambtemp', 'drugs']
    genotypes = ['WT', 'APP']
    
    # Create figure with 3 panels
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    
    for ax_idx, time_block in enumerate(time_blocks):
        ax = axes[ax_idx]
        
        df_block = df_animal[df_animal['time_block'] == time_block].copy()
        
        if len(df_block) == 0:
            ax.text(0.5, 0.5, f'No data for {time_block}',
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=14)
            ax.set_title(f'{time_block}', fontsize=16, fontweight='bold')
            continue
        
        # Plot each genotype×condition combination
        legend_labels = []
        
        for genotype in genotypes:
            for condition in conditions:
                df_subset = df_block[
                    (df_block['genotype'] == genotype) & 
                    (df_block['condition'] == condition)
                ]
                
                if len(df_subset) == 0:
                    continue
                
                # Get PSDs
                psds = np.array([p for p in df_subset['mean_psd'].values])
                freqs = df_subset['freqs'].iloc[0]
                
                # ✅ CONVERT TO dB FIRST, THEN compute mean/sem
                psds_db = convert_to_db_for_plot(psds)
                mean_psd = np.mean(psds_db, axis=0)
                sem_psd = stats.sem(psds_db, axis=0)
                
                # Get color
                color = COLORS_CONDITION[genotype][condition]
                
                # Line style: solid for WT, dashed for APP
                linestyle = '-' if genotype == 'WT' else '--'
                linewidth = 2.5 if genotype == 'WT' else 2.0
                
                # Plot
                label = f'{genotype} {condition.capitalize()} (n={len(df_subset)})'
                ax.plot(freqs, mean_psd, label=label,
                       color=color, linewidth=linewidth, 
                       linestyle=linestyle, alpha=0.9)
                ax.fill_between(freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                               color=color, alpha=0.15)
                
                legend_labels.append(label)
        
        # Format axes
        ax.set_xlabel('Frequency (Hz)', fontweight='bold', fontsize=13)
        ax.set_ylabel('Power (dB)', fontweight='bold', fontsize=13)
        ax.set_title(f'{time_block}', fontsize=16, fontweight='bold')
        ax.set_xlim([0, 100])
        ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add legend
        if legend_labels:
            ax.legend(loc='upper right', frameon=True, fontsize=9,
                     framealpha=0.9, edgecolor='gray')
    
    # Overall title
    fig.suptitle(f'{state} State - {signal_type} Power Spectral Density\n'
                f'WT vs APP across All Conditions and Time Blocks',
                fontsize=18, fontweight='bold', y=1.02)
    
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_comprehensive_PSD_by_timeblock.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()


# ===========================
# PLOT 5: Comprehensive Band Power Comparison
# ===========================

def plot_comprehensive_bandpower_summary(df_animal, state='NREM', signal_type='EEG'):
    """
    Create a comprehensive summary of band powers showing all conditions, genotypes,
    and time blocks in one figure (similar to baseline plot but extended).
    
    Left panel: PSD curves for all genotype×condition (averaged across time blocks)
    Right panel: Band powers grouped by frequency band with RM-ANOVA results
    """
    print(f"\n=== COMPREHENSIVE BAND POWER SUMMARY ===")
    
    conditions = ['baseline', 'ambtemp', 'drugs']
    genotypes = ['WT', 'APP']
    
    # Create figure
    fig = plt.figure(figsize=(22, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=0.3)
    
    # LEFT PANEL: PSD Curves (averaged across all time blocks)
    ax_psd = fig.add_subplot(gs[0])
    
    for genotype in genotypes:
        for condition in conditions:
            df_subset = df_animal[
                (df_animal['genotype'] == genotype) & 
                (df_animal['condition'] == condition)
            ]
            
            if len(df_subset) == 0:
                continue
            
            # Get PSDs
            psds = np.array([p for p in df_subset['mean_psd'].values])
            freqs = df_subset['freqs'].iloc[0]
            
            # ✅ CONVERT TO dB FIRST, THEN compute mean/sem
            psds_db = convert_to_db_for_plot(psds)
            mean_psd = np.mean(psds_db, axis=0)
            sem_psd = stats.sem(psds_db, axis=0)
            
            # Get color
            color = COLORS_CONDITION[genotype][condition]
            
            # Line style
            linestyle = '-' if genotype == 'WT' else '--'
            linewidth = 2.8 if genotype == 'WT' else 2.3
            
            # Plot
            label = f'{genotype} {condition.capitalize()} (n={len(df_subset)})'
            ax_psd.plot(freqs, mean_psd, label=label,
                       color=color, linewidth=linewidth,
                       linestyle=linestyle, alpha=0.9)
            ax_psd.fill_between(freqs, mean_psd - sem_psd, mean_psd + sem_psd,
                               color=color, alpha=0.15)
    
    ax_psd.set_xlabel('Frequency (Hz)', fontweight='bold', fontsize=13)
    ax_psd.set_ylabel('Power (dB)', fontweight='bold', fontsize=13)
    ax_psd.set_title(f'{signal_type} Power Spectral Density\n(Averaged Across Time Blocks)', 
                     fontsize=14, fontweight='bold')
    ax_psd.legend(loc='upper right', frameon=True, fontsize=10,
                 framealpha=0.9, edgecolor='gray', ncol=2)
    ax_psd.set_xlim([0, 100])
    ax_psd.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax_psd.spines['top'].set_visible(False)
    ax_psd.spines['right'].set_visible(False)
    
    # RIGHT PANEL: Band Powers with RM-ANOVA
    ax_bands = fig.add_subplot(gs[1])
    
    band_names = list(FREQ_BANDS.keys())
    n_bands = len(band_names)
    n_conditions = len(conditions)
    
    # Setup bar positions
    bar_width = 0.13
    group_spacing = 0.4
    
    x_base = np.arange(n_bands) * (n_conditions * 2 * bar_width + group_spacing)
    
    # Prepare data for RM-ANOVA across all bands
    all_anova_results = {}
    
    for band in band_names:
        power_col = f'{band}_power'
        
        # Prepare data for RM-ANOVA
        anova_data = []
        for idx, row in df_animal.iterrows():
            if power_col in row and not pd.isna(row[power_col]):
                power_db = convert_power_to_db(row[power_col])
                anova_data.append({
                    'mouse_id': row['mouse_id'],
                    'genotype': row['genotype'],
                    'condition': row['condition'],
                    'time_block': row['time_block'],
                    'power_db': power_db
                })
        
        df_anova = pd.DataFrame(anova_data)
        
        if len(df_anova) > 0:
            print(f"\n  Performing RM-ANOVA for {band} band (comprehensive summary)...")
            anova_results = perform_rm_anova(
                df_anova, 
                value_col='power_db',
                within_factors=['condition', 'time_block'],
                between_factor='genotype',
                subject_col='mouse_id'
            )
            
            if 'error' not in anova_results:
                all_anova_results[band] = anova_results
                print(f"    ✓ RM-ANOVA completed for {band}")
                if 'pvalues' in anova_results:
                    sig_effects = [(k, v) for k, v in results['pvalues'].items() if v < 0.05]
                    if sig_effects:
                        print(f"    Significant effects:")
                        for effect, pval in sig_effects:
                            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                            print(f"      {effect}: p={pval:.4f} {sig}")
    
    # Plot bars for each genotype×condition
    for cond_idx, condition in enumerate(conditions):
        for geno_idx, genotype in enumerate(genotypes):
            df_subset = df_animal[
                (df_animal['genotype'] == genotype) & 
                (df_animal['condition'] == condition)
            ]
            
            if len(df_subset) == 0:
                continue
            
            # Get band powers (convert to dB for display)
            means_db = [convert_power_to_db(df_subset[f'{band}_power'].mean()) for band in band_names]
            # Calculate SEM in dB scale
            sems_db = []
            for band in band_names:
                vals = df_subset[f'{band}_power'].values
                if len(vals) > 0:
                    vals_db = convert_power_to_db(vals)
                    sems_db.append(stats.sem(vals_db))
                else:
                    sems_db.append(0)
            
            # Calculate x positions
            offset = (cond_idx * 2 + geno_idx) * bar_width
            x_positions = x_base + offset
            
            color = COLORS_CONDITION[genotype][condition]
            
            # Plot bars
            ax_bands.bar(x_positions, means_db, bar_width,
                        yerr=sems_db, color=color, alpha=0.85,
                        capsize=3, edgecolor='black', linewidth=1,
                        label=f'{genotype} {condition.capitalize()}')
    
    # Add significance markers based on RM-ANOVA results
    # Show if there's a significant genotype × condition interaction or main effect
    for band_idx, band in enumerate(band_names):
        if band in all_anova_results:
            anova_res = all_anova_results[band]
            if 'pvalues' in anova_res:
                # Check for interaction or main effects
                significant_effects = []
                for effect, pval in anova_res['pvalues'].items():
                    if pval < 0.05 and ('genotype' in effect.lower() or 'condition' in effect.lower()):
                        significant_effects.append((effect, pval))
                
                if significant_effects:
                    # Add text annotation above the bars
                    x_pos = x_base[band_idx] + (n_conditions * 2 * bar_width) / 2 - bar_width/2
                    y_max = ax_bands.get_ylim()[1]
                    
                    # Find most significant effect
                    most_sig = min(significant_effects, key=lambda x: x[1])
                    sig_marker = '***' if most_sig[1] < 0.001 else '**' if most_sig[1] < 0.01 else '*'
                    
                    ax_bands.text(x_pos, y_max * 0.95, sig_marker,
                                ha='center', va='top', fontsize=10, 
                                fontweight='bold', color='red')
    
    # Format band power panel
    ax_bands.set_xticks(x_base + (n_conditions * 2 * bar_width) / 2 - bar_width/2)
    ax_bands.set_xticklabels([b.capitalize() for b in band_names], 
                             rotation=20, ha='right', fontsize=11)
    ax_bands.set_ylabel('Power (dB)', fontweight='bold', fontsize=13)
    ax_bands.set_title('Band Power Comparison\n(* = significant RM-ANOVA effect)', 
                      fontsize=14, fontweight='bold')
    ax_bands.spines['top'].set_visible(False)
    ax_bands.spines['right'].set_visible(False)
    ax_bands.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Add compact legend
    handles, labels = ax_bands.get_legend_handles_labels()
    # Remove duplicates
    by_label = dict(zip(labels, handles))
    ax_bands.legend(by_label.values(), by_label.keys(),
                   loc='upper right', frameon=True, fontsize=8,
                   framealpha=0.9, edgecolor='gray', ncol=1)
    
    # Overall title
    fig.suptitle(f'{state} State - {signal_type} Comprehensive Power Spectrum Analysis\n'
                f'WT vs APP: All Conditions Combined (with RM-ANOVA)',
                fontsize=18, fontweight='bold')
    
    fig.tight_layout()
    
    out_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_comprehensive_summary.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {out_path}")
    plt.close()
    
    # Save RM-ANOVA results for all bands
    if all_anova_results:
        anova_summary_path = os.path.join(OUTPUT_DIR, f'{state}_{signal_type}_comprehensive_RMANOVA_summary.txt')
        with open(anova_summary_path, 'w') as f:
            f.write(f"Repeated Measures ANOVA Results - Comprehensive Summary\n")
            f.write(f"State: {state}\n")
            f.write(f"Signal: {signal_type}\n")
            f.write(f"Design: Genotype (between) × Condition × Time Block (within)\n")
            f.write(f"\n{'='*80}\n\n")
    
            
            for band, anova_results in all_anova_results.items():  # CHANGED: results -> anova_results
                band_range = FREQ_BANDS[band]
                f.write(f"\n{'='*80}\n")
                f.write(f"BAND: {band.upper()} ({band_range[0]}-{band_range[1]} Hz)\n")
                f.write(f"State: {state}\n")
                f.write(f"Signal: {signal_type}\n")
                f.write(f"\n{anova_results['summary']}\n")  # CHANGED: results -> anova_results
                f.write(f"\n")
                
                if 'pvalues' in anova_results:  # CHANGED: results -> anova_results
                    f.write(f"\nSignificant effects (p < 0.05):\n")
                    sig_effects = [(k, v) for k, v in anova_results['pvalues'].items() if v < 0.05]  # CHANGED
                    if sig_effects:
                        for effect, pval in sig_effects:
                            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                            f.write(f"  • {effect}: p={pval:.6f} {sig}\n")
                    else:
                        f.write(f"  None\n")
                f.write(f"\n")
        
        print(f"✓ Saved comprehensive RM-ANOVA summary: {anova_summary_path}")


# ===========================
# Main (Updated)
# ===========================

def main():
    print("="*70)
    print("PSD AND BAND POWER ANALYSIS")
    print("Focus: WT vs APP Comparison")
    print("="*70)
    
    # Choose state
    print("\nAvailable states: NREM, REM, Wake")
    state_input = input("Enter state (default NREM): ").strip() or "NREM"
    
    state_map = {'0': 'Wake', '1': 'NREM', '2': 'REM',
                 'wake': 'Wake', 'nrem': 'NREM', 'rem': 'REM'}
    state = state_map.get(state_input.lower(), state_input)
    
    # Load signal data
    print(f"\nLoading {state} signal data...")
    df_signals = load_signal_data(state=state)
    
    if df_signals is None or len(df_signals) == 0:
        print("⚠️  No signal data found!")
        return
    
    # REMOVED: plot_notch_filter_comparison(df_signals, state=state)
    
    # Compute PSDs for EEG
    print("\n=== Computing PSDs for EEG ===")
    df_eeg_psd = process_all_bouts_psd(df_signals, signal_col='eeg_signal')
    df_eeg_animal = compute_per_animal_psd(df_eeg_psd)
    
    # Generate plots
    print(f"\n{'='*70}")
    print("GENERATING PLOTS")
    print(f"{'='*70}\n")
    
    print("1. Baseline PSD comparison...")
    plot_baseline_psd_comparison(df_eeg_animal, state=state, signal_type='EEG')
    
    print("\n2. Band powers by condition and time block...")
    plot_band_powers_by_condition(df_eeg_animal, state=state, signal_type='EEG')
    
    print("\n3. PSD curves by condition...")
    plot_psd_curves_by_condition(df_eeg_animal, state=state, signal_type='EEG')
    
    print("\n4. Comprehensive PSD by time block...")
    plot_comprehensive_psd_by_timeblock(df_eeg_animal, state=state, signal_type='EEG')
    
    print("\n5. Comprehensive band power summary...")
    plot_comprehensive_bandpower_summary(df_eeg_animal, state=state, signal_type='EEG')
    
    # If ACh signal available, repeat for ACh
    if df_signals['ach_signal'].notna().any():
        print("\n=== Computing PSDs for ACh ===")
        df_ach_psd = process_all_bouts_psd(df_signals, signal_col='ach_signal')
        df_ach_animal = compute_per_animal_psd(df_ach_psd)
        
               
        print("\n6. ACh signal analysis...")
        plot_baseline_psd_comparison(df_ach_animal, state=state, signal_type='ACh')
        plot_band_powers_by_condition(df_ach_animal, state=state, signal_type='ACh')
        plot_psd_curves_by_condition(df_ach_animal, state=state, signal_type='ACh')
        plot_comprehensive_psd_by_timeblock(df_ach_animal, state=state, signal_type='ACh')
        plot_comprehensive_bandpower_summary(df_ach_animal, state=state, signal_type='ACh')
    
    print(f"\n{'='*70}")
    print(f"✓ All PSD plots saved to: {OUTPUT_DIR}")
    print(f"{'='*70}")
    
    # Summary
    print("\nGenerated files:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if f.endswith('.png') or f.endswith('.csv'):
            print(f"  • {f}")


if __name__ == "__main__":
    main()