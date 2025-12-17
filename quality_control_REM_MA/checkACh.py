#!/usr/bin/env python3
"""
Script to check ACh signal frequency content and determine appropriate frequency band.

This will help you verify if PLV_ACH_BAND = (0.01, 0.5) Hz is correct for your data.

Usage:
    python check_ach_frequency_band.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.signal import welch, spectrogram
from scipy import stats

# ===========================
# Configuration
# ===========================

# Same directories as your main analysis
SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"

# How many files to analyze (set to None to analyze all)
MAX_FILES = 5  # Start with 5 files to get quick results


# ===========================
# Functions
# ===========================

def load_ach_signal(mat_path):
    """
    Load ACh signal and its sampling rate from .mat file.
    
    Returns
    -------
    ach : array or None
        ACh signal
    fs_ach : float or None
        ACh sampling rate
    """
    try:
        data = sio.loadmat(mat_path, squeeze_me=True)
        
        # Try different possible ACh field names
        ach = None
        ach_key_found = None
        ach_keys = ("ne", "ach", "ACh", "Ach", "acetylcholine", "NE")
        
        for key in ach_keys:
            if key in data:
                ach = np.asarray(data[key], dtype=float).squeeze()
                ach_key_found = key
                break
        
        if ach is None:
            return None, None
        
        # Try to find sampling rate
        fs_ach = None
        freq_candidates = [
            f"{ach_key_found}_frequency",
            f"{ach_key_found}_fs",
            "ne_frequency",
            "ach_frequency",
            "fs_ne",
            "fs_ach",
            "frequency",
            "fs",
        ]
        
        for fk in freq_candidates:
            if fk in data:
                fs_ach = float(np.squeeze(data[fk]))
                break
        
        # If still no fs found, try to get EEG fs as backup
        if fs_ach is None:
            if "eeg_frequency" in data:
                fs_ach = float(np.squeeze(data["eeg_frequency"]))
                print(f"  ⚠️  Using EEG frequency ({fs_ach} Hz) for ACh (no ACh fs found)")
        
        return ach, fs_ach
        
    except Exception as e:
        print(f"  Error loading {mat_path}: {e}")
        return None, None


def analyze_ach_spectrum(ach, fs_ach, recording_name):
    """
    Compute and plot power spectrum of ACh signal.
    
    Returns
    -------
    peak_freq : float
        Frequency with maximum power
    freq_range_90pct : tuple
        Frequency range containing 90% of power
    """
    # Remove DC offset
    ach = ach - np.mean(ach)
    
    # Compute power spectral density using Welch's method
    nperseg = min(len(ach), int(4 * fs_ach))  # 4-second windows
    f, psd = welch(ach, fs=fs_ach, nperseg=nperseg, scaling='density')
    
    # Focus on low frequencies (0-10 Hz)
    freq_mask = f <= 10.0
    f = f[freq_mask]
    psd = psd[freq_mask]
    
    # Find peak frequency
    peak_idx = np.argmax(psd)
    peak_freq = f[peak_idx]
    
    # Find frequency range containing 90% of power
    cumsum = np.cumsum(psd)
    cumsum = cumsum / cumsum[-1]  # Normalize to 0-1
    
    freq_5pct = f[np.argmin(np.abs(cumsum - 0.05))]
    freq_95pct = f[np.argmin(np.abs(cumsum - 0.95))]
    
    return peak_freq, (freq_5pct, freq_95pct), f, psd


def plot_ach_analysis(recordings_data, output_path):
    """
    Create comprehensive visualization of ACh frequency content.
    
    Parameters
    ----------
    recordings_data : list of dicts
        Each dict contains: name, ach, fs_ach, peak_freq, freq_range, f, psd
    """
    n_recordings = len(recordings_data)
    
    # Create figure with multiple panels
    fig = plt.figure(figsize=(16, 12))
    
    # Panel 1: Individual power spectra (linear scale)
    ax1 = plt.subplot(3, 2, 1)
    colors = plt.cm.viridis(np.linspace(0, 1, n_recordings))
    
    for i, rec in enumerate(recordings_data):
        ax1.plot(rec['f'], rec['psd'], alpha=0.7, color=colors[i], 
                label=f"{rec['name'][:20]}... (peak: {rec['peak_freq']:.2f} Hz)")
    
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel('Power Spectral Density')
    ax1.set_title('ACh Power Spectra (Linear Scale)')
    ax1.legend(fontsize=7, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 5)  # Focus on 0-5 Hz
    
    # Panel 2: Power spectra (log scale)
    ax2 = plt.subplot(3, 2, 2)
    
    for i, rec in enumerate(recordings_data):
        ax2.semilogy(rec['f'], rec['psd'], alpha=0.7, color=colors[i])
    
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Power Spectral Density (log scale)')
    ax2.set_title('ACh Power Spectra (Log Scale)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 5)
    
    # Panel 3: Peak frequencies distribution
    ax3 = plt.subplot(3, 2, 3)
    peak_freqs = [rec['peak_freq'] for rec in recordings_data]
    
    ax3.hist(peak_freqs, bins=15, edgecolor='black', alpha=0.7, color='steelblue')
    ax3.axvline(np.mean(peak_freqs), color='red', linestyle='--', linewidth=2,
               label=f'Mean: {np.mean(peak_freqs):.2f} Hz')
    ax3.axvline(np.median(peak_freqs), color='orange', linestyle='--', linewidth=2,
               label=f'Median: {np.median(peak_freqs):.2f} Hz')
    ax3.set_xlabel('Peak Frequency (Hz)')
    ax3.set_ylabel('Count')
    ax3.set_title('Distribution of Peak Frequencies Across Recordings')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Frequency ranges (5-95 percentile)
    ax4 = plt.subplot(3, 2, 4)
    
    freq_ranges = [rec['freq_range'] for rec in recordings_data]
    low_freqs = [fr[0] for fr in freq_ranges]
    high_freqs = [fr[1] for fr in freq_ranges]
    
    x_pos = np.arange(len(recordings_data))
    
    for i in x_pos:
        ax4.plot([i, i], [low_freqs[i], high_freqs[i]], 'o-', linewidth=2, markersize=8)
    
    ax4.axhline(np.mean(low_freqs), color='blue', linestyle='--', alpha=0.5,
               label=f'Mean low: {np.mean(low_freqs):.2f} Hz')
    ax4.axhline(np.mean(high_freqs), color='red', linestyle='--', alpha=0.5,
               label=f'Mean high: {np.mean(high_freqs):.2f} Hz')
    
    # Highlight current band
    ax4.axhspan(0.01, 0.5, alpha=0.2, color='yellow', label='Current band (0.01-0.5 Hz)')
    
    ax4.set_xlabel('Recording')
    ax4.set_ylabel('Frequency (Hz)')
    ax4.set_title('Frequency Range (5-95% Power) per Recording')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([f"Rec {i+1}" for i in x_pos], rotation=45)
    
    # Panel 5: Average spectrum across all recordings
    ax5 = plt.subplot(3, 2, 5)
    
    # Interpolate all PSDs to common frequency grid
    f_common = recordings_data[0]['f']
    psd_all = []
    
    for rec in recordings_data:
        psd_interp = np.interp(f_common, rec['f'], rec['psd'])
        psd_all.append(psd_interp)
    
    psd_mean = np.mean(psd_all, axis=0)
    psd_std = np.std(psd_all, axis=0)
    
    ax5.plot(f_common, psd_mean, 'b-', linewidth=2, label='Mean')
    ax5.fill_between(f_common, psd_mean - psd_std, psd_mean + psd_std, 
                     alpha=0.3, color='blue', label='±1 SD')
    
    # Mark current band
    ax5.axvspan(0.01, 0.5, alpha=0.2, color='yellow', label='Current band')
    
    ax5.set_xlabel('Frequency (Hz)')
    ax5.set_ylabel('Power Spectral Density')
    ax5.set_title('Average ACh Power Spectrum (Mean ± SD)')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim(0, 5)
    
    # Panel 6: Summary statistics and recommendations
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('off')
    
    # Compute summary statistics
    mean_peak = np.mean(peak_freqs)
    median_peak = np.median(peak_freqs)
    mean_low = np.mean(low_freqs)
    mean_high = np.mean(high_freqs)
    
    # Recommendation logic
    current_band = (0.01, 0.5)
    
    if mean_peak < 0.5:
        recommendation = "✅ CURRENT BAND OK"
        rec_band = current_band
        color = 'green'
    elif mean_peak < 1.0:
        recommendation = "⚠️ CONSIDER ADJUSTING"
        rec_band = (0.1, 2.0)
        color = 'orange'
    elif mean_peak < 2.0:
        recommendation = "⚠️ SHOULD ADJUST"
        rec_band = (0.5, 3.0)
        color = 'orange'
    else:
        recommendation = "❌ DEFINITELY ADJUST"
        rec_band = (1.0, 5.0)
        color = 'red'
    
    summary_text = f"""
    SUMMARY STATISTICS
    {'='*40}
    
    Number of recordings: {n_recordings}
    
    Peak Frequencies:
      Mean:   {mean_peak:.3f} Hz
      Median: {median_peak:.3f} Hz
      Range:  {min(peak_freqs):.3f} - {max(peak_freqs):.3f} Hz
    
    90% Power Range:
      Lower bound: {mean_low:.3f} Hz (mean)
      Upper bound: {mean_high:.3f} Hz (mean)
    
    {'='*40}
    CURRENT SETTING
    {'='*40}
    
    PLV_ACH_BAND = ({current_band[0]}, {current_band[1]}) Hz
    
    {'='*40}
    RECOMMENDATION
    {'='*40}
    
    {recommendation}
    
    Suggested band: ({rec_band[0]}, {rec_band[1]}) Hz
    
    {'='*40}
    INTERPRETATION
    {'='*40}
    
    • If peak < 0.5 Hz:
      Current band (0.01-0.5 Hz) is appropriate.
      ACh shows very slow fluctuations.
    
    • If peak 0.5-2 Hz:
      Consider (0.1-2 Hz) or (0.5-3 Hz).
      ACh has faster dynamics than expected.
    
    • If peak > 2 Hz:
      Check if ACh signal is correctly loaded.
      May need (1-5 Hz) range.
    
    {'='*40}
    NEXT STEPS
    {'='*40}
    
    1. Review power spectra above
    2. If recommendation suggests change:
       - Edit sleep_wavelet_analysis.py
       - Line ~66: PLV_ACH_BAND = ...
       - Set to recommended band
       - Re-run batch analysis
    
    3. Compare results before/after change
    """
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
            fontsize=9, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle('ACh Signal Frequency Analysis - Verification of PLV_ACH_BAND Setting',
                fontsize=14, fontweight='bold')
    fig.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved analysis to: {output_path}")
    
    plt.show()
    
    return mean_peak, rec_band, recommendation


# ===========================
# Main Analysis
# ===========================

def main():
    print("="*70)
    print("ACh FREQUENCY BAND VERIFICATION")
    print("="*70)
    print("\nThis script will analyze your ACh signal to determine if the")
    print("current frequency band (0.01-0.5 Hz) is appropriate.")
    print()
    
    # Find .mat files
    mat_files = sorted(glob.glob(os.path.join(SIGNAL_DIR, "*.mat")))
    
    if not mat_files:
        print(f"❌ No .mat files found in {SIGNAL_DIR}")
        return
    
    print(f"Found {len(mat_files)} .mat files")
    
    if MAX_FILES is not None and len(mat_files) > MAX_FILES:
        print(f"Analyzing first {MAX_FILES} files (set MAX_FILES=None to analyze all)")
        mat_files = mat_files[:MAX_FILES]
    
    print()
    
    # Analyze each file
    recordings_data = []
    
    for i, mat_path in enumerate(mat_files, 1):
        base = os.path.basename(mat_path)
        print(f"[{i}/{len(mat_files)}] Loading {base}...", end=" ")
        
        ach, fs_ach = load_ach_signal(mat_path)
        
        if ach is None:
            print("❌ No ACh data found")
            continue
        
        if fs_ach is None:
            print("❌ No sampling rate found")
            continue
        
        print(f"✓ ACh signal loaded (fs={fs_ach} Hz, length={len(ach)/fs_ach:.1f}s)")
        
        # Analyze spectrum
        peak_freq, freq_range, f, psd = analyze_ach_spectrum(ach, fs_ach, base)
        
        recordings_data.append({
            'name': base,
            'ach': ach,
            'fs_ach': fs_ach,
            'peak_freq': peak_freq,
            'freq_range': freq_range,
            'f': f,
            'psd': psd
        })
        
        print(f"     Peak frequency: {peak_freq:.3f} Hz")
        print(f"     90% power range: {freq_range[0]:.3f} - {freq_range[1]:.3f} Hz")
    
    if not recordings_data:
        print("\n❌ No ACh data could be loaded from any files!")
        print("\nPossible reasons:")
        print("  1. ACh field name is different (check .mat file contents)")
        print("  2. ACh data not present in these recordings")
        print("  3. Files are corrupted or in wrong format")
        return
    
    print(f"\n{'='*70}")
    print(f"Successfully analyzed {len(recordings_data)} recordings")
    print(f"{'='*70}\n")
    
    # Create comprehensive visualization
    output_path = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/ach_frequency_analysis.png"
    
    mean_peak, recommended_band, recommendation = plot_ach_analysis(recordings_data, output_path)
    
    # Print summary to console
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"\nAnalyzed {len(recordings_data)} recordings")
    print(f"Average peak frequency: {mean_peak:.3f} Hz")
    print(f"\nCurrent setting: PLV_ACH_BAND = (0.01, 0.5) Hz")
    print(f"\n{recommendation}")
    print(f"Recommended band: ({recommended_band[0]}, {recommended_band[1]}) Hz")
    
    if "OK" in recommendation:
        print(f"\n✅ Your current frequency band is appropriate!")
        print(f"   No changes needed.")
    elif "CONSIDER" in recommendation or "SHOULD" in recommendation:
        print(f"\n⚠️  Consider updating your frequency band:")
        print(f"\n   1. Open: sleep_wavelet_analysis.py")
        print(f"   2. Find line ~66: PLV_ACH_BAND = (0.01, 0.5)")
        print(f"   3. Change to: PLV_ACH_BAND = {recommended_band}")
        print(f"   4. Re-run batch analysis")
        print(f"\n   This will improve PLV accuracy by matching ACh dynamics.")
    else:
        print(f"\n❌ Your current frequency band needs adjustment!")
        print(f"\n   IMPORTANT: Update your analysis:")
        print(f"   1. Open: sleep_wavelet_analysis.py")
        print(f"   2. Find line ~66: PLV_ACH_BAND = (0.01, 0.5)")
        print(f"   3. Change to: PLV_ACH_BAND = {recommended_band}")
        print(f"   4. Re-run batch analysis")
    
    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()