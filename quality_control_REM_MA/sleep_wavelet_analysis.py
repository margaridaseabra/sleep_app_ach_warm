#!/usr/bin/env python3
"""
Interactive bout viewer with PLV/PAC analysis for ACh-EEG coupling.

New features:
- Press 's' to SELECT current bout for analysis
- Press 'u' to UNSELECT current bout
- Press 'c' to COMPUTE PLV/PAC on all selected bouts
- Press 'x' to EXPORT selected bout indices to CSV
- Results saved per recording

Navigation:
    → / ↓ / n : next bout
    ← / ↑ / p : previous bout
    s         : select current bout
    u         : unselect current bout
    c         : compute PLV/PAC on selected bouts
    x         : export selected bout list
    q / Esc   : quit viewer & script
"""

import os
import glob
import math
import numpy as np
import pandas as pd
import scipy
import scipy.io as sio
import scipy.fftpack
import scipy.signal
from scipy.signal import welch, ShortTimeFFT, hilbert, butter, filtfilt, resample
from scipy.signal.windows import hamming
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import pywt
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.widgets import SpanSelector

# ===========================
# Configuration
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"

PRE_SEC = 20.0
POST_SEC = 30.0
MIN_BOUT_LEN_SEC = 5.0

DELTA_BAND = (1.0, 4.0)
THETA_BAND = (4.0, 8.0)

HF_BAND = (40.0, 80.0)
HF_PERCENTILE = 95.0
HF_MIN_BURST_DUR = 0.050

CWT_MIN_FREQ = 1
CWT_MAX_FREQ = 12
CWT_N_FREQS = 40
CWT_WAVELET = "cmor0.5-1.0"

STATE_COLORS = {
    0: (0.90, 0.80, 0.95, 0.30),  # Wake
    1: (1.00, 0.60, 0.60, 0.30),  # NREM
    2: (0.60, 1.00, 0.60, 0.30),  # REM
    15: (0.70, 0.70, 1.00, 0.30), # MA
}
STATE_NAMES = {0: "Wake", 1: "NREM", 2: "REM", 15: "MA"}

BAND_PRESETS = {
    "1": ("0–5 Hz band",   0.2, 5.0),
    "2": ("0–15 Hz band",  0.2, 15.0),
    "3": ("0–30 Hz band",  0.2, 30.0),
    "4": ("theta 6–10 Hz",   6.0, 10.0),
    "5": ("spindle 9–16 Hz", 9.0, 16.0),
    "6": ("beta 15–30 Hz",   15.0, 30.0),
    "7": ("gamma 30–80 Hz",  30.0, 80.0),
}
DEFAULT_BAND_KEY = "3"

# ===========================
# PLV/PAC Configuration
# ===========================

# PLV: Phase-locking between EEG theta and ACh signal
PLV_EEG_BAND = (4.0, 8.0)   # theta in EEG
PLV_ACH_BAND = (0.01, 0.5)  # slow fluctuations in ACh (adjust based on your data)

# PAC: Phase-amplitude coupling
PAC_PHASE_BAND = (4.0, 8.0)    # theta phase (low freq)
PAC_AMP_BAND = (30.0, 80.0)    # gamma amplitude (high freq)

# Output directory for analysis results
RESULTS_DIR = "/Users/margaridaseabra/sleep_app_ach_warm/quality_control_REM_MA/plv_pac_results"
os.makedirs(RESULTS_DIR, exist_ok=True)


# ===========================
# Signal Processing Helpers
# ===========================

def bandpass_filter(sig, fs, low, high, order=4):
    """Butterworth bandpass filter."""
    nyq = fs / 2.0
    if low >= nyq or high >= nyq:
        raise ValueError(f"Filter frequencies {low}-{high} Hz exceed Nyquist {nyq} Hz")
    if low <= 0:
        low = 0.01  # avoid zero frequency
    b, a = butter(order, [low / nyq, high / nyq], btype='band')
    return filtfilt(b, a, sig)


def compute_plv(sig1, sig2, fs, band1=None, band2=None):
    """
    Compute Phase-Locking Value between two signals.
    
    Parameters
    ----------
    sig1, sig2 : array
        Input signals (must be same length)
    fs : float
        Sampling rate
    band1, band2 : tuple (low, high) or None
        If provided, bandpass filter signals first
        
    Returns
    -------
    plv : float
        Phase-locking value [0, 1]
    phase_diff : array
        Instantaneous phase difference (radians)
    """
    if len(sig1) != len(sig2):
        raise ValueError("Signals must have same length for PLV")
    
    # Bandpass filter if requested
    if band1 is not None:
        sig1 = bandpass_filter(sig1, fs, *band1)
    if band2 is not None:
        sig2 = bandpass_filter(sig2, fs, *band2)
    
    # Analytic signals via Hilbert transform
    analytic1 = hilbert(sig1)
    analytic2 = hilbert(sig2)
    
    # Instantaneous phases
    phase1 = np.angle(analytic1)
    phase2 = np.angle(analytic2)
    
    # Phase difference
    phase_diff = phase1 - phase2
    
    # PLV = |mean(exp(i * phase_diff))|
    plv = np.abs(np.mean(np.exp(1j * phase_diff)))
    
    return plv, phase_diff


def compute_pac_mi(sig, fs, phase_band, amp_band, n_bins=18):
    """
    Compute Phase-Amplitude Coupling using Modulation Index (Tort et al. 2010).
    
    Measures how much the amplitude of a high-frequency band is modulated
    by the phase of a low-frequency band.
    
    Parameters
    ----------
    sig : array
        Input signal (EEG)
    fs : float
        Sampling rate
    phase_band : tuple (low, high)
        Frequency band for phase (e.g., theta 4-8 Hz)
    amp_band : tuple (low, high)
        Frequency band for amplitude (e.g., gamma 30-80 Hz)
    n_bins : int
        Number of phase bins for histogram
        
    Returns
    -------
    mi : float
        Modulation index [0, 1] (0=no coupling, 1=max coupling)
    mean_amp_per_phase : array
        Mean amplitude in each phase bin
    phase_bins : array
        Phase bin centers (radians)
    """
    # Extract phase from low-freq band
    sig_phase = bandpass_filter(sig, fs, *phase_band)
    analytic_phase = hilbert(sig_phase)
    phase = np.angle(analytic_phase)
    
    # Extract amplitude envelope from high-freq band
    sig_amp = bandpass_filter(sig, fs, *amp_band)
    analytic_amp = hilbert(sig_amp)
    amplitude = np.abs(analytic_amp)
    
    # Bin phases into n_bins
    phase_bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    phase_centers = (phase_bins[:-1] + phase_bins[1:]) / 2
    
    # Compute mean amplitude per phase bin
    mean_amp = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (phase >= phase_bins[i]) & (phase < phase_bins[i + 1])
        if np.any(mask):
            mean_amp[i] = amplitude[mask].mean()
        else:
            mean_amp[i] = 0.0
    
    # Normalize to probability distribution
    p = mean_amp / (mean_amp.sum() + 1e-12)
    
    # Modulation Index = KL divergence from uniform distribution
    uniform = np.ones(n_bins) / n_bins
    kl = np.sum(p * np.log((p + 1e-12) / (uniform + 1e-12)))
    mi = kl / np.log(n_bins)  # normalize to [0, 1]
    
    return mi, mean_amp, phase_centers


def compute_pac_cross(eeg, ach, fs, phase_band_eeg, amp_band_ach, n_bins=18):
    """
    Cross-signal PAC: phase of EEG modulates amplitude of ACh.
    
    Returns
    -------
    mi : float
        Modulation index
    mean_amp : array
        Mean ACh amplitude per EEG phase bin
    phase_centers : array
        Phase bin centers
    """
    # Phase from EEG
    sig_phase = bandpass_filter(eeg, fs, *phase_band_eeg)
    analytic_phase = hilbert(sig_phase)
    phase = np.angle(analytic_phase)
    
    # Amplitude from ACh
    sig_amp = bandpass_filter(ach, fs, *amp_band_ach)
    analytic_amp = hilbert(sig_amp)
    amplitude = np.abs(analytic_amp)
    
    # Bin
    phase_bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    phase_centers = (phase_bins[:-1] + phase_bins[1:]) / 2
    
    mean_amp = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (phase >= phase_bins[i]) & (phase < phase_bins[i + 1])
        if np.any(mask):
            mean_amp[i] = amplitude[mask].mean()
    
    p = mean_amp / (mean_amp.sum() + 1e-12)
    uniform = np.ones(n_bins) / n_bins
    kl = np.sum(p * np.log((p + 1e-12) / (uniform + 1e-12)))
    mi = kl / np.log(n_bins)
    
    return mi, mean_amp, phase_centers


# ===========================
# Keep all your existing utility functions
# ===========================

def load_mat_signal(mat_path):
    """Load EEG, EMG, sampling rate, and ACh (ne) + its sampling rate."""
    data = sio.loadmat(mat_path, squeeze_me=True)

    eeg = np.asarray(data["eeg"], dtype=float)
    emg = np.asarray(data["emg"], dtype=float)
    fs = float(np.squeeze(data["eeg_frequency"]))

    # ACh / NE trace
    ach = None
    ach_key_found = None
    ach_keys_tried = ("ne", "ach", "ACh", "Ach", "acetylcholine")
    for key in ach_keys_tried:
        if key in data:
            ach = np.asarray(data[key], dtype=float).squeeze()
            ach_key_found = key
            break

    fs_ach = None
    if ach is not None:
        freq_candidates = [
            f"{ach_key_found}_frequency",
            f"{ach_key_found}_fs",
            "ne_frequency",
            "ach_frequency",
            "fs_ne",
            "fs_ach",
        ]
        for fk in freq_candidates:
            if fk in data:
                fs_ach = float(np.squeeze(data[fk]))
                break

        if fs_ach is None:
            print(f"Found ACh trace '{ach_key_found}' but no fs. Using EEG fs ({fs} Hz)")
            fs_ach = fs
        else:
            print(f"Found ACh trace '{ach_key_found}' with fs = {fs_ach} Hz")
    else:
        print(f"No ACh/NE trace found (tried: {', '.join(ach_keys_tried)})")

    return eeg, emg, fs, ach, fs_ach


def load_scores_1hz(csv_path):
    """Load per-second scoring from CSV."""
    df = pd.read_csv(csv_path)
    num_cols = list(df.select_dtypes(include="number").columns)
    if not num_cols:
        raise ValueError(f"No numeric columns in {csv_path}")

    score_col = None
    for c in num_cols:
        if "score" in c.lower() or "state" in c.lower():
            score_col = c
            break
    if score_col is None:
        score_col = min(num_cols, key=lambda c: df[c].nunique())

    states = df[score_col].to_numpy().astype(int)

    time_col = None
    for c in num_cols:
        if c == score_col:
            continue
        if "time" in c.lower() or "sec" in c.lower():
            time_col = c
            break

    if time_col is not None:
        time_s = df[time_col].to_numpy().astype(float)
    else:
        time_s = np.arange(len(states), dtype=float)

    if np.any(np.diff(time_s) <= 0):
        time_s = np.arange(len(states), dtype=float)

    return time_s, states


def find_state_bouts(states, min_len_sec=1):
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


def compute_psd(sig, fs):
    """Welch PSD."""
    nperseg = min(len(sig), int(4 * fs))
    f, Pxx = welch(sig, fs=fs, nperseg=nperseg)
    return f, Pxx


def generate_wavelet_fourier(len_wavelet, f_start, f_stop, deltafreq,
                             sample_rate, f0, normalisation):
    """Compute wavelet coefficients and Fourier transform."""
    freqs = np.arange(f_start, f_stop, deltafreq)
    scales = f0 / freqs * sample_rate
    xi = np.arange(-len_wavelet / 2., len_wavelet / 2.)
    xsd = xi[:, np.newaxis] / scales
    wavelet_coefs = np.exp(1j * 2. * np.pi * f0 * xsd) * np.exp(-xsd**2 / 2.)
    weighting_function = lambda x: x ** (-(1.0 + normalisation))
    weighted_wavelet_coefs = wavelet_coefs * weighting_function(scales[np.newaxis, :])
    wf = scipy.fftpack.fft(weighted_wavelet_coefs, axis=0)
    wf = wf.conj()
    return wf


def compute_scalogram_ephyviewer(data, min_freq, max_freq, freq_resolution,
                                 fs, f0=1, exp_corr=0,
                                 time_smooth=0.5, wanted_size=5.):
    """Ephyviewer-style Morlet scalogram."""
    n_samples = len(data)
    len_wavelet = int(2 ** np.ceil(np.log(wanted_size * fs) / np.log(2)))
    sig_chunk_size = wanted_size * fs
    downsample_ratio = int(np.ceil(sig_chunk_size / len_wavelet))
    sig_chunk_size = downsample_ratio * len_wavelet
    sub_sample_rate = fs / downsample_ratio

    wavelet_fourrier = generate_wavelet_fourier(
        len_wavelet, min_freq, max_freq, freq_resolution,
        sub_sample_rate, f0, exp_corr
    )

    if downsample_ratio > 1:
        n = 8
        q = downsample_ratio
        filter_sos = scipy.signal.cheby1(n, 0.05, 0.8 / q, output="sos")
    else:
        filter_sos = None

    i_start = 0
    if downsample_ratio > 1:
        i_start = i_start - (i_start % downsample_ratio)
    i_start = max(0, min(i_start, n_samples))
    if downsample_ratio > 1:
        i_start = i_start - (i_start % downsample_ratio)

    i_stop = min(i_start + sig_chunk_size, n_samples)
    sigs_chunk = data[i_start:i_stop]
    if sigs_chunk.dtype != "float32":
        sigs_chunk = sigs_chunk.astype("float32")

    if downsample_ratio > 1:
        small_sig = scipy.signal.sosfiltfilt(filter_sos, sigs_chunk)
        small_sig = small_sig[::downsample_ratio].copy()
    else:
        small_sig = sigs_chunk.copy()

    left_pad = 0
    if small_sig.shape[0] != wavelet_fourrier.shape[0]:
        z = np.zeros(wavelet_fourrier.shape[0], dtype=small_sig.dtype)
        left_pad = wavelet_fourrier.shape[0] - small_sig.shape[0]
        z[:small_sig.shape[0]] = small_sig
        small_sig = z

    small_sig -= small_sig.mean()
    small_sig_f = scipy.fftpack.fft(small_sig)
    wt_tmp = scipy.fftpack.ifft(
        small_sig_f[:, np.newaxis] * wavelet_fourrier, axis=0
    )
    wt = scipy.fftpack.fftshift(wt_tmp, axes=[0])
    wt = np.abs(wt).astype("float32")
    wt = 10 * np.log10(wt + 1e-20)

    if time_smooth > 0:
        sigma = time_smooth * sub_sample_rate
        wt = gaussian_filter(wt, sigma=[sigma, 0], mode="reflect")

    if left_pad > 0:
        wt = wt[:-left_pad]

    return wt.T


def compute_wavelet_scalogram_matlab_like(
    sig, fs, f_min=1.0, f_max=12.0, df=0.5,
    f0=1.0, exp_corr=0.0, time_smooth=0.0,
    baseline_mode="window"):
    """Wavelet scalogram similar to MATLAB code."""
    duration_sec = len(sig) / fs
    scalogram_db = compute_scalogram_ephyviewer(
        sig, min_freq=f_min, max_freq=f_max,
        freq_resolution=df, fs=fs, f0=f0,
        exp_corr=exp_corr, time_smooth=time_smooth,
        wanted_size=duration_sec
    )
    freqs = np.arange(f_min, f_max, df)

    if baseline_mode == "window":
        baseline_db = scalogram_db.mean(axis=1, keepdims=True)
        scalogram_db = scalogram_db - baseline_db

    return freqs, scalogram_db


def compute_wavelet_scalogram_pretty(
    sig, fs, f_min, f_max, df=0.25,
    time_smooth=0.15, return_linear=False):
    """Pretty scalogram for detail panels."""
    duration_sec = len(sig) / fs
    scalogram_db = compute_scalogram_ephyviewer(
        sig, min_freq=f_min, max_freq=f_max,
        freq_resolution=df, fs=fs, f0=1.0,
        exp_corr=0.0, time_smooth=time_smooth,
        wanted_size=duration_sec
    )
    freqs = np.arange(f_min, f_max, df)
    power_lin = 10 ** (scalogram_db / 10.0)
    low = np.percentile(power_lin, 5, axis=1, keepdims=True)
    high = np.percentile(power_lin, 95, axis=1, keepdims=True)
    power_norm = (power_lin - low) / (high - low + 1e-12)
    power_norm = np.clip(power_norm, 0.0, 1.0)

    if return_linear:
        return freqs, power_norm, power_lin
    else:
        return freqs, power_norm


def compute_spectrogram_and_ratio(eeg_win, fs, win_start_sec, window_duration=5.0, mfft=None):
    """Short-time FFT spectrogram and theta/delta ratio."""
    nperseg = round(fs * window_duration)
    hop = round(nperseg / 2)
    window = hamming(nperseg)

    SFT = ShortTimeFFT(
        window, hop=hop, fs=fs, fft_mode="onesided",
        mfft=mfft, scale_to="psd"
    )

    Sx = SFT.spectrogram(eeg_win)
    n_frames = Sx.shape[1]
    win_duration = len(eeg_win) / fs
    time = np.linspace(win_start_sec, win_start_sec + win_duration, n_frames, endpoint=False)

    frequencies = SFT.f
    freq_mask = frequencies <= 30.0
    frequencies = frequencies[freq_mask]
    Sx = Sx[freq_mask, :]
    Sx_db = 10 * np.log10(Sx + 1e-20)

    delta_mask = np.where((frequencies > 1) & (frequencies <= 4))[0]
    theta_mask = np.where((frequencies > 4) & (frequencies <= 8))[0]

    delta_power = np.mean(Sx_db[delta_mask, :], axis=0)
    theta_power = np.mean(Sx_db[theta_mask, :], axis=0)
    theta_delta_ratio = delta_power / (theta_power + 1e-12)
    theta_delta_ratio = gaussian_filter1d(theta_delta_ratio, 4)
    Sx_db = gaussian_filter(Sx_db, sigma=4)

    return frequencies, time, Sx_db, theta_delta_ratio


def bout_core_samples(start_s, end_s, fs):
    """Middle 50% of bout in samples."""
    duration = end_s - start_s
    if duration <= 0:
        return 0, 0
    cs = start_s + 0.25 * duration
    ce = start_s + 0.75 * duration
    return int(cs * fs), int(ce * fs)


def shade_states(ax, time_s, states, win_start_sec, win_end_sec):
    """Shade background by state."""
    start_idx = int(max(0, np.searchsorted(time_s, win_start_sec, side="left")))
    end_idx = int(min(len(time_s), np.searchsorted(time_s, win_end_sec, side="right")))

    for i in range(start_idx, end_idx):
        s = int(states[i])
        color = STATE_COLORS.get(s)
        if color is None:
            continue
        t0 = time_s[i]
        t1 = time_s[i + 1] if i + 1 < len(time_s) else t0 + 1.0
        left = max(t0, win_start_sec)
        right = min(t1, win_end_sec)
        if right > left:
            ax.axvspan(left, right, color=color, zorder=0)


# ===========================
# Modified Viewer Class with Selection
# ===========================

class BoutViewerWithSelection:
    """Interactive viewer with bout selection and PLV/PAC computation."""
    
    def __init__(self, base, eeg, emg, fs, time_s, states, bouts_idx, state_label, ach=None, fs_ach=None):
        self.base = base
        self.eeg = eeg
        self.emg = emg
        self.fs = fs
        self.time_s = time_s
        self.states = states
        self.bouts_idx = bouts_idx
        self.state_label = state_label
        self.idx = 0
        self.stop_all = False
        self.ach = ach
        self.fs_ach = fs_ach if fs_ach is not None else fs

        # detail-wavelet settings
        self.detail_band_key = DEFAULT_BAND_KEY  # key into BAND_PRESETS
        self.span_selector = None  # will be created on first update

        # Selection state
        self.selected_bouts = []  # indices of selected bouts

        self.fig = plt.figure(figsize=(12, 11))
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.update_plot()

    def on_key(self, event):
        """Handle keyboard events."""
        # Navigation
        if event.key in ["right", "down", "n"]:
            if self.idx < len(self.bouts_idx) - 1:
                self.idx += 1
                self.update_plot()
        elif event.key in ["left", "up", "p"]:
            if self.idx > 0:
                self.idx -= 1
                self.update_plot()
        elif event.key in ["q", "escape"]:
            self.stop_all = True
            plt.close(self.fig)
        
        # Detail-band shortcuts
        elif event.key in ("l", "h", "a") or event.key in BAND_PRESETS:
            if event.key == "l":
                key = "1"   # 0–5 Hz
            elif event.key == "h":
                key = "5"   # spindle
            elif event.key == "a":
                key = "3"   # 0–30 Hz
            else:
                key = event.key

            if key in BAND_PRESETS:
                self.detail_band_key = key
                band_label, f_min, f_max = BAND_PRESETS[key]
                print(f"Detail band set to {band_label} [{f_min}-{f_max} Hz]")

        # Selection controls
        elif event.key == 's':
            # Select current bout
            if self.idx not in self.selected_bouts:
                self.selected_bouts.append(self.idx)
                print(f"✓ Bout {self.idx + 1} SELECTED (total: {len(self.selected_bouts)})")
            else:
                print(f"  Bout {self.idx + 1} already selected")
            self.update_plot()
        
        elif event.key == 'u':
            # Unselect current bout
            if self.idx in self.selected_bouts:
                self.selected_bouts.remove(self.idx)
                print(f"✗ Bout {self.idx + 1} UNSELECTED (total: {len(self.selected_bouts)})")
            else:
                print(f"  Bout {self.idx + 1} not selected")
            self.update_plot()
        
        elif event.key == 'c':
            # Compute PLV/PAC
            if self.ach is None:
                print("⚠️  No ACh data available for analysis")
                return
            if not self.selected_bouts:
                print("⚠️  No bouts selected. Press 's' to select current bout.")
                return
            self.compute_plv_pac_analysis()
        
        elif event.key == 'x':
            # Export selected bout list
            self.export_selected_bouts()

    def on_select_span(self, tmin, tmax):
        """Callback when user drags a horizontal span on the EEG axis."""
        if tmax < tmin:
            tmin, tmax = tmax, tmin
        if tmax - tmin < 0.3:  # ignore very tiny selections
            return
        self.open_detail_wavelet(tmin, tmax)

    def open_detail_wavelet(self, t0, t1):
        """Create a new figure with a pretty Morlet scalogram for [t0, t1]."""
        fs = self.fs
        n_samples = len(self.eeg)

        i0 = max(0, int(t0 * fs))
        i1 = min(n_samples, int(t1 * fs))
        if i1 <= i0 + 10:
            return

        seg = self.eeg[i0:i1]

        # choose frequency band from presets
        band_label, f_min, f_max = BAND_PRESETS.get(
            self.detail_band_key, BAND_PRESETS[DEFAULT_BAND_KEY]
        )

        # power_norm: 0–1, power_lin: absolute (linear)
        freqs, power_norm, power_lin = compute_wavelet_scalogram_pretty(
            seg, fs, f_min=f_min, f_max=f_max,
            df=0.25, time_smooth=0.15, return_linear=True
        )

        n_freqs, n_times = power_norm.shape
        t_edges = np.linspace(t0, t1, n_times + 1)

        # frequency edges for pcolormesh
        if len(freqs) > 1:
            df_eff = freqs[1] - freqs[0]
        else:
            df_eff = 1.0
        f_edges = np.concatenate([freqs - df_eff / 2, [freqs[-1] + df_eff / 2]])

        # Figure with 3 panels: scalogram + mean spectrum + histogram
        fig, (ax_wv, ax_prof, ax_hist) = plt.subplots(
            3, 1, figsize=(7, 6),
            gridspec_kw={"height_ratios": [3, 1, 1]},
        )

        # --- top: scalogram (0–1) ---
        im = ax_wv.pcolormesh(
            t_edges,
            f_edges,
            power_norm,
            shading="auto",
            cmap="turbo",
            vmin=0.0,
            vmax=1.0,
        )
        ax_wv.set_ylabel("Frequency (Hz)")
        ax_wv.set_ylim(f_min, f_max)
        ax_wv.set_xlim(t0, t1)

        cbar = fig.colorbar(im, ax=ax_wv, pad=0.02)
        cbar.set_label("Relative power (0–1) per freq")

        # --- middle: mean power vs frequency in dB ---
        mean_power_lin = power_lin.mean(axis=1)  # average over time
        mean_power_db = 10 * np.log10(mean_power_lin + 1e-20)

        ax_prof.plot(freqs, mean_power_db, "-")
        ax_prof.set_xlabel("Frequency (Hz)")
        ax_prof.set_ylabel("Mean power (dB)")
        ax_prof.set_xlim(f_min, f_max)

        # mark dominant frequency
        dom_idx = np.argmax(mean_power_db)
        dom_freq = freqs[dom_idx]
        ax_prof.axvline(dom_freq, color="red", linestyle="--", alpha=0.7,
                        label=f"Peak ≈ {dom_freq:.2f} Hz")
        ax_prof.legend(fontsize=8, loc="best")

        # --- bottom: histogram of peak frequencies ---
        peak_idx = np.argmax(power_lin, axis=0)
        peak_freqs = freqs[peak_idx]

        n_bins = 12
        bins = np.linspace(f_min, f_max, n_bins + 1)
        counts, edges = np.histogram(peak_freqs, bins=bins)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        prop = counts / counts.sum() if counts.sum() > 0 else counts

        ax_hist.bar(
            bin_centers,
            prop,
            width=np.diff(edges),
            align="center",
            edgecolor="black",
        )

        mean_peak = peak_freqs.mean()
        med_peak = np.median(peak_freqs)
        ax_hist.axvline(mean_peak, color="red", linestyle="--",
                        label=f"Mean {mean_peak:.2f} Hz")
        ax_hist.axvline(med_peak, color="black", linestyle=":",
                        label=f"Median {med_peak:.2f} Hz")

        ax_hist.set_xlabel("Peak frequency per time bin (Hz)")
        ax_hist.set_ylabel("Proportion of time")
        ax_hist.set_xlim(f_min, f_max)
        ax_hist.legend(fontsize=8, loc="best")

        fig.suptitle(f"Detail wavelet {band_label}\n{t0:.2f}–{t1:.2f} s", fontsize=11)
        fig.tight_layout()
        fig.show()

    def update_plot(self):
        """Update the plot for the current bout."""
        self.fig.clear()
        n_bouts = len(self.bouts_idx)

        # 5 panels: EEG, EMG, spectrogram, wavelet scalogram, PSD
        ax1 = self.fig.add_subplot(5, 1, 1)
        ax2 = self.fig.add_subplot(5, 1, 2, sharex=ax1)
        ax3 = self.fig.add_subplot(5, 1, 3, sharex=ax1)
        ax4 = self.fig.add_subplot(5, 1, 4, sharex=ax1)
        ax5 = self.fig.add_subplot(5, 1, 5)

        # SpanSelector for interactive detail-wavelet popup
        if self.span_selector is None:
            self.span_selector = SpanSelector(
                ax1,
                self.on_select_span,
                "horizontal",
                useblit=True,
                props=dict(alpha=0.3, facecolor="yellow"),
                interactive=True,
            )
        else:
            self.span_selector.ax = ax1

        start_idx, end_idx = self.bouts_idx[self.idx]
        start_s = float(self.time_s[start_idx])
        end_s = float(self.time_s[end_idx - 1] + 1.0)

        fs = self.fs
        eeg = self.eeg
        emg = self.emg

        n_samples = len(eeg)
        total_sec = n_samples / fs

        win_start_sec = max(0.0, start_s - PRE_SEC)
        win_end_sec = min(total_sec, end_s + POST_SEC)
        win_start_sample = int(win_start_sec * fs)
        win_end_sample = int(win_end_sec * fs)

        eeg_win = eeg[win_start_sample:win_end_sample]
        emg_win = emg[win_start_sample:win_end_sample]
        t_win = np.arange(win_start_sample, win_end_sample) / fs

        for ax in [ax1, ax2, ax3, ax4]:
            shade_states(ax, self.time_s, self.states, win_start_sec, win_end_sec)
            ax.set_xlim(win_start_sec, win_end_sec)

        ax1.plot(t_win, eeg_win, linewidth=0.5, color="black")
        ax1.set_ylabel("EEG (a.u.)")

        legend_patches = []
        for s in sorted(np.unique(self.states)):
            name = STATE_NAMES.get(int(s), str(s))
            color = STATE_COLORS.get(int(s))
            if color is not None:
                legend_patches.append(Patch(facecolor=color, edgecolor="none", label=name))
        if legend_patches:
            ax1.legend(handles=legend_patches, loc="upper right", fontsize=8)

        ax2.plot(t_win, emg_win, linewidth=0.5, color="black")
        ax2.set_ylabel("EMG (a.u.)")

        # Spectrogram + theta/delta ratio
        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )
        im = ax3.pcolormesh(t_spec, freqs_spec, Sx_db, shading="auto", cmap="viridis")
        ax3.set_ylabel("Freq (Hz)")
        ax3.set_ylim(0, 30)
        win_end_sec_actual = win_start_sec + len(eeg_win) / fs
        ax3.set_xlim(win_start_sec, win_end_sec_actual)

        cbar = self.fig.colorbar(im, ax=ax3)
        cbar.set_label("Power (dB)")

        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.4, linewidth=1.0)
        ax3b.set_ylabel("Theta/Delta", color="white")
        ax3b.tick_params(axis="y", colors="white")

        # Wavelet scalogram over whole window (1–30 Hz)
        f_min, f_max, df = 1.0, 30.0, 0.5
        freqs_wv, scalogram_db = compute_wavelet_scalogram_matlab_like(
            eeg_win, fs, f_min=f_min, f_max=f_max, df=df,
            f0=1.0, exp_corr=0.0, time_smooth=0.0, baseline_mode="window",
        )
        n_freqs, n_times = scalogram_db.shape

        t_edges = np.linspace(win_start_sec, win_end_sec_actual, n_times + 1)
        if len(freqs_wv) > 1:
            df_eff = freqs_wv[1] - freqs_wv[0]
        else:
            df_eff = 1.0
        f_edges = np.concatenate(
            [freqs_wv - df_eff / 2, [freqs_wv[-1] + df_eff / 2]]
        )

        vlim = 6.0

        im2 = ax4.pcolormesh(
            t_edges,
            f_edges,
            scalogram_db,
            shading="auto",
            cmap="coolwarm",
            vmin=-vlim,
            vmax=vlim,
        )
        ax4.set_ylabel("Frequency (Hz)")
        ax4.set_xlabel("Time (s)")
        ax4.set_ylim(0, 30)
        ax4.set_title("Morlet scalogram (dB rel. window mean)")

        cbar2 = self.fig.colorbar(im2, ax=ax4)
        cbar2.set_label("Power (dB rel. mean)")

        # ACh trace (labeled as ACh even if variable is 'ne')
        if self.ach is not None and len(self.ach) > 1:
            ach = np.asarray(self.ach, dtype=float).squeeze()
            fs_ach = float(self.fs_ach)
            t_ach_full = np.arange(len(ach)) / fs_ach

            if (win_start_sec >= t_ach_full[-1]) or (win_end_sec_actual <= t_ach_full[0]):
                ax5.text(0.5, 0.5, "ACh trace does not cover this time window.",
                        ha="center", va="center", transform=ax5.transAxes)
                ax5.set_axis_off()
            else:
                ach_win = np.interp(t_win, t_ach_full, ach)
                ax5.plot(t_win, ach_win, linewidth=0.6, color="black")
                ax5.set_ylabel("ACh (a.u.)")
                ax5.set_xlabel("Time (s)")
                ax5.set_xlim(win_start_sec, win_end_sec_actual)
                ax5.set_ylim(-7, 7)
                ax5.set_title("ACh trace, aligned to EEG window")

        # Build title with selection indicator
        is_selected = self.idx in self.selected_bouts
        marker = "★ SELECTED" if is_selected else "☆ not selected"
        
        title_text = (
            f"{self.base} | {self.state_label} bout {self.idx+1}/{n_bouts}  "
            f"{start_s:.1f}–{end_s:.1f} s  "
            f"(window {win_start_sec:.1f}–{win_end_sec_actual:.1f} s)\n"
            f"[{marker}] {len(self.selected_bouts)} bouts selected. "
            "s=select, u=unselect, c=compute PLV/PAC, x=export list\n"
            "Use ←/→ or n/p to move, q/Esc to quit. "
            "Drag on EEG for detail wavelet. Bands: "
            "1=0–5, 2=0–15, 3=0–30, 4=theta, 5=spindle, 6=beta, 7=gamma (l/h/a→1/5/3)."
        )
        
        self.fig.suptitle(title_text, fontsize=11)

        self.fig.tight_layout(rect=[0, 0.08, 1, 0.93])
        self.fig.canvas.draw_idle()

    def export_selected_bouts(self):
        """Export list of selected bout indices to CSV."""
        if not self.selected_bouts:
            print("⚠️  No bouts selected to export")
            return
        
        records = []
        for bout_idx in sorted(self.selected_bouts):
            start_idx, end_idx = self.bouts_idx[bout_idx]
            start_s = float(self.time_s[start_idx])
            end_s = float(self.time_s[end_idx - 1] + 1.0)
            
            records.append({
                "bout_index": bout_idx,
                "start_s": start_s,
                "end_s": end_s,
                "duration_s": end_s - start_s,
            })
        
        df = pd.DataFrame(records)
        out_csv = os.path.join(RESULTS_DIR, f"{self.base}_selected_bouts_list.csv")
        df.to_csv(out_csv, index=False)
        print(f"✓ Exported {len(records)} selected bout indices to: {out_csv}")
    
    def compute_plv_pac_analysis(self):
        """Compute PLV and PAC for all selected bouts."""
        print(f"\n{'='*70}")
        print(f"Computing PLV/PAC for {len(self.selected_bouts)} selected {self.state_label} bouts...")
        print(f"{'='*70}")
        
        results = []
        
        for bout_idx in sorted(self.selected_bouts):
            start_idx, end_idx = self.bouts_idx[bout_idx]
            start_s = float(self.time_s[start_idx])
            end_s = float(self.time_s[end_idx - 1] + 1.0)
            
            # Extract EEG bout
            i_start = int(start_s * self.fs)
            i_stop = int(end_s * self.fs)
            eeg_bout = self.eeg[i_start:i_stop]
            
            # Extract ACh bout (resample if needed)
            t_eeg = np.arange(len(eeg_bout)) / self.fs + start_s
            t_ach_full = np.arange(len(self.ach)) / self.fs_ach
            ach_bout = np.interp(t_eeg, t_ach_full, self.ach)
            
            if len(eeg_bout) < 100 or len(ach_bout) < 100:
                print(f"  ⚠️  Bout {bout_idx + 1} too short ({len(eeg_bout)} samples), skipping")
                continue
            
            try:
                # 1. PLV: EEG theta vs ACh slow fluctuations
                plv_val, _ = compute_plv(
                    eeg_bout, ach_bout, self.fs,
                    band1=PLV_EEG_BAND,
                    band2=PLV_ACH_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PLV failed for bout {bout_idx + 1}: {e}")
                plv_val = np.nan
            
            try:
                # 2. PAC within EEG: theta phase -> gamma amplitude
                pac_eeg, _, _ = compute_pac_mi(
                    eeg_bout, self.fs,
                    phase_band=PAC_PHASE_BAND,
                    amp_band=PAC_AMP_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PAC (EEG) failed for bout {bout_idx + 1}: {e}")
                pac_eeg = np.nan
            
            try:
                # 3. Cross-signal PAC: EEG theta phase -> ACh amplitude
                pac_cross, _, _ = compute_pac_cross(
                    eeg_bout, ach_bout, self.fs,
                    phase_band_eeg=PAC_PHASE_BAND,
                    amp_band_ach=PLV_ACH_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PAC (cross) failed for bout {bout_idx + 1}: {e}")
                pac_cross = np.nan
            
            results.append({
                "file": self.base,
                "state": self.state_label,
                "bout_index": bout_idx,
                "start_s": start_s,
                "end_s": end_s,
                "duration_s": end_s - start_s,
                "plv_eeg_ach": plv_val,
                "pac_eeg_theta_gamma": pac_eeg,
                "pac_eeg_phase_ach_amp": pac_cross,
            })
            
            print(f"  Bout {bout_idx + 1}: PLV={plv_val:.3f}, "
                  f"PAC_EEG={pac_eeg:.3f}, PAC_cross={pac_cross:.3f}")
        
        # Save results
        df = pd.DataFrame(results)
        out_csv = os.path.join(RESULTS_DIR, f"{self.base}_plv_pac_results.csv")
        df.to_csv(out_csv, index=False)
        print(f"\n✓ Saved {len(results)} bout analyses to: {out_csv}\n")
        
        # Plot summary
        self.plot_analysis_summary(df)
    
    def plot_analysis_summary(self, df):
        """Plot PLV and PAC distributions."""
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        
        # PLV
        plv_vals = df["plv_eeg_ach"].dropna()
        if len(plv_vals) > 0:
            axes[0].hist(plv_vals, bins=12, edgecolor="black", alpha=0.7, color="steelblue")
            axes[0].axvline(plv_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={plv_vals.mean():.3f}")
            axes[0].set_xlabel("PLV (EEG theta - ACh)")
            axes[0].set_ylabel("Count")
            axes[0].set_title(f"Phase-Locking Value\n"
                             f"EEG {PLV_EEG_BAND[0]}-{PLV_EEG_BAND[1]} Hz vs "
                             f"ACh {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
            axes[0].legend()
        
        # PAC within EEG
        pac_vals = df["pac_eeg_theta_gamma"].dropna()
        if len(pac_vals) > 0:
            axes[1].hist(pac_vals, bins=12, edgecolor="black", alpha=0.7, color="orange")
            axes[1].axvline(pac_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={pac_vals.mean():.3f}")
            axes[1].set_xlabel("Modulation Index")
            axes[1].set_ylabel("Count")
            axes[1].set_title(f"PAC within EEG\n"
                             f"Phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz, "
                             f"Amp {PAC_AMP_BAND[0]}-{PAC_AMP_BAND[1]} Hz")
            axes[1].legend()
        
        # Cross-signal PAC
        cross_vals = df["pac_eeg_phase_ach_amp"].dropna()
        if len(cross_vals) > 0:
            axes[2].hist(cross_vals, bins=12, edgecolor="black", alpha=0.7, color="green")
            axes[2].axvline(cross_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={cross_vals.mean():.3f}")
            axes[2].set_xlabel("Modulation Index")
            axes[2].set_ylabel("Count")
            axes[2].set_title(f"Cross-signal PAC\n"
                             f"EEG phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz → "
                             f"ACh amp {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
            axes[2].legend()
        
        fig.suptitle(f"{self.base} - {len(df)} selected {self.state_label} bouts",
                    fontsize=13, fontweight="bold")
        fig.tight_layout()
        plt.show()


# Keep the original BoutViewer class definition for backwards compatibility
class BoutViewer:
    def __init__(self, base, eeg, emg, fs, time_s, states, bouts_idx, state_label, ach=None, fs_ach=None):
        self.base = base
        self.eeg = eeg
        self.emg = emg
        self.fs = fs
        self.time_s = time_s
        self.states = states
        self.bouts_idx = bouts_idx
        self.state_label = state_label
        self.idx = 0
        self.stop_all = False
        self.ach = ach
        self.fs_ach = fs_ach if fs_ach is not None else fs
        self.detail_band_key = DEFAULT_BAND_KEY
        self.span_selector = None

        self.fig = plt.figure(figsize=(12, 11))
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.update_plot()

    def on_key(self, event):
        if event.key in ["right", "down", "n"]:
            if self.idx < len(self.bouts_idx) - 1:
                self.idx += 1
                self.update_plot()
        elif event.key in ["left", "up", "p"]:
            if self.idx > 0:
                self.idx -= 1
                self.update_plot()
        elif event.key in ["q", "escape"]:
            self.stop_all = True
            plt.close(self.fig)
        elif event.key in ("l", "h", "a") or event.key in BAND_PRESETS:
            if event.key == "l":
                key = "1"
            elif event.key == "h":
                key = "5"
            elif event.key == "a":
                key = "3"
            else:
                key = event.key
            if key in BAND_PRESETS:
                self.detail_band_key = key
                band_label, f_min, f_max = BAND_PRESETS[key]
                print(f"Detail band set to {band_label} [{f_min}-{f_max} Hz]")

    def on_select_span(self, tmin, tmax):
        if tmax < tmin:
            tmin, tmax = tmax, tmin
        if tmax - tmin < 0.3:
            return
        self.open_detail_wavelet(tmin, tmax)

    def open_detail_wavelet(self, t0, t1):
        fs = self.fs
        n_samples = len(self.eeg)

        i0 = max(0, int(t0 * fs))
        i1 = min(n_samples, int(t1 * fs))
        if i1 <= i0 + 10:
            return

        seg = self.eeg[i0:i1]

        band_label, f_min, f_max = BAND_PRESETS.get(
            self.detail_band_key, BAND_PRESETS[DEFAULT_BAND_KEY]
        )

        freqs, power_norm, power_lin = compute_wavelet_scalogram_pretty(
            seg, fs, f_min=f_min, f_max=f_max,
            df=0.25, time_smooth=0.15, return_linear=True
        )

        n_freqs, n_times = power_norm.shape
        t_edges = np.linspace(t0, t1, n_times + 1)

        if len(freqs) > 1:
            df_eff = freqs[1] - freqs[0]
        else:
            df_eff = 1.0
        f_edges = np.concatenate([freqs - df_eff / 2, [freqs[-1] + df_eff / 2]])

        fig, (ax_wv, ax_prof, ax_hist) = plt.subplots(
            3, 1, figsize=(7, 6),
            gridspec_kw={"height_ratios": [3, 1, 1]},
        )

        im = ax_wv.pcolormesh(
            t_edges,
            f_edges,
            power_norm,
            shading="auto",
            cmap="turbo",
            vmin=0.0,
            vmax=1.0,
        )
        ax_wv.set_ylabel("Frequency (Hz)")
        ax_wv.set_ylim(f_min, f_max)
        ax_wv.set_xlim(t0, t1)

        cbar = fig.colorbar(im, ax=ax_wv, pad=0.02)
        cbar.set_label("Relative power (0–1) per freq")

        mean_power_lin = power_lin.mean(axis=1)  # average over time
        mean_power_db = 10 * np.log10(mean_power_lin + 1e-20)

        ax_prof.plot(freqs, mean_power_db, "-")
        ax_prof.set_xlabel("Frequency (Hz)")
        ax_prof.set_ylabel("Mean power (dB)")
        ax_prof.set_xlim(f_min, f_max)

        dom_idx = np.argmax(mean_power_db)
        dom_freq = freqs[dom_idx]
        ax_prof.axvline(dom_freq, color="red", linestyle="--", alpha=0.7,
                        label=f"Peak ≈ {dom_freq:.2f} Hz")
        ax_prof.legend(fontsize=8, loc="best")

        peak_idx = np.argmax(power_lin, axis=0)
        peak_freqs = freqs[peak_idx]

        n_bins = 12
        bins = np.linspace(f_min, f_max, n_bins + 1)
        counts, edges = np.histogram(peak_freqs, bins=bins)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        prop = counts / counts.sum() if counts.sum() > 0 else counts

        ax_hist.bar(
            bin_centers,
            prop,
            width=np.diff(edges),
            align="center",
            edgecolor="black",
        )

        mean_peak = peak_freqs.mean()
        med_peak = np.median(peak_freqs)
        ax_hist.axvline(mean_peak, color="red", linestyle="--",
                        label=f"Mean {mean_peak:.2f} Hz")
        ax_hist.axvline(med_peak, color="black", linestyle=":",
                        label=f"Median {med_peak:.2f} Hz")

        ax_hist.set_xlabel("Peak frequency per time bin (Hz)")
        ax_hist.set_ylabel("Proportion of time")
        ax_hist.set_xlim(f_min, f_max)
        ax_hist.legend(fontsize=8, loc="best")

        fig.suptitle(f"Detail wavelet {band_label}\n{t0:.2f}–{t1:.2f} s", fontsize=11)
        fig.tight_layout()
        fig.show()

    def update_plot(self):
        """Update the plot for the current bout."""
        self.fig.clear()
        n_bouts = len(self.bouts_idx)

        # 5 panels: EEG, EMG, spectrogram, wavelet scalogram, PSD
        ax1 = self.fig.add_subplot(5, 1, 1)
        ax2 = self.fig.add_subplot(5, 1, 2, sharex=ax1)
        ax3 = self.fig.add_subplot(5, 1, 3, sharex=ax1)
        ax4 = self.fig.add_subplot(5, 1, 4, sharex=ax1)
        ax5 = self.fig.add_subplot(5, 1, 5)

        # SpanSelector for interactive detail-wavelet popup
        if self.span_selector is None:
            self.span_selector = SpanSelector(
                ax1,
                self.on_select_span,
                "horizontal",
                useblit=True,
                props=dict(alpha=0.3, facecolor="yellow"),
                interactive=True,
            )
        else:
            self.span_selector.ax = ax1

        start_idx, end_idx = self.bouts_idx[self.idx]
        start_s = float(self.time_s[start_idx])
        end_s = float(self.time_s[end_idx - 1] + 1.0)

        fs = self.fs
        eeg = self.eeg
        emg = self.emg

        n_samples = len(eeg)
        total_sec = n_samples / fs

        win_start_sec = max(0.0, start_s - PRE_SEC)
        win_end_sec = min(total_sec, end_s + POST_SEC)
        win_start_sample = int(win_start_sec * fs)
        win_end_sample = int(win_end_sec * fs)

        eeg_win = eeg[win_start_sample:win_end_sample]
        emg_win = emg[win_start_sample:win_end_sample]
        t_win = np.arange(win_start_sample, win_end_sample) / fs

        for ax in [ax1, ax2, ax3, ax4]:
            shade_states(ax, self.time_s, self.states, win_start_sec, win_end_sec)
            ax.set_xlim(win_start_sec, win_end_sec)

        ax1.plot(t_win, eeg_win, linewidth=0.5, color="black")
        ax1.set_ylabel("EEG (a.u.)")

        legend_patches = []
        for s in sorted(np.unique(self.states)):
            name = STATE_NAMES.get(int(s), str(s))
            color = STATE_COLORS.get(int(s))
            if color is not None:
                legend_patches.append(Patch(facecolor=color, edgecolor="none", label=name))
        if legend_patches:
            ax1.legend(handles=legend_patches, loc="upper right", fontsize=8)

        ax2.plot(t_win, emg_win, linewidth=0.5, color="black")
        ax2.set_ylabel("EMG (a.u.)")

        # Spectrogram + theta/delta ratio
        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )
        im = ax3.pcolormesh(t_spec, freqs_spec, Sx_db, shading="auto", cmap="viridis")
        ax3.set_ylabel("Freq (Hz)")
        ax3.set_ylim(0, 30)
        win_end_sec_actual = win_start_sec + len(eeg_win) / fs
        ax3.set_xlim(win_start_sec, win_end_sec_actual)

        cbar = self.fig.colorbar(im, ax=ax3)
        cbar.set_label("Power (dB)")

        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.4, linewidth=1.0)
        ax3b.set_ylabel("Theta/Delta", color="white")
        ax3b.tick_params(axis="y", colors="white")

        # Wavelet scalogram over whole window (1–30 Hz)
        f_min, f_max, df = 1.0, 30.0, 0.5
        freqs_wv, scalogram_db = compute_wavelet_scalogram_matlab_like(
            eeg_win, fs, f_min=f_min, f_max=f_max, df=df,
            f0=1.0, exp_corr=0.0, time_smooth=0.0, baseline_mode="window",
        )
        n_freqs, n_times = scalogram_db.shape

        t_edges = np.linspace(win_start_sec, win_end_sec_actual, n_times + 1)
        if len(freqs_wv) > 1:
            df_eff = freqs_wv[1] - freqs_wv[0]
        else:
            df_eff = 1.0
        f_edges = np.concatenate(
            [freqs_wv - df_eff / 2, [freqs_wv[-1] + df_eff / 2]]
        )

        vlim = 6.0

        im2 = ax4.pcolormesh(
            t_edges,
            f_edges,
            scalogram_db,
            shading="auto",
            cmap="coolwarm",
            vmin=-vlim,
            vmax=vlim,
        )
        ax4.set_ylabel("Frequency (Hz)")
        ax4.set_xlabel("Time (s)")
        ax4.set_ylim(0, 30)
        ax4.set_title("Morlet scalogram (dB rel. window mean)")

        cbar2 = self.fig.colorbar(im2, ax=ax4)
        cbar2.set_label("Power (dB rel. mean)")

        # ACh trace (labeled as ACh even if variable is 'ne')
        if self.ach is not None and len(self.ach) > 1:
            ach = np.asarray(self.ach, dtype=float).squeeze()
            fs_ach = float(self.fs_ach)
            t_ach_full = np.arange(len(ach)) / fs_ach

            if (win_start_sec >= t_ach_full[-1]) or (win_end_sec_actual <= t_ach_full[0]):
                ax5.text(0.5, 0.5, "ACh trace does not cover this time window.",
                        ha="center", va="center", transform=ax5.transAxes)
                ax5.set_axis_off()
            else:
                ach_win = np.interp(t_win, t_ach_full, ach)
                ax5.plot(t_win, ach_win, linewidth=0.6, color="black")
                ax5.set_ylabel("ACh (a.u.)")
                ax5.set_xlabel("Time (s)")
                ax5.set_xlim(win_start_sec, win_end_sec_actual)
                ax5.set_ylim(-7, 7)
                ax5.set_title("ACh trace, aligned to EEG window")

        # Build title with selection indicator
        is_selected = self.idx in self.selected_bouts
        marker = "★ SELECTED" if is_selected else "☆ not selected"
        
        title_text = (
            f"{self.base} | {self.state_label} bout {self.idx+1}/{n_bouts}  "
            f"{start_s:.1f}–{end_s:.1f} s  "
            f"(window {win_start_sec:.1f}–{win_end_sec_actual:.1f} s)\n"
            f"[{marker}] {len(self.selected_bouts)} bouts selected. "
            "s=select, u=unselect, c=compute PLV/PAC, x=export list\n"
            "Use ←/→ or n/p to move, q/Esc to quit. "
            "Drag on EEG for detail wavelet. Bands: "
            "1=0–5, 2=0–15, 3=0–30, 4=theta, 5=spindle, 6=beta, 7=gamma (l/h/a→1/5/3)."
        )
        
        self.fig.suptitle(title_text, fontsize=11)

        self.fig.tight_layout(rect=[0, 0.08, 1, 0.93])
        self.fig.canvas.draw_idle()

    def export_selected_bouts(self):
        """Export list of selected bout indices to CSV."""
        if not self.selected_bouts:
            print("⚠️  No bouts selected to export")
            return
        
        records = []
        for bout_idx in sorted(self.selected_bouts):
            start_idx, end_idx = self.bouts_idx[bout_idx]
            start_s = float(self.time_s[start_idx])
            end_s = float(self.time_s[end_idx - 1] + 1.0)
            
            records.append({
                "bout_index": bout_idx,
                "start_s": start_s,
                "end_s": end_s,
                "duration_s": end_s - start_s,
            })
        
        df = pd.DataFrame(records)
        out_csv = os.path.join(RESULTS_DIR, f"{self.base}_selected_bouts_list.csv")
        df.to_csv(out_csv, index=False)
        print(f"✓ Exported {len(records)} selected bout indices to: {out_csv}")
    
    def compute_plv_pac_analysis(self):
        """Compute PLV and PAC for all selected bouts."""
        print(f"\n{'='*70}")
        print(f"Computing PLV/PAC for {len(self.selected_bouts)} selected {self.state_label} bouts...")
        print(f"{'='*70}")
        
        results = []
        
        for bout_idx in sorted(self.selected_bouts):
            start_idx, end_idx = self.bouts_idx[bout_idx]
            start_s = float(self.time_s[start_idx])
            end_s = float(self.time_s[end_idx - 1] + 1.0)
            
            # Extract EEG bout
            i_start = int(start_s * self.fs)
            i_stop = int(end_s * self.fs)
            eeg_bout = self.eeg[i_start:i_stop]
            
            # Extract ACh bout (resample if needed)
            t_eeg = np.arange(len(eeg_bout)) / self.fs + start_s
            t_ach_full = np.arange(len(self.ach)) / self.fs_ach
            ach_bout = np.interp(t_eeg, t_ach_full, self.ach)
            
            if len(eeg_bout) < 100 or len(ach_bout) < 100:
                print(f"  ⚠️  Bout {bout_idx + 1} too short ({len(eeg_bout)} samples), skipping")
                continue
            
            try:
                # 1. PLV: EEG theta vs ACh slow fluctuations
                plv_val, _ = compute_plv(
                    eeg_bout, ach_bout, self.fs,
                    band1=PLV_EEG_BAND,
                    band2=PLV_ACH_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PLV failed for bout {bout_idx + 1}: {e}")
                plv_val = np.nan
            
            try:
                # 2. PAC within EEG: theta phase -> gamma amplitude
                pac_eeg, _, _ = compute_pac_mi(
                    eeg_bout, self.fs,
                    phase_band=PAC_PHASE_BAND,
                    amp_band=PAC_AMP_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PAC (EEG) failed for bout {bout_idx + 1}: {e}")
                pac_eeg = np.nan
            
            try:
                # 3. Cross-signal PAC: EEG theta phase -> ACh amplitude
                pac_cross, _, _ = compute_pac_cross(
                    eeg_bout, ach_bout, self.fs,
                    phase_band_eeg=PAC_PHASE_BAND,
                    amp_band_ach=PLV_ACH_BAND
                )
            except Exception as e:
                print(f"  ⚠️  PAC (cross) failed for bout {bout_idx + 1}: {e}")
                pac_cross = np.nan
            
            results.append({
                "file": self.base,
                "state": self.state_label,
                "bout_index": bout_idx,
                "start_s": start_s,
                "end_s": end_s,
                "duration_s": end_s - start_s,
                "plv_eeg_ach": plv_val,
                "pac_eeg_theta_gamma": pac_eeg,
                "pac_eeg_phase_ach_amp": pac_cross,
            })
            
            print(f"  Bout {bout_idx + 1}: PLV={plv_val:.3f}, "
                  f"PAC_EEG={pac_eeg:.3f}, PAC_cross={pac_cross:.3f}")
        
        # Save results
        df = pd.DataFrame(results)
        out_csv = os.path.join(RESULTS_DIR, f"{self.base}_plv_pac_results.csv")
        df.to_csv(out_csv, index=False)
        print(f"\n✓ Saved {len(results)} bout analyses to: {out_csv}\n")
        
        # Plot summary
        self.plot_analysis_summary(df)
    
    def plot_analysis_summary(self, df):
        """Plot PLV and PAC distributions."""
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        
        # PLV
        plv_vals = df["plv_eeg_ach"].dropna()
        if len(plv_vals) > 0:
            axes[0].hist(plv_vals, bins=12, edgecolor="black", alpha=0.7, color="steelblue")
            axes[0].axvline(plv_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={plv_vals.mean():.3f}")
            axes[0].set_xlabel("PLV (EEG theta - ACh)")
            axes[0].set_ylabel("Count")
            axes[0].set_title(f"Phase-Locking Value\n"
                             f"EEG {PLV_EEG_BAND[0]}-{PLV_EEG_BAND[1]} Hz vs "
                             f"ACh {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
            axes[0].legend()
        
        # PAC within EEG
        pac_vals = df["pac_eeg_theta_gamma"].dropna()
        if len(pac_vals) > 0:
            axes[1].hist(pac_vals, bins=12, edgecolor="black", alpha=0.7, color="orange")
            axes[1].axvline(pac_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={pac_vals.mean():.3f}")
            axes[1].set_xlabel("Modulation Index")
            axes[1].set_ylabel("Count")
            axes[1].set_title(f"PAC within EEG\n"
                             f"Phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz, "
                             f"Amp {PAC_AMP_BAND[0]}-{PAC_AMP_BAND[1]} Hz")
            axes[1].legend()
        
        # Cross-signal PAC
        cross_vals = df["pac_eeg_phase_ach_amp"].dropna()
        if len(cross_vals) > 0:
            axes[2].hist(cross_vals, bins=12, edgecolor="black", alpha=0.7, color="green")
            axes[2].axvline(cross_vals.mean(), color="red", linestyle="--",
                           label=f"Mean={cross_vals.mean():.3f}")
            axes[2].set_xlabel("Modulation Index")
            axes[2].set_ylabel("Count")
            axes[2].set_title(f"Cross-signal PAC\n"
                             f"EEG phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz → "
                             f"ACh amp {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
            axes[2].legend()
        
        fig.suptitle(f"{self.base} - {len(df)} selected {self.state_label} bouts",
                    fontsize=13, fontweight="bold")
        fig.tight_layout()
        plt.show()


# ===========================
# Main
# ===========================

def main():
    mat_files = sorted(glob.glob(os.path.join(SIGNAL_DIR, "*.mat")))
    if not mat_files:
        print(f"No .mat files found in {SIGNAL_DIR}")
        return

    print(f"Found {len(mat_files)} signal files.")

    state_map = {0: "Wake", 1: "NREM", 2: "REM"}
    print("\nWhich state do you want to explore?")
    print("0 = Wake, 1 = NREM, 2 = REM")
    choice = input("Enter 0, 1, or 2 (default 2=REM): ").strip()
    
    if choice == "0":
        target_state = 0
    elif choice == "1":
        target_state = 1
    else:
        target_state = 2
    
    state_label = state_map[target_state]

    # NEW: Ask for analysis mode
    print("\nAnalysis mode:")
    print("1 = Interactive viewer with manual bout selection (press 's' to select)")
    print("2 = Batch mode: analyze ALL bouts automatically (no viewer)")
    mode_choice = input("Enter 1 or 2 (default 1=Interactive): ").strip()
    batch_mode = (mode_choice == "2")

    print(f"\n{'='*70}")
    print(f"PLV/PAC Analysis Configuration:")
    print(f"  PLV: EEG {PLV_EEG_BAND} Hz ↔ ACh {PLV_ACH_BAND} Hz")
    print(f"  PAC (EEG): Phase {PAC_PHASE_BAND} Hz → Amp {PAC_AMP_BAND} Hz")
    print(f"  Results will be saved to: {RESULTS_DIR}")
    if batch_mode:
        print(f"  MODE: Batch analysis (all bouts, no viewer)")
    else:
        print(f"  MODE: Interactive selection")
    print(f"{'='*70}\n")

    for mat_path in mat_files:
        base = os.path.basename(mat_path).replace(".mat", "")
        csv_path = os.path.join(SCORES_DIR, base + "_scored_scores_1Hz.csv")

        if not os.path.exists(csv_path):
            print(f"⚠️  No scores CSV for {base}, skipping.")
            continue

        print(f"\n=== File: {base} ===")
        
        eeg, emg, fs, ach, fs_ach = load_mat_signal(mat_path)
        time_s, states = load_scores_1hz(csv_path)

        total_sec = len(eeg) / fs
        valid_mask = time_s < total_sec
        time_s = time_s[valid_mask]
        states = states[valid_mask]

        bouts_by_state = find_state_bouts(states, min_len_sec=1)
        bouts_idx = bouts_by_state.get(target_state, [])
        bouts_idx = [(s, e) for (s, e) in bouts_idx if (e - s) >= MIN_BOUT_LEN_SEC]

        if not bouts_idx:
            print(f"  No {state_label} bouts ≥ {MIN_BOUT_LEN_SEC}s, skipping.")
            continue

        print(f"  {len(bouts_idx)} {state_label} bouts found.")
        
        if batch_mode:
            # Batch mode: analyze all bouts without viewer
            if ach is None:
                print("  ⚠️  No ACh data - skipping PLV/PAC analysis")
                continue
            
            print(f"  Running batch PLV/PAC analysis on all {len(bouts_idx)} bouts...")
            batch_analyze_bouts(base, eeg, ach, fs, fs_ach, time_s, bouts_idx, state_label)
        else:
            # Interactive mode
            if ach is None:
                print("  ⚠️  No ACh data - using basic viewer (no PLV/PAC)")
                viewer = BoutViewer(base, eeg, emg, fs, time_s, states, bouts_idx, state_label, ach=ach, fs_ach=fs_ach)
            else:
                print(f"  ✓ ACh data available - using PLV/PAC viewer")
                viewer = BoutViewerWithSelection(base, eeg, emg, fs, time_s, states, bouts_idx, state_label, ach=ach, fs_ach=fs_ach )
            
            plt.show()

            if viewer.stop_all:
                print("\nExiting at user request.")
                return

    print("\nDone.")


def batch_analyze_bouts(base, eeg, ach, fs, fs_ach, time_s, bouts_idx, state_label):
    """
    Analyze all bouts without interactive viewer.
    
    Parameters
    ----------
    base : str
        Recording name
    eeg : array
        EEG signal
    ach : array
        ACh signal
    fs : float
        EEG sampling rate
    fs_ach : float
        ACh sampling rate
    time_s : array
        Time vector for scoring
    bouts_idx : list of tuples
        List of (start_idx, end_idx) for each bout
    state_label : str
        State name (Wake/NREM/REM)
    """
    results = []
    
    for bout_num, (start_idx, end_idx) in enumerate(bouts_idx, 1):
        start_s = float(time_s[start_idx])
        end_s = float(time_s[end_idx - 1] + 1.0)
        
        # Extract EEG bout
        i_start = int(start_s * fs)
        i_stop = int(end_s * fs)
        eeg_bout = eeg[i_start:i_stop]
        
        # Extract ACh bout (resample if needed)
        t_eeg = np.arange(len(eeg_bout)) / fs + start_s
        t_ach_full = np.arange(len(ach)) / fs_ach
        ach_bout = np.interp(t_eeg, t_ach_full, ach)
        
        if len(eeg_bout) < 100 or len(ach_bout) < 100:
            print(f"  ⚠️  Bout {bout_num}/{len(bouts_idx)} too short ({len(eeg_bout)} samples), skipping")
            continue
        
        try:
            # 1. PLV: EEG theta vs ACh slow fluctuations
            plv_val, _ = compute_plv(
                eeg_bout, ach_bout, fs,
                band1=PLV_EEG_BAND,
                band2=PLV_ACH_BAND
            )
        except Exception as e:
            print(f"  ⚠️  Bout {bout_num}/{len(bouts_idx)}: PLV failed - {e}")
            plv_val = np.nan
        
        try:
            # 2. PAC within EEG: theta phase -> gamma amplitude
            pac_eeg, _, _ = compute_pac_mi(
                eeg_bout, fs,
                phase_band=PAC_PHASE_BAND,
                amp_band=PAC_AMP_BAND
            )
        except Exception as e:
            print(f"  ⚠️  Bout {bout_num}/{len(bouts_idx)}: PAC (EEG) failed - {e}")
            pac_eeg = np.nan
        
        try:
            # 3. Cross-signal PAC: EEG theta phase -> ACh amplitude
            pac_cross, _, _ = compute_pac_cross(
                eeg_bout, ach_bout, fs,
                phase_band_eeg=PAC_PHASE_BAND,
                amp_band_ach=PLV_ACH_BAND
            )
        except Exception as e:
            print(f"  ⚠️  Bout {bout_num}/{len(bouts_idx)}: PAC (cross) failed - {e}")
            pac_cross = np.nan
        
        results.append({
            "file": base,
            "state": state_label,
            "bout_index": bout_num - 1,  # 0-indexed
            "start_s": start_s,
            "end_s": end_s,
            "duration_s": end_s - start_s,
            "plv_eeg_ach": plv_val,
            "pac_eeg_theta_gamma": pac_eeg,
            "pac_eeg_phase_ach_amp": pac_cross,
        })
        
        # Progress indicator
        if bout_num % 10 == 0 or bout_num == len(bouts_idx):
            print(f"    Processed {bout_num}/{len(bouts_idx)} bouts...")
    
    # Save results
    df = pd.DataFrame(results)
    out_csv = os.path.join(RESULTS_DIR, f"{base}_{state_label}_all_bouts_plv_pac.csv")
    df.to_csv(out_csv, index=False)
    print(f"  ✓ Saved {len(results)} bout analyses to: {out_csv}")
    
    # Print summary statistics
    print(f"\n  Summary statistics for {base} ({state_label}):")
    print(f"    PLV (EEG-ACh):        mean={df['plv_eeg_ach'].mean():.3f}, std={df['plv_eeg_ach'].std():.3f}")
    print(f"    PAC (EEG θ→γ):        mean={df['pac_eeg_theta_gamma'].mean():.3f}, std={df['pac_eeg_theta_gamma'].std():.3f}")
    print(f"    PAC (EEG θ→ACh amp):  mean={df['pac_eeg_phase_ach_amp'].mean():.3f}, std={df['pac_eeg_phase_ach_amp'].std():.3f}")
    
    # Optional: create summary plot
    plot_batch_summary(df, base, state_label)


def plot_batch_summary(df, base, state_label):
    """Create summary plots for batch analysis."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    
    # Row 1: Histograms
    # PLV
    plv_vals = df["plv_eeg_ach"].dropna()
    if len(plv_vals) > 0:
        axes[0, 0].hist(plv_vals, bins=20, edgecolor="black", alpha=0.7, color="steelblue")
        axes[0, 0].axvline(plv_vals.mean(), color="red", linestyle="--",
                          label=f"Mean={plv_vals.mean():.3f}")
        axes[0, 0].set_xlabel("PLV (EEG theta - ACh)")
        axes[0, 0].set_ylabel("Count")
        axes[0, 0].set_title(f"Phase-Locking Value\n"
                            f"EEG {PLV_EEG_BAND[0]}-{PLV_EEG_BAND[1]} Hz vs "
                            f"ACh {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
        axes[0, 0].legend()
    
    # PAC within EEG
    pac_vals = df["pac_eeg_theta_gamma"].dropna()
    if len(pac_vals) > 0:
        axes[0, 1].hist(pac_vals, bins=20, edgecolor="black", alpha=0.7, color="orange")
        axes[0, 1].axvline(pac_vals.mean(), color="red", linestyle="--",
                          label=f"Mean={pac_vals.mean():.3f}")
        axes[0, 1].set_xlabel("Modulation Index")
        axes[0, 1].set_ylabel("Count")
        axes[0, 1].set_title(f"PAC within EEG\n"
                            f"Phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz, "
                            f"Amp {PAC_AMP_BAND[0]}-{PAC_AMP_BAND[1]} Hz")
        axes[0, 1].legend()
    
    # Cross-signal PAC
    cross_vals = df["pac_eeg_phase_ach_amp"].dropna()
    if len(cross_vals) > 0:
        axes[0, 2].hist(cross_vals, bins=20, edgecolor="black", alpha=0.7, color="green")
        axes[0, 2].axvline(cross_vals.mean(), color="red", linestyle="--",
                          label=f"Mean={cross_vals.mean():.3f}")
        axes[0, 2].set_xlabel("Modulation Index")
        axes[0, 2].set_ylabel("Count")
        axes[0, 2].set_title(f"Cross-signal PAC\n"
                            f"EEG phase {PAC_PHASE_BAND[0]}-{PAC_PHASE_BAND[1]} Hz → "
                            f"ACh amp {PLV_ACH_BAND[0]}-{PLV_ACH_BAND[1]} Hz")
        axes[0, 2].legend()
    
    # Row 2: Time series
    bout_nums = df["bout_index"].values
    
    axes[1, 0].plot(bout_nums, df["plv_eeg_ach"], 'o-', markersize=3, alpha=0.6)
    axes[1, 0].axhline(plv_vals.mean(), color="red", linestyle="--", alpha=0.7)
    axes[1, 0].set_xlabel("Bout number")
    axes[1, 0].set_ylabel("PLV")
    axes[1, 0].set_title("PLV across bouts")
    
    axes[1, 1].plot(bout_nums, df["pac_eeg_theta_gamma"], 'o-', markersize=3, alpha=0.6, color="orange")
    axes[1, 1].axhline(pac_vals.mean(), color="red", linestyle="--", alpha=0.7)
    axes[1, 1].set_xlabel("Bout number")
    axes[1, 1].set_ylabel("PAC (EEG)")
    axes[1, 1].set_title("PAC (EEG) across bouts")
    
    axes[1, 2].plot(bout_nums, df["pac_eeg_phase_ach_amp"], 'o-', markersize=3, alpha=0.6, color="green")
    axes[1, 2].axhline(cross_vals.mean(), color="red", linestyle="--", alpha=0.7)
    axes[1, 2].set_xlabel("Bout number")
    axes[1, 2].set_ylabel("PAC (cross)")
    axes[1, 2].set_title("PAC (cross) across bouts")
    
    fig.suptitle(f"{base} - {len(df)} {state_label} bouts (batch analysis)",
                fontsize=14, fontweight="bold")
    fig.tight_layout()
    
    # Save figure
    out_fig = os.path.join(RESULTS_DIR, f"{base}_{state_label}_batch_summary.png")
    fig.savefig(out_fig, dpi=150, bbox_inches="tight")
    print(f"  ✓ Saved summary figure to: {out_fig}")
    plt.close(fig)


if __name__ == "__main__":
    main()
