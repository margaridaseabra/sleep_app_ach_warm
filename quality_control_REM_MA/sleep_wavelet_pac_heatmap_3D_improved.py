#!/usr/bin/env python3
"""
IMPROVED Interactive bout viewer (REM or NREM) using per-second scoring (1 Hz).

KEY IMPROVEMENTS:
1. Better plot alignment using GridSpec
2. Enhanced ACh visualization with change detection
3. More interpretable PAC with z-scoring
4. Clearer labels and consistent styling
5. Added smoothing options for ACh trace

Navigation (when the Matplotlib window has focus):
    → / ↓ / n : next bout
    ← / ↑ / p : previous bout
    q / Esc   : quit viewer & script
    d         : open 3D PAC explorer

Detail-wavelet controls:
    1 : 0–5 Hz        l : shortcut → 1 (0–5 Hz)
    2 : 0–15 Hz       h : shortcut → 5 (spindle 9–16 Hz)
    3 : 0–30 Hz       a : shortcut → 3 (0–30 Hz)
    4 : theta 6–10 Hz
    5 : spindle 9–16 Hz
    6 : beta 15–30 Hz
    7 : gamma 30–80 Hz

    Drag on EEG (panel 1) to select a time window → popup wavelet figure.
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
from scipy.signal import welch, ShortTimeFFT
from scipy.signal.windows import hamming
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import pywt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.widgets import SpanSelector, Slider

# ===========================
# Configuration
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11_cropped_ambtemp"
SCORES_DIR = "/Users/margaridaseabra/24.11_cropped_ambtemp"

# context around bout (seconds)
PRE_SEC = 20.0
POST_SEC = 30.0

# minimum bout length (in seconds) to consider
MIN_BOUT_LEN_SEC = 5.0

# bands (Hz) for ratio 
DELTA_BAND = (1.0, 4.0)
THETA_BAND = (4.0, 8.0)

# High-frequency band 
HF_BAND = (40.0, 80.0)

# Within-bout burst detection parameters
HF_PERCENTILE = 95.0
HF_MIN_BURST_DUR = 0.050

# CWT settings
CWT_MIN_FREQ = 1
CWT_MAX_FREQ = 12
CWT_N_FREQS = 40
CWT_WAVELET = "cmor0.5-1.0"

STATE_COLORS = {
    0: (0.90, 0.80, 0.95, 0.30),  # Wake  -> purple
    1: (1.00, 0.60, 0.60, 0.30),  # NREM  -> red
    2: (0.60, 1.00, 0.60, 0.30),  # REM   -> green
    15: (0.70, 0.70, 1.00, 0.30), # MA    -> blue
}
STATE_NAMES = {0: "Wake", 1: "NREM", 2: "REM", 15: "MA"}

# Preset frequency bands for the detail wavelet popup
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
# PAC configuration (IMPROVED)
# ===========================
PAC_ENABLED = True
PAC_MODE = "heatmap"  # "line" or "heatmap"

# 3D Heatmap PAC parameters
PAC3D_PHASE_CENTERS = np.arange(4.0, 12.1, 1.0)
PAC3D_AMP_CENTERS   = np.arange(30.0, 80.1, 5.0)
PAC3D_PHASE_BW = 2.0
PAC3D_AMP_BW   = 10.0
PAC3D_WIN_SEC  = 4.0
PAC3D_STEP_SEC = 1.0
PAC3D_DOWNSAMPLE_FS = 250.0

# Heatmap amplitude bands (overlapping for better resolution)
PAC_AMP_BANDS = [(30,40),(35,45),(40,50),(45,55),(50,60),(55,65),(60,70),(65,75),(70,80)]

# PAC downsampling
PAC_DOWNSAMPLE_FS = 250.0

# Default: theta-phase -> gamma-amplitude
PAC_PHASE_BAND = (6.0, 10.0)
PAC_AMP_BAND   = (30.0, 80.0)

# Sliding-window PAC(t)
PAC_WIN_SEC  = 4.0
PAC_STEP_SEC = 1.0

# Surrogate testing (set to 0 for speed)
PAC_SURROGATES_N = 0
PAC_SURR_BLOCK_SEC = 0.5
PAC_SURR_PCT = 95.0
PAC_SURR_SEED = 0

# ===========================
# ACh visualization settings (NEW)
# ===========================
ACH_SMOOTH_WINDOW = 10.0  # seconds for smoothing (increase for more smoothing)
ACH_SHOW_DERIVATIVE = False  # show rate of change (set True if desired)
ACH_SHOW_CHANGE_HIGHLIGHTING = False  # show green/red shading (set True if desired)
ACH_CHANGE_THRESHOLD = 0.5  # standard deviations for marking significant changes (if highlighting enabled)

# ===========================
# Utility functions
# ===========================

def load_mat_signal(mat_path):
    """Load EEG, EMG, sampling rate, and (optionally) ACh + its own sampling rate."""
    data = sio.loadmat(mat_path, squeeze_me=True)

    eeg = np.asarray(data["eeg"], dtype=float)
    emg = np.asarray(data["emg"], dtype=float)
    fs = float(np.squeeze(data["eeg_frequency"]))

    # ---- ACh / NE trace ----
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
            print(
                f"Found ACh trace '{ach_key_found}' but no specific fs for it.\n"
                f"   Using EEG fs ({fs} Hz) as a fallback."
            )
            fs_ach = fs
        else:
            print(
                f"Found ACh trace '{ach_key_found}' with its own fs = {fs_ach} Hz"
            )
    else:
        print(
            f"No ACh/NE trace found in {mat_path} "
            f"(tried keys: {', '.join(ach_keys_tried)})"
        )

    return eeg, emg, fs, ach, fs_ach


def load_scores_1hz(csv_path):
    """Load per-second scoring from CSV."""
    df = pd.read_csv(csv_path)
    num_cols = list(df.select_dtypes(include="number").columns)
    if not num_cols:
        raise ValueError(f"No numeric columns found in {csv_path}")

    score_col = None
    for c in num_cols:
        lc = c.lower()
        if "score" in lc or "state" in lc:
            score_col = c
            break
    if score_col is None:
        score_col = min(num_cols, key=lambda c: df[c].nunique())

    states = df[score_col].to_numpy().astype(int)

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


def find_state_bouts(states, min_len_sec=1):
    """Find contiguous bouts for each state in 'states' (1 Hz)."""
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


# ========== Wavelet helpers ==========

def generate_wavelet_fourier(len_wavelet, f_start, f_stop, deltafreq,
                             sample_rate, f0, normalisation):
    """Compute the wavelet coefficients at all scales and compute its Fourier transform."""
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
        len_wavelet,
        min_freq,
        max_freq,
        freq_resolution,
        sub_sample_rate,
        f0,
        exp_corr,
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

    if left_pad > 0:
        wt = wt[:-left_pad]

    if time_smooth > 0:
        sigma = time_smooth * sub_sample_rate
        wt = gaussian_filter(wt, sigma=[sigma, 0], mode="reflect")

    return wt.T


def compute_wavelet_scalogram_matlab_like(
    sig,
    fs,
    f_min=1.0,
    f_max=12.0,
    df=0.5,
    f0=1.0,
    exp_corr=0.0,
    time_smooth=0.0,
    baseline_mode="window",
):
    """Wavelet scalogram similar to MATLAB code."""
    duration_sec = len(sig) / fs

    scalogram_db = compute_scalogram_ephyviewer(
        sig,
        min_freq=f_min,
        max_freq=f_max,
        freq_resolution=df,
        fs=fs,
        f0=f0,
        exp_corr=exp_corr,
        time_smooth=time_smooth,
        wanted_size=duration_sec,
    )

    freqs = np.arange(f_min, f_max, df)

    if baseline_mode == "window":
        baseline_db = scalogram_db.mean(axis=1, keepdims=True)
        scalogram_db = scalogram_db - baseline_db

    return freqs, scalogram_db


def compute_wavelet_scalogram_pretty(
    sig,
    fs,
    f_min,
    f_max,
    df=0.25,
    time_smooth=0.15,
    return_linear=False,
):
    """Convenience wrapper for "pretty" panels."""
    duration_sec = len(sig) / fs

    scalogram_db = compute_scalogram_ephyviewer(
        sig,
        min_freq=f_min,
        max_freq=f_max,
        freq_resolution=df,
        fs=fs,
        f0=1.0,
        exp_corr=0.0,
        time_smooth=time_smooth,
        wanted_size=duration_sec,
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


# ========== Spectrogram / ratio ==========

def compute_spectrogram_and_ratio(eeg_win, fs, win_start_sec, window_duration=5.0, mfft=None):
    """Short-time FFT spectrogram and theta/delta-like ratio."""
    nperseg = round(fs * window_duration)
    hop = round(nperseg / 2)
    window = hamming(nperseg)

    SFT = ShortTimeFFT(
        window,
        hop=hop,
        fs=fs,
        fft_mode="onesided",
        mfft=mfft,
        scale_to="psd",
    )

    Sx = SFT.spectrogram(eeg_win)
    n_frames = Sx.shape[1]

    win_duration = len(eeg_win) / fs

    time = np.linspace(
        win_start_sec,
        win_start_sec + win_duration,
        n_frames,
        endpoint=False,
    )

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


def shade_states(ax, time_s, states, win_start_sec, win_end_sec):
    """Shade background according to per-second states."""
    start_idx = int(max(0, np.searchsorted(time_s, win_start_sec, side="left")))
    end_idx = int(min(len(time_s), np.searchsorted(time_s, win_end_sec, side="right")))

    for i in range(start_idx, end_idx):
        s = int(states[i])
        color = STATE_COLORS.get(s)
        if color is None:
            continue

        t0 = time_s[i]
        if i + 1 < len(time_s):
            t1 = time_s[i + 1]
        else:
            t1 = t0 + 1.0
        left = max(t0, win_start_sec)
        right = min(t1, win_end_sec)
        if right > left:
            ax.axvspan(left, right, color=color, zorder=0)


# ========== PAC helpers ==========

def _bandpass_sos(sig, fs, f_lo, f_hi, order=4):
    """Zero-phase bandpass using SOS."""
    nyq = 0.5 * fs
    f_lo = max(0.001, float(f_lo))
    f_hi = min(float(f_hi), nyq * 0.999)
    if f_hi <= f_lo:
        raise ValueError(f"Invalid bandpass: [{f_lo}, {f_hi}] Hz with nyq={nyq}")

    sos = scipy.signal.butter(order, [f_lo / nyq, f_hi / nyq],
                              btype="bandpass", output="sos")
    return scipy.signal.sosfiltfilt(sos, sig)


def pac_inputs_hilbert(eeg_win, fs, phase_band, amp_band, filt_order=4):
    """Compute phase and amplitude envelope using Hilbert transform."""
    x_phase = _bandpass_sos(eeg_win, fs, phase_band[0], phase_band[1], order=filt_order)
    x_amp   = _bandpass_sos(eeg_win, fs, amp_band[0], amp_band[1], order=filt_order)

    phase = np.angle(scipy.signal.hilbert(x_phase))
    amp   = np.abs(scipy.signal.hilbert(x_amp))
    return phase, amp


def pac_mvl_from_phase_amp(phase, amp, normalize=True):
    """Mean Vector Length PAC."""
    z = amp * np.exp(1j * phase)
    mvl = float(np.abs(np.mean(z)))
    if normalize:
        mvl /= (np.mean(amp) + 1e-12)
    return mvl


def _resample_for_pac(x, fs, target_fs):
    """Resample x from fs -> target_fs using rational resample_poly."""
    if target_fs is None or fs <= target_fs:
        return x, fs
    
    if target_fs <= 2.2 * PAC_AMP_BANDS[-1][1]:
        return x, fs

    fs = float(fs)
    target_fs = float(target_fs)
    
    ratio = target_fs / fs
    
    for denom in [1, 2, 3, 4, 5, 8, 10, 16, 20, 25, 32, 40, 50, 64, 80, 100, 125, 160, 200]:
        num = round(ratio * denom)
        actual_ratio = num / denom
        if abs(actual_ratio - ratio) < 0.01:
            up, down = num, denom
            break
    else:
        try:
            from fractions import Fraction
            frac = Fraction(target_fs).limit_denominator(200) / Fraction(fs).limit_denominator(200)
            frac = frac.limit_denominator(200)
            up, down = frac.numerator, frac.denominator
        except:
            up = int(round(target_fs))
            down = int(round(fs))
            from math import gcd
            g = gcd(up, down)
            up //= g
            down //= g
    
    up = int(up)
    down = int(down)
    
    if up == down:
        return x, fs
    
    y = scipy.signal.resample_poly(x, up, down)
    actual_fs = fs * up / down
    
    return y, actual_fs


def pac_mvl_heatmap(eeg_win, fs, win_start_sec,
                    phase_band, amp_bands,
                    win_sec, step_sec,
                    filt_order=4,
                    downsample_fs=None):
    """Heatmap PAC: fixed phase band; MVL computed for multiple amplitude bands over time."""
    x, fs2 = _resample_for_pac(eeg_win, fs, downsample_fs)

    phase, _ = pac_inputs_hilbert(x, fs2, phase_band=phase_band, amp_band=amp_bands[0], filt_order=filt_order)

    amp_envs = []
    for (f_lo, f_hi) in amp_bands:
        _, amp = pac_inputs_hilbert(x, fs2, phase_band=phase_band, amp_band=(f_lo, f_hi), filt_order=filt_order)
        amp_envs.append(amp)
    amp_envs = np.asarray(amp_envs)

    n = len(phase)
    win_n  = int(round(win_sec * fs2))
    step_n = int(round(step_sec * fs2))
    if n < win_n or win_n < 3:
        return np.array([]), np.array([]), np.zeros((len(amp_bands), 0), dtype=float)

    i0s = np.arange(0, n - win_n + 1, step_n, dtype=int)
    n_win = len(i0s)
    n_b = len(amp_bands)
    pac_map = np.zeros((n_b, n_win), dtype=float)

    for w, i0 in enumerate(i0s):
        i1 = i0 + win_n
        ph = phase[i0:i1]
        for bi in range(n_b):
            pac_map[bi, w] = pac_mvl_from_phase_amp(ph, amp_envs[bi, i0:i1])

    t_edges = win_start_sec + np.concatenate([i0s / fs2, [(i0s[-1] + win_n) / fs2]])

    y_edges = np.array([amp_bands[0][0]] + [b[1] for b in amp_bands], dtype=float)
    return t_edges, y_edges, pac_map


def _moving_mean(x, win_n):
    """Fast moving mean using cumulative sum."""
    x = np.asarray(x)
    c = np.cumsum(np.concatenate([[0], x]))
    return (c[win_n:] - c[:-win_n]) / win_n


def pac_cube_mvl(
    eeg_win, fs, win_start_sec,
    phase_centers, amp_centers,
    phase_bw, amp_bw,
    win_sec, step_sec,
    downsample_fs=250.0,
    filt_order=4,
):
    """Compute 3D PAC cube (phase × amplitude × time)."""
    x, fs2 = _resample_for_pac(eeg_win, fs, downsample_fs)

    win_n  = int(round(win_sec * fs2))
    step_n = int(round(step_sec * fs2))
    n = len(x)
    if n < win_n or win_n < 3:
        return np.array([]), phase_centers, amp_centers, np.zeros((len(phase_centers), len(amp_centers), 0))

    valid_len = n - win_n + 1
    i0s = np.arange(0, valid_len, step_n, dtype=int)
    t_centers = win_start_sec + (i0s + win_n/2) / fs2

    amp_envs = []
    for fa in amp_centers:
        a_lo = fa - amp_bw/2
        a_hi = fa + amp_bw/2
        _, amp = pac_inputs_hilbert(x, fs2, phase_band=(phase_centers[0]-1, phase_centers[0]+1),
                                    amp_band=(a_lo, a_hi), filt_order=filt_order)
        amp_envs.append(amp)
    amp_envs = np.asarray(amp_envs)

    mean_amp = np.stack([_moving_mean(amp_envs[ai], win_n) for ai in range(len(amp_centers))], axis=0)
    mean_amp = mean_amp[:, i0s]

    cube = np.zeros((len(phase_centers), len(amp_centers), len(i0s)), dtype=float)

    for pi, fp in enumerate(phase_centers):
        p_lo = fp - phase_bw/2
        p_hi = fp + phase_bw/2
        phase, _ = pac_inputs_hilbert(x, fs2, phase_band=(p_lo, p_hi),
                                      amp_band=(amp_centers[0]-amp_bw/2, amp_centers[0]+amp_bw/2),
                                      filt_order=filt_order)
        e = np.exp(1j * phase)

        for ai in range(len(amp_centers)):
            z = amp_envs[ai] * e
            mean_z = _moving_mean(z, win_n)
            mean_z = mean_z[i0s]
            cube[pi, ai, :] = np.abs(mean_z) / (mean_amp[ai, :] + 1e-12)

    return t_centers, phase_centers, amp_centers, cube


# ========== ACh processing helpers (NEW) ==========

def process_ach_trace(ach_raw, fs_ach, smooth_window_sec=5.0):
    """
    Process ACh trace: smooth and compute derivative.
    
    Returns:
        ach_smooth: smoothed trace
        ach_deriv: rate of change (derivative)
        ach_z: z-scored smoothed trace
    """
    # Smooth the trace
    smooth_samples = int(smooth_window_sec * fs_ach)
    if smooth_samples >= 3:
        ach_smooth = gaussian_filter1d(ach_raw, sigma=smooth_samples / 4)
    else:
        ach_smooth = ach_raw.copy()
    
    # Compute derivative (rate of change)
    ach_deriv = np.gradient(ach_smooth) * fs_ach  # units per second
    
    # Z-score the smoothed trace
    ach_mean = np.mean(ach_smooth)
    ach_std = np.std(ach_smooth)
    if ach_std > 0:
        ach_z = (ach_smooth - ach_mean) / ach_std
    else:
        ach_z = np.zeros_like(ach_smooth)
    
    return ach_smooth, ach_deriv, ach_z


def detect_ach_changes(ach_z, threshold=0.5):
    """
    Detect significant increases and decreases in ACh trace.
    
    Returns:
        increases: boolean mask of significant increases
        decreases: boolean mask of significant decreases
    """
    increases = ach_z > threshold
    decreases = ach_z < -threshold
    return increases, decreases


# ===========================
# Viewer class (per bout)
# ===========================

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
        self.last_span = None

        # detail-wavelet settings
        self.detail_band_key = DEFAULT_BAND_KEY
        self.span_selector = None

        self.fig = plt.figure(figsize=(14, 12))
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
        elif event.key == "d":
            self.open_pac3d()
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
        """Callback when user drags a horizontal span on the EEG axis."""
        if tmax < tmin:
            tmin, tmax = tmax, tmin
        if tmax - tmin < 0.3:
            return
        self.last_span = (tmin, tmax)
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
            3, 1, figsize=(8, 7),
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
        ax_wv.set_ylabel("Frequency (Hz)", fontsize=10)
        ax_wv.set_ylim(f_min, f_max)
        ax_wv.set_xlim(t0, t1)

        cbar = fig.colorbar(im, ax=ax_wv, pad=0.02)
        cbar.set_label("Relative power (0–1) per freq", fontsize=9)

        mean_power_lin = power_lin.mean(axis=1)
        mean_power_db = 10 * np.log10(mean_power_lin + 1e-20)

        ax_prof.plot(freqs, mean_power_db, "-", linewidth=1.5)
        ax_prof.set_xlabel("Frequency (Hz)", fontsize=10)
        ax_prof.set_ylabel("Mean power (dB)", fontsize=10)
        ax_prof.set_xlim(f_min, f_max)
        ax_prof.grid(True, alpha=0.3)

        dom_idx = np.argmax(mean_power_db)
        dom_freq = freqs[dom_idx]
        ax_prof.axvline(dom_freq, color="red", linestyle="--", alpha=0.7,
                        label=f"Peak ≈ {dom_freq:.2f} Hz")
        ax_prof.legend(fontsize=8, loc="best")

        peak_idx = np.argmax(power_lin, axis=0)
        peak_freqs = freqs[peak_idx]

        n_bins = 15
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
            alpha=0.7,
        )

        mean_peak = peak_freqs.mean()
        med_peak = np.median(peak_freqs)
        ax_hist.axvline(mean_peak, color="red", linestyle="--",
                        label=f"Mean {mean_peak:.2f} Hz")
        ax_hist.axvline(med_peak, color="black", linestyle=":",
                        label=f"Median {med_peak:.2f} Hz")

        ax_hist.set_xlabel("Peak frequency per time bin (Hz)", fontsize=10)
        ax_hist.set_ylabel("Proportion of time", fontsize=10)
        ax_hist.set_xlim(f_min, f_max)
        ax_hist.legend(fontsize=8, loc="best")
        ax_hist.grid(True, alpha=0.3)

        fig.suptitle(
            f"Detail wavelet {band_label}\n{t0:.2f}–{t1:.2f} s",
            fontsize=12,
            fontweight='bold'
        )

        fig.tight_layout()
        fig.show()

    def update_plot(self):
        self.fig.clear()
        n_bouts = len(self.bouts_idx)

        # Use GridSpec for better alignment
        gs = gridspec.GridSpec(7, 1, figure=self.fig, hspace=0.3,
                              height_ratios=[1, 1, 1.2, 1.2, 1.2, 1.2, 0.05])
        
        ax1 = self.fig.add_subplot(gs[0])  # EEG
        ax2 = self.fig.add_subplot(gs[1], sharex=ax1)  # EMG
        ax3 = self.fig.add_subplot(gs[2], sharex=ax1)  # Spectrogram
        ax4 = self.fig.add_subplot(gs[3], sharex=ax1)  # Wavelet
        ax_pac = self.fig.add_subplot(gs[4], sharex=ax1)  # PAC
        ax5 = self.fig.add_subplot(gs[5], sharex=ax1)  # ACh
        
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

        # Shade all axes with state colors
        for ax in [ax1, ax2, ax3, ax4, ax_pac, ax5]:
            shade_states(ax, self.time_s, self.states, win_start_sec, win_end_sec)
            ax.set_xlim(win_start_sec, win_end_sec)

        # Panel 1: EEG
        ax1.plot(t_win, eeg_win, linewidth=0.5, color="black")
        ax1.set_ylabel("EEG (μV)", fontsize=10, fontweight='bold')
        ax1.tick_params(labelbottom=False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # Add legend
        legend_patches = []
        for s in sorted(np.unique(self.states)):
            name = STATE_NAMES.get(int(s), str(s))
            color = STATE_COLORS.get(int(s))
            if color is not None:
                legend_patches.append(Patch(facecolor=color, edgecolor="none", label=name))
        if legend_patches:
            ax1.legend(handles=legend_patches, loc="upper right", fontsize=8, framealpha=0.9)

        # Panel 2: EMG
        ax2.plot(t_win, emg_win, linewidth=0.5, color="black")
        ax2.set_ylabel("EMG (μV)", fontsize=10, fontweight='bold')
        ax2.tick_params(labelbottom=False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        # Panel 3: Spectrogram + theta/delta ratio
        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )
        im = ax3.pcolormesh(
            t_spec, freqs_spec, Sx_db,
            shading="auto", cmap="viridis"
        )
        ax3.set_ylabel("Frequency (Hz)", fontsize=10, fontweight='bold')
        ax3.set_ylim(0, 30)
        ax3.tick_params(labelbottom=False)
        ax3.spines['top'].set_visible(False)

        cbar = self.fig.colorbar(im, ax=ax3, pad=0.01, fraction=0.046)
        cbar.set_label("Power (dB)", fontsize=8)

        # Theta/Delta ratio overlay
        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.5, linewidth=1.5, label="Δ/θ ratio")
        ax3b.set_ylabel("Δ/θ ratio", color="white", fontsize=9)
        ax3b.tick_params(axis="y", colors="white", labelsize=8)
        ax3b.spines['top'].set_visible(False)
        ax3b.legend(loc="upper left", fontsize=7, framealpha=0.7)

        # Panel 4: Wavelet scalogram
        f_min, f_max, df = 1.0, 30.0, 0.5
        freqs_wv, scalogram_db = compute_wavelet_scalogram_matlab_like(
            eeg_win, fs, f_min=f_min, f_max=f_max, df=df,
            f0=1.0, exp_corr=0.0, time_smooth=0.0, baseline_mode="window",
        )
        n_freqs, n_times = scalogram_db.shape

        t_edges = np.linspace(win_start_sec, win_end_sec, n_times + 1)
        if len(freqs_wv) > 1:
            df_eff = freqs_wv[1] - freqs_wv[0]
        else:
            df_eff = 1.0
        f_edges = np.concatenate([freqs_wv - df_eff / 2, [freqs_wv[-1] + df_eff / 2]])

        vlim = 6.0
        im2 = ax4.pcolormesh(
            t_edges, f_edges, scalogram_db,
            shading="auto", cmap="coolwarm", vmin=-vlim, vmax=vlim,
        )
        ax4.set_ylabel("Frequency (Hz)", fontsize=10, fontweight='bold')
        ax4.set_ylim(0, 30)
        ax4.tick_params(labelbottom=False)
        ax4.spines['top'].set_visible(False)

        cbar2 = self.fig.colorbar(im2, ax=ax4, pad=0.01, fraction=0.046)
        cbar2.set_label("Power (dB, norm.)", fontsize=8)

        # Panel 5: PAC
        if PAC_ENABLED:
            try:
                if PAC_MODE == "heatmap":
                    t_edges, y_edges, pac_map = pac_mvl_heatmap(
                        eeg_win, fs, win_start_sec=win_start_sec,
                        phase_band=PAC_PHASE_BAND, amp_bands=PAC_AMP_BANDS,
                        win_sec=PAC_WIN_SEC, step_sec=PAC_STEP_SEC,
                        filt_order=4, downsample_fs=PAC_DOWNSAMPLE_FS
                    )

                    if pac_map.shape[1] > 0:
                        # Z-score normalize for better visualization
                        pac_mean = pac_map.mean()
                        pac_std = pac_map.std()
                        if pac_std > 0:
                            pac_map_z = (pac_map - pac_mean) / pac_std
                            vmin, vmax = -2, 3  # Show 2 SD below to 3 SD above mean
                        else:
                            pac_map_z = pac_map
                            vmin, vmax = pac_map.min(), pac_map.max()
                        
                        im_pac = ax_pac.pcolormesh(
                            t_edges, y_edges, pac_map_z,
                            shading="auto", cmap="magma", vmin=vmin, vmax=vmax
                        )
                        ax_pac.set_ylabel("Amp (Hz)", fontsize=10, fontweight='bold')
                        ax_pac.tick_params(labelbottom=False)
                        ax_pac.spines['top'].set_visible(False)
                        
                        cbar_pac = self.fig.colorbar(im_pac, ax=ax_pac, pad=0.01, fraction=0.046)
                        cbar_pac.set_label("PAC (z-score)", fontsize=8)
                        
                        # Add title with interpretation
                        phase_str = f"{PAC_PHASE_BAND[0]:.0f}-{PAC_PHASE_BAND[1]:.0f} Hz"
                        ax_pac.set_title(
                            f"Phase-Amplitude Coupling: {phase_str} phase modulates amplitude\n"
                            f"(Warmer colors = stronger coupling)",
                            fontsize=9, pad=5
                        )
                    else:
                        ax_pac.text(0.5, 0.5, "PAC: window too short", ha="center", va="center",
                                   transform=ax_pac.transAxes, fontsize=10)
                        ax_pac.set_axis_off()
            except Exception as e:
                ax_pac.text(0.5, 0.5, f"PAC error:\n{str(e)[:50]}", ha="center", va="center",
                           transform=ax_pac.transAxes, fontsize=9)
                ax_pac.set_axis_off()
        else:
            ax_pac.text(0.5, 0.5, "PAC disabled", ha="center", va="center",
                       transform=ax_pac.transAxes, fontsize=10)
            ax_pac.set_axis_off()

        # Panel 6: ACh trace (CLEAN VERSION - like reference image)
        if self.ach is not None and len(self.ach) > 1:
            ach = np.asarray(self.ach, dtype=float).squeeze()
            fs_ach = float(self.fs_ach)
            t_ach_full = np.arange(len(ach)) / fs_ach

            if (win_start_sec >= t_ach_full[-1]) or (win_end_sec <= t_ach_full[0]):
                ax5.text(0.5, 0.5, "ACh trace does not cover this time window.",
                        ha="center", va="center", transform=ax5.transAxes, fontsize=10)
                ax5.set_axis_off()
            else:
                # Interpolate to window
                ach_win = np.interp(t_win, t_ach_full, ach)
                
                # Process ACh trace (smoothing)
                ach_smooth, ach_deriv, ach_z = process_ach_trace(
                    ach_win, fs, smooth_window_sec=ACH_SMOOTH_WINDOW
                )
                
                # Main plot: Clean black line (like reference image)
                ax5.plot(t_win, ach_smooth, linewidth=1.5, color="black", zorder=3)
                
                # Optional: Highlight changes if enabled
                if ACH_SHOW_CHANGE_HIGHLIGHTING:
                    increases, decreases = detect_ach_changes(ach_z, threshold=ACH_CHANGE_THRESHOLD)
                    ax5.fill_between(t_win, ach_smooth.min(), ach_smooth.max(),
                                    where=increases, alpha=0.15, color="green", zorder=1)
                    ax5.fill_between(t_win, ach_smooth.min(), ach_smooth.max(),
                                    where=decreases, alpha=0.15, color="red", zorder=1)
                
                # Labels and styling (clean, like reference)
                ax5.set_ylabel("ACh (a.u.)", fontsize=10, fontweight='bold')
                ax5.set_xlabel("Time (s)", fontsize=10, fontweight='bold')
                ax5.set_xlim(win_start_sec, win_end_sec)
                
                # Adaptive y-axis based on data, with some padding
                y_min, y_max = ach_smooth.min(), ach_smooth.max()
                y_range = y_max - y_min
                y_pad = y_range * 0.15  # 15% padding
                ax5.set_ylim(y_min - y_pad, y_max + y_pad)
                
                ax5.spines['top'].set_visible(False)
                ax5.spines['right'].set_visible(False)
                # Clean look - no grid (like reference image)
                ax5.grid(False)
                
                # Optional: Add derivative subplot if enabled
                if ACH_SHOW_DERIVATIVE:
                    ax5b = ax5.twinx()
                    ax5b.plot(t_win, ach_deriv, color="blue", alpha=0.3, linewidth=0.8)
                    ax5b.set_ylabel("dACh/dt (a.u./s)", fontsize=8, color="blue")
                    ax5b.tick_params(axis="y", labelcolor="blue", labelsize=7)
                    ax5b.axhline(0, color="blue", linestyle=":", alpha=0.3, linewidth=0.5)
                    ax5b.spines['top'].set_visible(False)
        else:
            ax5.text(0.5, 0.5, "No ACh data available",
                    ha="center", va="center", transform=ax5.transAxes, fontsize=10)
            ax5.set_axis_off()

        # Main title
        self.fig.suptitle(
            f"{self.base} | {self.state_label} bout {self.idx+1}/{n_bouts}  "
            f"{start_s:.1f}–{end_s:.1f} s  (window {win_start_sec:.1f}–{win_end_sec:.1f} s)\n"
            "Navigation: ←/→ or n/p | Quit: q/Esc | Detail wavelet: drag on EEG | 3D PAC: d key\n"
            "Wavelet bands: 1=0–5, 2=0–15, 3=0–30, 4=theta, 5=spindle, 6=beta, 7=gamma (shortcuts: l/h/a)",
            fontsize=11, fontweight='bold'
        )

        self.fig.tight_layout(rect=[0, 0.02, 1, 0.95])
        self.fig.canvas.draw_idle()

    def open_pac3d(self):
        """Open 3D PAC explorer for current or selected window."""
        if self.last_span is not None:
            t0, t1 = self.last_span
        else:
            start_idx, end_idx = self.bouts_idx[self.idx]
            start_s = float(self.time_s[start_idx])
            end_s = float(self.time_s[end_idx - 1] + 1.0)
            t0 = max(0.0, start_s - PRE_SEC)
            t1 = min(len(self.eeg)/self.fs, end_s + POST_SEC)

        i0 = max(0, int(t0 * self.fs))
        i1 = min(len(self.eeg), int(t1 * self.fs))
        if i1 <= i0 + int(self.fs * PAC3D_WIN_SEC):
            print("PAC3D: selection too short.")
            return

        eeg_seg = self.eeg[i0:i1]
        win_start_sec = t0

        print(f"Computing 3D PAC cube for {t0:.1f}-{t1:.1f}s...")
        t_centers, pC, aC, cube = pac_cube_mvl(
            eeg_seg, self.fs, win_start_sec,
            phase_centers=PAC3D_PHASE_CENTERS,
            amp_centers=PAC3D_AMP_CENTERS,
            phase_bw=PAC3D_PHASE_BW,
            amp_bw=PAC3D_AMP_BW,
            win_sec=PAC3D_WIN_SEC,
            step_sec=PAC3D_STEP_SEC,
            downsample_fs=PAC3D_DOWNSAMPLE_FS
        )

        print(f"Opening 3D explorer: {len(pC)} phase × {len(aC)} amp × {len(t_centers)} time bins")
        explorer = PAC3DExplorer(t_centers, pC, aC, cube)
        plt.show()


class PAC3DExplorer:
    """Interactive 3D PAC comodulogram viewer."""
    def __init__(self, t, phase_c, amp_c, cube):
        self.t = t
        self.phase_c = phase_c
        self.amp_c = amp_c
        self.cube = cube
        self.ti = 0
        self.sel_pi = int(len(phase_c)//2)
        self.sel_ai = int(len(amp_c)//2)

        self.fig = plt.figure(figsize=(12, 6))
        self.ax_map = self.fig.add_subplot(1, 2, 1)
        self.ax_ts  = self.fig.add_subplot(1, 2, 2)

        ax_sl = self.fig.add_axes([0.12, 0.06, 0.76, 0.03])
        self.slider = Slider(ax_sl, "time idx", 0, max(len(t)-1, 0), valinit=0, valstep=1)
        self.slider.on_changed(self._on_slider)

        self.cid = self.fig.canvas.mpl_connect("button_press_event", self._on_click)

        self.im = None
        self.line, = self.ax_ts.plot([], [])
        self._draw()

    def _on_slider(self, val):
        self.ti = int(val)
        self._draw()

    def _on_click(self, event):
        if event.inaxes != self.ax_map:
            return
        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            return
        self.sel_ai = int(np.argmin(np.abs(self.amp_c - x)))
        self.sel_pi = int(np.argmin(np.abs(self.phase_c - y)))
        self._draw()

    def _draw(self):
        self.ax_map.clear()
        self.ax_ts.clear()

        if self.cube.shape[2] == 0:
            self.ax_map.text(0.5, 0.5, "Window too short for PAC cube",
                           ha="center", va="center", transform=self.ax_map.transAxes)
            self.fig.canvas.draw_idle()
            return

        C = self.cube[:, :, self.ti]

        # Robust scaling
        vmin, vmax = np.percentile(self.cube, 5), np.percentile(self.cube, 95)

        self.im = self.ax_map.imshow(
            C, origin="lower", aspect="auto",
            extent=[self.amp_c[0], self.amp_c[-1], self.phase_c[0], self.phase_c[-1]],
            vmin=vmin, vmax=vmax, cmap="magma"
        )
        self.ax_map.set_xlabel("Amplitude frequency (Hz)", fontsize=11)
        self.ax_map.set_ylabel("Phase frequency (Hz)", fontsize=11)
        self.ax_map.set_title(f"PAC Comodulogram @ t={self.t[self.ti]:.1f}s (norm MVL)", fontsize=12)

        # Crosshair
        self.ax_map.axvline(self.amp_c[self.sel_ai], color="white", linestyle="--", alpha=0.8, linewidth=1.5)
        self.ax_map.axhline(self.phase_c[self.sel_pi], color="white", linestyle="--", alpha=0.8, linewidth=1.5)

        # Time series
        y = self.cube[self.sel_pi, self.sel_ai, :]
        self.ax_ts.plot(self.t, y, linewidth=1.5)
        self.ax_ts.set_xlabel("Time (s)", fontsize=11)
        self.ax_ts.set_ylabel("PAC (norm MVL)", fontsize=11)
        self.ax_ts.set_title(
            f"PAC(t) at phase={self.phase_c[self.sel_pi]:.1f} Hz, amp={self.amp_c[self.sel_ai]:.1f} Hz",
            fontsize=12
        )
        self.ax_ts.grid(True, alpha=0.3)

        # Colorbar
        if hasattr(self, 'cbar'):
            self.cbar.remove()
        self.cbar = self.fig.colorbar(self.im, ax=self.ax_map, fraction=0.046, pad=0.04)
        self.cbar.set_label("PAC (norm MVL)", fontsize=10)

        self.fig.tight_layout(rect=[0, 0.1, 1, 1])
        self.fig.canvas.draw_idle()


def main():
    mat_files = sorted(glob.glob(os.path.join(SIGNAL_DIR, "*.mat")))
    if not mat_files:
        print(f"No .mat files found in {SIGNAL_DIR}")
        return

    print(f"Found {len(mat_files)} signal files:")
    for i, f in enumerate(mat_files):
        print(f"  {i+1}. {os.path.basename(f)}")
    
    choice = input("\nEnter file number to analyze (or press Enter for all): ").strip()
    
    if choice.isdigit():
        idx = int(choice) - 1
        if 0 <= idx < len(mat_files):
            mat_files = [mat_files[idx]]
        else:
            print(f"Invalid choice. Analyzing all files.")
    
    state_map = {1: "NREM", 2: "REM"}
    print("\nWhich state do you want to explore?")
    print("1 = NREM, 2 = REM")
    choice = input("Enter 1 or 2 (default 2=REM): ").strip()
    if choice == "1":
        target_state = 1
    else:
        target_state = 2
    state_label = state_map[target_state]

    for mat_path in mat_files:
        base = os.path.basename(mat_path).replace(".mat", "")
        csv_path = os.path.join(SCORES_DIR, base + "_scored_scores_1Hz.csv")

        if not os.path.exists(csv_path):
            print(f"⚠️  No scores CSV for {base}, skipping.")
            continue

        print(f"\n=== File: {base} ===")
        print(f"  MAT  : {mat_path}")
        print(f"  CSV  : {csv_path}")

        eeg, emg, fs, ach, fs_ach = load_mat_signal(mat_path)
        time_s, states = load_scores_1hz(csv_path)

        total_sec = len(eeg) / fs
        valid_mask = time_s < total_sec
        time_s = time_s[valid_mask]
        states = states[valid_mask]

        print(f"  Unique states in scoring: {sorted(np.unique(states))}")

        bouts_by_state = find_state_bouts(states, min_len_sec=1)
        bouts_idx = bouts_by_state.get(target_state, [])
        bouts_idx = [(s, e) for (s, e) in bouts_idx if (e - s) >= MIN_BOUT_LEN_SEC]

        if not bouts_idx:
            print(f"  No {state_label} bouts ≥ {MIN_BOUT_LEN_SEC}s, skipping file.")
            continue

        print(f"  {len(bouts_idx)} {state_label} bouts to inspect.")

        viewer = BoutViewer(base, eeg, emg, fs, time_s, states, bouts_idx, state_label, ach=ach, fs_ach=fs_ach)
        plt.show()

        if viewer.stop_all:
            break


if __name__ == "__main__":
    main()