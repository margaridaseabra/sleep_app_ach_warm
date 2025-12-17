#!/usr/bin/env python3
"""
Interactive bout viewer (REM or NREM) using per-second scoring (1 Hz).

Inputs per recording:
  - MAT: eeg, emg, eeg_frequency
  - CSV:  ..._scored_scores_1Hz.csv
          (per-second scores: 0=wake,1=NREM,2=REM,15=MA)
          Can be either:
            * one numeric column with the scores, or
            * two numeric columns (time_s, score). We auto-detect.

What the script does:
  1. Builds bouts = contiguous runs of the same score (0/1/2/15) from the
     per-second scores.
  2. You choose a target state (NREM or REM).
  3. For each bout of that state:
     - defines a window [start_s-20, end_s+30]
     - plots EEG, EMG, spectrogram+Δ/Θ ratio, Morlet wavelet scalogram,
       and PSD of bout core
     - colours the background according to *per-second* states.

Navigation (when the Matplotlib window has focus):
    → / ↓ / n : next bout
    ← / ↑ / p : previous bout
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
from scipy.signal import welch, ShortTimeFFT
from scipy.signal.windows import hamming
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import pywt
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ===========================
# Config – adjust paths here
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"

# context around bout (seconds)
PRE_SEC = 20.0   # before bout start
POST_SEC = 30.0  # after bout end

# minimum bout length (in seconds) to consider
MIN_BOUT_LEN_SEC = 5.0

# bands (Hz) for ratio (following Yizhao's code)
DELTA_BAND = (1.0, 4.0)
THETA_BAND = (4.0, 8.0)

# High-frequency band where you suspect "extra" oscillations (tune as needed)
HF_BAND = (40.0, 80.0)      # Hz

# Within-bout burst detection parameters (no external baseline)
HF_PERCENTILE = 95.0        # bursts = above this percentile of the bout's HF envelope
HF_MIN_BURST_DUR = 0.050    # seconds, min burst duration (e.g. 50 ms)

# CWT settings (for the optional PyWavelets CWT; not used in main panel now)
CWT_MIN_FREQ = 1
CWT_MAX_FREQ = 30
CWT_N_FREQS = 40
CWT_WAVELET = "cmor0.5-1.0"

STATE_COLORS = {
    0: (0.90, 0.80, 0.95, 0.30),  # Wake  -> purple-ish
    1: (1.00, 0.60, 0.60, 0.30),  # NREM  -> red-ish
    2: (0.60, 1.00, 0.60, 0.30),  # REM   -> green-ish
    15: (0.70, 0.70, 1.00, 0.30), # MA    -> bluish
}
STATE_NAMES = {0: "Wake", 1: "NREM", 2: "REM", 15: "MA"}


# ===========================
# Utility functions
# ===========================

def load_mat_signal(mat_path):
    """Load EEG, EMG and sampling rate from .mat file."""
    data = sio.loadmat(mat_path, squeeze_me=True)
    eeg = np.asarray(data["eeg"], dtype=float)
    emg = np.asarray(data["emg"], dtype=float)
    fs = float(np.squeeze(data["eeg_frequency"]))
    return eeg, emg, fs


def load_scores_1hz(csv_path):
    """
    Load per-second scoring from CSV.

    Accepts:
      - one numeric column: assumed to be the scores (0/1/2/15) at 1 Hz,
      - OR two numeric columns: time_s and scores. We auto-detect the score column.

    Returns:
      time_s : float array of seconds
      states : int array of scores (0/1/2/15)
    """
    df = pd.read_csv(csv_path)
    num_cols = list(df.select_dtypes(include="number").columns)
    if not num_cols:
        raise ValueError(f"No numeric columns found in {csv_path}")

    # pick score column
    score_col = None
    for c in num_cols:
        lc = c.lower()
        if "score" in lc or "state" in lc:
            score_col = c
            break
    if score_col is None:
        # heuristic: scores column has few unique values (< 20)
        score_col = min(num_cols, key=lambda c: df[c].nunique())

    states = df[score_col].to_numpy().astype(int)

    # pick time column if present
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
        # assume 1 Hz from 0, 1, 2, ...
        time_s = np.arange(len(states), dtype=float)

    # make sure time_s is strictly increasing
    if np.any(np.diff(time_s) <= 0):
        # fallback: regenerate from indices
        time_s = np.arange(len(states), dtype=float)

    return time_s, states


def find_state_bouts(states, min_len_sec=1):
    """
    Find contiguous bouts for each state in 'states' (1 Hz).

    Returns:
      dict: state -> list of (start_idx, end_idx) with end_idx exclusive
    """
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
                end = i  # first index of new state
                length = end - start
                if length >= min_len_sec:
                    bouts_by_state[current_state].append((start, end))
                # start new bout
                start = i
                current_state = s

    # final bout
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


# ========== Wavelet helpers (ephyviewer-style) ==========

def generate_wavelet_fourier(len_wavelet, f_start, f_stop, deltafreq,
                             sample_rate, f0, normalisation):
    """
    Compute the wavelet coefficients at all scales and compute its Fourier transform.

    Returns
    -------
    wf : array, shape (len_wavelet, n_freqs)
        Fourier transform of the wavelet coefficients (after weighting).
        Axis 0 is time; axis 1 is frequency.
    """
    # frequencies for this family
    freqs = np.arange(f_start, f_stop, deltafreq)
    scales = f0 / freqs * sample_rate

    # time vector for the wavelet window
    xi = np.arange(-len_wavelet / 2., len_wavelet / 2.)
    xsd = xi[:, np.newaxis] / scales

    # Morlet wavelet in time (for each scale)
    wavelet_coefs = np.exp(1j * 2. * np.pi * f0 * xsd) * np.exp(-xsd**2 / 2.)

    # frequency-dependent weighting
    weighting_function = lambda x: x ** (-(1.0 + normalisation))
    weighted_wavelet_coefs = wavelet_coefs * weighting_function(scales[np.newaxis, :])

    # Fourier transform of wavelets along time axis
    wf = scipy.fftpack.fft(weighted_wavelet_coefs, axis=0)
    wf = wf.conj()
    return wf


def compute_scalogram_ephyviewer(data, min_freq, max_freq, freq_resolution,
                                 fs, f0=1, exp_corr=0,
                                 time_smooth=0.5, wanted_size=5.):
    """
    Ephyviewer-style Morlet scalogram:
    - Linear frequency grid [min_freq:max_freq) with step freq_resolution.
    - Returns wt.T of shape (n_freqs, n_samples) in dB.
    """
    n_samples = len(data)

    # choose a wavelet length (power of 2, at least wanted_size seconds)
    len_wavelet = int(2 ** np.ceil(np.log(wanted_size * fs) / np.log(2)))
    sig_chunk_size = wanted_size * fs
    downsample_ratio = int(np.ceil(sig_chunk_size / len_wavelet))

    sig_chunk_size = downsample_ratio * len_wavelet
    sub_sample_rate = fs / downsample_ratio

    # Fourier-domain wavelet family
    wavelet_fourrier = generate_wavelet_fourier(
        len_wavelet,
        min_freq,
        max_freq,
        freq_resolution,
        sub_sample_rate,
        f0,
        exp_corr,
    )

    # Optional anti-alias filter for downsampling
    if downsample_ratio > 1:
        n = 8
        q = downsample_ratio
        filter_sos = scipy.signal.cheby1(n, 0.05, 0.8 / q, output="sos")
    else:
        filter_sos = None

    # Single chunk starting at 0 for now
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
        # pad signal to match wavelet length
        z = np.zeros(wavelet_fourrier.shape[0], dtype=small_sig.dtype)
        left_pad = wavelet_fourrier.shape[0] - small_sig.shape[0]
        z[:small_sig.shape[0]] = small_sig
        small_sig = z

    # avoid border effect
    small_sig -= small_sig.mean()

    small_sig_f = scipy.fftpack.fft(small_sig)
    if small_sig_f.shape[0] != wavelet_fourrier.shape[0]:
        print("oulala", small_sig_f.shape, wavelet_fourrier.shape)

    wt_tmp = scipy.fftpack.ifft(
        small_sig_f[:, np.newaxis] * wavelet_fourrier, axis=0
    )
    wt = scipy.fftpack.fftshift(wt_tmp, axes=[0])
    wt = np.abs(wt).astype("float32")

    # Convert to dB
    wt = 10 * np.log10(wt + 1e-20)

    # remove padding
    if left_pad > 0:
        wt = wt[:-left_pad]

    # Smooth in time (first axis = time, second = freq)
    if time_smooth > 0:
        sigma = time_smooth * sub_sample_rate
        wt = gaussian_filter(wt, sigma=[sigma, 0], mode="reflect")

    return wt.T  # shape (n_freqs, n_samples)



def compute_wavelet_scalogram_matlab_like(
    sig,
    fs,
    f_min=1.0,
    f_max=30.0,
    df=0.5,
    f0=1.0,
    exp_corr=0.0,
    time_smooth=0.0,
    baseline_mode="window",  # "window" or None
):
    """
    Wavelet scalogram similar to your MATLAB code:

    - Uses the FULL length of 'sig' (no 5 s truncation).
    - Linear frequency grid [f_min:f_max) with step df.
    - Morlet wavelets (ephyviewer implementation).
    - Output in dB, optionally baseline-normalised per frequency.
    """

    duration_sec = len(sig) / fs  # 👈 full window length in seconds

    # 1) Compute dB scalogram (freq x time) in absolute dB
    scalogram_db = compute_scalogram_ephyviewer(
        sig,
        min_freq=f_min,
        max_freq=f_max,
        freq_resolution=df,
        fs=fs,
        f0=f0,
        exp_corr=exp_corr,
        time_smooth=time_smooth,
        wanted_size=duration_sec,   # 👈 KEY CHANGE: no more fixed 5 s
    )

    # 2) Frequency vector
    freqs = np.arange(f_min, f_max, df)

    # 3) Optional baseline normalisation (subtract mean over time per freq)
    if baseline_mode == "window":
        baseline_db = scalogram_db.mean(axis=1, keepdims=True)
        scalogram_db = scalogram_db - baseline_db
    elif baseline_mode is None:
        pass
    else:
        raise ValueError(f"Unknown baseline_mode: {baseline_mode}")

    return freqs, scalogram_db


# ========== Optional: PyWavelets CWT helpers (not used in main viewer panel now) ==========

CWT_BW_OCT = 0.5        # wavelet bandwidth in octaves (~Q ≈ 7, like 6–7 cycles)
CWT_DELTA_OCT = None    # if None -> CWT_BW_OCT / 4 (4 freqs per bandwidth)
CWT_FREQ_SHIFT_FACTOR = 1.0
CWT_USE_OCTAVE_GRID = True  # set False to go back to plain log/linear spacing


def define_cwt_frequencies(
    foi_start: float,
    foi_end: float,
    bw_oct: float = 0.5,
    delta_oct: float | None = None,
    freq_shift_factor: float = 1.0,
):
    """
    Construct log2-spaced frequencies in octaves, similar to the MNE wavelet code.
    """
    from math import sqrt, log, log2, pi

    # Convert between bandwidth in octaves and Morlet Q (characteristic parameter)
    def bw2qt(bw: float) -> float:
        L = math.sqrt(2 * math.log(2))
        qt_val = ((2**bw + 2**(-bw) + 2) / (2**bw - 2**(-bw)) * L)
        return qt_val

    if delta_oct is None:
        delta_oct = bw_oct / 4.0

    foi = 2 ** np.arange(
        math.log2(foi_start),
        math.log2(foi_end + 1e-5),
        delta_oct,
    )
    foi *= freq_shift_factor

    foi_min = 2 * foi / (2**bw_oct + 1)
    foi_max = 2 * foi / (2**-bw_oct + 1)
    sigma_freq = (foi_max - foi_min) / (2 * math.sqrt(2 * math.log(2)))
    sigma_time = 1.0 / (2 * math.pi * sigma_freq)
    qt = bw2qt(bw_oct)
    return foi, sigma_time, sigma_freq, bw_oct, qt


def compute_cwt_morlet(
    sig,
    fs,
    f_min,
    f_max,
    n_freqs,
    wavelet_name=CWT_WAVELET,
    log_spaced=True,
):
    """
    Compute Morlet-like CWT using PyWavelets (not used in main viewer panel now).
    """
    dt = 1.0 / fs

    if CWT_USE_OCTAVE_GRID:
        freqs, sigma_time, sigma_freq, bw_oct_eff, qt_eff = define_cwt_frequencies(
            foi_start=f_min,
            foi_end=f_max,
            bw_oct=CWT_BW_OCT,
            delta_oct=CWT_DELTA_OCT,
            freq_shift_factor=CWT_FREQ_SHIFT_FACTOR,
        )
    else:
        if log_spaced:
            if f_min <= 0:
                raise ValueError("f_min must be > 0 for log-spaced frequencies.")
            freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freqs)
        else:
            freqs = np.linspace(f_min, f_max, n_freqs)

    cf = pywt.central_frequency(wavelet_name)
    scales = cf / (freqs * dt)

    coef, _ = pywt.cwt(sig, scales, wavelet_name, sampling_period=dt)
    power = np.abs(coef) ** 2
    return freqs, power


# ========== Spectrogram / ratio and other helpers ==========

def compute_spectrogram_and_ratio(eeg_win, fs, win_start_sec, window_duration=5.0, mfft=None):
    """
    Short-time FFT spectrogram and theta/delta-like ratio.
    """
    nperseg = round(fs * window_duration)
    hop = round(nperseg / 2)  # overlap 50%
    window = hamming(nperseg)

    SFT = ShortTimeFFT(
        window,
        hop=hop,
        fs=fs,
        fft_mode="onesided",
        mfft=mfft,
        scale_to="psd",
    )

    Sx = SFT.spectrogram(eeg_win)  # shape: (n_freqs, n_frames)
    time = SFT.t(len(eeg_win)) + win_start_sec

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
    """
    Middle 50% of a bout in seconds, converted to samples.
    """
    duration = end_s - start_s
    if duration <= 0:
        return 0, 0
    cs = start_s + 0.25 * duration
    ce = start_s + 0.75 * duration
    return int(cs * fs), int(ce * fs)


def shade_states(ax, time_s, states, win_start_sec, win_end_sec):
    """
    Shade background according to per-second states inside [win_start_sec, win_end_sec].
    """
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


# (HF envelope / burst helpers kept for future use – not used in plotting right now)

def hf_envelope_within_bout(power_cwt, freqs_cwt, t_cwt,
                            hf_band=HF_BAND):
    f_low, f_high = hf_band
    idx = (freqs_cwt >= f_low) & (freqs_cwt <= f_high)
    if not np.any(idx):
        raise ValueError(f"HF band {hf_band} does not overlap with CWT freqs.")
    env = power_cwt[idx, :].mean(axis=0)
    env_log = np.log10(env + 1e-20)
    median = np.median(env_log)
    mad = np.median(np.abs(env_log - median)) + 1e-12
    env_z = (env_log - median) / (1.4826 * mad)
    return t_cwt, env_z


def detect_bursts_from_env(t_env, env_z,
                           percentile=HF_PERCENTILE,
                           min_dur=HF_MIN_BURST_DUR):
    thr = np.percentile(env_z, percentile)
    above = env_z >= thr

    bursts = []
    mask = np.zeros_like(above, dtype=bool)

    in_burst = False
    start_idx = 0

    for i, flag in enumerate(above):
        if flag and not in_burst:
            in_burst = True
            start_idx = i
        elif not flag and in_burst:
            end_idx = i
            dur = t_env[end_idx - 1] - t_env[start_idx]
            if dur >= min_dur:
                bursts.append((t_env[start_idx], t_env[end_idx - 1]))
                mask[start_idx:end_idx] = True
            in_burst = False

    if in_burst:
        end_idx = len(above)
        dur = t_env[end_idx - 1] - t_env[start_idx]
        if dur >= min_dur:
            bursts.append((t_env[start_idx], t_env[end_idx - 1]))
            mask[start_idx:end_idx] = True

    return bursts, mask, thr


def summarise_bursts(bursts, bout_start, bout_end):
    if not bursts:
        return 0.0, 0.0, 0.0
    total_burst_time = sum(be - bs for bs, be in bursts)
    duration = max(bout_end - bout_start, 1e-6)
    duty_cycle = total_burst_time / duration
    rate_per_sec = len(bursts) / duration
    mean_dur = total_burst_time / len(bursts)
    return duty_cycle, rate_per_sec, mean_dur


# ===========================
# Viewer class (per bout)
# ===========================

class BoutViewer:
    def __init__(self, base, eeg, emg, fs, time_s, states, bouts_idx, state_label):
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

    def update_plot(self):
        self.fig.clear()
        n_bouts = len(self.bouts_idx)

        # 5 panels: EEG, EMG, spectrogram, wavelet scalogram, PSD
        ax1 = self.fig.add_subplot(5, 1, 1)
        ax2 = self.fig.add_subplot(5, 1, 2, sharex=ax1)
        ax3 = self.fig.add_subplot(5, 1, 3, sharex=ax1)
        ax4 = self.fig.add_subplot(5, 1, 4, sharex=ax1)
        ax5 = self.fig.add_subplot(5, 1, 5)

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
                legend_patches.append(Patch(facecolor=color, edgecolor="none",
                                            label=name))
        if legend_patches:
            ax1.legend(handles=legend_patches, loc="upper right", fontsize=8)

        ax2.plot(t_win, emg_win, linewidth=0.5, color="black")
        ax2.set_ylabel("EMG (a.u.)")

        # Spectrogram + theta/delta ratio
        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )
        im = ax3.pcolormesh(
            t_spec, freqs_spec, Sx_db,
            shading="auto", cmap="viridis"
        )
        ax3.set_ylabel("Freq (Hz)")
        ax3.set_ylim(0, 30)
        win_end_sec = win_start_sec + len(eeg_win) / fs
        ax3.set_xlim(win_start_sec, win_end_sec)

        cbar = self.fig.colorbar(im, ax=ax3)
        cbar.set_label("Power (dB)")

        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.4, linewidth=1.0)
        ax3b.set_ylabel("Theta/Delta", color="white")
        ax3b.tick_params(axis="y", colors="white")

        # MATLAB-style Morlet wavelet scalogram (ephyviewer-like)
                # --- Morlet scalogram ---
        f_min, f_max, df = 1.0, 30.0, 0.5
        freqs_wv, scalogram_db = compute_wavelet_scalogram_matlab_like(
            eeg_win,
            fs,
            f_min=f_min,
            f_max=f_max,
            df=df,
            f0=1.0,
            exp_corr=0.0,
            time_smooth=0.0,
            baseline_mode="window",
        )
        n_freqs, n_times = scalogram_db.shape

        # Duration that the CWT actually covers
        dur_wv = n_times / fs  # in seconds
        wv_start = win_start_sec
        wv_end   = min(win_start_sec + dur_wv, win_end_sec)

        # Time edges for pcolormesh: only [wv_start, wv_end]
        t_edges = np.linspace(wv_start, wv_end, n_times + 1)

        # Frequency edges as before
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
        ax4.set_ylim(f_min, f_max)
        ax4.set_title(
            f"Morlet scalogram (covers {wv_start:.1f}–{wv_end:.1f} s)"
        )

        cbar2 = self.fig.colorbar(im2, ax=ax4)
        cbar2.set_label("Power (dB rel. mean)")

         
        # PSD of bout core (middle 50%)
        cs, ce = bout_core_samples(start_s, end_s, fs)
        cs = max(cs, 0)
        ce = min(ce, n_samples)
        if ce > cs:
            eeg_core = eeg[cs:ce]
            f_psd, Pxx = compute_psd(eeg_core, fs)
            ax5.semilogy(f_psd, Pxx, label="PSD (bout core)", color="tab:blue")
            ax5.set_xlim(0, 80)
            ax5.set_xlabel("Freq (Hz)")
            ax5.set_ylabel("PSD")
            ax5.legend(loc="upper right", fontsize=8)
        else:
            ax5.text(0.5, 0.5, "Bout too short for core PSD",
                     ha="center", va="center", transform=ax5.transAxes)
            ax5.set_axis_off()

        self.fig.suptitle(
            f"{self.base} | {self.state_label} bout {self.idx+1}/{n_bouts}  "
            f"{start_s:.1f}–{end_s:.1f} s  "
            f"(window {win_start_sec:.1f}–{win_end_sec:.1f} s)\n"
            "Use ←/→ or n/p to move between bouts, q/Esc to quit.",
            fontsize=13,
        )

        self.fig.tight_layout(rect=[0, 0.06, 1, 0.93])
        self.fig.canvas.draw_idle()


# ===========================
# Main
# ===========================

def main():
    mat_files = sorted(glob.glob(os.path.join(SIGNAL_DIR, "*.mat")))
    if not mat_files:
        print(f"No .mat files found in {SIGNAL_DIR}")
        return

    print(f"Found {len(mat_files)} signal files.")

    state_map = {1: "NREM", 2: "REM"}
    print("Which state do you want to explore?")
    print("1 = NREM, 2 = REM")
    choice = input("Enter 1 or 2 (default 2=REM): ").strip()
    if choice == "1":
        target_state = 1
    else:
        target_state = 2
    state_label = state_map[target_state]

    for mat_path in mat_files:
        base = os.path.basename(mat_path).replace(".mat", "")
        csv_path = os.path.join(
            SCORES_DIR, base + "_scored_scores_1Hz.csv"
        )

        if not os.path.exists(csv_path):
            print(f"⚠️  No scores CSV for {base}, skipping.")
            continue

        print(f"\n=== File: {base} ===")
        print(f"  MAT  : {mat_path}")
        print(f"  CSV  : {csv_path}")

        eeg, emg, fs = load_mat_signal(mat_path)
        time_s, states = load_scores_1hz(csv_path)

        # trim states if longer than signal
        total_sec = len(eeg) / fs
        valid_mask = time_s < total_sec
        time_s = time_s[valid_mask]
        states = states[valid_mask]

        print(f"  Unique states in scoring: {sorted(np.unique(states))}")

        # build bouts for each state
        bouts_by_state = find_state_bouts(states, min_len_sec=1)
        bouts_idx = bouts_by_state.get(target_state, [])
        # filter by duration in seconds
        bouts_idx = [
            (s, e) for (s, e) in bouts_idx
            if (e - s) >= MIN_BOUT_LEN_SEC
        ]

        if not bouts_idx:
            print(f"  No {state_label} bouts ≥ {MIN_BOUT_LEN_SEC}s, skipping file.")
            continue

        print(f"  {len(bouts_idx)} {state_label} bouts to inspect.")

        viewer = BoutViewer(base, eeg, emg, fs, time_s, states, bouts_idx, state_label)
        plt.show()

        if viewer.stop_all:
            print("Exiting at user request.")
            return

    print("Done.")


if __name__ == "__main__":
    main()
