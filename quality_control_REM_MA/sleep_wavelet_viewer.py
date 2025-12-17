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
     - plots EEG, EMG, spectrogram+Δ/Θ ratio, Morlet CWT, and PSD of bout core
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
import scipy.io as sio
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


# CWT settings (to look at high-frequency bursts)
CWT_MIN_FREQ = 1
CWT_MAX_FREQ = 12
CWT_N_FREQS = 40
#CWT_WAVELET = "cmor1.5-1.0"
#CWT_WAVELET = "cmor3.0-1.0"  # more cycles (better freq, worse time)
CWT_WAVELET = "cmor0.5-1.0"  # fewer cycles (better time, worse freq)

# Colour for each *state* (per-second scores):
# 0 = wake, 1 = NREM, 2 = REM, 15 = MA
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


# Log-spaced / Linear Morlet CWT using PyWavelets (true or false)
# New CWT config (you can tweak these)
CWT_BW_OCT = 0.5        # wavelet bandwidth in octaves (~Q ≈ 7, like 6–7 cycles)
CWT_DELTA_OCT = None    # if None -> CWT_BW_OCT / 4 (4 freqs per bandwidth)
CWT_FREQ_SHIFT_FACTOR = 1.0
CWT_USE_OCTAVE_GRID = True  # set False to go back to plain log/linear spacing


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
    Compute Morlet-like CWT using PyWavelets.

    If CWT_USE_OCTAVE_GRID is True, ignore n_freqs and build an
    octave/log2-spaced frequency grid using define_cwt_frequencies().
    Otherwise, use either log- or linearly-spaced freqs between f_min/f_max.

    Returns
    -------
    freqs : array (n_freqs,)
        Frequencies in Hz.
    power : array (n_freqs, n_times)
        Wavelet power (|coef|^2).
    """
    dt = 1.0 / fs

    # --- Choose frequency grid ---
    if CWT_USE_OCTAVE_GRID:
        freqs, sigma_time, sigma_freq, bw_oct_eff, qt_eff = define_cwt_frequencies(
            foi_start=f_min,
            foi_end=f_max,
            bw_oct=CWT_BW_OCT,
            delta_oct=CWT_DELTA_OCT,
            freq_shift_factor=CWT_FREQ_SHIFT_FACTOR,
        )
        # (sigma_time, sigma_freq, bw_oct_eff, qt_eff are available if you want to inspect)
    else:
        if log_spaced:
            if f_min <= 0:
                raise ValueError("f_min must be > 0 for log-spaced frequencies.")
            freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freqs)
        else:
            freqs = np.linspace(f_min, f_max, n_freqs)

    # --- Map frequencies -> scales for this wavelet ---
    cf = pywt.central_frequency(wavelet_name)
    scales = cf / (freqs * dt)

    # --- CWT via PyWavelets ---
    coef, _ = pywt.cwt(sig, scales, wavelet_name, sampling_period=dt)
    power = np.abs(coef) ** 2

    return freqs, power




def compute_spectrogram_and_ratio(eeg_win, fs, win_start_sec, window_duration=5.0, mfft=None):
    """
    DIRECT port of Yizhao's get_fft_plots(), but for a windowed EEG segment.

    Inputs
    ------
    eeg_win      : 1D array, EEG for full window (20s + bout + 30s)
    fs           : sampling frequency
    win_start_sec: absolute start time (in seconds) of eeg_win in the recording
    window_duration : STFT window length in seconds (default 5)
    mfft         : passed through to ShortTimeFFT (usually None)

    Returns
    -------
    freqs        : frequencies (<= 30 Hz)
    time         : absolute time centres (same as Plotly version)
    Sx_db_smooth : spectrogram in dB, Gaussian-smoothed
    theta_delta  : theta/delta ratio (delta_power / theta_power), smoothed
    """
    nperseg = round(fs * window_duration)
    hop = round(nperseg / 2)  # this is "noverlap" in the original code
    window = hamming(nperseg)

    SFT = ShortTimeFFT(
        window,
        hop=hop,
        fs=fs,
        fft_mode="onesided",
        mfft=mfft,
        scale_to="psd",
    )

    # STFT on the *windowed* EEG
    Sx = SFT.spectrogram(eeg_win)  # shape: (n_freqs, n_frames)

    # time centres in seconds, shifted by win_start_sec (exactly like start_time)
    time = SFT.t(len(eeg_win)) + win_start_sec

    frequencies = SFT.f
    freq_mask = frequencies <= 30.0
    frequencies = frequencies[freq_mask]
    Sx = Sx[freq_mask, :]

    # dB
    Sx_db = 10 * np.log10(Sx + 1e-20)

    # delta/theta masks exactly as in the original function
    delta_mask = np.where((frequencies > 1) & (frequencies <= 4))[0]
    theta_mask = np.where((frequencies > 4) & (frequencies <= 8))[0]

    delta_power = np.mean(Sx_db[delta_mask, :], axis=0)
    theta_power = np.mean(Sx_db[theta_mask, :], axis=0)

    theta_delta_ratio = delta_power / (theta_power + 1e-12)
    theta_delta_ratio = gaussian_filter1d(theta_delta_ratio, 4)

    # smooth spectrogram
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
    # find indices overlapping the window
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
            t1 = t0 + 1.0  # assume 1 s for last one
        left = max(t0, win_start_sec)
        right = min(t1, win_end_sec)
        if right > left:
            ax.axvspan(left, right, color=color, zorder=0)

def hf_envelope_within_bout(power_cwt, freqs_cwt, t_cwt,
                            hf_band=HF_BAND):
    """
    Compute a broadband high-frequency envelope from the CWT and
    normalise it *within that bout* (robust z, using median+MAD).

    power_cwt : 2D array (freq x time) from Morlet CWT (linear power).
    freqs_cwt: 1D array of frequencies.
    t_cwt    : 1D array of times (absolute seconds).
    hf_band  : (f_low, f_high) in Hz.

    Returns:
        t_env   : same as t_cwt
        env_z   : robust z-scored HF envelope
    """
    f_low, f_high = hf_band
    idx = (freqs_cwt >= f_low) & (freqs_cwt <= f_high)
    if not np.any(idx):
        raise ValueError(f"HF band {hf_band} does not overlap with CWT freqs.")

    # broadband HF envelope: average power over band
    env = power_cwt[idx, :].mean(axis=0)        # shape (time,)

    # log to compress dynamic range
    env_log = np.log10(env + 1e-20)

    # robust normalisation within this bout: median + MAD
    median = np.median(env_log)
    mad = np.median(np.abs(env_log - median)) + 1e-12
    env_z = (env_log - median) / (1.4826 * mad)   # ≈ z-score if Gaussian

    return t_cwt, env_z


def detect_bursts_from_env(t_env, env_z,
                           percentile=HF_PERCENTILE,
                           min_dur=HF_MIN_BURST_DUR):
    """
    Turn a within-bout HF envelope (env_z) into bursts:

    - threshold = given percentile of env_z (e.g. 95th).
    - bursts = contiguous segments above threshold
               with duration >= min_dur.

    Returns:
        bursts     : list of (t_start, t_end)
        mask       : boolean array, True inside bursts
        thr_value  : threshold in env_z units
    """
    # high percentile threshold within this bout
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

    # burst extending to the end
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

def define_cwt_frequencies(
    foi_start: float,
    foi_end: float,
    bw_oct: float = 0.5,
    delta_oct: float | None = None,
    freq_shift_factor: float = 1.0,
):
    """
    Construct log2-spaced frequencies in octaves, similar to the MNE wavelet code.

    Parameters
    ----------
    foi_start : float
        Lowest frequency of interest (Hz).
    foi_end : float
        Highest frequency of interest (Hz).
    bw_oct : float
        Wavelet bandwidth in octaves. Larger -> broader wavelet in freq (more smoothing).
        Roughly corresponds to Morlet 'Q' (number of cycles).
    delta_oct : float | None
        Step size in octaves between consecutive frequencies.
        If None, uses bw_oct / 4 (about 4 wavelets per bandwidth).
    freq_shift_factor : float
        Factor to shift the entire frequency grid in log2-space (usually keep at 1.0).

    Returns
    -------
    foi : np.ndarray
        Centre frequencies in Hz (1D array).
    sigma_time : np.ndarray
        Approximate temporal SD of each wavelet (seconds).
    sigma_freq : np.ndarray
        Approximate spectral SD of each wavelet (Hz).
    bw_oct : float
        Effective bandwidth in octaves (returned for convenience).
    qt : float
        Approximate Morlet Q (characteristic "number of cycles"-like parameter).
    """
    from math import sqrt, log, log2, pi

    # Convert between bandwidth in octaves and Morlet Q (characteristic parameter)
    def bw2qt(bw: float) -> float:
        L = sqrt(2 * log(2))
        qt_val = ((2**bw + 2**(-bw) + 2) / (2**bw - 2**(-bw)) * L)
        return qt_val

    if delta_oct is None:
        delta_oct = bw_oct / 4.0

    # log2-spaced freqs: 2^(k * delta_oct) * foi_start
    foi = 2 ** np.arange(
        math.log2(foi_start),
        math.log2(foi_end + 1e-5),
        delta_oct,
    )
    foi *= freq_shift_factor

    # Approximate min/max of each wavelet in freq, then sigma_freq, sigma_time
    foi_min = 2 * foi / (2**bw_oct + 1)    # arithmetic mean approx
    foi_max = 2 * foi / (2**-bw_oct + 1)
    sigma_freq = (foi_max - foi_min) / (2 * sqrt(2 * log(2)))
    sigma_time = 1.0 / (2 * pi * sigma_freq)

    qt = bw2qt(bw_oct)

    return foi, sigma_time, sigma_freq, bw_oct, qt


# ===========================
# Viewer class (per bout)
# ===========================

class BoutViewer:
    def __init__(self, base, eeg, emg, fs, time_s, states, bouts_idx, state_label):
        """
        bouts_idx : list of (start_idx, end_idx) index pairs for the chosen state.
        """
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
        # print("key:", event.key)  # uncomment to debug on your backend
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

        # 5 panels: EEG, EMG, spectrogram, CWT, PSD
        ax1 = self.fig.add_subplot(5, 1, 1)
        ax2 = self.fig.add_subplot(5, 1, 2, sharex=ax1)
        ax3 = self.fig.add_subplot(5, 1, 3, sharex=ax1)
        ax4 = self.fig.add_subplot(5, 1, 4, sharex=ax1)
        ax5 = self.fig.add_subplot(5, 1, 5)

        start_idx, end_idx = self.bouts_idx[self.idx]
        # per-second times: bout spans [time_s[start_idx], time_s[end_idx-1] + 1)
        start_s = float(self.time_s[start_idx])
        end_s = float(self.time_s[end_idx - 1] + 1.0)

        fs = self.fs
        eeg = self.eeg
        emg = self.emg

        n_samples = len(eeg)
        total_sec = n_samples / fs

        # ---- analysis window around current bout ----
        win_start_sec = max(0.0, start_s - PRE_SEC)
        win_end_sec = min(total_sec, end_s + POST_SEC)
        win_start_sample = int(win_start_sec * fs)
        win_end_sample = int(win_end_sec * fs)

        eeg_win = eeg[win_start_sample:win_end_sample]
        emg_win = emg[win_start_sample:win_end_sample]
        t_win = np.arange(win_start_sample, win_end_sample) / fs  # absolute time

        # ---- coloured background according to per-second states ----
        for ax in [ax1, ax2, ax3, ax4]:
            shade_states(ax, self.time_s, self.states, win_start_sec, win_end_sec)
            ax.set_xlim(win_start_sec, win_end_sec)

        # ---- raw EEG & EMG ----
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

        # ---- spectrogram + Δ/Θ-like ratio ----
        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )

        # pcolormesh with time/freq *centres* (like Plotly) – let Matplotlib handle shading
        im = ax3.pcolormesh(
            t_spec, freqs_spec, Sx_db,
            shading="auto", cmap="viridis"
        )
        ax3.set_ylabel("Freq (Hz)")
        ax3.set_ylim(0, 30)

        # make sure x-limits match the full window (EEG/EMG)
        win_end_sec = win_start_sec + len(eeg_win) / fs
        ax3.set_xlim(win_start_sec, win_end_sec)

        cbar = self.fig.colorbar(im, ax=ax3)
        cbar.set_label("Power (dB)")

        # theta/delta ratio line – exactly like Plotly version
        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.4, linewidth=1.0)
        ax3b.set_ylabel("Theta/Delta", color="white")
        ax3b.tick_params(axis="y", colors="white")


        # ---- Morlet CWT ----
       
        # Compute CWT on EEG window
        # Compute CWT using the existing function
        freqs_cwt, power_cwt = compute_cwt_morlet(
            eeg_win, fs, f_min=1.0, f_max=12, n_freqs=60
        )

        # Convert to dB
        cwt_db = 10 * np.log10(power_cwt + 1e-20)

        # Per-frequency z-score across time to highlight transients
        mean_f = np.mean(cwt_db, axis=1, keepdims=True)
        std_f = np.std(cwt_db, axis=1, keepdims=True) + 1e-12
        cwt_z = (cwt_db - mean_f) / std_f
    
        # Time and frequency edges for pcolormesh
        n_freqs, n_times = cwt_z.shape
        t_edges = np.linspace(win_start_sec, win_end_sec, n_times + 1)
        f_edges = np.linspace(freqs_cwt[0], freqs_cwt[-1], n_freqs + 1)

        # Plot z-scored CWT with a symmetric colour scale
        vmin, vmax = -3.0, 3.0  # ±3 s.d. highlights bursts
        im2 = ax4.pcolormesh(
            t_edges, f_edges, cwt_z,
            shading="auto", cmap="viridis",
            vmin=vmin, vmax=vmax
        )
        ax4.set_ylabel("Frequency (Hz)")
        ax4.set_xlabel("Time (s)")
        ax4.set_title("Morlet CWT (1–50 Hz, z‑scored)")
        self.fig.colorbar(im2, ax=ax4, label="Power (z)")


        # ---- PSD of bout core (middle 50% of this bout) ----
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

        # ---- title ----
        n_bouts = len(self.bouts_idx)
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

    # choose state
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
