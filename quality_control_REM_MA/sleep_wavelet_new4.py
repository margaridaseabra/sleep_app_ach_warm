#!/usr/bin/env python3
"""
Interactive bout viewer (REM or NREM) using per-second scoring (1 Hz).
Enhanced with delta/theta wavelet panels and delta/theta power time courses.
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
from matplotlib.widgets import SpanSelector

# ===========================
# Config
# ===========================

SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"

PRE_SEC = 20.0   # before bout start
POST_SEC = 30.0  # after bout end
MIN_BOUT_LEN_SEC = 5.0

DELTA_BAND = (1.0, 4.0)
THETA_BAND = (4.0, 8.0)

HF_BAND = (40.0, 80.0)
HF_PERCENTILE = 95.0
HF_MIN_BURST_DUR = 0.050  # seconds

# Wavelet band presets (including delta/theta)
BAND_PRESETS = {
    "1": ("0–5 Hz band",   0.2, 5.0),
    "2": ("0–15 Hz band",  0.2, 15.0),
    "3": ("0–30 Hz band",  0.2, 30.0),
    "4": ("theta 6–10 Hz",   6.0, 10.0),
    "5": ("spindle 9–16 Hz", 9.0, 16.0),
    "6": ("beta 15–30 Hz",   15.0, 30.0),
    "7": ("gamma 30–80 Hz",  30.0, 80.0),
    "d": ("delta 1–4 Hz", 1.0, 4.0),
    "t": ("theta 4–8 Hz", 4.0, 8.0),
}
DEFAULT_BAND_KEY = "3"

STATE_COLORS = {
    0: (0.90, 0.80, 0.95, 0.30),
    1: (1.00, 0.60, 0.60, 0.30),
    2: (0.60, 1.00, 0.60, 0.30),
    15: (0.70, 0.70, 1.00, 0.30),
}
STATE_NAMES = {0: "Wake", 1: "NREM", 2: "REM", 15: "MA"}


# ===========================
# Utility functions
# ===========================

def load_mat_signal(mat_path):
    data = sio.loadmat(mat_path, squeeze_me=True)
    eeg = np.asarray(data["eeg"], dtype=float)
    emg = np.asarray(data["emg"], dtype=float)
    fs = float(np.squeeze(data["eeg_frequency"]))
    return eeg, emg, fs

def load_scores_1hz(csv_path):
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
        if c == score_col:
            continue
        lc = c.lower()
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
    nperseg = min(len(sig), int(4 * fs))
    f, Pxx = welch(sig, fs=fs, nperseg=nperseg)
    return f, Pxx

def generate_wavelet_fourier(len_wavelet, f_start, f_stop, deltafreq,
                             sample_rate, f0, normalisation):
    freqs = np.arange(f_start, f_stop, deltafreq)
    scales = f0 / freqs * sample_rate

    xi = np.arange(-len_wavelet / 2., len_wavelet / 2.)
    xsd = xi[:, np.newaxis] / scales

    wavelet_coefs = np.exp(1j * 2. * np.pi * f0 * xsd) * np.exp(-xsd**2 / 2.0)
    weighting_function = lambda x: x ** (-(1.0 + normalisation))
    weighted = wavelet_coefs * weighting_function(scales[np.newaxis, :])

    wf = scipy.fftpack.fft(weighted, axis=0).conj()
    return wf

def compute_scalogram_ephyviewer(data, min_freq, max_freq, freq_resolution,
                                 fs, f0=1, exp_corr=0,
                                 time_smooth=0.5, wanted_size=5.0):
    n_samples = len(data)
    len_wavelet = int(2 ** np.ceil(np.log(wanted_size * fs) / np.log(2)))
    sig_chunk_size = wanted_size * fs
    downsample_ratio = int(np.ceil(sig_chunk_size / len_wavelet))
    sig_chunk_size = downsample_ratio * len_wavelet
    sub_fs = fs / downsample_ratio

    wf = generate_wavelet_fourier(len_wavelet, min_freq, max_freq,
                                  freq_resolution, sub_fs, f0, exp_corr)

    if downsample_ratio > 1:
        sos = scipy.signal.cheby1(8, 0.05, 0.8 / downsample_ratio, output="sos")
    else:
        sos = None

    sig_chunk = data[:sig_chunk_size]
    if sig_chunk.dtype != np.float32:
        sig_chunk = sig_chunk.astype("float32")

    if downsample_ratio > 1:
        small = scipy.signal.sosfiltfilt(sos, sig_chunk)
        small = small[::downsample_ratio].copy()
    else:
        small = sig_chunk.copy()

    left_pad = 0
    if small.shape[0] != wf.shape[0]:
        z = np.zeros(wf.shape[0], dtype=small.dtype)
        left_pad = wf.shape[0] - small.shape[0]
        z[:small.shape[0]] = small
        small = z

    small -= small.mean()
    small_f = scipy.fftpack.fft(small)
    wt_tmp = scipy.fftpack.ifft(small_f[:, np.newaxis] * wf, axis=0)
    wt = scipy.fftpack.fftshift(wt_tmp, axes=[0])
    wt = np.abs(wt).astype("float32")

    wt_db = 10 * np.log10(wt + 1e-20)
    if left_pad > 0:
        wt_db = wt_db[:-left_pad]

    if time_smooth > 0:
        sigma = time_smooth * sub_fs
        wt_db = gaussian_filter(wt_db, sigma=[sigma, 0], mode="reflect")

    return wt_db.T

def compute_wavelet_scalogram_matlab_like(sig, fs, f_min=1.0, f_max=12.0,
                                          df=0.5, f0=1.0, exp_corr=0.0,
                                          time_smooth=0.0, baseline_mode="window"):
    duration_sec = len(sig) / fs
    scalogram_db = compute_scalogram_ephyviewer(
        sig, f_min, f_max, df, fs,
        f0=f0, exp_corr=exp_corr,
        time_smooth=time_smooth, wanted_size=duration_sec
    )
    freqs = np.arange(f_min, f_max, df)
    if baseline_mode == "window":
        baseline = scalogram_db.mean(axis=1, keepdims=True)
        scalogram_db = scalogram_db - baseline
    elif baseline_mode is None:
        pass
    else:
        raise ValueError(f"Unknown baseline_mode: {baseline_mode}")
    return freqs, scalogram_db

def compute_wavelet_scalogram_pretty(sig, fs, f_min, f_max,
                                     df=0.25, time_smooth=0.15, return_linear=False):
    scalogram_db = compute_scalogram_ephyviewer(
        sig, f_min, f_max, df, fs,
        f0=1.0, exp_corr=0.0,
        time_smooth=time_smooth,
        wanted_size=len(sig)/fs
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

def bout_core_samples(start_s, end_s, fs):
    duration = end_s - start_s
    cs = start_s + 0.25 * duration
    ce = start_s + 0.75 * duration
    return int(cs * fs), int(ce * fs)

def shade_states(ax, time_s, states, win_start_sec, win_end_sec):
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

def compute_spectrogram_and_ratio(sig, fs, t0=0.0):
    """
    Computes spectrogram and theta/delta ratio from EEG signal.

    Returns:
        freqs (Hz), times (s), power (dB), theta/delta ratio
    """
    nperseg = int(4.0 * fs)
    noverlap = int(3.0 * fs)
    freqs, times, Sxx = scipy.signal.spectrogram(
        sig, fs=fs, window="hann", nperseg=nperseg,
        noverlap=noverlap, scaling="density", mode="psd"
    )
    Sx_db = 10 * np.log10(Sxx + 1e-20)
    times += t0

    # Compute theta/delta ratio
    delta_idx = np.logical_and(freqs >= 1.0, freqs < 4.0)
    theta_idx = np.logical_and(freqs >= 4.0, freqs < 8.0)
    delta_power = Sxx[delta_idx].sum(axis=0)
    theta_power = Sxx[theta_idx].sum(axis=0)
    ratio = theta_power / (delta_power + 1e-20)

    return freqs, times, Sx_db, ratio

# ========== Viewer class (per bout) ==========

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

        self.detail_band_key = DEFAULT_BAND_KEY
        self.span_selector = None

        self.fig = plt.figure(figsize=(14, 10))
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
        elif event.key in ("d", "t"):
            self.detail_band_key = event.key
            band_label, f_min, f_max = BAND_PRESETS[event.key]
            print(f"Detail band set to {band_label} [{f_min}-{f_max} Hz]")
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

        mean_power_lin = power_lin.mean(axis=1)
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

        fig.suptitle(
            f"Detail wavelet {band_label}\n{t0:.2f}–{t1:.2f} s",
            fontsize=11,
        )

        fig.tight_layout()
        fig.show()

    def update_plot(self):
        self.fig.clear()
        n_bouts = len(self.bouts_idx)

        ax1 = self.fig.add_subplot(7, 1, 1)  # EEG
        ax2 = self.fig.add_subplot(7, 1, 2, sharex=ax1)  # EMG
        ax3 = self.fig.add_subplot(7, 1, 3, sharex=ax1)  # Spectrogram
        ax4_delta = self.fig.add_subplot(7, 2, 7, sharex=ax1)  # Delta scalogram
        ax4_theta = self.fig.add_subplot(7, 2, 8, sharex=ax1)  # Theta scalogram
        ax5 = self.fig.add_subplot(7, 1, 6, sharex=ax1)  # Power & ratio
        ax6 = self.fig.add_subplot(7, 1, 7, sharex=ax1)  # PSD

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

        for ax in (ax1, ax2):
            shade_states(ax, self.time_s, self.states, win_start_sec, win_end_sec)
            ax.set_xlim(win_start_sec, win_end_sec)

        ax1.plot(t_win, eeg_win, linewidth=0.5, color="black")
        ax1.set_ylabel("EEG (a.u.)")

        legend_patches = []
        for s in sorted(np.unique(self.states)):
            name = STATE_NAMES.get(int(s), str(s))
            color = STATE_COLORS.get(s)
            if color is not None:
                legend_patches.append(Patch(facecolor=color, edgecolor="none", label=name))
        if legend_patches:
            ax1.legend(handles=legend_patches, loc="upper right", fontsize=8)

        ax2.plot(t_win, emg_win, linewidth=0.5, color="black")
        ax2.set_ylabel("EMG (a.u.)")

        freqs_spec, t_spec, Sx_db, ratio = compute_spectrogram_and_ratio(
            eeg_win, fs, win_start_sec
        )
        im = ax3.pcolormesh(
            t_spec, freqs_spec, Sx_db,
            shading="auto", cmap="viridis"
        )
        ax3.set_ylabel("Freq (Hz)")
        ax3.set_ylim(0, 30)
        ax3.set_xlim(win_start_sec, win_end_sec)
        cbar = self.fig.colorbar(im, ax=ax3)
        cbar.set_label("Power (dB)")

        ax3b = ax3.twinx()
        ax3b.plot(t_spec, ratio, color="white", alpha=0.4, linewidth=1.0)
        ax3b.set_ylabel("Θ/Δ (STFT)", color="white")
        ax3b.tick_params(axis="y", colors="white")

        band_defs = [
            (DELTA_BAND[0], DELTA_BAND[1], ax4_delta, "Delta 1–4 Hz"),
            (THETA_BAND[0], THETA_BAND[1], ax4_theta, "Theta 4–8 Hz")
        ]
        power_bands = {}
        for f_min, f_max, ax_band, label in band_defs:
            freqs, power_norm, power_lin = compute_wavelet_scalogram_pretty(
                eeg_win, fs, f_min=f_min, f_max=f_max,
                df=0.25, time_smooth=0.15, return_linear=True
            )
            t_edges = np.linspace(win_start_sec, win_end_sec, power_norm.shape[1] + 1)
            if len(freqs) > 1:
                df_eff = freqs[1] - freqs[0]
            else:
                df_eff = 1.0
            f_edges = np.concatenate([freqs - df_eff / 2, [freqs[-1] + df_eff / 2]])

            im2 = ax_band.pcolormesh(
                t_edges, f_edges, power_norm,
                shading="auto", cmap="turbo", vmin=0.0, vmax=1.0
            )
            ax_band.set_ylabel("Hz")
            ax_band.set_ylim(f_min, f_max)
            ax_band.set_title(label)

            band_power = power_lin.mean(axis=0)
            power_bands[label] = (t_edges[:-1], band_power)
            ax_band.plot(t_edges[:-1], band_power, color='white', alpha=0.6, lw=1.0)

        t_delta, delta_pow = power_bands["Delta 1–4 Hz"]
        t_theta, theta_pow = power_bands["Theta 4–8 Hz"]
        ax5.plot(t_delta, delta_pow, label='Delta (1–4 Hz)', color='blue')
        ax5.plot(t_theta, theta_pow, label='Theta (4–8 Hz)', color='green')
        ratio_dt = delta_pow / (theta_pow + 1e-12)
        ax5.plot(t_delta, ratio_dt, label='Δ/Θ ratio', color='purple', linestyle='--')
        ax5.set_ylabel("Power / Ratio")
        ax5.legend(loc="upper right", fontsize=8)

        cs, ce = bout_core_samples(start_s, end_s, fs)
        cs = max(cs, 0)
        ce = min(ce, n_samples)
        if ce > cs + 10:
            eeg_core = eeg[cs:ce]
            f_psd, Pxx = compute_psd(eeg_core, fs)
            ax6.semilogy(f_psd, Pxx, label="PSD (bout core)", color="tab:blue")
            ax6.set_xlim(0, 80)
            ax6.set_xlabel("Freq (Hz)")
            ax6.set_ylabel("PSD")
            ax6.legend(loc="upper right", fontsize=8)
        else:
            ax6.text(0.5, 0.5, "Bout too short for core PSD",
                     ha="center", va="center", transform=ax6.transAxes)
            ax6.set_axis_off()

        self.fig.suptitle(
            f"{self.base} | {self.state_label} bout {self.idx+1}/{n_bouts}  "
            f"{start_s:.1f}–{end_s:.1f} s  "
            f"(window {win_start_sec:.1f}–{win_end_sec:.1f} s)\n"
            "Use ←/→ or n/p to move between bouts, q/Esc to quit. "
            "Drag on EEG to open detail wavelet. Bands: "
            "1=0–5, 2=0–15, 3=0–30, 4=theta, 5=spindle, 6=beta, 7=gamma, d=delta, t=theta.",
            fontsize=12,
        )

        self.fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        self.fig.canvas.draw_idle()


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
        csv_path = os.path.join(SCORES_DIR, base + "_scored_scores_1Hz.csv")
        if not os.path.exists(csv_path):
            print(f"⚠️  No scores CSV for {base}, skipping.")
            continue

        print(f"\n=== File: {base} ===")
        print(f"  MAT  : {mat_path}")
        print(f"  CSV  : {csv_path}")

        eeg, emg, fs = load_mat_signal(mat_path)
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

        viewer = BoutViewer(base, eeg, emg, fs, time_s, states, bouts_idx, state_label)
        plt.show()
        if viewer.stop_all:
            print("Exiting at user request.")
            return

    print("Done.")


if __name__ == "__main__":
    main()
