import os, glob, re
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.signal
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# ===========================
# CONFIG
# ===========================
SIGNAL_DIR = "/Users/margaridaseabra/24.11 signalnotscored"
SCORES_DIR = "/Users/margaridaseabra/24.11scores"

# xlsx containing ambient warming start/end per recording
WARMING_XLSX = "/Users/margaridaseabra/sleep_app_ach_warm/Ambtemp prep-3h-selection/ambtemp-times-all-mice.xlsx"

# bands of interest
DELTA_BAND = (0.5, 5.0)
SIGMA_BAND = (10.0, 15.0)

# power timebase to align with scoring (seconds)
POWER_BIN_SEC = 1.0

# smoothing in seconds for the bandpower traces (set 0 to disable)
POWER_SMOOTH_SEC = 5.0

# plot EEG downsample 
EEG_MAX_PLOT_POINTS = 250_000

STATE_COLORS = {
    0: (0.90, 0.80, 0.95, 0.35),  # Wake
    1: (1.00, 0.60, 0.60, 0.35),  # NREM
    2: (0.60, 1.00, 0.60, 0.35),  # REM
    15:(0.70, 0.70, 1.00, 0.35),  # MA
}
STATE_NAMES = {0: "Wake", 1: "NREM", 2: "REM", 15: "MA"}

# Files to exclude from notch filtering
NOTCH_EXCLUDE_FILES = [
    "20251001_baseline_mouse1_APP.mat",
    "20251002_baseline_mouse2_WT.mat",
    "20251003_ambtemp_mouse1_APP.mat",
    "20251005_baseline_mouse8_WT.mat",
    "20251006_baseline_mouse4_WT.mat"
]

# Notch filter settings (50 Hz powerline and harmonics)
NOTCH_FREQ = 50.0  # Hz
NOTCH_Q = 30.0     # Quality factor (higher = narrower notch)
NOTCH_HARMONICS = [1, 2, 3, 4, 5]  # Apply at 50, 100, 150, 200, 250 Hz

# ===========================
# Functions
# ===========================
def apply_notch_filter(signal, fs, freq=50.0, Q=30.0, harmonics=[1, 2, 3, 4, 5]):
    """
    Apply notch filter at specified frequency and its harmonics.
    
    Parameters:
    -----------
    signal : array
        Input signal
    fs : float
        Sampling frequency
    freq : float
        Fundamental frequency to notch (e.g., 50 Hz)
    Q : float
        Quality factor (higher = narrower notch)
    harmonics : list
        List of harmonic multipliers to apply (e.g., [1, 2, 3] for 50, 100, 150 Hz)
    
    Returns:
    --------
    filtered_signal : array
        Signal with notch filter applied
    """
    filtered = signal.copy()
    nyq = fs / 2.0
    
    for h in harmonics:
        notch_f = freq * h
        if notch_f >= nyq:
            print(f"Skipping {notch_f} Hz (above Nyquist frequency {nyq} Hz)")
            continue
        
        # Design notch filter
        b, a = scipy.signal.iirnotch(notch_f, Q, fs)
        
        # Apply filter
        filtered = scipy.signal.filtfilt(b, a, filtered)
        print(f"  ✓ Applied notch at {notch_f} Hz")
    
    return filtered

def load_mat_signal(mat_path):
    data = sio.loadmat(mat_path, squeeze_me=True)
    eeg = np.asarray(data["eeg"], dtype=float).squeeze()
    emg = np.asarray(data["emg"], dtype=float).squeeze()
    fs  = float(np.squeeze(data["eeg_frequency"]))
    
    # Check if this file should be notch filtered
    basename = os.path.basename(mat_path)
    if basename not in NOTCH_EXCLUDE_FILES:
        print(f"  Applying 50 Hz notch filter and harmonics...")
        eeg = apply_notch_filter(eeg, fs, NOTCH_FREQ, NOTCH_Q, NOTCH_HARMONICS)
    else:
        print(f"Skipping notch filter (file in exclusion list)")
    
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

# ===========================
# Ambient warming intervals from XLSX
# ===========================
def _to_seconds(x):
    """Accept numeric seconds, hh:mm:ss strings, pandas Timestamps, or timedeltas."""
    if pd.isna(x):
        return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    if isinstance(x, pd.Timedelta):
        return x.total_seconds()
    if isinstance(x, pd.Timestamp):
        # if you store clock times only, pandas may parse as Timestamp.
        # interpret as time-of-day and convert to seconds since midnight:
        return x.hour * 3600 + x.minute * 60 + x.second + x.microsecond * 1e-6
    s = str(x).strip()
    # hh:mm:ss
    m = re.match(r"^(\d+):(\d+):(\d+(?:\.\d+)?)$", s)
    if m:
        hh, mm, ss = m.groups()
        return int(hh) * 3600 + int(mm) * 60 + float(ss)
    # mm:ss
    m = re.match(r"^(\d+):(\d+(?:\.\d+)?)$", s)
    if m:
        mm, ss = m.groups()
        return int(mm) * 60 + float(ss)
    # fallback numeric-in-string
    try:
        return float(s)
    except:
        return None

def load_warming_intervals(xlsx_path, base):
    """
    Returns list of (start_sec, end_sec) for this recording base.
    Expected columns:
      - Column 0: Mouse identifier (e.g., "1 WT", "2 WT", "3 APP")
      - Column 1: Genotype (WT, APP)
      - Columns with 'start'/'started' and 'finish'/'end' for timing
    """
    if not xlsx_path or not os.path.exists(xlsx_path):
        print(f"Warming XLSX not found: {xlsx_path}")
        return []

    # Only load warming intervals for ambtemp recordings
    if "ambtemp" not in base.lower():
        print(f"  No warming intervals (not an ambtemp recording)")
        return []

    df = pd.read_excel(xlsx_path)

    # find start/end columns (more flexible matching)
    cols = {c.lower(): c for c in df.columns}
    start_col = next((cols[k] for k in cols if "start" in k), None)
    end_col   = next((cols[k] for k in cols if "finish" in k or ("end" in k and "start" not in k)), None)
    
    if start_col is None or end_col is None:
        print(f"⚠️  XLSX missing start/end columns: {list(df.columns)}")
        return []

    # Extract mouse number from filename (e.g., "mouse1" from "20251003_ambtemp_mouse1_APP")
    mouse_match = re.search(r'mouse\s*(\d+)', base, re.IGNORECASE)
    if not mouse_match:
        print(f"Could not extract mouse number from {base}")
        return []
    
    mouse_num = mouse_match.group(1)  # just the number
    
    # First column contains mouse identifier - match the number
    first_col = df.columns[0]
    mask = df[first_col].astype(str).str.contains(rf'\b{mouse_num}\b', case=False, regex=True, na=False)
    sub = df[mask].copy()
    
    if sub.empty:
        print(f"No warming intervals found for mouse {mouse_num} in {base}")
        return []

    intervals = []
    for _, row in sub.iterrows():
        s = _to_seconds(row[start_col])
        e = _to_seconds(row[end_col])
        if s is None or e is None:
            continue
        if e < s:
            s, e = e, s
        intervals.append((float(s), float(e)))

    # de-dup and sort
    intervals = sorted(set((round(a, 6), round(b, 6)) for a, b in intervals))
    print(f"  Found {len(intervals)} warming interval(s)")
    return [(a, b) for a, b in intervals]

# ===========================
# Signal processing
# ===========================
def bandpower_1hz_hilbert(eeg, fs, band, bin_sec=1.0, smooth_sec=0.0):
    """Bandpass -> Hilbert amp -> power -> average in 1s bins -> (optional) smooth."""
    nyq = 0.5 * fs
    lo, hi = band
    lo = max(0.001, float(lo))
    hi = min(float(hi), nyq * 0.999)
    sos = scipy.signal.butter(4, [lo/nyq, hi/nyq], btype="bandpass", output="sos")
    x = scipy.signal.sosfiltfilt(sos, eeg)

    amp = np.abs(scipy.signal.hilbert(x))
    p = amp**2

    total_sec = len(p) / fs
    n_bins = int(np.floor(total_sec / bin_sec))
    if n_bins <= 1:
        return np.array([]), np.array([])

    edges_s = np.arange(0, n_bins + 1) * bin_sec
    edges_i = np.clip((edges_s * fs).astype(int), 0, len(p))
    # ensure strictly increasing
    edges_i = np.unique(edges_i)
    if len(edges_i) < 2:
        return np.array([]), np.array([])

    sums = np.add.reduceat(p, edges_i[:-1])
    counts = np.diff(edges_i).astype(float)
    mean_p = sums / np.maximum(counts, 1.0)

    t = (edges_i[:-1] / fs)  # seconds
    # log power looks nicer across hours
    mean_db = 10.0 * np.log10(mean_p + 1e-20)

    if smooth_sec and smooth_sec > 0:
        sigma = smooth_sec / bin_sec
        mean_db = gaussian_filter1d(mean_db, sigma=sigma, mode="nearest")

    return t, mean_db

def compress_state_segments(time_s, states):
    """Return list of (t0, t1, state) segments for fast shading."""
    time_s = np.asarray(time_s, float)
    states = np.asarray(states, int)
    if len(time_s) == 0:
        return []
    segs = []
    s0 = states[0]
    t0 = time_s[0]
    for i in range(1, len(states)):
        if states[i] != s0:
            t1 = time_s[i]
            segs.append((t0, t1, int(s0)))
            s0 = states[i]
            t0 = time_s[i]
    # last
    segs.append((t0, time_s[-1] + 1.0, int(s0)))
    return segs

# ===========================
# Plotting
# ===========================
def plot_full_recording(mat_path):
    base = os.path.basename(mat_path).replace(".mat", "")
    csv_path = os.path.join(SCORES_DIR, base + "_scored_scores_1Hz.csv")
    if not os.path.exists(csv_path):
        print(f"⚠️  No scores CSV for {base}, skipping.")
        return

    print(f"Loading data for {base}...")
    eeg, emg, fs = load_mat_signal(mat_path)
    print(f"  EEG: {len(eeg)} samples, fs={fs} Hz")
    
    time_s, states = load_scores_1hz(csv_path)
    print(f"  Scores: {len(time_s)} epochs")

    total_sec = len(eeg) / fs
    valid = time_s < total_sec
    time_s = time_s[valid]
    states = states[valid]
    segs = compress_state_segments(time_s, states)
    print(f"  State segments: {len(segs)}")

    warm_intervals = load_warming_intervals(WARMING_XLSX, base)

    # power traces (1 Hz)
    print("Computing bandpower...")
    tD, delta_db = bandpower_1hz_hilbert(eeg, fs, DELTA_BAND, bin_sec=POWER_BIN_SEC, smooth_sec=POWER_SMOOTH_SEC)
    tS, sigma_db = bandpower_1hz_hilbert(eeg, fs, SIGMA_BAND, bin_sec=POWER_BIN_SEC, smooth_sec=POWER_SMOOTH_SEC)
    print(f"  Delta: {len(tD)} points, Sigma: {len(tS)} points")

    # downsample EEG for plotting
    n = len(eeg)
    step = max(1, int(np.ceil(n / EEG_MAX_PLOT_POINTS)))
    t_eeg = (np.arange(0, n, step) / fs)
    eeg_p = eeg[::step]
    print(f"  EEG downsampled to {len(eeg_p)} points")

    print("Creating plot...")
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(4, 1, height_ratios=[1.3, 1.0, 1.0, 0.6], hspace=0.15)

    ax_eeg   = fig.add_subplot(gs[0])
    ax_delta = fig.add_subplot(gs[1], sharex=ax_eeg)
    ax_sigma = fig.add_subplot(gs[2], sharex=ax_eeg)
    ax_hyp   = fig.add_subplot(gs[3], sharex=ax_eeg)

    # shade states on EEG/delta/sigma
    for ax in (ax_eeg, ax_delta, ax_sigma):
        for (a, b, st) in segs:
            col = STATE_COLORS.get(st, None)
            if col is not None:
                ax.axvspan(a, b, color=col, zorder=0)

    # shade warming intervals across all panels
    if len(warm_intervals) > 0:
        print(f"  Adding {len(warm_intervals)} warming interval(s)...")
        for (ws, we) in warm_intervals:
            for ax in (ax_eeg, ax_delta, ax_sigma, ax_hyp):
                ax.axvspan(ws, we, color=(1.0, 0.85, 0.20, 0.20), zorder=1)
            ax_eeg.text(ws, 0.98, "Warming ON", transform=ax_eeg.get_xaxis_transform(),
                        va="top", ha="left", fontsize=9)

    # EEG panel
    ax_eeg.plot(t_eeg, eeg_p, color="black", lw=0.4)
    ax_eeg.set_ylabel("EEG (μV)")
    
    # Add notch filter indicator to title
    basename = os.path.basename(mat_path)
    notch_status = "" if basename in NOTCH_EXCLUDE_FILES else " [50Hz notch applied]"
    ax_eeg.set_title(f"{base}  |  fs={fs:.1f} Hz  |  duration={total_sec/3600:.2f} h{notch_status}")

    # Delta panel
    if len(tD) > 0:
        ax_delta.plot(tD, delta_db, color="black", lw=1.0)
    ax_delta.set_ylabel("Delta power\n(0.5–5 Hz, dB)")

    # Sigma panel
    if len(tS) > 0:
        ax_sigma.plot(tS, sigma_db, color="black", lw=1.0)
    ax_sigma.set_ylabel("Sigma power\n(10–15 Hz, dB)")

    # Hypnogram 
    ax_hyp.set_ylim(0, 1)
    ax_hyp.set_yticks([])
    
    # Fill entire panel with light grey background
    ax_hyp.axhspan(0, 1, color=(0.95, 0.95, 0.95, 1.0), zorder=0)
    
    # Overlay red during warming periods
    if len(warm_intervals) > 0:
        for (ws, we) in warm_intervals:
            ax_hyp.axvspan(ws, we, color=(1.0, 0.4, 0.4, 0.6), zorder=1)
    
    ax_hyp.set_xlabel("Time (s)")
    ax_hyp.set_ylabel("Warming")

    # legends
    patches = []
    for st in sorted(set(states)):
        st = int(st)
        if st in STATE_COLORS:
            patches.append(Patch(facecolor=STATE_COLORS[st], edgecolor="none", label=STATE_NAMES.get(st, str(st))))
    if len(warm_intervals) > 0:
        patches.append(Patch(facecolor=(1.0, 0.85, 0.20, 0.20), edgecolor="none", label="Ambient warming window"))
    ax_eeg.legend(handles=patches, loc="upper right", fontsize=9, framealpha=0.95)

    # "making it pretty"
    for ax in (ax_eeg, ax_delta, ax_sigma):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax_hyp.spines["top"].set_visible(False)
    ax_hyp.spines["right"].set_visible(False)
    ax_hyp.spines["left"].set_visible(False)

    plt.tight_layout()
    print("✓ Plot ready!")
    plt.show()

def main():
    mat_files = sorted(glob.glob(os.path.join(SIGNAL_DIR, "*.mat")))
    if not mat_files:
        print(f"No .mat files found in {SIGNAL_DIR}")
        return

    print(f"Found {len(mat_files)} .mat files.")
    
    for i, f in enumerate(mat_files):
        basename = os.path.basename(f)
        notch_marker = "" if basename in NOTCH_EXCLUDE_FILES else " [notch]"
        ambtemp_marker = " [ambtemp]" if "ambtemp" in basename.lower() else ""
        print(f"{i+1:3d}. {basename}{notch_marker}{ambtemp_marker}")
    
    idx = int(input("\nFile number: ").strip()) - 1
    if 0 <= idx < len(mat_files):
        plot_full_recording(mat_files[idx])
    else:
        print("Invalid index.")

if __name__ == "__main__":
    main()
