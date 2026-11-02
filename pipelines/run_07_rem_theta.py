#!/usr/bin/env python3
"""REM theta stability analysis for the sleep_app_ach_warm session manifest.

Outputs are descriptive and exploratory: no claim that theta dominance is
equivalent to Karvat LAVI, phasic REM, or a validated REM-fragility biomarker.
"""
from __future__ import annotations

import argparse
import os
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, signal

FS_TARGET = 200.0
THETA = (6.0, 10.0)
WINDOW_S = 4.0
STEP_S = 1.0
MIN_REM_S = 20.0


def mat_vars(path: Path, keys: tuple[str, ...]) -> dict:
    try:
        return io.loadmat(str(path), variable_names=list(keys), squeeze_me=True)
    except (NotImplementedError, ValueError) as exc:
        # MATLAB v7.3/HDF5 fallback. Do not hide unrelated load errors.
        try:
            import h5py
            if not h5py.is_hdf5(path):
                raise exc
            with h5py.File(path) as f:
                return {k: np.asarray(f[k]).squeeze() for k in keys if k in f}
        except ImportError as err:
            raise RuntimeError("MATLAB v7.3 files require h5py") from err


def score_vector(path: Path) -> np.ndarray:
    if path.suffix.lower() == ".mat":
        data = mat_vars(path, ("sleep_scores",))
        if "sleep_scores" not in data:
            raise ValueError(f"No sleep_scores in {path}")
        return np.asarray(data["sleep_scores"], dtype=float).reshape(-1)
    df = pd.read_csv(path)
    tcol = "time_s" if "time_s" in df.columns else "t0_s"
    scol = "score" if "score" in df.columns else "final_code"
    if tcol not in df or scol not in df:
        raise ValueError(f"Score CSV requires time_s/score or t0_s/final_code: {path}")
    times = pd.to_numeric(df[tcol], errors="coerce").to_numpy(dtype=float)
    values = pd.to_numeric(df[scol], errors="coerce").to_numpy(dtype=float)
    if len(times) == 0 or not np.isfinite(times).all():
        raise ValueError("Invalid score time axis")
    times = times - times[0]
    if np.any(np.diff(times) <= 0) or np.max(np.abs(times - np.round(times))) > 0.15:
        raise ValueError("Expected distinct, integer-second score epochs")
    idx = np.round(times).astype(int)
    scores = np.full(int(idx.max()) + 1, -1.0)
    scores[idx] = values
    return scores


def reclassify_short_wake(scores: np.ndarray, maximum: int = 15) -> np.ndarray:
    out = scores.copy()
    starts = np.r_[0, np.where(np.diff(out) != 0)[0] + 1]
    ends = np.r_[starts[1:], len(out)]
    for a, b in zip(starts, ends):
        if out[a] == 0 and b - a <= maximum:
            out[a:b] = 15
    return out


def rem_bouts(scores: np.ndarray) -> list[tuple[int, int, int | None, bool]]:
    """Return REM [start,end) in seconds; next state and end censoring."""
    starts = np.r_[0, np.where(np.diff(scores) != 0)[0] + 1]
    ends = np.r_[starts[1:], len(scores)]
    bouts = []
    for a, b in zip(starts, ends):
        if scores[a] != 2:
            continue
        censored = b == len(scores)
        next_state = None if censored else int(scores[b])
        bouts.append((int(a), int(b), next_state, bool(censored)))
    return bouts


def period_label(second: float) -> str:
    if second < 10800:
        return "ZT3_6"
    if second < 21600:
        return "ZT6_9"
    if second < 32400:
        return "ZT9_12"
    return "OTHER"


def _band_power(f: np.ndarray, pxx: np.ndarray, lo: float, hi: float) -> np.ndarray:
    select = (f >= lo) & (f <= hi)
    if np.sum(select) < 2:
        return np.full(pxx.shape[1], np.nan)
    return np.trapezoid(pxx[select, :], x=f[select], axis=0)


def theta_features(eeg: np.ndarray, fs: float, start_s: int, end_s: int):
    """One row per 1-s step, using centered 4-s EEG windows within REM.

    No smoothing across the scored REM boundary; no data after a feature
    window's end enter that feature. All spectral features are power-based.
    """
    duration = end_s - start_s
    if duration < max(MIN_REM_S, WINDOW_S):
        return pd.DataFrame()
    i0 = int(round(start_s * fs))
    i1 = int(round(end_s * fs))
    raw = np.asarray(eeg[i0:i1], dtype=np.float64)
    if len(raw) < WINDOW_S * fs:
        return pd.DataFrame()
    if not np.isfinite(raw).all():
        return pd.DataFrame()
    ratio = Fraction(FS_TARGET / float(fs)).limit_denominator(20000)
    x = signal.resample_poly(raw, ratio.numerator, ratio.denominator)
    nper = round(WINDOW_S * FS_TARGET)
    hop = round(STEP_S * FS_TARGET)
    if len(x) < nper:
        return pd.DataFrame()
    f, t, pxx = signal.spectrogram(
        x, fs=FS_TARGET, window="hann", nperseg=nper,
        noverlap=nper - hop, nfft=nper, detrend="constant",
        scaling="density", mode="psd",
    )
    theta = _band_power(f, pxx, *THETA)
    delta = _band_power(f, pxx, 1, 4)
    broad = _band_power(f, pxx, 1, 20)
    theta_rel = theta / np.maximum(broad, np.finfo(float).tiny)
    theta_delta = theta / np.maximum(delta, np.finfo(float).tiny)
    band = (f >= THETA[0]) & (f <= THETA[1])
    th_idx = np.flatnonzero(band)
    peaks = np.argmax(pxx[band, :], axis=0)
    ridge = f[th_idx[peaks]]
    peak_power = pxx[th_idx[peaks], np.arange(len(t))]
    flank = ((f >= 4) & (f < 6)) | ((f > 10) & (f <= 12))
    flank_ref = np.median(pxx[flank, :], axis=0)
    prominence = peak_power / np.maximum(flank_ref, np.finfo(float).tiny)
    # Explicitly a provisional, power-based theta-dominance threshold.
    dominant = (theta_rel >= 0.30) & (prominence >= 1.5)
    window_end = t + WINDOW_S / 2
    return pd.DataFrame({
        "bout_time_s": t,
        "observed_until_s": window_end,
        "seconds_to_end_s": duration - window_end,
        "theta_power": theta,
        "delta_power": delta,
        "theta_relative_1_20": theta_rel,
        "theta_delta_ratio": theta_delta,
        "theta_peak_hz": ridge,
        "theta_peak_prominence": prominence,
        "theta_dominant": dominant.astype(int),
    })


def _streak_metrics(flags: np.ndarray) -> tuple[float, float]:
    flags = np.asarray(flags, dtype=bool)
    edges = np.diff(np.r_[False, flags, False].astype(int))
    lengths = np.where(edges == -1)[0] - np.where(edges == 1)[0]
    return (float(np.max(lengths)) if len(lengths) else 0.0,
            float(np.median(lengths)) if len(lengths) else 0.0)


def process_session(recording: Path, score_path: Path, meta: dict):
    m = mat_vars(recording, ("eeg", "eeg_frequency"))
    if "eeg" not in m or "eeg_frequency" not in m:
        raise ValueError(f"Missing eeg or eeg_frequency in {recording}")
    eeg = np.asarray(m["eeg"]).reshape(-1)
    fs = float(np.asarray(m["eeg_frequency"]).squeeze())
    if not np.isfinite(fs) or fs <= 25:
        raise ValueError(f"Invalid eeg_frequency={fs}")
    scores = reclassify_short_wake(score_vector(score_path))
    available_s = len(eeg) / fs
    if abs(available_s - len(scores)) > 5:
        raise ValueError(
            f"EEG/score duration mismatch: {available_s:.1f} versus {len(scores)} s. "
            "Check score-to-EEG time alignment before analysis."
        )
    bout_rows, window_frames = [], []
    for number, (a, b, next_state, right_cens) in enumerate(rem_bouts(scores), start=1):
        if b > available_s + 0.5:
            continue
        windows = theta_features(eeg, fs, a, b)
        if windows.empty:
            continue
        windows["session_id"] = meta["session_id"]
        windows["bout_id"] = number
        windows["start_s"] = a
        windows["duration_s"] = b - a
        windows["next_state"] = "CENSORED" if next_state is None else {0: "WK", 1: "NREM", 15: "MA"}.get(next_state, str(next_state))
        windows["right_censored"] = right_cens
        windows["time_bin"] = period_label(a)
        windows["time_s"] = a + windows["observed_until_s"]
        # For genuine exits the outcome is known; for session-end censoring,
        # only windows with >=10 s available future qualify as negatives.
        if right_cens:
            windows["exit_within_10s"] = np.where(windows["seconds_to_end_s"] >= 10, 0.0, np.nan)
        else:
            windows["exit_within_10s"] = (windows["seconds_to_end_s"] <= 10).astype(float)
        flags = windows["theta_dominant"].to_numpy(dtype=bool)
        longest, typical = _streak_metrics(flags)
        middle = windows[(windows["bout_time_s"] >= 10) & (windows["seconds_to_end_s"] >= 20)]
        final = windows[windows["seconds_to_end_s"] <= 10]
        row = {
            **meta, "bout_id": number, "start_s": a, "end_s": b,
            "duration_s": b-a, "time_bin": period_label(a),
            "next_state": windows["next_state"].iloc[0],
            "right_censored": right_cens,
            "n_windows": len(windows),
            "theta_peak_median_hz": float(windows["theta_peak_hz"].median()),
            "theta_peak_sd_hz": float(windows["theta_peak_hz"].std(ddof=0)),
            "theta_ridge_step_median_hz": float(np.median(np.abs(np.diff(windows["theta_peak_hz"])))) if len(windows) > 1 else np.nan,
            "theta_relative_median": float(windows["theta_relative_1_20"].median()),
            "theta_delta_median": float(windows["theta_delta_ratio"].median()),
            "theta_dominant_fraction": float(flags.mean()),
            "theta_dominant_longest_streak_s": longest,
            "theta_dominant_median_streak_s": typical,
            "theta_relative_middle": float(middle["theta_relative_1_20"].mean()) if len(middle) else np.nan,
            "theta_relative_last10": float(final["theta_relative_1_20"].mean()) if len(final) and not right_cens else np.nan,
            "theta_peak_sd_middle_hz": float(middle["theta_peak_hz"].std(ddof=0)) if len(middle) > 1 else np.nan,
            "theta_peak_sd_last10_hz": float(final["theta_peak_hz"].std(ddof=0)) if len(final) > 1 and not right_cens else np.nan,
        }
        bout_rows.append(row)
        window_frames.append(windows)
    if not bout_rows:
        return pd.DataFrame(), pd.DataFrame()
    bouts = pd.DataFrame(bout_rows)
    frames = pd.concat(window_frames, ignore_index=True)
    for key, value in meta.items():
        frames[key] = value
    return bouts, frames


def summarize(bouts: pd.DataFrame, windows: pd.DataFrame):
    grouping = ["mouse_id", "genotype", "cohort", "condition", "time_bin"]
    numeric = [c for c in bouts.columns if c.startswith("theta_") and c not in {"theta_relative_last10"}]
    numeric.append("theta_relative_last10")
    mouse = bouts.groupby(grouping, dropna=False).agg(
        n_rem_bouts=("bout_id", "count"),
        median_duration_s=("duration_s", "median"),
        **{c: (c, "median") for c in numeric},
    ).reset_index()
    # Time-to-exit averaged FIRST within mouse to prevent treating windows
    # or bouts as independent biological replicates.
    w = windows.loc[~windows["right_censored"]].copy()
    w = w[(w["seconds_to_end_s"] >= 0) & (w["seconds_to_end_s"] <= 30)]
    w["seconds_before_end_bin"] = (5 * np.floor(w["seconds_to_end_s"] / 5)).astype(int)
    w = w.groupby(grouping + ["seconds_before_end_bin"], dropna=False).agg(
        n_windows=("theta_relative_1_20", "size"),
        theta_relative_mean=("theta_relative_1_20", "mean"),
        theta_peak_mean_hz=("theta_peak_hz", "mean"),
        theta_dominant_fraction=("theta_dominant", "mean"),
    ).reset_index()
    return mouse, w


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--data-root", type=Path, default=Path(os.environ.get(
        "AMBTEMP_DATA_ROOT", str(Path.home() / "Desktop" / "AmbTemp_proj"))))
    ap.add_argument("--output", type=Path)
    ap.add_argument("--cohort", choices=["all", "thesis", "new"], default="all")
    ap.add_argument("--condition", choices=["all", "baseline", "warming", "donepezil"], default="all")
    ap.add_argument("--session-id", default="")
    ap.add_argument("--single-mat", type=Path, help="Test a scored MAT directly (sleep_scores inside MAT)")
    ap.add_argument("--single-score", type=Path, help="Optional alternative score CSV for --single-mat")
    ap.add_argument("--mouse-id", default="mouse1", help="Metadata for --single-mat")
    ap.add_argument("--genotype", default="APP", choices=["APP", "WT"])
    args = ap.parse_args()
    out = args.output or args.repo / "outputs" / "eeg" / "rem_theta"
    out.mkdir(parents=True, exist_ok=True)

    sessions = []
    if args.single_mat:
        sessions = [(args.single_mat, args.single_score or args.single_mat, {
            "session_id": args.single_mat.stem, "mouse_id": args.mouse_id,
            "genotype": args.genotype, "cohort": "test", "condition": "baseline",
        })]
    else:
        manifest = pd.read_csv(args.repo / "config" / "sessions.csv").fillna("")
        mice = pd.read_csv(args.repo / "config" / "mice.csv").fillna("")
        genotypes = dict(zip(mice["mouse_id"], mice["genotype"]))
        for _, s in manifest.iterrows():
            if str(s["include_sleep"]).lower() in {"0", "false", "no"}:
                continue
            if args.cohort != "all" and s["cohort"] != args.cohort:
                continue
            if args.condition != "all" and s["condition"] != args.condition:
                continue
            if args.session_id and s["session_id"] != args.session_id:
                continue
            mat = args.data_root / s["recording_file"]
            score = args.data_root / s["score_file"]
            standard = args.repo / "outputs" / "sleep_architecture" / "all" / "standardized_scores" / f"{s['session_id']}.csv"
            if standard.exists():
                score = standard
            sessions.append((mat, score, {
                "session_id": s["session_id"], "mouse_id": s["mouse_id"],
                "genotype": genotypes.get(s["mouse_id"], "UNKNOWN"),
                "cohort": s["cohort"], "condition": s["condition"],
            }))

    all_bouts, all_windows, qc = [], [], []
    for i, (rec, sc, meta) in enumerate(sessions, start=1):
        try:
            if not rec.is_file() or not sc.is_file():
                raise FileNotFoundError(f"Missing recording or scoring: {rec} / {sc}")
            print(f"[{i}/{len(sessions)}] {meta['session_id']}", flush=True)
            bouts, windows = process_session(rec, sc, meta)
            if not bouts.empty:
                all_bouts.append(bouts)
                all_windows.append(windows)
            qc.append({**meta, "status": "OK", "n_bouts": len(bouts), "message": ""})
        except Exception as e:
            warnings.warn(f"{meta['session_id']}: {e}")
            qc.append({**meta, "status": "ERROR", "n_bouts": 0, "message": str(e)})
    pd.DataFrame(qc).to_csv(out / "qc_sessions.csv", index=False)
    if not all_bouts:
        print(f"No eligible bouts. See {out / 'qc_sessions.csv'}")
        return
    bouts = pd.concat(all_bouts, ignore_index=True)
    windows = pd.concat(all_windows, ignore_index=True)
    mice, preexit = summarize(bouts, windows)
    bouts.to_csv(out / "rem_theta_per_bout.csv", index=False)
    windows.to_csv(out / "rem_theta_per_window.csv", index=False)
    mice.to_csv(out / "rem_theta_per_mouse.csv", index=False)
    preexit.to_csv(out / "rem_theta_preexit_per_mouse.csv", index=False)
    print(f"Saved {len(bouts)} REM bouts from {len(mice)} mouse/condition/bin groups to {out}")


if __name__ == "__main__":
    main()
