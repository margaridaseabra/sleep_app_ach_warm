# Spindle analysis

This module performs one spindle detection analysis across baseline,
ambient-warming and donepezil sessions.

## Reference method

The detector is based on the wavelet-driven spindle detector described by
Bandarabadi et al. (Nature Communications, 2020), which was also the basis of
the MSc-thesis spindle analysis.

## Locked primary detector

- EEG resampled to 200 Hz.
- 50-Hz notch filtering.
- complex frequency B-spline wavelet `fbsp2-1-2`.
- wavelet-energy band: 9–16 Hz in 0.5-Hz steps.
- no extra frequency weighting in the paper-primary detector.
- 200-ms Hann smoothing.
- detection threshold: NREM-like mean + 3 SD.
- onset/offset threshold: NREM-like mean + 1 SD.
- maximum rejection threshold: mean + 20 SD.
- values > mean + 10 SD are excluded while estimating thresholds.
- duration: 0.4–2.0 s.
- ≥80% overlap with NREM-like time.
- 5–30 spindle-band oscillatory cycles.
- spindle-band PSD peak must exceed the peak in both neighboring
  6–8.5-Hz and 16.5–20-Hz bands.

Short Wake bouts ≤15 s are first reclassified as microarousals, matching the
sleep-architecture pipeline. Microarousals are assimilated into NREM when
defining NREM-like macro-bouts.

## Primary transition endpoint

The pre-specified primary transition endpoint is spindle rate during the final
25 s before:

- NREM → REM
- NREM → Wake

The pipeline additionally exports a 60-s, 1-Hz spindle-occupancy trajectory
before each transition for descriptive plots.

Bout-level observations are summarized within mouse before group inference.

## Analysis windows

The same detector is used for:

- baseline `BASELINE_ZT3_9`;
- warming `WARM_PRE` and the actual registered `WARM_ACTIVE` period;
- donepezil `DRUG_PRE` and the actual `DRUG_0_3H` period;
- mouse-specific matched baseline windows for warming and donepezil.

This permits both simple within-session change and matched
difference-in-differences analyses.

## Run order

```matlab
cd('/path/to/sleep_app_ach_warm')

run('pipelines/run_10_spindle_analysis.m')
run('pipelines/run_12_export_spindle_prism.m')
```

Raw outputs are written under:

```text
outputs/eeg/spindles/
```

Prism-ready outputs are written under:

```text
results/prism/spindles/
```

## Main Figure 1 files

- `baseline_spindle_density_ZT3_9.csv`
- `baseline_preREM_vs_preWake_spindle_rate.csv`
- `baseline_transition_timecourse.csv`

## Main Figure 3 / intervention files

- `warming_spindle_density_response.csv`
- `warming_transition_spindle_metrics.csv`
- `donepezil_spindle_density_response.csv`
- `donepezil_transition_spindle_metrics.csv`

The PAC pipeline should reuse the same NREM transition definitions rather than
re-deriving a separate set of transitions.
