# REM theta stability: exploratory extension

This step complements the existing PAC, spindle, and survival work without
changing its code or existing output tables.

**Primary question:** Is the organization/continuity of REM theta altered in
APP/PS1 mice, and do theta features change before REM termination?

## Run from the repository

From `~/sleep_app_ach_warm` on the feature branch:

```bash
python -m pip install numpy scipy pandas h5py
export AMBTEMP_DATA_ROOT="$HOME/Desktop/AmbTemp_proj"
python pipelines/run_07_rem_theta.py --cohort thesis --condition baseline
python pipelines/run_07_rem_theta.py --cohort all --condition warming
python pipelines/run_07_rem_theta.py --cohort all --condition donepezil
```

To run only one manifest session:

```bash
python pipelines/run_07_rem_theta.py --session-id 20251001_mouse1_baseline
```

To perform an isolated smoke test on a *scored* MAT file containing EEG,
`eeg_frequency`, and `sleep_scores`:

```bash
python pipelines/run_07_rem_theta.py \
  --single-mat /path/to/20251001_mouse1_baseline_scored.mat \
  --mouse-id mouse1 --genotype APP \
  --output /tmp/rem_theta_smoketest
```

Tests (do not need raw experimental recordings):

```bash
python -m unittest discover -s tests -p 'test_rem_theta.py' -v
```

The script uses `config/sessions.csv` and `config/mice.csv` as sources of
metadata. For the new cohort it can read scored CSVs with `t0_s` and
`final_code`; for the historical cohort it reads `sleep_scores` from MAT.
It prefers `outputs/sleep_architecture/all/standardized_scores/<session_id>.csv`
if available, following the existing `run_01_sleep_architecture.m` pipeline.

**Important:** Recording MAT must contain `eeg` and `eeg_frequency`.
Matching files must share a zero-second origin. A recording/score duration
discrepancy greater than five seconds is rejected. Resolve alignment or
different recording origins explicitly; never silently shift traces to match.

## Outputs

By default under `outputs/eeg/rem_theta/` (kept out of git):

- `rem_theta_per_bout.csv`: one row per complete or censored REM episode;
  duration, next state, theta peak frequency variability, theta-to-delta,
  median relative theta, and theta-dominant fraction/streak measures.
- `rem_theta_per_window.csv`: a 4-s periodogram every 1 s within REM,
  time since REM onset, remaining time to termination, and
  `exit_within_10s` for subsequent prospective analysis.
- `rem_theta_per_mouse.csv`: bout medians summarized **within mouse**
  × condition × ZT3–6/6–9/9–12 bin.
- `rem_theta_preexit_per_mouse.csv`: five-second pre-exit bins summarized
  within mouse before cross-genotype comparison.
- `qc_sessions.csv`: included/failed sessions, number of qualifying bouts.

All recordings are analyzed in their entirety, with bouts attributed to the
ZT bin in which they begin. For warming and donepezil sessions, ZT3–6 is
**pre-intervention**; ZT6–9 is intervention; ZT9–12 is cool-off or later drug
exposure. A baseline-session ZT6–9 bin is **not** a warming exposure.

## Definitions and scope

- Original state codes: Wake=0, NREM=1, REM=2, MA=15, undefined=-1.
  Wake runs <=15 s are reclassified as MA. REM remains strictly contiguous
  across any intervening non-REM epoch; no merging of interrupted REM.
- Eligible REM bouts last >=20 s. Bouts truncated by missing/undefined
  labels or recording end are flagged right-censored.
- EEG resampled separately within each REM bout to 200 Hz with a polyphase
  anti-aliasing resampler. A 4-s Hann spectrogram with 1-s step gives
  0.25-Hz bins. Only EEG within the scored REM bout is used.
- Theta 6–10 Hz, delta 1–4 Hz, relative theta = theta power / 1–20 Hz power.
  Theta peak frequency is tracked separately every 1 s.
- **Provisional theta dominance**: relative theta >=0.30 and theta peak
  prominence >=1.5 versus 4–6/10–12 Hz flanks. This is a thresholded,
  power-based operational measure. It is **not** Karvat's LAVI, an
  independently validated sustained-rhythm measure, or a phasic/tonic REM label.
  Test sensitivity to thresholds before treating its streak durations as robust.
- The onset/exit trajectory fields must not be used to predict an exit from
  post-event signals: each window is labeled by its **window end**.
  `exit_within_10s` is only a label, never a predictor.
- Short bouts and 4-s windowing affect frequency reliability. The script
  does not automatically run genotype hypothesis tests or Cox regression.

## Next scientific steps

1. QC 5–10 bouts visually in each cohort and review threshold sensitivity;
   check mouse-by-mouse effects before group statistics.
2. Compare theta-ridge stability and relative power in baseline APP/PS1
   versus WT, and warming versus matched baseline ZT bins.
3. Integrate Somnotate posterior probabilities as additional pre-exit
   predictors only after validating the classifier and preventing future-data
   leakage (standard HMM smoothing is not causal).
4. Add LAVI/ABBA **separately**, after testing reliability in short REM bouts.
5. Validate REM ACh photometry before interpreting within-REM ACh dynamics;
   pre-REM NREM ACh can be evaluated with the already established QC pipeline.

This initial step is descriptive, not a causal demonstration of REM fragility.
The expanded mouse cohort should be evaluated using mouse-level summaries
or mouse-clustered uncertainty, not treating 1-s windows as independent mice.
