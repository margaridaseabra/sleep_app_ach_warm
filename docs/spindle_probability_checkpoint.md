# Spindle / probability checkpoint (2026-10-22)

## Purpose

Freeze the **actual tested MATLAB scripts from the local `~/sleep_app_ach_warm` repository**, rather than overwriting them with earlier ZIP versions. The following detector run completed on 47 sessions: 34,400 detected spindle events and 3,891 NREM transitions. The original raw-EEG processing should not need to run again for probability-only analyses.

## Order of operations

1. Run `prepare_github_checkpoint.sh` from Terminal. It checks Git origin, fetches references, creates a new **local branch**, and stages source code and a small set of mouse-level output CSVs.
2. Review `git diff --cached` and especially inspect any added `config/` paths for machine-specific content or confidential information.
3. Commit and push only once satisfied. Confirm the pushed GitHub branch and its contents in a browser.
4. Keep a separate local or institutional backup of **raw recordings**, `outputs/eeg/spindles/` event-level CSVs, QC figures, and any licensed/private data. These are deliberately excluded from Git.

## Essential files to preserve

- `pipelines/run_10_spindle_analysis.m`: source spindle detector, including the tested compatibility fix for the 50-Hz notch (`iirnotch` or Butterworth fallback).
- `pipelines/run_12_export_spindle_prism.m`: figure-source export.
- `pipelines/run_13_plot_spindle_results.m` (if present): spindle density / transition exploratory plots.
- `pipelines/run_14_spindle_probability.m`: 1-second pre-transition spindle occupancy and mouse-level outcomes.
- `pipelines/run_15_spindle_probability_stats.m`: cohort-stratified tests and baseline-controlled intervention contrasts. Ensure it includes the MATLAB-compatible `[char(name) '.png']` style filename fix.
- `config/mice.csv`, `config/sessions.csv`, `config/paths_config.m` (if present): genotype and session mapping.
- `docs/spindle_analysis.md` and relevant `tools/` utilities (if present).

## Analysis convention

- Biological units for genotype inference: **mice**, not bouts, detected spindles or overlapping epochs.
- Baseline matching period: ZT3–9, with October/June cohort considered in comparisons.
- Analyze NREM→REM and NREM→Wake separately and compare their within-mouse contrast.
- Time-resolved measure: spindle occupancy for each second before the transition.
- Primary late-transition measures: mean occupancy over [-25, -15) seconds, and the fraction of transitions with any spindle in the final 25 seconds.
- Keep the full-window eligibility rules (at least 25 s NREM; optionally 60 s for extended curves) and report how many bouts and mice contribute.
- For warming and donepezil, preserve within-mouse pre/post windows and corresponding baseline-matched changes (difference-in-differences). A non-significant genotype-by-intervention contrast is **not** evidence of no effect within either genotype.

## Interpretation of exploratory completed analyses

In the expanded 16-mouse sample (8 WT / 8 APP/PS1), baseline pre-REM spindle occupancy in the [-25, -15) s interval was lower in APP/PS1, with a negative effect in both cohorts. A similar reduction appeared in the probability of any spindle in the final 25 s. These observations should be presented with individual-mouse results, cohort stratification, and detector QC before strong biological interpretation. Warming/donepezil probability results did not provide evidence for genotype-specific rescue.

## Next module: PAC

Keep the existing spindle detector frozen. The PAC module should reuse the **same scored bouts and transition timestamps**, while operating on raw EEG for time-frequency estimates. Validate low-frequency phase extraction, time-window duration, nonstationarity, surrogate nulls, and mouse-level statistics on a small pilot before computing all sessions. Do not mix a new PAC implementation with the thesis results without labeling the methodological differences.
