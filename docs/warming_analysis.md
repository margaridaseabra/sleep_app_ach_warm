# Ambient warming sleep analysis

This folder documents the reproducible ambient-warming sleep analysis for the expanded APP/PS1 and WT cohort.

## Experimental structure

The warming recording is divided into:

- `ZT3_6`: pre-warming at baseline room temperature
- `ZT6_9`: ambient warming
- `ZT9_12`: cool-down

Fixed windows require at least 95% valid scored coverage. Wake bouts <=15 s are reclassified as micro-arousals (MA), matching the baseline pipeline.

## Run order

From the MATLAB repository root:

```matlab
run('pipelines/run_00_prepare_data.m')
run('pipelines/run_01_sleep_architecture.m')
run('pipelines/run_warming_complete.m')
run('pipelines/run_07_export_warming_prism.m')
```


## Actual warming timing

The source-of-truth timing registry is:

```text
config/warming_times.csv
```

It contains both the October and June batches. The preferred machine-readable
fields are:

```text
warming_start_elapsed_s
warming_end_elapsed_s
```

and the same values are also retained as human-readable:

```text
warming_start_elapsed_hms
warming_end_elapsed_hms
```

The October rows were converted directly from the existing timing workbook, so
sub-minute precision is preserved. The June rows were transcribed from the
handwritten lab notebook and are marked as `nearest_minute`.

The analysis prefers the numeric seconds columns and falls back to `HH:MM:SS`
when needed.

For each mouse the analysis windows are:

- `ZT3_6` (legacy phase label): up to 3 h immediately before actual warming onset.
  If warming starts before 3 h elapsed, all available data from recording start
  to warming onset are used.
- `ZT6_9` (legacy phase label): actual registered warming onset to actual
  registered warming end.
- `ZT9_12` (legacy phase label): 3 h immediately following warming end.

These names are retained only for compatibility with the thesis/Prism workflow.
The exact numerical boundaries for every mouse are written to
`warming_raw_metrics_long.csv` as `window_start_s` and `window_end_s`, together
with `timing_status`.

The same mouse-specific boundaries are applied to the separate baseline
recording for the matched baseline-vs-warming and difference-in-differences
analyses.

The timing registry also contains provenance fields (`batch`, `timing_source`,
`source_file`, `precision`, `notes`) and `include_analysis`.

A second registry is stored at:

```text
config/drug_times.csv
```

For June it stores both the saline and donepezil timestamps. The October source
workbook contains only the drug timestamp, so the October saline fields are
left blank rather than inferred.

## Analysis hierarchy

### 1. Thesis-compatible within-session response

For each mouse:

```text
delta_bin = warming_session_bin - warming_session_ZT3_6
```

The pre-warming bin is therefore zero by definition.

This reproduces the response framing used in the thesis.

### 2. Matched baseline-vs-warming comparison

Each mouse also has a separate baseline recording. The pipeline compares baseline and warming sessions at the same ZT window.

The key direct manipulation model at `ZT6_9` is:

```text
Value ~ Genotype * Condition + Cohort + (1|Mouse)
```

for the combined cohort.

The genotype x condition term asks whether warming affects APP/PS1 differently from WT.

### 3. Difference-in-differences sensitivity analysis

To remove ordinary time-of-day change more explicitly:

```text
DiD =
  (warming_target - warming_pre)
  -
  (baseline_target - baseline_pre)
```

This is calculated for `ZT6_9` and `ZT9_12`.

### 4. Survival

NREM and REM bout survival are analyzed with:

- mouse-balanced Kaplan-Meier curves
- Cox models with mouse-clustered robust SE
- warming-only genotype comparisons
- baseline-vs-warming genotype x condition models
- leave-one-mouse-out sensitivity for the condition interaction

## Outputs

Regenerable analysis outputs:

```text
outputs/sleep_architecture/warming_complete/
```

Version-controlled Prism source data:

```text
results/prism/warming/
```

### Prism folders

- `source/`: one mouse per row with raw warming, delta, matched baseline, and DiD values.
- `grouped/`: two-header-row Prism grouped tables with WT and APP mouse replicates already arranged.
- `survival/`: exact mouse-balanced KM coordinates.
- `stats/`: MATLAB statistical outputs.

## Primary interpretation

The primary warming question is not simply whether APP/PS1 and WT differ during the warming session. It is whether the **response to warming differs by genotype**, with particular emphasis on REM sleep expression and REM bout maintenance.

For final inference prioritize:

1. `ZT6_9` genotype x condition effect (baseline vs warming)
2. within-session delta from pre
3. DiD sensitivity analysis
4. `ZT9_12` as persistence/cool-down
5. survival interaction and its leave-one-mouse-out robustness

## Git policy

Commit:

- `pipelines/run_warming_complete.m`
- `pipelines/run_07_export_warming_prism.m`
- `docs/warming_analysis.md`
- `results/prism/warming/`

Do not commit `outputs/` or standardized per-second score files.


## Validation safeguards

The v2 pipeline explicitly enforces one row per mouse x metric x window for the delta, matched-condition, and difference-in-differences tables. It also reads figure metadata directly from `config/mice.csv`. This prevents missing sex values in the historical cohort from creating duplicate analysis rows.

Mixed-model ANOVA uses a Satterthwaite degrees-of-freedom method when supported by the installed MATLAB version and falls back to the default `anova(lme)` implementation otherwise.
