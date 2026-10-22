# Donepezil sleep-analysis plan

## Why this pipeline differs from the thesis version

The thesis described the drug session using nominal ZT3–6 pre-manipulation,
ZT6–9 acute drug and ZT9–12 later/washout bins. The expanded dataset now has
recorded injection times for every active donepezil session, and those injection
times vary substantially across animals.

The publication pipeline therefore uses injection-aligned windows:

- `PRE`: up to 3 h immediately before the recorded donepezil injection.
- `DRUG_0_3H`: injection to +3 h. This is the **primary acute drug window**.
- `DRUG_3_6H`: +3 h to +6 h. This is a secondary late-response window.

The same elapsed-time boundaries are applied to each mouse's separate baseline
recording.

## Main analyses

For state occupancy, bout rate, bout duration and fragmentation:

1. Raw donepezil-session values.
2. Within-session change from `PRE`.
3. Matched baseline vs donepezil at the same elapsed time.
4. Difference-in-differences:
   `(drug target - drug PRE) - (baseline target - baseline PRE)`.
5. Mixed-effects models with mouse as random intercept and cohort adjustment
   for the combined cohort.

For NREM and REM bout survival:

1. Donepezil-only genotype comparison.
2. **Primary acute-response survival test:** within donepezil session,
   `PRE` vs `DRUG_0_3H`, testing genotype × phase.
3. Matched baseline-vs-donepezil genotype × condition sensitivity analysis.
4. Leave-one-mouse-out sensitivity.

## Saline

June sessions contain recorded saline-injection times; the October notebook
does not. Saline timing is therefore not imputed for the primary combined
analysis. A June-only saline/injection-control sensitivity analysis can be
added separately and must be labelled as such.

## Cohorts

Combined models adjust for cohort (`original/thesis` vs `new`). Original and
new cohorts are also exported separately for replication/sensitivity checks.

## Run order

```matlab
run('pipelines/run_00_prepare_data.m')
run('pipelines/run_01_sleep_architecture.m')
run('pipelines/run_donepezil_complete.m')
run('pipelines/run_08_export_donepezil_prism.m')
```

Outputs:

- `outputs/sleep_architecture/donepezil_complete/`
- `results/prism/donepezil/`
