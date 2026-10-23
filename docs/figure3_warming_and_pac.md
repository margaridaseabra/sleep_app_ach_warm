# Figure 3 and Figure 1 supplementary PAC — source and interpretation

## Figure 3 — effects of ambient warming

This is a descriptive data figure generated from existing `warming_complete` tables. It makes **no new statistical claims**. Data plotted as difference-in-differences are mouse-level values, using matched time windows:

`(warming active - warming pre) - (baseline matched active - baseline matched pre)`.

Figure 3 panels show a protocol schematic (individual warm start/end times vary), REM percent time, REM mean bout duration, NREM percent time, microarousal percent time, and REM bout termination hazard ratios by cohort. Main active period: ZT6–9. Hazard estimates and p-values are reproduced from the original cluster-aware survival-model reference table; mouse-level figure panels are descriptive, with exploratory original Welch p-values for genotype contrasts.

**Do not label the full panel 'rescue' as if rescue were established.** There are different directions/magnitudes by phenotype and cohort, and no directly measured core temperature or targeted preoptic intervention in the supplied dataset. The October ZT6–9 REM Cox genotype × warming p=.0037 did not replicate in June p=.638; pooled p=.709. The ZT6–9 NREM percent-time genotype DiD is p=.0269 uncorrected, Sidak p=.0530 for 2 phases. NREM percent-time DiD is WT +16.64 pp, APP -0.38 pp, so this is not evidence of an APP-specific NREM rescue.

Survival curves: descriptive KM at ZT6–9, separated by October/June and genotype with baseline/warming conditions. Bouts are not independent mice; don't compute inferential p-values by naïvely treating bouts as independent.

## Extended Data: PAC

`Fig1S_PAC_mouse_level` shows mouse-level baseline Morlet ROI PAC:
- NREM slow-frequency phase 1–4 Hz; sigma amplitude 9–16 Hz, before NREM→REM and NREM→Wake
- REM theta phase 4–10 Hz; gamma amplitude 30–80 Hz

Full EEG-derived Morlet comodulogram PNG/PDF figures are copied from `real_comodulograms_v22_fine/figures` **only if already present**. No maps are invented or interpolated in this script. The group heatmaps are descriptive, not pixelwise inferential tests. Baseline Morlet preREM−preWake genotype contrast p~.017, BH q~.069 (four outcomes), exploratory. No confirmed generalized REM theta–gamma genotype deficit.

`Fig3S_warming_PAC_mouse_level`: Hilbert estimator from the prior intervention analysis, baseline-matched individual mouse PAC Z responses (only if local results exist). Warming PAC does not establish a genotype-specific pre-REM rescue; no new significance testing.

## Data provenance

All plots are derived directly from the CSVs under `outputs/sleep_architecture/warming_complete` and `outputs/eeg/pac`. The `outputs/` directory stays ignored in public Git. Prism CSVs are derived summaries, not raw EEG or individual animal identity beyond mouse IDs. Check institutional data-sharing policies before pushing.
