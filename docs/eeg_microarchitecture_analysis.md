# EEG microarchitecture analysis: spindle and PAC specification

This document locks the analysis choices **before** inspecting the expanded
cohort results.

## 1. Spindles

Reference method: Bandarabadi et al., Nature Communications (2020),
"A role for spindles in the onset of rapid eye movement sleep."

### Keep from the thesis

- EEG resampled to 200 Hz.
- 50-Hz notch filtering.
- Complex frequency B-spline wavelet detector.
- Spindle band 9–16 Hz.
- Wavelet-energy trace smoothed with a 200-ms Hann window.
- High threshold = mean + 3 SD.
- Low threshold = mean + 1 SD for onset/offset refinement.
- Duration 0.4–2.0 s.
- 5–30 oscillatory cycles.
- Reject events whose spindle-band power does not exceed both flanking
  6–8.5-Hz and 16.5–20-Hz bands.
- Event features: duration, center frequency, cycles, symmetry and
  peak-to-peak amplitude.
- Spindle density normalized to NREM duration.

### Changes for the paper

1. **Use the full 9–16-Hz wavelet-energy band.**
   The thesis text described 10–14-Hz wavelet centers even though the
   reference detector defines spindle energy over 9–16 Hz.

2. **Do not multiply wavelet power by frequency (no additional 1/f weighting)
   in the primary detector.**
   This is not part of the Bandarabadi detection method. If desired it can be
   retained as a sensitivity analysis, not as the primary detector.

3. **Primary transition endpoint: spindle rate in the final 25 s before
   NREM→REM**, with NREM→Wake as the comparator.
   The reference paper found the maximal pre-REM increase at approximately
   25 s. A 60-s aligned time course can still be shown descriptively.

4. **Mouse is the unit of inference.**
   Transition/bout values are summarized within mouse before genotype or
   condition statistics. Do not test every 1-s time point as independent.

5. **Blinded detector QC.**
   Visually inspect a pre-specified random subset of detections and misses
   across genotypes/cohorts/conditions. Record sensitivity/FDR-style QC rather
   than tuning thresholds after seeing group outcomes.

6. **Sex/cohort sensitivity.**
   The new WT cohort is all female while new APP contains both sexes, and the
   reference spindle paper notes possible sex effects. Main inference should
   not pretend genotype and sex can be separated in the new cohort. Report
   original/new separately and use female-only new-cohort plots as descriptive
   sensitivity if useful.

## 2. PAC

Reference method: Munia & Aviyente, Scientific Reports (2019),
"Time-Frequency Based Phase-Amplitude Coupling Measure For Neuronal
Oscillations."

### Keep from the thesis

- RID-Rihaczek complex time-frequency representation.
- tf-MVL as the coupling metric.
- NREM slow-phase / sigma-amplitude analysis.
- REM theta-phase / gamma-amplitude analysis.
- Surrogate-normalized maps and per-mouse ROI summaries.
- NREM microarousals assimilated into NREM macro-bouts.

Primary a-priori ROIs remain:

- NREM: phase 1–4 Hz, amplitude 9–16 Hz.
- REM: phase 4–10 Hz, amplitude 30–80 Hz.

### Changes for the paper

1. **Use the surrogate method exactly as described in the PAC paper.**
   Split the amplitude envelope at a random time point and swap the two
   blocks (equivalent to a random circular shift), repeat 200 times, and use
   the 95th percentile. Do not use permutation of many fixed 0.5-s blocks as
   the primary surrogate.

2. **Use surrogate-normalized ROI mean (Z) as the primary scalar.**
   Thresholded hotspot maps are visualization. Statistical inference should
   not be based on how visually large a hotspot looks.

3. **No random window subsampling as the primary analysis.**
   Use deterministic windows and hierarchical averaging:
   windows → bout → mouse. This gives equal mouse-level inference while
   preventing animals with more sleep from dominating.

4. **Avoid the thesis's ≥180-s NREM-bout requirement as the only analysis.**
   Requiring very long NREM bouts can preferentially select the least
   fragmented APP bouts. For transition PAC, require only enough continuous
   NREM-like data for the pre-specified transition window and report the
   number of contributing bouts.

5. **Longer NREM windows.**
   Use 10-s windows as the primary NREM PAC estimate (with a 5-s-window
   sensitivity analysis if direct thesis comparability is needed). This
   provides more cycles for the slowest 1-Hz phase component.

6. **Pre-specify transition context.**
   Primary NREM PAC:
   - final 60 s before NREM→REM
   - final 60 s before NREM→Wake
   summarized per bout and mouse.
   The contrast `preREM - preWake` is secondary but pre-specified.

7. **PAC confound QC.**
   Check state-resolved PSD/oscillatory power and reject obvious EEG/EMG
   artifacts. Interpret theta–gamma PAC cautiously because broadband,
   nonsinusoidal or muscle-related activity can produce apparent coupling.

## 3. Condition comparisons

For the paper:

- Figure 1: baseline phenotype first.
- Figure 3: baseline-matched vs active ambient-warming window.
- Donepezil: mechanistic/pharmacological comparator, preferably using the
  actual 0–3 h post-injection window.

Do not choose parameters or exclude mice based on whether the expanded cohort
reproduces the thesis result.
