# PAC comodulograms — real EEG data, v2.3 terminal workflow

This module reproduces full frequency-resolved Morlet comodulograms from the 16 baseline recordings. It uses the 6,696 Morlet-v2.0 accepted windows, the already-frozen spindle transitions and bandwidth/FWHM calibration, with 200 circular-shift surrogates and 200 Hz EEG preprocessing. **No synthetic heatmaps are used for paper figures.**

- **NREM**, before NREM→REM or NREM→Wake: phase 1–4 Hz; sigma amplitude 9–16 Hz.
- **REM**, during REM bouts: theta phase 4–10 Hz; gamma amplitude 30–80 Hz. This is not NREM slow–sigma coupling.
- `locked` is the original frequency grid and checks each mouse/context ROI against saved v2.0 output (1e-3 Z tolerance).
- `fine` calculates additional actual EEG frequencies, NOT image interpolation. Fine-grid results are exploratory until reviewed; run locked first.
- PAC Z is mean vector length minus mean surrogate divided by surrogate SD at each frequency cell.
- Aggregation: window → bout → mouse → genotype within cohort → equal-weight October/June group average.
- Blue/white/red colormap uses zero as white; negative PAC Z can mean below surrogate expectation; **not** negative physiological coupling.
- APP−WT and preREM−preWake maps are **differences of mean Z**, not statistical Z maps or pixelwise significance tests.
- Group heatmaps are descriptive; inferential tests use mouse-level ROI values; baseline PAC results did not show a corrected genotype deficit.
- High REM gamma near 80 Hz requires artifact/harmonic checks. Do not imply a REM PAC deficit based on heatmap color.
- This package creates **baseline** maps only. A real donepezil or warming **Morlet** map needs the corresponding Morlet intervention EEG processing; never render an illustrative treatment map as measured.

## Terminal workflow

```bash
bash ~/Downloads/pac_comodulograms_terminal_v23/tools/install_pac_figures.sh
cd ~/sleep_app_ach_warm
bash tools/run_pac_figures.sh locked
bash tools/run_pac_figures.sh fine  # optional, if locked QC passes
bash tools/stage_pac_figures.sh locked
# Inspect staged files, particularly data; repo is public!
git diff --cached --name-only
git diff --cached --check
git commit -m "Add real-data NREM and REM PAC comodulograms"
git push -u origin HEAD
```

The run generates PNG/PDF/editable `.fig` figures under ignored `outputs/eeg/pac/real_comodulograms_v22_{locked|fine}/figures/`. The stage script copies only selected derived maps and figure files to `results/prism/pac/...` and `results/figures/pac/...`. Raw EEG, recorded scores and bulky per-bout analysis outputs are not staged.

Original scalar point plot (no EEG recomputation):

```bash
conda run -n pac_review python tools/plot_saved_pac_roi.py
```

## Publication caution

The code was inspected, but MATLAB execution must verify that v2.0 ROI reproduction passes. If locked-grid checks fail, **do not stage or commit figures**. Preserve the console message and investigate signal preprocessing/phase support. Keep all spindle scripts and prior PAC estimators unchanged.
