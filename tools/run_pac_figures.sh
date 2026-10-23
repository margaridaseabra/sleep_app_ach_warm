#!/usr/bin/env bash
# Generate REAL Morlet NREM and REM PAC maps from EEG, using MATLAB headless.
set -euo pipefail
MODE="${1:-locked}"
[[ "$MODE" == locked || "$MODE" == fine ]] || { echo "Usage: bash tools/run_pac_figures.sh [locked|fine]" >&2; exit 2; }
REPO="${PAC_REPO:-$HOME/sleep_app_ach_warm}"
[[ -f "$REPO/pipelines/run_22_pac_real_comodulograms.m" ]] || { echo "Install package first." >&2; exit 1; }
[[ -f "$REPO/outputs/eeg/pac/baseline_morlet_bandwidth_v20/paired_pac_windows.csv" ]] || { echo "Missing v2.0 Morlet baseline output." >&2; exit 1; }
[[ -f "$REPO/outputs/eeg/spindles/spindle_transition_bouts.csv" ]] || { echo "Missing frozen spindle transitions." >&2; exit 1; }
if [[ -n "${MATLAB_BIN:-}" ]]; then
  M="$MATLAB_BIN"
elif command -v matlab >/dev/null 2>&1; then
  M="$(command -v matlab)"
elif [[ -x /Applications/MATLAB_R2026a.app/bin/matlab ]]; then
  M=/Applications/MATLAB_R2026a.app/bin/matlab
else
  echo "MATLAB not found. Set MATLAB_BIN=/Applications/MATLAB_R2026a.app/bin/matlab" >&2
  exit 1
fi
LOGS="$REPO/outputs/eeg/pac/map_generation_logs"
mkdir -p "$LOGS"
LOG="$LOGS/real_comodulograms_${MODE}.log"
# Paths are local filesystem paths; MATLAB -batch makes this fully Terminal-run.
echo "Running real EEG PAC maps: $MODE"
echo "MATLAB: $M"
echo "Log: $LOG"
"$M" -batch "cd('$REPO'); setenv('PAC_MAP_GRID','$MODE'); run('pipelines/run_22_pac_real_comodulograms.m');" 2>&1 | tee "$LOG"
echo "Figures: $REPO/outputs/eeg/pac/real_comodulograms_v22_${MODE}/figures/"
