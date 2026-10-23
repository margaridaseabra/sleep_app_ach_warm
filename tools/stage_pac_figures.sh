#!/usr/bin/env bash
# QC check and STAGE only source + selected derived figures/maps for GitHub.
# Never commits, never pushes, never stages continuous EEG or raw recordings.
set -euo pipefail
MODE="${1:-locked}"
[[ "$MODE" == locked || "$MODE" == fine ]] || { echo "Usage: bash tools/stage_pac_figures.sh [locked|fine]" >&2; exit 2; }
REPO="${PAC_REPO:-$HOME/sleep_app_ach_warm}"
OUT="$REPO/outputs/eeg/pac/real_comodulograms_v22_${MODE}"
[[ -d "$REPO/.git" ]] || { echo "Missing Git repository" >&2; exit 1; }
[[ -f "$OUT/roi_reference_comparison.csv" ]] || { echo "First run: bash tools/run_pac_figures.sh $MODE" >&2; exit 1; }
cd "$REPO"
[[ "$(git branch --show-current)" == paper/pac-* ]] || { echo "Refusing to stage PAC outside paper/pac-* branch" >&2; exit 1; }
if ! git diff --cached --quiet; then
   echo "Other files are already staged. Review or commit those separately first." >&2
   exit 1
fi
# Validate the scientific deliverables before copying them to version control.
python3 - "$OUT" "$MODE" <<'PY'
import csv, sys, math
from pathlib import Path
out, mode = Path(sys.argv[1]), sys.argv[2]
def rows(filename):
    p = out/filename
    if not p.is_file(): raise SystemExit(f'Missing {p}')
    with p.open(newline='') as f: return list(csv.DictReader(f))
qc = rows('pac_map_qc.csv')
roi = rows('roi_reference_comparison.csv')
mouse = rows('mouse_comodulograms_long.csv')
group = rows('group_comodulograms_long.csv')
if len(qc) != 16: raise SystemExit(f'QC: need 16 sessions, got {len(qc)}')
if len(roi) != 48: raise SystemExit(f'QC: need 48 mouse/context rows, got {len(roi)}')
if len({r['mouse_id'] for r in roi}) != 16: raise SystemExit('QC: expected 16 mice')
if {r['context'] for r in roi} != {'NREM_to_REM','NREM_to_WAKE','REM'}:
    raise SystemExit('QC: missing NREM or REM context')
if len({r['session_id'] for r in qc}) != 16: raise SystemExit('QC: sessions repeated')
if not all(int(float(r['n_cohorts'])) == 2 for r in group): raise SystemExit('QC: missing cohort in group maps')
if not all(math.isfinite(float(r['pac_z'])) for r in mouse): raise SystemExit('QC: nonfinite PAC in mouse maps')
if mode=='locked':
    d=max(abs(float(r['difference_z'])) for r in roi)
    if d>=0.001: raise SystemExit(f'QC: locked-grid reference failed, max discrepancy {d}')
    print(f'Locked-grid reference passed: max ROI difference {d:.6g}')
print('QC passed: 16 sessions, 16 mice, 48 mouse/context ROI rows; REM theta-gamma included')
PY
for ext in png pdf fig; do
 for name in nrem_group_comodulograms nrem_genotype_difference_comodulograms rem_theta_gamma_group_comodulograms mouse_level_roi_pac; do
  [[ -s "$OUT/figures/$name.$ext" ]] || { echo "Missing figure $name.$ext" >&2; exit 1; }
 done
done
MAP_DST="results/prism/pac/baseline/comodulograms_${MODE}"
FIG_DST="results/figures/pac/baseline/comodulograms_${MODE}"
mkdir -p "$MAP_DST" "$FIG_DST"
for name in group_comodulograms_long cohort_comodulograms_long mouse_comodulograms_long roi_reference_comparison pac_map_qc; do
  cp "$OUT/$name.csv" "$MAP_DST/"
done
for ext in png pdf fig; do
  for name in nrem_group_comodulograms nrem_genotype_difference_comodulograms rem_theta_gamma_group_comodulograms mouse_level_roi_pac; do
    cp "$OUT/figures/$name.$ext" "$FIG_DST/"
  done
done
# Standardize CSV line endings ONLY in copied Git files; never mutate outputs.
python3 - "$MAP_DST" <<'PY'
from pathlib import Path
import sys
for p in Path(sys.argv[1]).glob('*.csv'):
    b=p.read_bytes();p.write_bytes(b.replace(b'\r\n',b'\n'))
PY
git add -- \
  pipelines/run_22_pac_real_comodulograms.m \
  tools/plot_saved_pac_roi.py \
  tools/run_pac_figures.sh \
  tools/stage_pac_figures.sh \
  docs/pac_comodulogram_figures.md \
  "$MAP_DST" "$FIG_DST"
if ! git diff --cached --check; then
  echo "Staged whitespace check failed. Review before committing." >&2
  exit 1
fi
echo
echo '=== Files staged for review ==='
git diff --cached --stat
echo
echo 'No commit or push performed. Review with: git diff --cached --name-only'
echo 'Then: git commit -m "Add real-data NREM and REM PAC comodulograms" && git push -u origin HEAD'
