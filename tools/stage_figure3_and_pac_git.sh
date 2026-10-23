#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
# No raw EEG / full generated-output directories. Review public data before committing.
git add tools/build_figure3_and_pac.py tools/run_figure3_and_pac.sh tools/stage_figure3_and_pac_git.sh docs/figure3_warming_and_pac.md
for p in results/figures/paper/figure3_warming results/prism/paper/figure3_warming results/figures/paper/figure1_baseline/pac_extended results/prism/paper/figure1_baseline/pac_extended; do
    [ ! -d "$p" ] || git add "$p"
done
printf '\nFiles staged:\n';git diff --cached --stat
printf '\nReview with git diff --cached and git status before commit/push.\n'
