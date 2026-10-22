#!/usr/bin/env python3
"""
Populate recording_file for June/new-cohort rows in config/sessions.csv.

This is intentionally separate from the session-manifest builder: it only
fills raw EEG paths and does not modify session selection or score-file
choices.

Run from anywhere:
    python tools/update_new_recording_paths.py
"""

from pathlib import Path
import csv
import os
import re
import shutil

repo = Path(__file__).resolve().parents[1]

default_data_root = Path.home() / "Desktop" / "AmbTemp_proj"
data_root = Path(
    os.environ.get("AMBTEMP_DATA_ROOT", str(default_data_root))
).expanduser()

sessions_file = repo / "config" / "sessions.csv"
backup_file = repo / "config" / "sessions_before_recording_path_update.csv"

OVERRIDES = {
    "20260605_mouse63_warming_2":
        "scored_files_june26batch/050626-m63-ammbt-2.mat",
}

if not sessions_file.is_file():
    raise FileNotFoundError(sessions_file)

shutil.copy2(sessions_file, backup_file)

with sessions_file.open(newline="") as f:
    rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys())

found = 0
missing = []

for row in rows:

    if row.get("cohort") != "new":
        continue

    sid = row["session_id"]

    if sid in OVERRIDES:
        relative = Path(OVERRIDES[sid])
    else:
        score_rel = Path(row["score_file"])

        raw_name = re.sub(
            r"_final_scoring_dash\.csv$",
            "",
            score_rel.name,
            flags=re.IGNORECASE,
        )

        if not raw_name.lower().endswith(".mat"):
            raw_name += ".mat"

        relative = score_rel.parent / raw_name

    full = data_root / relative

    if full.is_file():
        row["recording_file"] = str(relative)
        found += 1
    else:
        missing.append(str(relative))

with sessions_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Resolved new-cohort raw recordings: {found}/24")

if missing:
    print("\nMissing:")
    for x in missing:
        print(" ", x)
else:
    print("All 24 new-cohort raw recordings found.")

print(f"\nBackup:\n{backup_file}")
