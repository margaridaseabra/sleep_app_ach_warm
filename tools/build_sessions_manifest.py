from pathlib import Path
import os
from datetime import datetime
import csv
import re

repo = Path.home() / "sleep_app_ach_warm"
default_data_root = Path.home() / "Desktop" / "AmbTemp_proj"
data_root = Path(
    os.environ.get("AMBTEMP_DATA_ROOT", str(default_data_root))
).expanduser()

old_dir = data_root / "scored_files_oct25batch"
new_dir = data_root / "scored_files_june26batch"

out_file = repo / "config" / "sessions.csv"

condition_map = {
    "baseline": "baseline",
    "bas": "baseline",
    "ambtemp": "warming",
    "ambt": "warming",
    "warm": "warming",
    "warming": "warming",
    "drug": "donepezil",
    "drugs": "donepezil",
    "donepezil": "donepezil",
}

rows = []
unparsed = []

def relative_path(p):
    return str(p.relative_to(data_root))

def add_old_file(p):
    name = p.name

    patterns = [
        re.compile(
            r"^(?P<date>\d{8})-"
            r"(?P<condition>baseline|ambtemp|ambt|drug|drugs)-"
            r"(?P<mouse>mouse\d+)-scored\.mat$",
            re.IGNORECASE,
        ),
        re.compile(
            r"^(?P<date>\d{8})-"
            r"(?P<mouse>mouse\d+)-"
            r"(?P<condition>baseline|ambtemp|ambt|drug|drugs)-"
            r"scored\.mat$",
            re.IGNORECASE,
        ),
    ]

    match = None
    for pattern in patterns:
        match = pattern.match(name)
        if match:
            break

    if not match:
        unparsed.append(str(p))
        return

    date = match.group("date")
    mouse = match.group("mouse").lower()
    raw_condition = match.group("condition").lower()
    condition = condition_map[raw_condition]

    session_id = f"{date}_{mouse}_{condition}"

    rows.append({
        "session_id": session_id,
        "mouse_id": mouse,
        "date": date,
        "condition": condition,
        "condition_raw": raw_condition,
        "cohort": "thesis",
        "session": "",
        "score_format": "legacy_mat",
        "score_file": relative_path(p),
        "recording_file": relative_path(p),
        "include_sleep": 1,
        "notes": "",
    })


def add_new_file(p):
    name = p.name

    pattern = re.compile(
        r"^(?P<date>\d{6})-"
        r"m(?P<mouse>\d+)-"
        r"(?P<condition>bas|baseline|ambt|ambtemp|drug|drugs)-"
        r"(?P<session>\d+)"
        r"(?:\.mat)?"
        r"_final_scoring_dash\.csv$",
        re.IGNORECASE,
    )

    match = pattern.match(name)

    if not match:
        unparsed.append(str(p))
        return

    raw_date = match.group("date")
    date = datetime.strptime(raw_date, "%d%m%y").strftime("%Y%m%d")

    mouse = "mouse" + match.group("mouse")
    raw_condition = match.group("condition").lower()
    condition = condition_map[raw_condition]
    session = match.group("session")

    session_id = f"{date}_{mouse}_{condition}_{session}"

    rows.append({
        "session_id": session_id,
        "mouse_id": mouse,
        "date": date,
        "condition": condition,
        "condition_raw": raw_condition,
        "cohort": "new",
        "session": session,
        "score_format": "dash_csv",
        "score_file": relative_path(p),
        "recording_file": "",
        "include_sleep": 1,
        "notes": "",
    })


for p in sorted(old_dir.glob("*.mat")):
    add_old_file(p)

for p in sorted(new_dir.glob("*_final_scoring_dash.csv")):
    add_new_file(p)


rows.sort(
    key=lambda r: (
        r["cohort"],
        r["mouse_id"],
        r["date"],
        r["condition"],
    )
)

fieldnames = [
    "session_id",
    "mouse_id",
    "date",
    "condition",
    "condition_raw",
    "cohort",
    "session",
    "score_format",
    "score_file",
    "recording_file",
    "include_sleep",
    "notes",
]

with out_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print()
print("========================================")
print("SESSION MANIFEST")
print("========================================")
print(f"Written: {out_file}")
print(f"Sessions: {len(rows)}")

print()
print("Counts by cohort / condition:")

for cohort in ["thesis", "new"]:
    for condition in ["baseline", "warming", "donepezil"]:
        n = sum(
            r["cohort"] == cohort and
            r["condition"] == condition
            for r in rows
        )
        print(f"{cohort:8s} {condition:10s}: {n}")

if unparsed:
    print()
    print("WARNING - files that could not be parsed:")
    for x in unparsed:
        print(x)
else:
    print()
    print("All scoring files parsed successfully.")
