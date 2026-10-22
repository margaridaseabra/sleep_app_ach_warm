# Manipulation timing registries

The experimental manipulation times are version-controlled in `config/`.

## Ambient warming

`config/warming_times.csv` contains actual heating onset/end times for the
October and June batches.

Columns:

- `session_id`
- `mouse_id`
- `cohort`
- `genotype`
- `batch`
- `warming_start_elapsed_s`
- `warming_start_elapsed_hms`
- `warming_end_elapsed_s`
- `warming_end_elapsed_hms`
- `timing_source`
- `source_file`
- `precision`
- `include_analysis`
- `notes`

October values come from `ambtemp-times-oct.xlsx`. June values come from the
handwritten notebook photograph `IMG_0416.jpg`.

The row for historical mouse6 from the October workbook is retained for
provenance but has `include_analysis=0` because that mouse is not present in the
current session manifest.

## Donepezil

`config/drug_times.csv` contains the manipulation times for the donepezil
sessions.

June notebook entries contain both saline and drug times. The October workbook
contains only the drug timing, so October saline columns are intentionally
blank.

Columns:

- `session_id`
- `mouse_id`
- `cohort`
- `genotype`
- `batch`
- `saline_elapsed_s`
- `saline_elapsed_hms`
- `drug_elapsed_s`
- `drug_elapsed_hms`
- `timing_source`
- `source_file`
- `precision`
- `include_analysis`
- `notes`

No missing intervention times should be inferred silently. If a time is later
recovered from another source, update the corresponding row and document the
source.
