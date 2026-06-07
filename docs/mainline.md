# Python asari to R mainline

This document keeps the conversion focused on the main processing path instead
of translating files in random order.

## Main Path

```text
command_line.py
  -> workflow.py
  -> samples.py / experiment.py
  -> constructors.py
  -> chromatograms.py / mass_functions.py
  -> peaks.py
  -> preferred_Feature_table.tsv
```

For the R prototype, the current translation path is:

```text
r-prototype/R/asari_process.R
  -> r-prototype/R/mzml_reader.R
  -> r-prototype/R/mass_tracks.R
  -> r-prototype/R/peak_detection.R
  -> preferred_Feature_table.tsv-style columns
```

## Python Modules and R Status

| Python module | Role in Python asari | R prototype status |
| --- | --- | --- |
| `command_line.py` | Parses user commands and parameters. | Partly represented by `asari_process.R`; no full CLI yet. |
| `workflow.py` | Orchestrates project processing, output folders, sample registry, mass track extraction, experiment processing, export. | Partly represented by `asari_process.R`; only project shell and output skeleton are implemented. |
| `samples.py` | Represents individual sample data. | Not ported yet. Current sample registry is a simple data frame. |
| `experiment.py` | Holds whole-project experiment state and coordinates LC/GC/DIMS processing. | Not ported yet. |
| `constructors.py` | Builds `MassGrid`, `CompositeMap`, aligned feature maps, and final feature structures. | Not ported yet. This is the next major missing mainline section. |
| `chromatograms.py` | Reads MS1 data and builds mass tracks / EICs. | Initial prototype in `mass_tracks.R`, following the broad binning and track-building structure. |
| `mass_functions.py` | Provides ppm matching, clustering, and mass alignment helpers. | Only minimal ppm helper logic exists in `mass_tracks.R`; full matching/clustering is not ported. |
| `peaks.py` | Detects and evaluates elution peaks on mass tracks. | Initial prototype in `peak_detection.R`; Gaussian fitting and full prominence logic are not ported. |
| output export | Writes feature tables and project metadata. | Skeleton output exists; detected features are not yet wired into `asari_process()`. |

## Current Implemented Chain

The R prototype currently proves this chain:

```text
.mzML
  -> MS1 points: scan_id, rtime, mz, intensity
  -> mass tracks: track_id, mz, scan_id, rtime, intensity
  -> detected feature columns:
       id_number, mz, rtime, rtime_left_base, rtime_right_base,
       parent_masstrack_id, peak_area, cSelectivity,
       goodness_fitting, snr, detection_counts
```

## Important Gaps

The next work should stay on the mainline and close these gaps in order:

1. Wire `mass_tracks.R` and `peak_detection.R` into `asari_process()`.
2. Add a single-sample output mode that writes a real non-empty feature table.
3. Port the relevant `mass_functions.py` matching helpers.
4. Start `constructors.py` equivalents for multi-sample mass grid alignment.
5. Compare R output against Python `preferred_Feature_table.tsv` on the same input data.

## Rule

Do not translate side modules such as dashboard, annotation, notebooks, or GC
annotation until the main path above can produce a comparable feature table.

