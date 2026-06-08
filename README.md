# asari Python-to-R

This repository is a working R port of
[asari](https://github.com/shuzhao-li/asari), a Python program for
high-resolution LC-MS metabolomics data processing.

The goal is to reproduce asari's scientific processing mainline in R without
mechanically copying the Python source tree. The current work is organized
around a focused prototype in [`r-prototype/`](r-prototype/) and a mainline map
in [`docs/mainline.md`](docs/mainline.md).

## Scope

The R implementation is being built around the core asari processing path:

1. Read explicit mzML inputs.
2. Register samples and run parameters.
3. Extract mass tracks / extracted ion chromatograms from MS1 scans.
4. Detect chromatographic peaks on mass tracks.
5. Write asari-style feature table output.
6. Extend from single-sample output to multi-sample alignment.

This repository is not trying to clone the Python package layout file by file.
When an R-native structure is clearer, the R port should use it while keeping
the scientific behavior aligned with Python asari.

## Current Status

The active code lives in [`r-prototype/`](r-prototype/).

| Area | Status |
| --- | --- |
| mzML reading | Initial reader implemented with `mzR`. |
| Parameter defaults | Prototype defaults translated from Python asari. |
| Mass tracks / EICs | Rough single-sample prototype implemented. |
| Peak detection | Rough prototype emits feature-table-style columns. |
| Single-sample process | Can run from an explicit mzML path and write a non-empty feature table. |
| Multi-sample alignment | Not implemented yet. |
| Full Python parity | Not claimed yet; comparison against Python output is the next checkpoint. |

For the detailed scientific route and missing pieces, see
[`docs/mainline.md`](docs/mainline.md).

## Quick Start

Run the prototype tests from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
Rscript r-prototype/tests/mzml_reader_test.R
Rscript r-prototype/tests/mass_tracks_test.R
Rscript r-prototype/tests/peak_detection_test.R
Rscript r-prototype/tests/single_sample_process_test.R
```

Run a single-sample prototype with an explicit mzML file:

```sh
Rscript r-prototype/scripts/process_single_sample.R "path/to/file.mzML" "r-output" 20
```

The optional third argument limits the number of scans for quick testing. Omit
it only when you are ready to process the full file.

Preview individual stages with explicit mzML input:

```sh
Rscript r-prototype/scripts/preview_mzml.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_mass_tracks.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_peaks.R "path/to/file.mzML"
```

No preview or processing script should rely on a hidden default mzML file.

## Repository Layout

```text
docs/mainline.md      Scientific mainline and porting checkpoints
r-prototype/R/        Current R prototype implementation
r-prototype/scripts/  Explicit-input preview and processing scripts
r-prototype/tests/    Prototype tests
```

## Development Priorities

1. Compare the current single-sample R output against Python asari's
   `preferred_Feature_table.tsv`.
2. Port the core mass matching and clustering helpers needed for alignment.
3. Implement multi-sample `MassGrid` / `CompositeMap` behavior in an R-native
   shape.

Until the mainline can produce comparable feature tables, dashboard,
annotation, notebooks, and GC annotation work should stay out of scope.
