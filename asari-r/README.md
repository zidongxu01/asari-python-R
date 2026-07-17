# asariR

This is a study-oriented R reimplementation of selected core parts of asari.

The goal is not to translate every Python line directly. The goal is to keep a
clear file-by-file correspondence for the core algorithm modules, so each R file
can be compared with the original asari Python file while still using natural R
structure.

## Current Study Scope

Main focus:

1. Construct mass tracks / EICs from MS1 data.
2. Understand m/z tolerance, binning, clustering, and merging.
3. Detect chromatographic peaks and filter noise.
4. Compare behavior with the original asari implementation.

Not the first priority:

- dashboard
- command line interface
- annotation
- GC-MS-specific workflow

## File Map

| Python asari file | R file | Purpose |
| --- | --- | --- |
| `chromatograms.py` | `R/chromatograms.R` | Mass track / EIC construction and RT alignment |
| `mass_functions.py` | `R/mass_functions.R` | m/z distance, clustering, seeds, close m/z helpers |
| `peaks.py` | `R/peaks.R` | Peak detection, noise estimation, peak evaluation |
| `default_parameters.py` | `R/default_parameters.R` | Default parameters |
| `constructors.py` | `R/constructors.R` | Higher-level data structures and mass grid construction |
| `workflow.py` | `R/workflow.R` | Pipeline orchestration |
| `samples.py` / `experiment.py` | `R/samples.R`, `R/experiment.R` | Sample and experiment-level organization |
| `utils.py` | `R/utils.R` | Shared helper functions |

## Suggested Order

1. Connect the completed `R/mass_functions.R` behavior to constructors and the MassGrid.
2. Complete `R/peaks.R`, especially real peak vs noise decisions.
3. Add module-level and end-to-end comparisons after each integration step.
