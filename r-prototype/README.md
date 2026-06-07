# asari R prototype

This folder is the first R prototype for reproducing Python asari behavior.

The first milestone is intentionally small:

1. Find `.mzML` files in an input folder.
2. Register samples with stable sample IDs.
3. Create an asari-style output folder.
4. Write a skeleton `preferred_Feature_table.tsv`.
5. Write a `project.json` that records the run parameters.

Default processing parameters live in `R/parameters.R`, translated from Python
asari's `default_parameters.py`.

The prototype uses `mzR` for the first mzML reader and can preview MS1 points
from real `.mzML` files.

Run the smoke test from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
```

Preview mzML points from a specific test file:

```sh
Rscript r-prototype/scripts/preview_mzml.R "实验数据/batch10_MT_20210804_001A.mzML"
```

Preview the first mass track prototype:

```sh
Rscript r-prototype/scripts/preview_mass_tracks.R "实验数据/batch10_MT_20210804_001A.mzML"
```

The mass track prototype follows the Python asari structure in
`asari/chromatograms.py`: filter MS1 points, group by thousandth-m/z bins,
screen by consecutive scans and peak height, then build an intensity vector over
the full retention-time range. It is not a complete line-by-line port yet.

Preview the first peak detection prototype:

```sh
Rscript r-prototype/scripts/preview_peaks.R "实验数据/batch10_MT_20210804_001A.mzML"
```

The peak detection prototype follows the broad structure of `asari/peaks.py`:
estimate track noise, split signal into ROIs, detect local apexes, estimate left
and right bases, and emit the fixed feature columns used by
`preferred_Feature_table.tsv`. Gaussian fitting and full SciPy-style prominence
logic are not ported yet.
