# asari R prototype

This folder is the first R prototype for reproducing Python asari behavior.

The first milestone is intentionally small:

1. Find `.mzML` files in an input folder.
2. Register samples with stable sample IDs.
3. Create an asari-style output folder.
4. Write a skeleton `preferred_Feature_table.tsv`.
5. Write a `project.json` that records the run parameters.

The prototype uses `mzR` for the first mzML reader and can preview MS1 points
from real `.mzML` files.

Run the smoke test from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
```

Preview mzML points from one test file:

```sh
Rscript r-prototype/scripts/preview_mzml.R
```

Preview the first mass track prototype:

```sh
Rscript r-prototype/scripts/preview_mass_tracks.R
```

The mass track prototype follows the Python asari structure in
`asari/chromatograms.py`: filter MS1 points, group by thousandth-m/z bins,
screen by consecutive scans and peak height, then build an intensity vector over
the full retention-time range. It is not a complete line-by-line port yet.
