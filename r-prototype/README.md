# asari R prototype

This folder is the first R prototype for reproducing Python asari behavior.

The current prototype covers the first single-sample processing path:

1. Find `.mzML` files in an input folder.
2. Register samples with stable sample IDs.
3. Create an asari-style output folder.
4. Extract mass tracks from MS1 data.
5. Build a single-sample mass grid and composite tracks.
6. Detect global elution peaks.
7. Write `full_Feature_table.tsv`, `preferred_Feature_table.tsv`, and
   `project.json`.

Default processing parameters live in `R/parameters.R`, translated from Python
asari's `default_parameters.py`.

The prototype uses `mzR` for the first mzML reader and can preview MS1 points
from real `.mzML` files.

Run the smoke test from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
```

Run a single-sample prototype that writes a non-empty feature table. You must
explicitly provide the mzML file:

```sh
Rscript r-prototype/scripts/process_single_sample.R "实验数据/batch10_MT_20210804_001A.mzML" "r-output" 20
```

The optional third argument limits the number of scans for quick testing. Omit
it only when you are ready to process the full file.

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
split wide mz clusters by mass seeds, merge very close tracks, then build an
intensity vector over the full retention-time range.

Preview the first peak detection prototype:

```sh
Rscript r-prototype/scripts/preview_peaks.R "实验数据/batch10_MT_20210804_001A.mzML"
```

The peak detection prototype follows the broad structure of `asari/peaks.py`:
audit each mass track, split signal into ROIs, detect local apexes with
prominence and width constraints, evaluate Gaussian shape, compute SNR and
cSelectivity, and emit the fixed feature columns used by asari feature tables.

Compare R checkpoints with Python asari output:

```sh
Rscript r-prototype/scripts/compare_with_asari.R "实验数据/MT02Dataset/batch14_MT_20210808_005.mzML"
```

For a quick partial check:

```sh
Rscript r-prototype/scripts/compare_with_asari.R "实验数据/MT02Dataset/batch14_MT_20210808_005.mzML" --max-scans 50
```

If you already ran Python asari and have its output folder, pass it in:

```sh
Rscript r-prototype/scripts/compare_with_asari.R "实验数据/MT02Dataset/batch14_MT_20210808_005.mzML" --python-output "path/to/python_asari_output"
```
