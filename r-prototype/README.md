# asari R prototype

This folder is the first R prototype for reproducing Python asari behavior.

The first milestone is intentionally small:

1. Find `.mzML` files in an input folder.
2. Register samples with stable sample IDs.
3. Create an asari-style output folder.
4. Write a skeleton `preferred_Feature_table.tsv`.
5. Write a `project.json` that records the run parameters.

The prototype does not parse mzML spectra yet. That will require an R mzML reader
such as `mzR`, `Spectra`, or `xcms`.
The current prototype uses `mzR` for the first mzML reader.

Run the smoke test from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
```

Preview mzML points from one test file:

```sh
Rscript r-prototype/scripts/preview_mzml.R
```
