source("r-prototype/R/asari_process.R")

tmp <- tempfile("asari-r-prototype-")
result <- asari_process(
  input_dir = "实验数据",
  output_dir = tmp,
  project_name = "smoke"
)

expected_files <- c(
  file.path(result$parameters$outdir, "preferred_Feature_table.tsv"),
  file.path(result$parameters$outdir, "project.json"),
  file.path(result$parameters$export_outdir, "full_Feature_table.tsv")
)

missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0L) {
  stop("Missing expected files: ", paste(missing_files, collapse = ", "))
}

expected_columns <- c(
  "id_number",
  "mz",
  "rtime",
  "rtime_left_base",
  "rtime_right_base",
  "parent_masstrack_id",
  "peak_area",
  "cSelectivity",
  "goodness_fitting",
  "snr",
  "detection_counts",
  result$samples$sample_name
)

actual_columns <- names(utils::read.delim(expected_files[1], nrows = 0, check.names = FALSE))
if (!identical(actual_columns, expected_columns)) {
  stop("Feature table columns do not match expected asari-style columns.")
}

cat("Smoke test passed.\n")

