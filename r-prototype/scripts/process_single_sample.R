source("r-prototype/R/asari_process.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Missing mzML file path.\n",
    "Usage: Rscript r-prototype/scripts/process_single_sample.R \"path/to/file.mzML\" [output_dir] [max_scans]",
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
input_dir <- dirname(input_file)
output_dir <- if (length(args) >= 2L) args[[2]] else "r-output"
max_scans <- if (length(args) >= 3L) as.integer(args[[3]]) else NULL

if (!is.null(max_scans) && is.na(max_scans)) {
  stop("max_scans must be an integer when provided.", call. = FALSE)
}

cat("Input mzML:\n")
cat(input_file, "\n")
cat("Output root:\n")
cat(output_dir, "\n")
if (!is.null(max_scans)) {
  cat("Max scans:", max_scans, "\n")
}
cat("\n")

result <- asari_process(
  input_dir = input_dir,
  output_dir = output_dir,
  project_name = "single_sample",
  process_features = TRUE,
  selected_file = input_file,
  max_scans = max_scans
)

cat("Created output folder:\n")
cat(result$parameters$outdir, "\n")
cat("Feature rows:", nrow(result$feature_table), "\n")

