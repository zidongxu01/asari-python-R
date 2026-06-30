source("r-prototype/R/asari_process.R")

input_file <- "实验数据/batch10_MT_20210804_001A.mzML"
tmp <- tempfile("asari-r-single-sample-")

result <- asari_process(
  input_dir = "实验数据",
  output_dir = tmp,
  project_name = "single",
  process_features = TRUE,
  selected_file = input_file,
  max_scans = 20
)

feature_path <- file.path(result$parameters$outdir, "preferred_Feature_table.tsv")
if (!file.exists(feature_path)) {
  stop("Expected preferred feature table to be written.")
}

feature_table <- utils::read.delim(feature_path, check.names = FALSE)
full_feature_path <- file.path(result$parameters$export_outdir, "full_Feature_table.tsv")
full_feature_table <- utils::read.delim(full_feature_path, check.names = FALSE)
if (nrow(full_feature_table) == 0L) {
  stop("Expected single-sample processing to write a non-empty full feature table.")
}

if (nrow(full_feature_table) < nrow(feature_table)) {
  stop("Expected full feature table to contain at least as many rows as preferred feature table.")
}

sample_name <- tools::file_path_sans_ext(basename(input_file))
if (!sample_name %in% names(feature_table)) {
  stop("Expected single-sample intensity column in feature table.")
}

if (any(full_feature_table[[sample_name]] < 0)) {
  stop("Expected single-sample peak areas to be non-negative.")
}

if (nrow(result$samples) != 1L || result$samples$sample_name[[1]] != sample_name) {
  stop("Expected single-sample processing to register only the selected mzML file.")
}

cat("single sample process test passed.\n")
