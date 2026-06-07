source("r-prototype/R/parameters.R")

parameters <- default_parameters()

required_names <- c(
  "project_name",
  "outdir",
  "database_mode",
  "mode",
  "mass_range",
  "mz_tolerance_ppm",
  "correction_tolerance_ppm",
  "min_timepoints",
  "signal_noise_ratio",
  "min_intensity_threshold",
  "min_peak_height",
  "min_peak_ratio",
  "wlen",
  "gaussian_shape",
  "peak_area",
  "workflow",
  "rt_align_method",
  "rtime_tolerance",
  "output_feature_table",
  "mass_grid_mapping"
)

missing_names <- setdiff(required_names, names(parameters))
if (length(missing_names) > 0L) {
  stop("Missing expected parameters: ", paste(missing_names, collapse = ", "))
}

if (!identical(parameters$mz_tolerance_ppm, 5)) {
  stop("Expected mz_tolerance_ppm to match Python asari default.")
}

if (!identical(parameters$min_peak_height, 100000)) {
  stop("Expected min_peak_height to match Python asari default.")
}

constants <- asari_constants()
if (!identical(constants$proton, 1.00727646677)) {
  stop("Expected proton constant to match Python asari default.")
}

cat("parameters test passed.\n")

