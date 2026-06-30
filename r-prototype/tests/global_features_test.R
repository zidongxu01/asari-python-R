source("r-prototype/R/global_features.R")
source("r-prototype/R/parameters.R")

make_track <- function(mz, id_number, intensity) {
  list(
    id_number = id_number,
    mz = mz,
    intensity = intensity
  )
}

parameters <- default_parameters()
parameters$mz_tolerance_ppm <- 5
parameters$project_sample_number_small <- 10
parameters$min_peak_height <- 10000
parameters$min_intensity_threshold <- 1000
parameters$min_timepoints <- 3
parameters$signal_noise_ratio <- 2
parameters$min_prominence_threshold <- 5000
parameters$gaussian_shape <- 0.5
parameters$wlen <- 5

sample_mass_tracks <- list(
  sample_a = list(
    make_track(100.00000, 0L, c(0, 1000, 5000, 12000, 25000, 12000, 5000, 1000, 0)),
    make_track(200.00000, 1L, c(0, 0, 0, 2000, 3000, 2000, 0, 0, 0))
  ),
  sample_b = list(
    make_track(100.00030, 0L, c(0, 1000, 4000, 10000, 20000, 10000, 4000, 1000, 0)),
    make_track(300.00000, 1L, c(0, 0, 0, 2000, 3000, 2000, 0, 0, 0))
  )
)

rt_table <- data.frame(
  scan_id = seq_len(9),
  scan_index = seq_len(9),
  rtime = seq(10, 90, by = 10)
)

result <- process_global_features(
  sample_mass_tracks,
  parameters = parameters,
  rt_table = rt_table
)

if (nrow(result$mass_grid) != 3L) {
  stop("Expected MassGrid to align one shared track and keep two unmatched tracks.")
}

if (!identical(result$composite_tracks[[1]]$intensity, c(0, 2000, 9000, 22000, 45000, 22000, 9000, 2000, 0))) {
  stop("Expected first composite track to sum shared sample tracks.")
}

if (nrow(result$feature_table) != 1L) {
  stop("Expected one global feature from the shared composite track.")
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
  "sample_a",
  "sample_b"
)

if (!identical(names(result$feature_table), expected_columns)) {
  stop("Expected feature table to include per-sample peak area columns.")
}

if (!identical(names(result$full_feature_table), expected_columns)) {
  stop("Expected full feature table to use asari export columns.")
}

if (nrow(result$preferred_feature_table) != 1L) {
  stop("Expected preferred feature table to keep the high-quality global feature.")
}

if (!identical(result$feature_table$id_number, "F0")) {
  stop("Expected first global feature id to be F0.")
}

if (result$feature_table$parent_masstrack_id[[1]] != 0L) {
  stop("Expected global feature to point to composite mass track 0.")
}

if (result$feature_table$rtime[[1]] != 50) {
  stop("Expected global feature apex retention time to be 50.")
}

if (result$feature_table$peak_area[[1]] <= 0) {
  stop("Expected positive global peak area.")
}

if (result$feature_table$sample_a[[1]] != sum(sample_mass_tracks$sample_a[[1]]$intensity[3:7])) {
  stop("Expected sample_a peak area to be extracted from sample track.")
}

if (result$feature_table$sample_b[[1]] != sum(sample_mass_tracks$sample_b[[1]]$intensity[3:7])) {
  stop("Expected sample_b peak area to be extracted from sample track.")
}

cat("global features test passed.\n")
