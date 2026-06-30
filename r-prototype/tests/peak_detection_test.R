source("r-prototype/R/peak_detection.R")
source("r-prototype/R/parameters.R")

rt_table <- data.frame(
  scan_id = 1:12,
  scan_index = 1:12,
  rtime = seq(10, 120, by = 10)
)

track <- list(
  mz = 100.001,
  intensity = c(0, 1000, 5000, 12000, 25000, 12000, 5000, 1000, 0, 0, 0, 0)
)

mass_track_result <- list(
  rt_table = rt_table,
  tracks = list(track),
  track_table = data.frame()
)

parameters <- default_parameters()
parameters$min_peak_height <- 10000
parameters$min_timepoints <- 3
parameters$min_intensity_threshold <- 1000
parameters$signal_noise_ratio <- 2
parameters$min_peak_ratio <- 0.001
parameters$min_prominence_threshold <- 5000
parameters$gaussian_shape <- 0.5

features <- detect_peaks(
  mass_track_result,
  parameters = parameters
)

if (nrow(features) != 1L) {
  stop("Expected one detected feature.")
}

if (!identical(features$id_number, "F0")) {
  stop("Expected first feature id to be F0.")
}

if (features$rtime[[1]] != 50) {
  stop("Expected peak apex retention time to be 50.")
}

if (features$peak_area[[1]] <= 0) {
  stop("Expected positive peak area.")
}

if (is.na(features$goodness_fitting[[1]]) || features$goodness_fitting[[1]] <= 0.5) {
  stop("Expected positive Gaussian goodness fitting score.")
}

flat_track <- list(
  mz = 101.001,
  intensity = c(0, 1000, 9000, 11000, 11500, 11000, 9000, 1000, 0, 0, 0, 0)
)

flat_result <- list(
  rt_table = rt_table,
  tracks = list(flat_track),
  track_table = data.frame()
)

flat_parameters <- parameters
flat_parameters$min_prominence_threshold <- 8000

flat_features <- detect_peaks(
  flat_result,
  parameters = flat_parameters
)

if (nrow(flat_features) != 0L) {
  stop("Expected flat low-prominence track to be filtered out.")
}

audited <- audit_mass_track(
  track$intensity,
  min_fwhm = 2,
  min_intensity_threshold = 1000,
  min_peak_height = 10000
)

if (!all(c("baseline", "noise", "intensity", "min_peak_height") %in% names(audited))) {
  stop("Expected audited mass track fields.")
}

gaussian_fit <- evaluate_gaussian_peak(
  values = c(0, 1, 4, 9, 4, 1, 0),
  apex_index = 4,
  left_index = 1,
  right_index = 7,
  min_fwhm = 2
)

if (gaussian_fit$goodness_fitting <= 0.7) {
  stop("Expected Gaussian-like peak to fit well.")
}

scipy_like <- find_peaks_asari(
  values = c(0, 1, 5, 1, 0, 1, 8, 1, 0),
  height = 4,
  distance = 3,
  prominence = 4,
  width = 1,
  wlen = 5
)

if (!identical(scipy_like$peaks, c(3L, 7L))) {
  stop("Expected scipy-like peak positions.")
}

if (any(scipy_like$properties$prominences < 4)) {
  stop("Expected peak prominences to satisfy threshold.")
}

if (any(scipy_like$properties$widths < 1)) {
  stop("Expected peak widths to satisfy threshold.")
}

cat("peak detection test passed.\n")
