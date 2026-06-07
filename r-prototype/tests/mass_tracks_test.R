source("r-prototype/R/mass_tracks.R")

points <- data.frame(
  scan_id = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
  rtime = c(10, 20, 30, 40, 50, 10, 20, 30, 40, 50),
  mz = c(100.00000, 100.00002, 100.00001, 100.00003, 100.00002,
         101.00000, 101.00002, 101.00001, 101.00003, 101.00002),
  intensity = c(1000, 2000, 3000, 2000, 1000, 500, 600, 700, 600, 500)
)

result <- extract_mass_tracks_from_points(
  points,
  mz_tolerance_ppm = 5,
  min_timepoints = 3,
  min_peak_height = 900
)

if (length(result$tracks) != 1L) {
  stop("Expected one mass track above min_peak_height.")
}

expected_columns <- c("track_id", "mz", "scan_id", "rtime", "intensity")
if (!identical(names(result$track_table), expected_columns)) {
  stop("Mass track table columns do not match expected columns.")
}

if (nrow(result$track_table) != 5L) {
  stop("Expected five nonzero scan points in the retained track.")
}

single_point <- data.frame(
  scan_id = 1,
  rtime = 10,
  mz = 100,
  intensity = 1000
)

single_result <- extract_mass_tracks_from_points(
  single_point,
  mz_tolerance_ppm = 5,
  min_timepoints = 1,
  min_peak_height = 100
)

if (length(single_result$tracks) != 1L) {
  stop("Expected single-point input to produce one track when thresholds allow it.")
}

cat("mass track test passed.\n")
