source("r-prototype/R/peak_detection.R")

rt_table <- data.frame(
  scan_id = 1:9,
  scan_index = 1:9,
  rtime = seq(10, 90, by = 10)
)

track <- list(
  mz = 100.001,
  intensity = c(0, 1000, 5000, 20000, 5000, 1000, 0, 0, 0)
)

mass_track_result <- list(
  rt_table = rt_table,
  tracks = list(track),
  track_table = data.frame()
)

features <- detect_peaks(
  mass_track_result,
  min_peak_height = 10000,
  min_timepoints = 3,
  min_intensity_threshold = 1000,
  signal_noise_ratio = 2
)

if (nrow(features) != 1L) {
  stop("Expected one detected feature.")
}

if (!identical(features$id_number, "F0")) {
  stop("Expected first feature id to be F0.")
}

if (features$rtime[[1]] != 40) {
  stop("Expected peak apex retention time to be 40.")
}

if (features$peak_area[[1]] <= 0) {
  stop("Expected positive peak area.")
}

cat("peak detection test passed.\n")

