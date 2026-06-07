source("r-prototype/R/peak_detection.R")
source("r-prototype/R/preview_utils.R")

file <- preview_input_file("r-prototype/scripts/preview_peaks.R")
print_input_file_summary(file)

mass_tracks <- extract_mass_tracks(
  file,
  mz_tolerance_ppm = 5,
  min_intensity = 1000,
  min_timepoints = 3,
  min_peak_height = 1000,
  max_scans = 20
)

features <- detect_peaks(
  mass_tracks,
  min_peak_height = 100000,
  min_timepoints = 3,
  min_intensity_threshold = 1000
)

cat("Features:", nrow(features), "\n")
cat("Feature table preview:\n")
print(head(features, 20))
