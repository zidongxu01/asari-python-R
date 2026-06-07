source("r-prototype/R/mass_tracks.R")
source("r-prototype/R/preview_utils.R")

file <- preview_input_file("实验数据/batch10_MT_20210804_001.mzML")
print_input_file_summary(file)

result <- extract_mass_tracks(
  file,
  mz_tolerance_ppm = 5,
  min_intensity = 1000,
  min_timepoints = 3,
  min_peak_height = 1000,
  max_scans = 20
)

cat("Tracks:", length(result$tracks), "\n")
cat("Track table preview:\n")
print(head(result$track_table, 20))
