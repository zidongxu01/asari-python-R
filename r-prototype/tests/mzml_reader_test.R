source("r-prototype/R/mzml_reader.R")

file <- "实验数据/batch10_MT_20210804_001.mzML"

scan_header <- read_mzml_header(file)
if (nrow(scan_header) == 0L) {
  stop("Expected at least one scan in mzML header.")
}

points <- read_mzml_points(file, max_scans = 2, min_intensity = 1000)
expected_columns <- c("scan_id", "rtime", "mz", "intensity")
if (!identical(names(points), expected_columns)) {
  stop("mzML point columns do not match expected columns.")
}
if (nrow(points) == 0L) {
  stop("Expected mzML point preview to contain rows.")
}
if (any(points$intensity <= 1000)) {
  stop("min_intensity filter did not apply correctly.")
}

cat("mzML reader test passed.\n")

