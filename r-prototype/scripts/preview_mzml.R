source("r-prototype/R/mzml_reader.R")
source("r-prototype/R/preview_utils.R")

file <- preview_input_file("实验数据/batch10_MT_20210804_001.mzML")
print_input_file_summary(file)

scan_header <- read_mzml_header(file)
cat("Header dimensions:\n")
print(dim(scan_header))
cat("First MS1 scans:\n")
print(head(scan_header[, c("seqNum", "msLevel", "retentionTime", "peaksCount")], 3))

points <- read_mzml_points(file, max_scans = 2, min_intensity = 1000)
cat("Point preview:\n")
print(head(points, 10))
cat("Rows kept:", nrow(points), "\n")
