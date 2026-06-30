source("r-prototype/R/feature_tables.R")

feature_table <- data.frame(
  id_number = c("F0", "F1", "F2"),
  mz = c(100.123456, 200.987654, 300.111111),
  rtime = c(10.123, 20.456, 30.789),
  rtime_left_base = c(9.111, 19.222, 29.333),
  rtime_right_base = c(11.111, 21.222, 31.333),
  parent_masstrack_id = c(0L, 1L, 2L),
  peak_area = c(1000, 2000, 3000),
  cSelectivity = c(0.95, 0.95, 0.2),
  goodness_fitting = c(0.95, 0.2, 0.95),
  snr = c(10, 10, 1),
  detection_counts = c(0L, 0L, 0L),
  sample_a = c(100, 0, 50),
  sample_b = c(200, 0, 0),
  stringsAsFactors = FALSE
)

tables <- prepare_feature_tables(
  feature_table,
  sample_names = c("sample_a", "sample_b")
)

if (!identical(tables$full$detection_counts, c(2, 0, 1))) {
  stop("Expected detection_counts to count sample areas greater than one.")
}

if (tables$full$mz[[1]] != 100.1235) {
  stop("Expected mz to be rounded to four decimals.")
}

if (tables$full$rtime[[1]] != 10.12) {
  stop("Expected rtime to be rounded to two decimals.")
}

if (nrow(tables$preferred) != 1L || tables$preferred$id_number[[1]] != "F0") {
  stop("Expected preferred table to keep only high-quality detected features.")
}

cat("feature tables test passed.\n")
