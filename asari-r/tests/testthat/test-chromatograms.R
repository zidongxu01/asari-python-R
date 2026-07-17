test_that("extract_single_track_full_rt_length builds one full-length track", {
  bin <- list(
    c(100.0, 0, 10),
    c(100.2, 1, 20),
    c(99.9, 1, 25)
  )

  track <- extract_single_track_full_rt_length(bin, rt_length = 4)

  expect_equal(track[[1]], 99.95)
  expect_equal(track[[2]], c(10, 25, 0, 0))
})

test_that("bin_to_mass_tracks builds one track for a narrow mz range", {
  bin <- list(
    c(100.0002, 1, 20),
    c(100.0000, 0, 10)
  )

  tracks <- bin_to_mass_tracks(bin, rt_length = 3, mz_tolerance_ppm = 5)

  expect_length(tracks, 1)
  expect_equal(tracks[[1]][[1]], 100.00015)
  expect_equal(tracks[[1]][[2]], c(10, 20, 0))
})

test_that("merge_two_mass_tracks averages mz and sums intensities", {
  track_1 <- list(100.0001, c(10, 20, 0, 5))
  track_2 <- list(100.0003, c(0, 15, 30, 5))

  merged <- merge_two_mass_tracks(track_1, track_2)

  expect_equal(merged[[1]], 100.0002)
  expect_equal(merged[[2]], c(10, 35, 30, 10))
})

test_that("get_thousandth_bins merges adjacent bins and filters noise", {
  make_points <- function(mz, scans, intensities) {
    Map(function(scan, intensity) c(mz, scan, intensity), scans, intensities)
  }

  mz_tree <- list(
    `100000` = make_points(100.0002, 10:14, c(200, 400, 1200, 500, 300)),
    `100001` = make_points(100.0012, 10:14, c(300, 500, 1300, 600, 400)),
    `100010` = make_points(100.0102, c(10, 20, 30, 40, 50), rep(2000, 5)),
    `100020` = make_points(100.0202, 10:14, rep(900, 5))
  )

  bins <- get_thousandth_bins(mz_tree)

  expect_length(bins, 1)
  expect_length(bins[[1]], 10)
})
