test_that("default parameters contain core mass-track settings", {
  params <- default_parameters()

  expect_equal(params$mz_tolerance_ppm, 5)
  expect_equal(params$min_intensity, 100)
  expect_equal(params$min_timepoints, 5)
  expect_equal(params$min_peak_height, 1000)
})
