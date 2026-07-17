test_that("ppm_to_mz_tolerance converts ppm to absolute m/z tolerance", {
  expect_equal(ppm_to_mz_tolerance(100, 5), 0.0005)
  expect_equal(ppm_to_mz_tolerance(800, 5), 0.004)
})
