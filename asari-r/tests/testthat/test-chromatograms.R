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
