test_that("flatten_tuplelist flattens tuple pairs and removes duplicates", {
  tuples <- list(c(1, 2), c(2, 3), c(3, 4))

  expect_equal(flatten_tuplelist(tuples), c(1, 2, 3, 4))
})

test_that("flatten_tuplelist validates its tuple input", {
  expect_null(flatten_tuplelist(list()))
  expect_error(flatten_tuplelist(c(1, 2)), "must be a list")
  expect_error(flatten_tuplelist(list(1, c(2, 3))), "at least two")
})

test_that("check_close_mzs returns one-based adjacent close pairs", {
  mzs <- c(100.0000, 100.0003, 100.0020, 100.0022)

  expect_equal(
    check_close_mzs(mzs, mz_tolerance_ppm = 5),
    list(c(2L, 1L), c(4L, 3L))
  )
})

test_that("check_close_mzs checks adjacent values rather than all pairs", {
  mzs <- c(100.0000, 100.0002, 100.0003)

  expect_equal(
    check_close_mzs(mzs, mz_tolerance_ppm = 5),
    list(c(2L, 1L), c(3L, 2L))
  )
})

test_that("check_close_mzs handles short input and validates sorting", {
  expect_equal(check_close_mzs(numeric()), list())
  expect_equal(check_close_mzs(100), list())
  expect_error(check_close_mzs(c(100.1, 100)), "ascending order")
})

test_that("check_close_mzs uses a strict tolerance comparison", {
  upper_mz <- 100
  exact_lower_boundary <- upper_mz * (1 - 5e-6)

  expect_equal(
    check_close_mzs(c(exact_lower_boundary, upper_mz), 5),
    list()
  )
  expect_equal(
    check_close_mzs(c(exact_lower_boundary + 1e-7, upper_mz), 5),
    list(c(2L, 1L))
  )
})

test_that("calculate_selectivity matches the Python four-value calculation", {
  mzs <- c(100.0000, 100.0005, 100.0010, 100.0015)
  expected <- c(
    0.54657191622213974,
    0.34549779899513522,
    0.34549632904743122,
    0.54656787982816335
  )

  expect_equal(calculate_selectivity(mzs, std_ppm = 5), expected, tolerance = 1e-14)
})

test_that("calculate_selectivity handles interior values and cutoff branches", {
  expect_equal(
    calculate_selectivity(c(100, 100.000001, 100.001, 100.02, 100.2)),
    c(0, 0, 0.74741079425419588, 1, 1),
    tolerance = 1e-14
  )
})

test_that("calculate_selectivity validates its documented preconditions", {
  expect_error(calculate_selectivity(c(100, 101, 102)), "more than three")
  expect_error(
    calculate_selectivity(c(100, 100.2, 100.1, 100.3)),
    "ascending order"
  )
  expect_error(calculate_selectivity(100:103, std_ppm = 0), "positive")
})
