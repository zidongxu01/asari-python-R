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

test_that("mass_paired_mapping matches unsorted lists and preserves positions", {
  result <- mass_paired_mapping(
    c(200.0000, 100.0000),
    c(200.0008, 100.0003),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(2L, 2L), c(1L, 1L)))
  expect_equal(
    result$ratio_deltas,
    c(2.999991000022228e-6, 3.999984000121733e-6),
    tolerance = 1e-16
  )
})

test_that("mass_paired_mapping leaves an ambiguous candidate unmatched", {
  result <- mass_paired_mapping(100, c(100.0002, 100.0003), std_ppm = 5)

  expect_equal(result$mapped, list())
  expect_equal(result$ratio_deltas, numeric())
})

test_that("mass_paired_mapping keeps the list2 minus list1 shift sign", {
  result <- mass_paired_mapping(100.0003, 100, std_ppm = 5)

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(
    result$ratio_deltas,
    -2.999991000022228e-6,
    tolerance = 1e-16
  )
})

test_that("mass_paired_mapping handles empty lists and validates inputs", {
  expect_equal(
    mass_paired_mapping(numeric(), numeric()),
    list(mapped = list(), ratio_deltas = numeric())
  )
  expect_error(mass_paired_mapping("100", 100), "numeric vectors")
  expect_error(mass_paired_mapping(100, 100, std_ppm = -1), "non-negative")
})

test_that("complete_mass_paired_mapping chooses the closest candidate", {
  result <- complete_mass_paired_mapping(
    c(100.0000, 150.0000, 200.0000),
    c(100.0002, 100.0003, 150.0004, 300.0000),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 3L)))
  expect_equal(result$list1_unmapped, 3L)
  expect_equal(result$list2_unmapped, c(2L, 4L))
})

test_that("complete_mass_paired_mapping resolves a shared-list2 conflict", {
  result <- complete_mass_paired_mapping(
    c(99.9998, 100.0001),
    100,
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(2L, 1L)))
  expect_equal(result$list1_unmapped, 1L)
  expect_equal(result$list2_unmapped, integer())
})

test_that("complete_mass_paired_mapping returns all unmatched positions", {
  expect_equal(
    complete_mass_paired_mapping(c(100, 200), c(300, 400)),
    list(
      mapped = list(),
      list1_unmapped = c(1L, 2L),
      list2_unmapped = c(1L, 2L)
    )
  )
  expect_equal(
    complete_mass_paired_mapping(numeric(), 100),
    list(
      mapped = list(),
      list1_unmapped = integer(),
      list2_unmapped = 1L
    )
  )
})

test_that("complete_mass_paired_mapping validates inputs", {
  expect_error(
    complete_mass_paired_mapping("100", 100),
    "numeric vectors"
  )
  expect_error(
    complete_mass_paired_mapping(100, 100, std_ppm = -1),
    "non-negative"
  )
})

test_that("all_mass_paired_mapping keeps every candidate within tolerance", {
  result <- all_mass_paired_mapping(
    c(100.0000, 100.0004, 200.0000),
    c(100.0002, 300.0000),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 1L)))
  expect_equal(result$list1_unmapped, 3L)
  expect_equal(result$list2_unmapped, 2L)
})

test_that("all_mass_paired_mapping searches adjacent Centurion buckets", {
  result <- all_mass_paired_mapping(99.9999, 100.0001, std_ppm = 5)

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(result$list1_unmapped, integer())
  expect_equal(result$list2_unmapped, integer())
})

test_that("all_mass_paired_mapping uses a strict ppm tolerance", {
  result <- all_mass_paired_mapping(
    c(100.00049, 100.00051),
    100,
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(result$list1_unmapped, 2L)
  expect_equal(result$list2_unmapped, integer())
})

test_that("all_mass_paired_mapping returns unmatched positions", {
  expect_equal(
    all_mass_paired_mapping(c(100, 200), c(300, 400)),
    list(
      mapped = list(),
      list1_unmapped = c(1L, 2L),
      list2_unmapped = c(1L, 2L)
    )
  )
  expect_equal(
    all_mass_paired_mapping(numeric(), numeric()),
    list(
      mapped = list(),
      list1_unmapped = integer(),
      list2_unmapped = integer()
    )
  )
})

test_that("all_mass_paired_mapping validates inputs", {
  expect_error(all_mass_paired_mapping("100", 100), "numeric vectors")
  expect_error(
    all_mass_paired_mapping(100, 100, std_ppm = -1),
    "non-negative"
  )
})

test_that("mass mapping correction can recover a new positive-shift match", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200, 300),
    c(100.0004, 200.0008, 300.0018),
    std_ppm = 5,
    correction_tolerance_ppm = 1
  )

  expect_equal(
    result$mapped,
    list(c(1L, 1L), c(2L, 2L), c(3L, 3L))
  )
  expect_equal(
    result$correction_ratio,
    3.999984000054678e-6,
    tolerance = 1e-16
  )
})

test_that("mass mapping correction leaves a sub-threshold shift unchanged", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200),
    c(100.00005, 200.0001),
    correction_tolerance_ppm = 1
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 2L)))
  expect_equal(
    result$correction_ratio,
    4.999997500167233e-7,
    tolerance = 1e-17
  )
})

test_that("mass mapping correction preserves Python negative-shift behavior", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200),
    c(99.9997, 199.9994),
    correction_tolerance_ppm = 1
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 2L)))
  expect_equal(
    result$correction_ratio,
    -2.999999999957481e-6,
    tolerance = 1e-16
  )
})

test_that("mass mapping correction returns NaN when there are no anchors", {
  result <- mass_paired_mapping_with_correction(100, 200)

  expect_equal(result$mapped, list())
  expect_true(is.nan(result$correction_ratio))
})

test_that("mass mapping correction validates inputs", {
  expect_error(
    mass_paired_mapping_with_correction("100", 100),
    "numeric vectors"
  )
  expect_error(
    mass_paired_mapping_with_correction(100, 100, std_ppm = -1),
    "std_ppm"
  )
  expect_error(
    mass_paired_mapping_with_correction(
      100,
      100,
      correction_tolerance_ppm = -1
    ),
    "correction_tolerance_ppm"
  )
})

test_that("landmark guided mapping corrects and extends the reference", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200, 300),
      c(1L, 2L),
      c(100.0004, 200.0008, 300.0018, 400.0016),
      c(1L, 2L, 4L),
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(
    result$new_reference_mzlist,
    c(
      100.00000000079999,
      200.00000000159997,
      300.0003000012,
      400.00000000639994
    ),
    tolerance = 1e-12
  )
  expect_equal(result$new_reference_map2, c(1L, 2L, 3L, 4L))
  expect_equal(result$REF_landmarks, c(1L, 2L, 4L))
  expect_equal(
    result$correction_ratio,
    3.999984000054678e-6,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping uses NA for an absent sample track", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200),
      1L,
      c(100.00005, 300),
      1L,
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(result$new_reference_mzlist, c(100.000025, 200, 300))
  expect_equal(result$new_reference_map2, c(1L, NA_integer_, 2L))
  expect_equal(result$REF_landmarks, 1L)
  expect_equal(
    result$correction_ratio,
    4.999997500167233e-7,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping corrects a negative mass shift", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200, 300),
      c(1L, 2L),
      c(99.9997, 199.9994, 299.9991),
      c(1L, 2L),
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(result$new_reference_map2, c(1L, 2L, 3L))
  expect_equal(
    result$correction_ratio,
    -2.999999999957481e-6,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping leaves correction unset with too few anchors", {
  result <- suppressMessages(
    landmark_guided_mapping(
      100 * seq_len(10),
      1L,
      c(100.0003, 1100),
      1L,
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_null(result$correction_ratio)
  expect_equal(result$new_reference_map2[[1L]], 1L)
  expect_equal(tail(result$new_reference_map2, 1L), 2L)
})

test_that("landmark guided mapping validates landmark positions", {
  expect_error(
    landmark_guided_mapping(100, 0, 100, 1L),
    "REF_landmarks"
  )
  expect_error(
    landmark_guided_mapping(100, 1L, 100, 2L),
    "SM_landmarks"
  )
})
