# All 19 top-level functions of peaks.py must exist to prevent missing functions when migrating the entire module.
test_that("peaks module contains every Python top-level function", {
  expected_functions <- c(
    "batch_deep_detect_elution_peaks",
    "iter_peak_detection_parameters",
    "stats_detect_elution_peaks",
    "compute_noise_by_flanks",
    "audit_mass_track",
    "detect_evaluate_peaks_on_roi",
    "evaluate_roi_peak_json_",
    "__peaks_cSelectivity_stats_",
    "gaussian_function__",
    "goodness_fitting__",
    "evaluate_gaussian_peak_on_intensity_list",
    "get_gaussian_peakarea_on_intensity_list",
    "lowess_smooth_track",
    "quick_detect_unique_elution_peak",
    "check_overlap_peaks",
    "_check_overlap",
    "_merge_peak_cluster",
    "cleanup_peak_cluster",
    "extend_ROI"
  )

  # Check function objects one by one, not just variables with the same name.
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# These reference values come directly from scipy.signal.find_peaks on the same input.
test_that("SciPy peak properties are reproduced", {
  values <- c(0, 1, 4, 1, 0, 2, 7, 2, 0)
  detected <- .find_peaks_scipy(
    values,
    height = 1,
    distance = 3,
    prominence = 0.5,
    width = 1,
    wlen = 5
  )

  # Peak position, prominence, base, FWHM interpolated position all use Python 0-based semantics.
  expect_equal(detected$peaks, c(2L, 6L))
  expect_equal(detected$properties$peak_heights, c(4, 7))
  expect_equal(detected$properties$prominences, c(4, 7))
  expect_equal(detected$properties$left_bases, c(0L, 4L))
  expect_equal(detected$properties$right_bases, c(4L, 8L))
  expect_equal(
    detected$properties$widths,
    c(1.333333333333, 1.4),
    tolerance = 1e-12
  )
})

# The trade-off between plateau peak center and equal peak distance is a situation where SciPy is prone to boundary differences.
test_that("plateau and distance rules match SciPy", {
  plateau <- .find_peaks_scipy(
    c(0, 1, 1, 1, 0),
    height = 1,
    distance = 3,
    prominence = 0.5,
    width = 1,
    wlen = 5
  )
  equal_peaks <- .find_peaks_scipy(
    c(0, 2, 0, 2, 0),
    height = 1,
    distance = 3,
    prominence = 0.5,
    width = 1,
    wlen = 5
  )

  # Even/odd platforms are rounded down to the SciPy midpoint; contour conflict peaks retain the peak on the right.
  expect_equal(plateau$peaks, 2L)
  expect_equal(plateau$properties$widths, 3)
  expect_equal(equal_peaks$peaks, 3L)
})

# The Python curve_fit reference result for a perfect Gaussian is R2=1, sigma=3.
test_that("Gaussian fitting and area match Python curve_fit", {
  scans <- 0:30
  intensity <- 1000 * exp(-((scans - 15)^2) / (2 * 3^2))
  fitted <- evaluate_gaussian_peak_on_intensity_list(
    intensity, 1000, 15, 0, 30
  )
  area <- get_gaussian_peakarea_on_intensity_list(intensity, 0, 31)

  # The SciPy result for the Gaussian full integration is 7519.884823893001.
  expect_equal(fitted[[1L]], 1, tolerance = 1e-10)
  expect_equal(fitted[[2L]], 3, tolerance = 1e-8)
  expect_equal(area, 7519.884823893001, tolerance = 1e-7)
})

# Verifying Python indexing semantics for local noise slicing and cSelectivity.
test_that("noise and cSelectivity calculations match Python", {
  scans <- 0:30
  intensity <- 1000 * exp(-((scans - 15)^2) / (2 * 3^2))
  noise <- compute_noise_by_flanks(
    list(apex = 15L), intensity, 0:30, 10, 9
  )
  selectivity <- `__peaks_cSelectivity_stats_`(
    intensity,
    list(list(left_base = 8L, right_base = 22L, height = 1000))
  )

  # The reference value is calculated by Python peaks.py on the same vector.
  expect_equal(noise, 310.40833071931445, tolerance = 1e-12)
  expect_equal(selectivity[[1L]], 1)
  expect_equal(selectivity[[2L]], 8:22)
})

# Verify the complete single track process, not just independent mathematical auxiliary functions.
test_that("full statistical peak detection matches Python output", {
  scans <- 0:100
  intensity <- as.integer(
    1000 + 200000 * exp(-((scans - 50)^2) / (2 * 5^2))
  )
  parameters <- list(
    min_peak_height = 100000,
    min_peak_ratio = 0.001,
    min_timepoints = 6,
    min_intensity_threshold = 1000,
    wlen = 25,
    signal_noise_ratio = 2,
    gaussian_shape = 0.5,
    min_prominence_threshold = 33000
  )
  track <- list(
    id_number = 7,
    mz = 123.456,
    rt_scan_numbers = scans,
    intensity = intensity
  )

  # At the same time, the three-layer interfaces of job generation, batch flattening and core detection are tested.
  jobs <- iter_peak_detection_parameters(list(track), 101, parameters)
  direct <- stats_detect_elution_peaks(jobs[[1L]])
  batched <- batch_deep_detect_elution_peaks(list(track), 101, parameters)
  peak <- direct[[1L]]

  expect_length(jobs, 1L)
  expect_equal(batched, direct)
  expect_equal(peak$apex, 50L)
  expect_equal(peak$peak_area, 2466782L)
  expect_equal(peak$height, 193264L)
  expect_equal(peak$left_base, 38L)
  expect_equal(peak$right_base, 62L)
  expect_equal(peak$goodness_fitting, 0.9999996292132238, tolerance = 1e-9)
  expect_equal(peak$cSelectivity, 1)
  expect_equal(peak$parent_masstrack_id, 7)
  expect_equal(peak$mz, 123.456)
  expect_equal(peak$snr, 193L)
})

# Verify discrete logic for overlapping peak boundaries correction, merging, and ROI expansion.
test_that("overlap cleanup and ROI extension match Python", {
  peak1 <- list(apex = 5L, left_base = 0L, right_base = 10L)
  peak2 <- list(apex = 12L, left_base = 6L, right_base = 16L)
  peak3 <- list(apex = 18L, left_base = 12L, right_base = 22L)

  # Two peaks overlap by 4 scans and should be identified as overlapping and re-split at four boundaries.
  expect_true(`_check_overlap`(peak1, peak2))
  cleaned <- cleanup_peak_cluster(list(peak1, peak2))
  expect_equal(cleaned[[1L]]$left_base, 0L)
  expect_equal(cleaned[[1L]]$right_base, 6L)
  expect_equal(cleaned[[2L]]$left_base, 10L)
  expect_equal(cleaned[[2L]]$right_base, 16L)

  # Trimodal clusters inherit the first maximum span peak after merging and extend to the entire cluster boundary.
  merged <- `_merge_peak_cluster`(list(peak1, peak2, peak3))
  expect_equal(merged$left_base, 0L)
  expect_equal(merged$right_base, 22L)
  expect_equal(extend_ROI(c(2L, 3L, 4L), 10L), 0:7)
})
