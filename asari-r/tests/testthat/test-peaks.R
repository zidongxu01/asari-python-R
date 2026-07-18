# peaks.py 的 19 个顶层函数必须全部存在，防止整模块迁移时漏函数。
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

  # 逐个检查函数对象，而不只检查同名变量。
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# 这些参考值直接来自相同输入下的 scipy.signal.find_peaks。
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

  # 峰位置、prominence、base、FWHM 插值位置都使用 Python 0-based 语义。
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

# 平台峰中心和等高峰 distance 取舍是 SciPy 容易产生边界差异的情况。
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

  # 偶数/奇数平台按 SciPy 中点向下取整；等高冲突峰保留右侧峰。
  expect_equal(plateau$peaks, 2L)
  expect_equal(plateau$properties$widths, 3)
  expect_equal(equal_peaks$peaks, 3L)
})

# 完美 Gaussian 的 Python curve_fit 参考结果是 R2=1、sigma=3。
test_that("Gaussian fitting and area match Python curve_fit", {
  scans <- 0:30
  intensity <- 1000 * exp(-((scans - 15)^2) / (2 * 3^2))
  fitted <- evaluate_gaussian_peak_on_intensity_list(
    intensity, 1000, 15, 0, 30
  )
  area <- get_gaussian_peakarea_on_intensity_list(intensity, 0, 31)

  # Gaussian 全积分的 SciPy 结果为 7519.884823893001。
  expect_equal(fitted[[1L]], 1, tolerance = 1e-10)
  expect_equal(fitted[[2L]], 3, tolerance = 1e-8)
  expect_equal(area, 7519.884823893001, tolerance = 1e-7)
})

# 验证局部噪声切片和 cSelectivity 的 Python 索引语义。
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

  # 参考值由 Python peaks.py 在同一向量上计算得到。
  expect_equal(noise, 310.40833071931445, tolerance = 1e-12)
  expect_equal(selectivity[[1L]], 1)
  expect_equal(selectivity[[2L]], 8:22)
})

# 验证完整单 track 流程，不只验证独立数学辅助函数。
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

  # 同时测试 job 生成、批量展平和核心检测三层接口。
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

# 验证重叠峰边界修正、合并和 ROI 扩展的离散逻辑。
test_that("overlap cleanup and ROI extension match Python", {
  peak1 <- list(apex = 5L, left_base = 0L, right_base = 10L)
  peak2 <- list(apex = 12L, left_base = 6L, right_base = 16L)
  peak3 <- list(apex = 18L, left_base = 12L, right_base = 22L)

  # 两峰重叠 4 scans，应被识别为重叠并重新拆分四个边界。
  expect_true(`_check_overlap`(peak1, peak2))
  cleaned <- cleanup_peak_cluster(list(peak1, peak2))
  expect_equal(cleaned[[1L]]$left_base, 0L)
  expect_equal(cleaned[[1L]]$right_base, 6L)
  expect_equal(cleaned[[2L]]$left_base, 10L)
  expect_equal(cleaned[[2L]]$right_base, 16L)

  # 三峰簇合并后继承第一个最大跨度峰，并扩展到整个簇边界。
  merged <- `_merge_peak_cluster`(list(peak1, peak2, peak3))
  expect_equal(merged$left_base, 0L)
  expect_equal(merged$right_base, 22L)
  expect_equal(extend_ROI(c(2L, 3L, 4L), 10L), 0:7)
})
