# Functions corresponding to asari/chromatograms.py.
#
# This file should contain mass track / EIC construction and RT alignment logic.
# The first implementation target is extract_mass_tracks().

# 从 mzML 文件中读取 MS1 和 MS2 扫描。
# MS1 数据经过强度过滤、m/z 分箱和合并后生成 mass tracks；MS2 谱图单独保存。
extract_mass_tracks <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  alldata <- list()
  rt_times <- numeric()
  ms2_spectra <- list()
  ii <- 0L

  ms_expt <- mzR::openMSfile(infile)
  on.exit(mzR::close(ms_expt), add = TRUE)
  scan_header <- mzR::header(ms_expt)

  timestamp <- tryCatch(
    as.numeric(as.POSIXct(
      mzR::runInfo(ms_expt)$startTimeStamp,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )),
    error = function(e) NULL
  )

  for (spec in seq_len(nrow(scan_header))) {
    if (scan_header$msLevel[[spec]] == 1L) {
      rt_times <- c(rt_times, scan_header$retentionTime[[spec]])
      spectrum <- mzR::peaks(ms_expt, spec)
      intensities <- trunc(spectrum[, "intensity"])
      good_positions <- intensities > min_intensity

      for (jj in which(good_positions)) {
        alldata[[length(alldata) + 1L]] <- c(
          spectrum[jj, "mz"], ii, intensities[[jj]]
        )
      }
      ii <- ii + 1L
    } else if (scan_header$msLevel[[spec]] == 2L) {
      precursor_mz <- scan_header$precursorMZ[[spec]]
      if (is.na(precursor_mz)) next

      spectrum <- mzR::peaks(ms_expt, spec)
      keep <- spectrum[, "intensity"] >= min_intensity &
        spectrum[, "mz"] < precursor_mz - 1
      peaks <- lapply(which(keep), function(jj) {
        c(spectrum[jj, "mz"], spectrum[jj, "intensity"])
      })

      if (length(peaks) > 0L) {
        ms2_spectra[[length(ms2_spectra) + 1L]] <- list(
          precursor_mz = precursor_mz,
          rtime = scan_header$retentionTime[[spec]],
          peaks = peaks
        )
      }
    }
  }

  mz_tree <- list()
  for (x in alldata) {
    key <- as.character(as.integer(x[[1L]] * 1000))
    mz_tree[[key]] <- c(mz_tree[[key]], list(x))
  }

  rm(alldata)
  rt_numbers <- seq_len(length(rt_times)) - 1L
  rt_length <- length(rt_numbers)

  good_bins <- get_thousandth_bins(
    mz_tree, mz_tolerance_ppm, min_timepoints, min_peak_height
  )
  tracks <- list()
  for (bin in good_bins) {
    tracks <- c(tracks, bin_to_mass_tracks(bin, rt_length, mz_tolerance_ppm))
  }

  tracks <- tracks[order(vapply(tracks, function(x) x[[1L]], numeric(1)))]
  merged <- list()
  to_remove <- integer()
  tracks_to_merge <- check_close_mzs(
    vapply(tracks, function(x) x[[1L]], numeric(1)),
    mz_tolerance_ppm
  )
  for (pair in tracks_to_merge) {
    a <- pair[[1L]]
    b <- pair[[2L]]
    merged[[length(merged) + 1L]] <- merge_two_mass_tracks(
      tracks[[a]], tracks[[b]]
    )
    to_remove <- c(to_remove, a, b)
  }
  updated_tracks <- c(
    tracks[!seq_along(tracks) %in% to_remove],
    merged
  )

  list(
    rt_numbers = rt_numbers,
    rt_times = rt_times,
    tracks = updated_tracks,
    ms2_spectra = ms2_spectra,
    acquisition_time = timestamp
  )
}

# 对应Python原名 extract_massTracks_；保留其大小写和尾部下划线接口。
extract_massTracks_ <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  extract_mass_tracks(
    infile, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height
  )
}

# 把有限 m/z 范围内的数据点整理成一条覆盖完整 RT 范围的 mass track。
# 代表 m/z 结合中位数和最高强度点计算；同一扫描存在多个点时保留最大强度。
extract_single_track_full_rt_length <- function(bin,
                                                rt_length,
                                                intensity_data_type = numeric) {
  mzs <- vapply(bin, function(x) x[[1L]], numeric(1))
  ints <- vapply(bin, function(x) x[[3L]], numeric(1))
  mz <- 0.5 * (stats::median(mzs) + mzs[[which.max(ints)]])
  intensity_track <- intensity_data_type(rt_length)

  for (r in bin) {
    scan_index <- r[[2L]] + 1L
    intensity_track[[scan_index]] <- max(
      r[[3L]], intensity_track[[scan_index]]
    )
  }

  list(mz, intensity_track)
}

# 对应Python原名 extract_single_track_fullrt_length。
extract_single_track_fullrt_length <- function(
    bin, rt_length, INTENSITY_DATA_TYPE = numeric) {
  extract_single_track_full_rt_length(bin, rt_length, INTENSITY_DATA_TYPE)
}

# 判断一个 m/z bin 应构成一条还是多条 mass tracks。
# m/z 范围较小时直接构建一条；范围较大时先按 m/z 聚类，再逐组构建。
bin_to_mass_tracks <- function(bin_data_tuples, rt_length, mz_tolerance_ppm = 5) {
  bin_data_tuples <- bin_data_tuples[
    order(vapply(bin_data_tuples, function(x) x[[1L]], numeric(1)))
  ]
  mz_range <- bin_data_tuples[[length(bin_data_tuples)]][[1L]] -
    bin_data_tuples[[1L]][[1L]]
  mz_tolerance <- bin_data_tuples[[1L]][[1L]] * 0.000001 * mz_tolerance_ppm

  if (mz_range < mz_tolerance * 2) {
    return(list(
      extract_single_track_full_rt_length(bin_data_tuples, rt_length)
    ))
  }

  rois <- build_chromatogram_by_mz_clustering(
    bin_data_tuples, rt_length, mz_tolerance
  )
  lapply(rois, function(roi) {
    extract_single_track_full_rt_length(roi, rt_length)
  })
}

# 对应 build_chromatogram_intensity_aware：从最高强度seed开始逐组吸收m/z容差内的点。
build_chromatogram_intensity_aware <- function(
    bin_data_tuples, rt_length, mz_tolerance) {
  invisible(rt_length)
  if (length(bin_data_tuples) == 0L) return(list())
  intensity <- vapply(bin_data_tuples, function(point) point[[3L]], numeric(1))
  remaining <- bin_data_tuples[order(intensity, decreasing = TRUE)]
  assigned <- list()
  while (length(remaining)) {
    seed <- list(remaining[[1L]])
    next_remaining <- list()
    if (length(remaining) > 1L) {
      for (point in remaining[-1L]) {
        if (abs(point[[1L]] - seed[[1L]][[1L]]) < mz_tolerance) {
          seed[[length(seed) + 1L]] <- point
        } else {
          next_remaining[[length(next_remaining) + 1L]] <- point
        }
      }
    }
    assigned[[length(assigned) + 1L]] <- seed
    remaining <- next_remaining
  }
  assigned
}

# 把需要拆分的 m/z 数据交给基于 m/z seeds 的最近邻聚类函数。
# rt_length 目前不参与计算，仅保留以对应 Python 函数签名。
build_chromatogram_by_mz_clustering <- function(bin_data_tuples,
                                                rt_length,
                                                mz_tolerance) {
  nn_cluster_by_mz_seeds(
    bin_data_tuples,
    mz_tolerance,
    presorted = FALSE
  )
}

# 合并两条属于同一离子信号的 mass tracks。
# 两条 track 的代表 m/z 取平均，完整 RT 范围内的强度逐位置相加。
merge_two_mass_tracks <- function(t1, t2) {
  list(
    0.5 * (t1[[1L]] + t2[[1L]]),
    t1[[2L]] + t2[[2L]]
  )
}

# 对应 get_thousandth_bins 内部 __rough_check_consecutive_scans__。
`__rough_check_consecutive_scans__` <- function(
    datatuples, gap_allowed = 2, min_timepoints = 5) {
  checked <- TRUE
  check_max_len <- 4 * min_timepoints
  if (length(datatuples) < check_max_len) {
    min_check_val <- gap_allowed + min_timepoints - 1
    rts <- sort(vapply(datatuples, function(point) point[[2L]], numeric(1)))
    stops <- seq.int(min_timepoints, length(rts))
    steps <- vapply(
      stops,
      function(index) rts[[index]] - rts[[index - min_timepoints + 1L]],
      numeric(1)
    )
    if (min(steps) > min_check_val) checked <- FALSE
  }
  checked
}

# 对应 get_thousandth_bins 内部 __check_min_peak_height__。
`__check_min_peak_height__` <- function(datatuples, min_peak_height) {
  max(vapply(datatuples, function(point) point[[3L]], numeric(1))) >= min_peak_height
}

# 合并相邻或 ppm 范围内的千分位 m/z bins，并过滤弱信号和零散信号。
# 返回的 good bins 会继续交给 bin_to_mass_tracks() 构建 mass tracks。
get_thousandth_bins <- function(mz_tree,
                                mz_tolerance_ppm = 5,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  tolerance <- 0.000001 * mz_tolerance_ppm
  keys <- sort(as.integer(
    names(mz_tree)[lengths(mz_tree) >= min_timepoints]
  ))
  bins_of_bins <- list()
  current_bin <- keys[[1L]]

  if (length(keys) > 1L) {
    for (ii in seq.int(2L, length(keys))) {
      delta <- keys[[ii]] - keys[[ii - 1L]]
      if (delta == 1L || delta < tolerance * keys[[ii]]) {
        current_bin <- c(current_bin, keys[[ii]])
      } else {
        bins_of_bins[[length(bins_of_bins) + 1L]] <- current_bin
        current_bin <- keys[[ii]]
      }
    }
  }

  bins_of_bins[[length(bins_of_bins) + 1L]] <- current_bin
  good_bins <- list()

  for (bin in bins_of_bins) {
    datatuples <- list()
    for (key in bin) {
      datatuples <- c(datatuples, mz_tree[[as.character(key)]])
    }

    if (`__check_min_peak_height__`(datatuples, min_peak_height) &&
        `__rough_check_consecutive_scans__`(
          datatuples, min_timepoints = min_timepoints
        )) {
      good_bins[[length(good_bins) + 1L]] <- datatuples
    }
  }

  good_bins
}

# 对应 SciPy interp1d(..., fill_value = "extrapolate") 的线性插值与外推。
linear_interpolate_with_extrapolation <- function(x, y, xout) {
  interpolation_order <- order(x)
  x <- x[interpolation_order]
  y <- y[interpolation_order]
  result <- stats::approx(x, y, xout = xout)$y

  left <- xout < x[[1L]]
  result[left] <- y[[1L]] + (xout[left] - x[[1L]]) *
    (y[[2L]] - y[[1L]]) / (x[[2L]] - x[[1L]])

  last <- length(x)
  right <- xout > x[[last]]
  result[right] <- y[[last]] + (xout[right] - x[[last]]) *
    (y[[last]] - y[[last - 1L]]) / (x[[last]] - x[[last - 1L]])

  result
}

# 使用匹配的 landmark peaks 拟合 LOWESS，将样本 scan 映射到参考样本 scan。
# 同时返回正向和反向映射，只保存发生变化且没有超出真实范围的 scan。
rt_lowess_calibration <- function(good_landmark_peaks,
                                  selected_reference_landmark_peaks,
                                  sample_rt_numbers,
                                  reference_rt_numbers,
                                  num_iterations,
                                  sample_name,
                                  outdir) {
  reference_rt_bound <- max(reference_rt_numbers)
  sample_rt_bound <- max(sample_rt_numbers)
  rt_rightend <- 1.1 * sample_rt_bound
  xx <- rep(-0.1 * sample_rt_bound, 3L)
  yy <- rep(-0.1 * sample_rt_bound, 3L)

  rt_cal <- clean_rt_calibration_points(Map(
    function(sample_peak, reference_peak) {
      c(sample_peak$apex, reference_peak$apex)
    },
    good_landmark_peaks,
    selected_reference_landmark_peaks
  ))

  xx <- c(
    xx,
    vapply(rt_cal, function(pair) pair[[1L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  yy <- c(
    yy,
    vapply(rt_cal, function(pair) pair[[2L]], numeric(1)),
    rep(rt_rightend, 3L)
  )

  frac <- 0.6 - 0.004 * (length(rt_cal) - 50)
  frac <- max(0.2, min(frac, 0.6))

  lowess_predicted <- hacked_lowess(
    yy,
    xx,
    frac = frac,
    it = num_iterations,
    xvals = sample_rt_numbers
  )

  ref_interpolated <- linear_interpolate_with_extrapolation(
    lowess_predicted,
    sample_rt_numbers,
    reference_rt_numbers
  )

  lowess_predicted <- as.integer(round(lowess_predicted))
  keep_forward <- sample_rt_numbers != lowess_predicted &
    lowess_predicted >= 0L &
    lowess_predicted <= reference_rt_bound
  rt_cal_dict <- stats::setNames(
    lowess_predicted[keep_forward],
    sample_rt_numbers[keep_forward]
  )

  ref_interpolated <- as.integer(round(ref_interpolated))
  keep_reverse <- reference_rt_numbers != ref_interpolated &
    ref_interpolated >= 0L &
    ref_interpolated <= sample_rt_bound
  reverse_rt_cal_dict <- stats::setNames(
    ref_interpolated[keep_reverse],
    reference_rt_numbers[keep_reverse]
  )

  list(rt_cal_dict, reverse_rt_cal_dict)
}

# 删除重复的 RT landmark 配对，并过滤偏差超出三倍总体标准差的异常点。
clean_rt_calibration_points <- function(rt_cal_pairs) {
  pair_keys <- vapply(
    rt_cal_pairs,
    function(pair) paste(pair, collapse = "\r"),
    character(1)
  )
  unique_pairs <- rt_cal_pairs[!duplicated(pair_keys)]
  deltas <- vapply(
    unique_pairs,
    function(pair) pair[[1L]] - pair[[2L]],
    numeric(1)
  )
  mean_delta <- mean(deltas)
  std3x <- sqrt(mean((deltas - mean_delta)^2)) * 3
  low <- mean_delta - std3x
  high <- mean_delta + std3x

  kept <- unique_pairs[deltas > low & deltas < high]
  kept[order(
    vapply(kept, function(pair) pair[[1L]], numeric(1)),
    vapply(kept, function(pair) pair[[2L]], numeric(1))
  )]
}

# 使用 R 的 LOWESS 拟合，并把拟合值线性插值到指定的 scan 坐标。
hacked_lowess <- function(yy, xx, frac, it, xvals) {
  fitted <- stats::lowess(
    xx,
    yy,
    f = frac,
    iter = it,
    delta = 0
  )
  stats::approx(
    fitted$x,
    fitted$y,
    xout = xvals,
    ties = "ordered"
  )$y
}

# 对应Python原名 __hacked_lowess__。
`__hacked_lowess__` <- function(yy, xx, frac, it, xvals) {
  hacked_lowess(yy, xx, frac, it, xvals)
}

# RT 校准的调试版本：返回相同映射，并输出 landmark 与校准曲线图。
rt_lowess_calibration_debug <- function(good_landmark_peaks,
                                         selected_reference_landmark_peaks,
                                         sample_rt_numbers,
                                         reference_rt_numbers,
                                         num_iterations,
                                         sample_name,
                                         outdir) {
  calibration <- rt_lowess_calibration(
    good_landmark_peaks,
    selected_reference_landmark_peaks,
    sample_rt_numbers,
    reference_rt_numbers,
    num_iterations,
    sample_name,
    outdir
  )

  sample_rt_bound <- max(sample_rt_numbers)
  rt_rightend <- 1.1 * sample_rt_bound
  rt_cal <- clean_rt_calibration_points(Map(
    function(sample_peak, reference_peak) {
      c(sample_peak$apex, reference_peak$apex)
    },
    good_landmark_peaks,
    selected_reference_landmark_peaks
  ))
  xx <- c(
    rep(-0.1 * sample_rt_bound, 3L),
    vapply(rt_cal, function(pair) pair[[1L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  yy <- c(
    rep(-0.1 * sample_rt_bound, 3L),
    vapply(rt_cal, function(pair) pair[[2L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  frac <- 0.6 - 0.004 * (length(rt_cal) - 50)
  frac <- max(0.2, min(frac, 0.6))
  lowess_predicted <- as.integer(round(hacked_lowess(
    yy,
    xx,
    frac = frac,
    it = num_iterations,
    xvals = sample_rt_numbers
  )))

  plot_file <- file.path(
    outdir,
    "export",
    paste0(sample_name, "_rtime_alignment_result.png")
  )
  grDevices::png(plot_file)
  graphics::plot(
    xx[4L:(length(xx) - 3L)],
    yy[4L:(length(yy) - 3L)],
    xlab = "Sample Retention Time (sec)",
    ylab = "Reference Retention Time (sec)",
    main = "Retention Time Alignment",
    col = "black",
    pch = 1
  )
  graphics::lines(sample_rt_numbers, lowess_predicted, col = "red")
  graphics::legend(
    "topleft",
    legend = c("Landmark peaks", "Alignment function"),
    col = c("black", "red"),
    pch = c(1, NA),
    lty = c(NA, 1)
  )
  grDevices::dev.off()

  calibration
}

# Python 中尚未实现的 Savitzky-Golay spline 占位函数。
savitzky_golay_spline <- function(good_landmark_peaks,
                                   selected_reference_landmark_peaks,
                                   sample_rt_numbers,
                                   reference_rt_numbers) {
  invisible(NULL)
}

# Python 中明确标记为尚未实现的 DWT RT 校准占位函数。
dwt_rt_calibrate <- function(good_landmark_peaks,
                             selected_reference_landmark_peaks,
                             sample_rt_numbers,
                             reference_rt_numbers) {
  stop("Not implemented", call. = FALSE)
}

# 根据 RT 校准映射，把 intensity 值复制到校准后的 scan 位置。
remap_intensity_track <- function(intensity_track, new, rt_cal_dict) {
  new[seq_along(intensity_track)] <- intensity_track

  for (ii in seq_along(rt_cal_dict)) {
    source_index <- as.integer(names(rt_cal_dict)[[ii]]) + 1L
    target_index <- as.integer(round(rt_cal_dict[[ii]])) + 1L
    if (source_index >= 1L && source_index <= length(intensity_track) &&
        target_index >= 1L && target_index <= length(new)) {
      new[[target_index]] <- intensity_track[[source_index]]
    }
  }

  new
}

# 使用最近边界值填充窗口，对 intensity track 做简单移动平均。
smooth_moving_average <- function(list_intensity, size = 9) {
  left_width <- floor(size / 2)
  right_width <- size - left_width - 1L

  vapply(seq_along(list_intensity), function(ii) {
    positions <- seq.int(ii - left_width, ii + right_width)
    positions <- pmin(length(list_intensity), pmax(1L, positions))
    mean(list_intensity[positions])
  }, numeric(1))
}

# 使用 LOWESS 平滑噪声较大的 intensity track。
smooth_lowess <- function(list_intensity, frac = 0.02) {
  scan_numbers <- seq_along(list_intensity) - 1L
  hacked_lowess(
    list_intensity,
    scan_numbers,
    frac = frac,
    it = 1,
    xvals = scan_numbers
  )
}
