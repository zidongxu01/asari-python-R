# 对应 Python asari/peaks.py：色谱峰检测、噪声估计和峰质量评价。
#
# 兼容性约定：
# 1. 对外暴露的 scan/apex/left_base/right_base 均保持 Python 的 0-based 编号。
# 2. mass track 使用具名 list，字段与 Python dict 一致。
# 3. 本文件只使用 base R 和 stats，避免为了 SciPy 的少量功能引入大型依赖。

# 把 Python 切片 [start:stop] 的行为复现在 R 中。
# start 和 stop 都使用 Python 的 0-based 边界；stop 所在位置不包含在结果中。
.python_slice <- function(values, start = NULL, stop = NULL) {
  n <- length(values)

  # Python 会把负边界解释为从向量末尾反向计数。
  normalize_bound <- function(bound, default) {
    if (is.null(bound)) return(default)
    bound <- as.integer(bound)
    if (bound < 0L) bound <- max(n + bound, 0L)
    min(max(bound, 0L), n)
  }

  start <- normalize_bound(start, 0L)
  stop <- normalize_bound(stop, n)
  if (stop <= start) return(values[0])

  # R 下标从 1 开始，因此把 Python 的 [start, stop) 转换成 R 区间。
  values[seq.int(start + 1L, stop)]
}

# 复现 Python int() 的向零截断，同时避免 R 32 位 integer 的溢出 NA。
.python_int <- function(value) {
  truncated <- trunc(value)
  if (length(truncated) == 1L && is.finite(truncated) &&
      truncated >= -.Machine$integer.max &&
      truncated <= .Machine$integer.max) {
    return(as.integer(truncated))
  }

  # 超出 R integer 范围时用无小数 double 保存同一个整数数值。
  truncated
}

# 复现 scipy.signal.detrend(..., type = "linear") 的线性去趋势。
.detrend_linear <- function(values) {
  values <- as.numeric(values)
  if (length(values) < 2L) return(values - mean(values))

  # SciPy 对样本序号拟合一次直线，然后返回原信号减拟合值。
  x <- seq_along(values) - 1L
  fit <- stats::lm.fit(cbind(1, x), values)
  values - as.vector(cbind(1, x) %*% fit$coefficients)
}

# 找出局部极大值及平台峰中心，行为对应 SciPy find_peaks 的第一阶段。
.local_maxima_1d <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  if (n < 3L) return(integer())

  peaks <- integer()
  ii <- 2L
  while (ii <= n - 1L) {
    # 峰的左侧必须严格上升；相等值由后面的平台逻辑统一处理。
    if (values[[ii - 1L]] < values[[ii]]) {
      right_edge <- ii
      while (right_edge < n && values[[right_edge + 1L]] == values[[ii]]) {
        right_edge <- right_edge + 1L
      }

      # 峰的右侧也必须严格下降；平台峰取中点并向下取整。
      if (right_edge < n && values[[right_edge]] > values[[right_edge + 1L]]) {
        peaks <- c(peaks, as.integer(floor((ii + right_edge) / 2)))
      }
      ii <- right_edge + 1L
    } else {
      ii <- ii + 1L
    }
  }

  # 私有峰查找器内部使用 R 的 1-based 下标。
  peaks
}

# 按 SciPy 的 distance 规则保留高优先级峰。
.select_by_peak_distance <- function(peaks, heights, distance) {
  if (length(peaks) < 2L) return(seq_along(peaks))
  distance <- as.integer(ceiling(distance))
  keep <- rep(TRUE, length(peaks))

  # SciPy 以峰高为优先级；同高峰按位置稳定排序后从右向左处理。
  priority_order <- order(heights, seq_along(heights))
  for (position in rev(priority_order)) {
    if (!keep[[position]]) next

    # 删除距离当前高优先级峰小于 distance 的左右邻峰。
    left <- position - 1L
    while (left >= 1L && peaks[[position]] - peaks[[left]] < distance) {
      keep[[left]] <- FALSE
      left <- left - 1L
    }
    right <- position + 1L
    while (right <= length(peaks) && peaks[[right]] - peaks[[position]] < distance) {
      keep[[right]] <- FALSE
      right <- right + 1L
    }
  }

  which(keep)
}

# 为单个峰计算 prominence、左右 base 以及半高宽。
.peak_prominence_and_width <- function(values, peak, wlen = NULL) {
  n <- length(values)
  peak_height <- values[[peak]]

  # SciPy 会把 wlen 向上取整为奇数，并以峰顶为窗口中心。
  if (is.null(wlen)) {
    left_limit <- 1L
    right_limit <- n
  } else {
    rounded_wlen <- as.integer(ceiling(wlen))
    if (rounded_wlen %% 2L == 0L) rounded_wlen <- rounded_wlen + 1L
    half_window <- rounded_wlen %/% 2L
    left_limit <- max(1L, peak - half_window)
    right_limit <- min(n, peak + half_window)
  }

  # 从峰顶向左扫描；遇到更高点停止，并记录最低点作为 left base。
  left_base <- peak
  left_min <- peak_height
  cursor <- peak
  while (cursor >= left_limit && values[[cursor]] <= peak_height) {
    if (values[[cursor]] < left_min) {
      left_min <- values[[cursor]]
      left_base <- cursor
    }
    cursor <- cursor - 1L
  }

  # 从峰顶向右执行同样的扫描，得到 right base。
  right_base <- peak
  right_min <- peak_height
  cursor <- peak
  while (cursor <= right_limit && values[[cursor]] <= peak_height) {
    if (values[[cursor]] < right_min) {
      right_min <- values[[cursor]]
      right_base <- cursor
    }
    cursor <- cursor + 1L
  }

  prominence <- peak_height - max(left_min, right_min)
  evaluation_height <- peak_height - 0.5 * prominence

  # 在 left base 与峰顶之间寻找半高交点，并进行线性插值。
  left_ip <- as.numeric(peak)
  cursor <- peak
  while (cursor > left_base && values[[cursor]] > evaluation_height) {
    cursor <- cursor - 1L
  }
  if (values[[cursor]] < evaluation_height && cursor < peak) {
    left_ip <- cursor + (evaluation_height - values[[cursor]]) /
      (values[[cursor + 1L]] - values[[cursor]])
  } else {
    left_ip <- as.numeric(cursor)
  }

  # 在峰顶与 right base 之间寻找半高交点，并进行线性插值。
  right_ip <- as.numeric(peak)
  cursor <- peak
  while (cursor < right_base && values[[cursor]] > evaluation_height) {
    cursor <- cursor + 1L
  }
  if (values[[cursor]] < evaluation_height && cursor > peak) {
    right_ip <- cursor - (evaluation_height - values[[cursor]]) /
      (values[[cursor - 1L]] - values[[cursor]])
  } else {
    right_ip <- as.numeric(cursor)
  }

  list(
    prominence = prominence,
    left_base = left_base,
    right_base = right_base,
    width = right_ip - left_ip,
    left_ip = left_ip,
    right_ip = right_ip
  )
}

# 实现 peaks.py 实际用到的 scipy.signal.find_peaks 参数子集。
.find_peaks_scipy <- function(values,
                              height = NULL,
                              distance = NULL,
                              prominence = NULL,
                              width = NULL,
                              wlen = NULL) {
  values <- as.numeric(values)

  # 第一步寻找所有局部极大值和平台峰中心。
  peaks <- .local_maxima_1d(values)
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(), prominences = numeric(),
      left_bases = integer(), right_bases = integer(), widths = numeric()
    )))
  }

  # 第二步应用最小峰高；此处仍保留 R 内部的 1-based 峰位置。
  if (!is.null(height)) {
    peaks <- peaks[values[peaks] >= height]
  }
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(), prominences = numeric(),
      left_bases = integer(), right_bases = integer(), widths = numeric()
    )))
  }

  # 第三步按峰高优先级执行最小峰间距筛选。
  if (!is.null(distance)) {
    selected <- .select_by_peak_distance(peaks, values[peaks], distance)
    peaks <- peaks[selected]
  }

  # 第四步计算 prominence、bases 和 FWHM width。
  metrics <- lapply(peaks, function(peak) {
    .peak_prominence_and_width(values, peak, wlen)
  })
  prominences <- vapply(metrics, `[[`, numeric(1), "prominence")
  widths <- vapply(metrics, `[[`, numeric(1), "width")

  # SciPy 先筛 prominence，再筛 width；每次都同步筛选全部属性。
  keep <- rep(TRUE, length(peaks))
  if (!is.null(prominence)) keep <- keep & prominences >= prominence
  if (!is.null(width)) keep <- keep & widths >= width
  peaks <- peaks[keep]
  metrics <- metrics[keep]
  prominences <- prominences[keep]
  widths <- widths[keep]

  # 对外转换回 Python 使用的 0-based 峰和 base 编号。
  properties <- list(
    peak_heights = values[peaks],
    prominences = prominences,
    left_bases = vapply(metrics, `[[`, integer(1), "left_base") - 1L,
    right_bases = vapply(metrics, `[[`, integer(1), "right_base") - 1L,
    widths = widths,
    left_ips = vapply(metrics, `[[`, numeric(1), "left_ip") - 1,
    right_ips = vapply(metrics, `[[`, numeric(1), "right_ip") - 1
  )

  list(peaks = as.integer(peaks - 1L), properties = properties)
}

# 使用非线性最小二乘目标拟合 Gaussian；失败时由调用方采用 Python 相同的回退值。
.fit_gaussian <- function(x, y, start) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  start <- as.numeric(start)

  # curve_fit 的目标是最小化 Gaussian 预测值与观测值的残差平方和。
  objective <- function(parameters) {
    predicted <- gaussian_function__(
      x, parameters[[1L]], parameters[[2L]], parameters[[3L]]
    )
    residual_sum <- sum((predicted - y)^2)
    if (is.finite(residual_sum)) residual_sum else .Machine$double.xmax
  }

  # BFGS 对完美 Gaussian 也能收敛；R 的 nls 在零残差附近会报奇异梯度。
  parameter_scale <- pmax(abs(start), c(1, 1, 1))
  fit <- stats::optim(
    start,
    objective,
    method = "BFGS",
    control = list(
      maxit = 10000L,
      reltol = 1e-12,
      parscale = parameter_scale
    )
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par))) {
    stop("Gaussian fitting failed to converge")
  }

  # 再从 BFGS 解做少量 Gauss-Newton 迭代，以贴近 curve_fit 的机器精度结果。
  refined <- tryCatch(
    suppressWarnings(stats::nls(
      y ~ a * exp(-((x - mu)^2) / (2 * sigma^2)),
      start = list(a = fit$par[[1L]], mu = fit$par[[2L]], sigma = fit$par[[3L]]),
      control = stats::nls.control(
        maxiter = 20L, tol = 1e-10, warnOnly = TRUE
      )
    )),
    error = function(e) NULL
  )
  if (!is.null(refined) && all(is.finite(stats::coef(refined)))) {
    return(stats::coef(refined))
  }

  stats::setNames(fit$par, c("a", "mu", "sigma"))
}

# 批量检测多条 mass tracks，并把每条 track 的结果展平成一个 peak list。
batch_deep_detect_elution_peaks <- function(list_mass_tracks,
                                            number_of_scans,
                                            parameters) {
  # 先生成与 Python multiprocessing 输入一致的 job 列表。
  jobs <- iter_peak_detection_parameters(
    list_mass_tracks, number_of_scans, parameters
  )

  # R 版当前顺序执行；顺序不会改变数值结果或输出排列。
  results <- lapply(jobs, stats_detect_elution_peaks)
  if (length(results) == 0L) return(list())
  unlist(results, recursive = FALSE)
}

# 为超过最低强度阈值的 mass tracks 生成峰检测任务。
iter_peak_detection_parameters <- function(list_mass_tracks,
                                           number_of_scans,
                                           parameters) {
  # Python np.any(...) 只要任意强度严格大于阈值便保留该 track。
  selected_tracks <- Filter(function(track) {
    any(track$intensity > parameters$min_intensity_threshold)
  }, list_mass_tracks)

  # 每个 job 保持 (mass_track, number_of_scans, parameters) 三元结构。
  lapply(selected_tracks, function(track) {
    list(track, number_of_scans, parameters)
  })
}

# 对一条 mass track 执行 asari 的统计引导峰检测主流程。
stats_detect_elution_peaks <- function(job) {
  # 解包 Python job 三元组，并读取全部检测参数。
  mass_track <- job[[1L]]
  number_of_scans <- as.integer(job[[2L]])
  parameters <- job[[3L]]
  min_peak_height <- parameters$min_peak_height
  min_peak_ratio <- parameters$min_peak_ratio
  min_fwhm <- round(0.5 * parameters$min_timepoints)
  min_intensity_threshold <- parameters$min_intensity_threshold
  wlen <- parameters$wlen
  snr <- parameters$signal_noise_ratio
  peakshape <- parameters$gaussian_shape
  min_prominence_ratio <- 0.02
  min_prominence_threshold <- parameters$min_prominence_threshold

  # 先清理、缩放并估计整条 mass track 的基线和噪声。
  audited <- audit_mass_track(
    mass_track$intensity,
    min_fwhm,
    min_intensity_threshold,
    min_peak_height,
    min_peak_ratio
  )
  baseline <- audited[[1L]]
  noise_level <- audited[[2L]]
  scaling_factor <- audited[[3L]]
  min_peak_height <- audited[[4L]]
  list_intensity <- audited[[5L]]
  list_scans <- seq_len(number_of_scans) - 1L
  list_json_peaks <- list()
  list_peaks <- list()

  # 用高于噪声水平的扫描点建立 ROI；扫描编号保持 0-based。
  selected_scans <- list_scans[list_intensity > noise_level]

  # Python `array.any()` 对仅包含 0 的数组返回 FALSE，这里明确复现该行为。
  if (length(selected_scans) > 0L && any(selected_scans != 0L)) {
    rois <- list()
    current_roi <- selected_scans[[1L]]

    # 相邻信号允许最多两个扫描的缺口，并把缺口补入 ROI。
    if (length(selected_scans) > 1L) {
      for (scan in selected_scans[-1L]) {
        if (scan - tail(current_roi, 1L) < 3L) {
          current_roi <- c(
            current_roi,
            seq.int(tail(current_roi, 1L) + 1L, scan)
          )
        } else {
          rois[[length(rois) + 1L]] <- current_roi
          current_roi <- scan
        }
      }
    }
    rois[[length(rois) + 1L]] <- current_roi

    # 删除过短 ROI，避免峰宽和 Gaussian 拟合没有足够点数。
    rois <- Filter(function(roi) length(roi) >= min_fwhm + 2L, rois)

    # 分别检测每个 ROI，并汇总 JSON 风格的峰记录。
    for (roi in rois) {
      if (length(roi) < 3 * min_fwhm) {
        roi <- extend_ROI(roi, number_of_scans)
      }
      intensity_roi <- list_intensity[roi + 1L]
      detected <- detect_evaluate_peaks_on_roi(
        intensity_roi,
        roi,
        min_peak_height,
        min_peak_ratio,
        min_fwhm,
        min_prominence_threshold,
        wlen,
        snr,
        peakshape,
        min_prominence_ratio,
        noise_level
      )
      list_json_peaks <- c(list_json_peaks, detected)
    }

    # 计算峰的 cSelectivity，并收集不属于任何峰的噪声扫描点。
    selectivity_stats <- `__peaks_cSelectivity_stats_`(
      list_intensity, list_json_peaks
    )
    list_cSelectivity <- selectivity_stats[[1L]]
    peak_datapoints <- selectivity_stats[[2L]]
    noise_data_points <- setdiff(list_scans, peak_datapoints)

    # 将选择性、母 mass track 信息、原始量级高度和局部 SNR 写回每个峰。
    for (ii in seq_along(list_json_peaks)) {
      peak <- list_json_peaks[[ii]]
      peak$cSelectivity <- list_cSelectivity[[ii]]
      peak$parent_masstrack_id <- mass_track$id_number
      peak$mz <- mass_track$mz
      peak$height <- scaling_factor * peak$height + baseline
      capped_height <- min(peak$height, 1e100)
      local_noise <- compute_noise_by_flanks(
        peak,
        mass_track$intensity,
        noise_data_points,
        min_intensity_threshold,
        noise_level
      )
      peak$snr <- .python_int(capped_height / local_noise)

      # 只有达到指定 SNR 的峰才进入最终结果；高度按 Python int() 截断。
      if (peak$snr >= snr) {
        peak$height <- .python_int(peak$height)
        list_peaks[[length(list_peaks) + 1L]] <- peak
      }
    }
  }

  list_peaks
}

# 用峰两侧的非峰数据点估计局部噪声。
compute_noise_by_flanks <- function(peak,
                                    list_intensity,
                                    noise_data_points,
                                    min_intensity_threshold,
                                    old_noise_level) {
  padding <- 30L
  number_of_points <- 100L

  # 加入 apex 后排序，并采用 Python list.index() 的第一个匹配位置。
  points <- sort(c(noise_data_points, peak$apex))
  center <- match(peak$apex, points) - 1L

  # 严格复现两个 Python 切片，包括负起点从列表尾部计数的语义。
  left_points <- .python_slice(
    points,
    center - padding - number_of_points,
    center - padding
  )
  right_points <- .python_slice(
    points,
    center + padding,
    center + number_of_points + padding
  )
  use_index <- c(left_points, right_points)

  # Python 在强度上加 1 后取均值，并以最低强度阈值作为下界。
  if (length(use_index) > 0L) {
    return(max(
      min_intensity_threshold,
      mean(list_intensity[use_index + 1L] + 1)
    ))
  }
  old_noise_level
}

# 清理 mass track，并返回基线、噪声、缩放因子、动态峰高阈值和清理后强度。
audit_mass_track <- function(list_intensity,
                             min_fwhm,
                             min_intensity_threshold,
                             min_peak_height,
                             min_peak_ratio) {
  # SciPy uniform_filter1d 会保留输入 dtype；先记录 Python ndarray 是否为整数型。
  working_is_integer <- is.integer(list_intensity)
  list_intensity <- as.numeric(list_intensity)
  scaling_factor <- 1
  low <- min_intensity_threshold
  high <- 1e8
  baseline <- low
  noise_level <- low
  max_intensity <- max(list_intensity)

  # 超过 1e8 的 track 先缩放，防止 prominence 计算被极大数破坏。
  if (max_intensity > high) {
    scaling_factor <- max_intensity / high
    list_intensity <- list_intensity / scaling_factor
    working_is_integer <- FALSE
  }

  # 只有中位数超过最低阈值时才重新估计基线和噪声。
  median_intensity <- stats::median(list_intensity)
  if (median_intensity > low) {
    above_low <- list_intensity[list_intensity > low]

    # 高基线 track 先执行与 scipy.signal.detrend 相同的线性去趋势。
    if (length(above_low) > length(list_intensity) * 0.5 &&
        median_intensity > 10 * min_peak_height) {
      list_intensity <- .detrend_linear(list_intensity)
      working_is_integer <- FALSE
    }

    # 取低于“下四分位数 + LOW”的点估计背景分布。
    bottom_cutoff <- low + as.numeric(stats::quantile(
      list_intensity, 0.25, names = FALSE, type = 7
    ))
    bottom_values <- list_intensity[list_intensity < bottom_cutoff]
    baseline <- mean(bottom_values)

    # NumPy std 默认 ddof=0，不能直接使用 R 默认 ddof=1 的 sd()。
    noise_level <- sqrt(mean((bottom_values - baseline)^2))
  }

  # 基线和噪声都不能低于预设的最低强度阈值。
  baseline <- max(baseline, low)
  noise_level <- max(noise_level, low)

  # 极高、经过缩放或相对噪声较大的 track 按 Python 规则做移动平均。
  if (scaling_factor > 1 || max_intensity > 1000 * min_peak_height) {
    list_intensity <- smooth_moving_average(
      list_intensity, size = min_fwhm * 3
    )
    if (working_is_integer) list_intensity <- trunc(list_intensity)
  } else if (100 * noise_level > max_intensity ||
             max_intensity < 10 * min_peak_height) {
    list_intensity <- smooth_moving_average(
      list_intensity, size = min_fwhm + 2
    )
    if (working_is_integer) list_intensity <- trunc(list_intensity)
  }

  # 最后扣除基线，并按全 track 最大强度动态调整最低峰高。
  list_intensity <- list_intensity - baseline
  min_peak_height <- max(
    min_peak_height,
    max_intensity * min_peak_ratio
  ) - baseline

  list(
    baseline,
    noise_level,
    scaling_factor,
    min_peak_height,
    list_intensity
  )
}

# 在一个 ROI 内检测候选峰、评价 Gaussian 峰形并清理重叠峰。
detect_evaluate_peaks_on_roi <- function(list_intensity_roi,
                                         rt_numbers_roi,
                                         min_peak_height,
                                         min_peak_ratio,
                                         min_fwhm,
                                         min_prominence_threshold,
                                         wlen,
                                         snr,
                                         peakshape,
                                         min_prominence_ratio,
                                         noise_level) {
  # min_peak_ratio 和 snr 在 Python 函数签名中保留，但本层不直接使用。
  invisible(min_peak_ratio)
  invisible(snr)
  list_peaks <- list()
  max_intensity <- max(list_intensity_roi)

  # prominence 阈值同时受固定阈值、ROI 高度、峰高阈值和噪声控制。
  roi_prominence <- max(
    min_prominence_threshold,
    max_intensity * min_prominence_ratio,
    0.33 * min_peak_height,
    2 * noise_level
  )

  # 首轮严格使用 wlen，参数集合与 Python scipy.find_peaks 调用一致。
  detected <- .find_peaks_scipy(
    list_intensity_roi,
    height = min_peak_height,
    distance = min_fwhm,
    prominence = roi_prominence,
    width = min_fwhm,
    wlen = wlen
  )

  # 若首轮无峰且 ROI 长于 wlen，则取消 prominence 窗口重新检测。
  if (length(detected$peaks) == 0L && length(list_intensity_roi) > wlen) {
    detected <- .find_peaks_scipy(
      list_intensity_roi,
      height = min_peak_height,
      distance = min_fwhm,
      prominence = roi_prominence,
      width = min_fwhm,
      wlen = NULL
    )
  }

  # 逐峰检查 base 跨度，再计算 Gaussian 拟合和最终 JSON 字段。
  for (ii in seq_along(detected$peaks)) {
    if (detected$properties$right_bases[[ii]] -
        detected$properties$left_bases[[ii]] >= min_fwhm + 2L) {
      evaluated <- evaluate_roi_peak_json_(
        ii - 1L,
        list_intensity_roi,
        rt_numbers_roi,
        detected$peaks,
        detected$properties,
        peakshape,
        min_fwhm
      )
      if (!is.null(evaluated)) {
        list_peaks[[length(list_peaks) + 1L]] <- evaluated
      }
    }
  }

  check_overlap_peaks(list_peaks)
}

# 把第 ii 个 ROI 候选峰转换为带边界和峰形分数的具名 list。
evaluate_roi_peak_json_ <- function(ii,
                                    list_intensity_roi,
                                    rt_numbers_roi,
                                    peaks,
                                    properties,
                                    peakshape,
                                    min_fwhm) {
  # ii 来自 Python 语义，所以访问 R 属性向量时先加 1。
  position <- as.integer(ii) + 1L
  left_index <- properties$left_bases[[position]]
  right_index <- properties$right_bases[[position]]

  # Gaussian 拟合的 apex 和边界都相对于当前 ROI，且保持 0-based。
  gaussian <- evaluate_gaussian_peak_on_intensity_list(
    list_intensity_roi,
    properties$peak_heights[[position]],
    peaks[[position]],
    left_index,
    right_index
  )
  goodness_fitting <- gaussian[[1L]]
  sigma <- gaussian[[2L]]

  # 只有峰形分数严格超过阈值时才保留该峰。
  if (goodness_fitting > peakshape) {
    halfwidth <- .python_int(max(abs(sigma), min_fwhm) * 6)
    proposed_left <- peaks[[position]] - halfwidth
    proposed_right <- peaks[[position]] + halfwidth
    left_index <- max(left_index, proposed_left)
    right_index <- min(right_index, proposed_right)

    # 把 ROI 内部位置映射回完整 mass track 的扫描编号。
    left_base <- rt_numbers_roi[[left_index + 1L]]
    right_base <- rt_numbers_roi[[right_index + 1L]]
    apex <- rt_numbers_roi[[peaks[[position]] + 1L]]
    peak_area <- .python_int(sum(
      list_intensity_roi[seq.int(left_index + 1L, right_index + 1L)]
    ))

    return(list(
      apex = apex,
      peak_area = peak_area,
      height = .python_int(properties$peak_heights[[position]]),
      left_base = left_base,
      right_base = right_base,
      goodness_fitting = goodness_fitting
    ))
  }

  NULL
}

# 计算一条 mass track 上每个峰的 cSelectivity 和全部峰覆盖的扫描点。
`__peaks_cSelectivity_stats_` <- function(list_intensity, jpeaks) {
  # 汇总每个峰左右 base 之间的闭区间；unique 对应 Python set。
  peak_datapoints <- unique(unlist(lapply(jpeaks, function(peak) {
    seq.int(peak$left_base, peak$right_base)
  }), use.names = FALSE))

  # 空峰列表必须产生空强度和空分数，而不是错误地下标。
  if (length(peak_datapoints) == 0L) {
    return(list(numeric(), integer()))
  }
  peak_intensities <- list_intensity[peak_datapoints + 1L]

  # 对每个峰，以一半峰高为阈值比较峰区点数和全 track 背景点数。
  selectivity <- vapply(jpeaks, function(peak) {
    threshold <- 0.5 * peak$height
    peak_level <- sum(peak_intensities > threshold)
    background_level <- sum(list_intensity > threshold)

    # 与 Python 一致：只有 0 < peak_level <= background_level 时才计算比例。
    if (background_level >= peak_level && peak_level > 0L) {
      peak_level / background_level
    } else {
      0
    }
  }, numeric(1))

  list(selectivity, sort(as.integer(peak_datapoints)))
}

# 标准 Gaussian 函数，a 为高度、mu 为中心、sigma 为标准差。
gaussian_function__ <- function(x, a, mu, sigma) {
  # 公式与 Python a*np.exp(-(x-mu)**2/(2*sigma**2)) 完全相同。
  a * exp(-((x - mu)^2) / (2 * sigma^2))
}

# 用决定系数 R^2 评价拟合质量。
goodness_fitting__ <- function(y_orignal, y_fitted) {
  # 保留 Python 参数名 y_orignal 的原始拼写，避免接口漂移。
  1 - sum((y_fitted - y_orignal)^2) /
    sum((y_orignal - mean(y_orignal))^2)
}

# 对给定边界内的强度点进行 Gaussian 拟合，返回 R^2 和 sigma。
evaluate_gaussian_peak_on_intensity_list <- function(intensity_list,
                                                     height,
                                                     apex,
                                                     left,
                                                     right) {
  # 反向边界在 Python 中产生空数组，拟合失败且初始 sigma 为 NaN。
  if (right < left) return(list(0, NaN))

  # Python arange(left, right+1) 包含左右边界，扫描位置保持 0-based。
  x <- seq.int(left, right)
  y <- intensity_list[x + 1L]
  initial_sigma <- sqrt(mean((x - mean(x))^2))
  goodness <- 0
  sigma <- initial_sigma

  # curve_fit 在观测点少于三个参数数目时会失败并返回本函数的回退值。
  if (length(x) < 3L) return(list(goodness, sigma))

  # 拟合失败时保留 goodness=0 和初始 sigma，与 Python except 分支一致。
  coefficients <- tryCatch(
    .fit_gaussian(x, y, c(height, apex, initial_sigma)),
    error = function(e) NULL
  )
  if (!is.null(coefficients)) {
    fitted <- gaussian_function__(
      x,
      coefficients[["a"]],
      coefficients[["mu"]],
      coefficients[["sigma"]]
    )
    goodness <- goodness_fitting__(y, fitted)
    sigma <- unname(coefficients[["sigma"]])
  }

  list(goodness, sigma)
}

# 用 Gaussian 拟合估算给定区间内的峰面积。
get_gaussian_peakarea_on_intensity_list <- function(intensity_list,
                                                    left,
                                                    right) {
  # Python 会把 right 限制到 length，并使用 [left, right) 半开区间。
  right <- min(right, length(intensity_list))
  if (right <= left) return(0)
  x <- seq.int(left, right - 1L)
  y <- intensity_list[x + 1L]
  area <- 0

  # 以区间最高点、位置均值和总体标准差作为 curve_fit 初值。
  if (length(x) >= 3L) {
    initial_sigma <- sqrt(mean((x - mean(x))^2))
    coefficients <- tryCatch(
      .fit_gaussian(x, y, c(max(y), mean(x), initial_sigma)),
      error = function(e) NULL
    )

    # Gaussian 的全积分为 a * sqrt(2*pi*sigma^2)。
    if (!is.null(coefficients)) {
      area <- coefficients[["a"]] * sqrt(
        2 * pi * coefficients[["sigma"]]^2
      )
    }
  }

  unname(area)
}

# 使用 chromatograms.R 的 LOWESS 实现平滑整条 mass track。
lowess_smooth_track <- function(list_intensity, number_of_scans) {
  fraction <- 0.05

  # 短于 200 scans 时扩大局部窗口，但上限保持 0.8。
  if (number_of_scans < 200) {
    fraction <- min(10 / number_of_scans, 0.8)
  }

  as.numeric(smooth_lowess(list_intensity, frac = fraction))
}

# 快速检测是否恰好存在一个高且突出的洗脱峰。
quick_detect_unique_elution_peak <- function(
    intensity_track,
    min_peak_height = 100000,
    min_fwhm = 3,
    min_prominence_threshold_ratio = 0.2) {
  max_intensity <- max(intensity_track)
  prominence <- min_prominence_threshold_ratio * max_intensity
  unique_peak <- NULL

  # Python 使用严格的大于号，刚好等于最低峰高时不进入检测。
  if (max_intensity > min_peak_height) {
    detected <- .find_peaks_scipy(
      intensity_track,
      height = min_peak_height,
      prominence = prominence,
      width = min_fwhm
    )

    # 只有检测到且仅检测到一个峰时才返回 apex 和 height。
    if (length(detected$peaks) == 1L) {
      unique_peak <- list(
        apex = detected$peaks[[1L]],
        height = detected$properties$peak_heights[[1L]]
      )
    }
  }

  unique_peak
}

# 将按 RT 排序的峰分组，并清理相互重叠的峰。
check_overlap_peaks <- function(list_peaks) {
  if (length(list_peaks) < 2L) return(list_peaks)

  # 连续峰只要与当前簇最后一个峰重叠超过 3 scans，就归入同一簇。
  clusters <- list()
  current_cluster <- list(list_peaks[[1L]])
  for (peak in list_peaks[-1L]) {
    if (`_check_overlap`(peak, tail(current_cluster, 1L)[[1L]])) {
      current_cluster[[length(current_cluster) + 1L]] <- peak
    } else {
      clusters[[length(clusters) + 1L]] <- current_cluster
      current_cluster <- list(peak)
    }
  }
  clusters[[length(clusters) + 1L]] <- current_cluster

  # 每个簇分别执行双峰边界拆分或多峰合并。
  unlist(lapply(clusters, cleanup_peak_cluster), recursive = FALSE)
}

# 判断两个峰的边界重叠是否严格超过 3 scans。
`_check_overlap` <- function(peak1, peak2) {
  interval1 <- c(peak1$left_base, peak1$right_base)
  interval2 <- c(peak2$left_base, peak2$right_base)
  left_peak <- interval1
  right_peak <- interval2

  # Python 按 right_base 大小区分左峰与右峰，而不是按 apex 排序。
  if (interval1[[2L]] > interval2[[2L]]) {
    left_peak <- interval2
    right_peak <- interval1
  }
  overlap <- max(0, left_peak[[2L]] - right_peak[[1L]])
  overlap > 3
}

# 把三个及以上的重叠峰合并成一个峰。
`_merge_peak_cluster` <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) return(cluster_peaks[[1L]])

  # 选择 base 跨度最大的峰继承其它属性；并列时保留第一个。
  peak_sizes <- vapply(cluster_peaks, function(peak) {
    peak$right_base - peak$left_base
  }, numeric(1))
  largest <- cluster_peaks[[which.max(peak_sizes)]]

  # 合并后边界扩展为整个簇的最左和最右边界。
  largest$left_base <- min(vapply(
    cluster_peaks, `[[`, numeric(1), "left_base"
  ))
  largest$right_base <- max(vapply(
    cluster_peaks, `[[`, numeric(1), "right_base"
  ))
  largest
}

# 修正一个重叠峰簇的边界。
cleanup_peak_cluster <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) return(cluster_peaks)

  if (length(cluster_peaks) == 2L) {
    peak1 <- cluster_peaks[[1L]]
    peak2 <- cluster_peaks[[2L]]

    # 两峰时排序四个 base，并把前两个和后两个分别分配给两个峰。
    bases <- sort(c(
      peak1$left_base, peak1$right_base,
      peak2$left_base, peak2$right_base
    ))
    peak1$left_base <- bases[[1L]]
    peak1$right_base <- bases[[2L]]
    peak2$left_base <- bases[[3L]]
    peak2$right_base <- bases[[4L]]

    # 若原 apex 不再位于新边界内部，就用边界中点并按 int() 截断。
    if (!(peak1$left_base < peak1$apex && peak1$apex < peak1$right_base)) {
      peak1$apex <- .python_int(0.5 * (
        peak1$left_base + peak1$right_base
      ))
    }
    if (!(peak2$left_base < peak2$apex && peak2$apex < peak2$right_base)) {
      peak2$apex <- .python_int(0.5 * (
        peak2$left_base + peak2$right_base
      ))
    }

    return(list(peak1, peak2))
  }

  # 三个及以上重叠峰按 Python 规则合并为一个峰。
  list(`_merge_peak_cluster`(cluster_peaks))
}

# 在过短 ROI 的左右两端各补最多三个扫描点。
extend_ROI <- function(ROI, number_of_scans) {
  # 左侧只保留大于等于 0 的合法 0-based 扫描编号。
  left <- c(ROI[[1L]] - 3L, ROI[[1L]] - 2L, ROI[[1L]] - 1L)
  left <- left[left >= 0L]

  # 右侧只保留严格小于总扫描数的编号。
  last_scan <- tail(ROI, 1L)
  right <- c(last_scan + 1L, last_scan + 2L, last_scan + 3L)
  right <- right[right < number_of_scans]

  c(left, ROI, right)
}
