source("r-prototype/R/mass_tracks.R")

smooth_moving_average <- function(values, size) {
  size <- max(1L, as.integer(size))
  if (size <= 1L || length(values) <= 2L) {
    return(values)
  }

  half <- floor(size / 2)
  smoothed <- numeric(length(values))
  for (i in seq_along(values)) {
    left <- max(1L, i - half)
    right <- min(length(values), i + half)
    smoothed[[i]] <- mean(values[left:right])
  }

  smoothed
}

estimate_track_noise <- function(intensity, min_intensity_threshold = 1000) {
  positive <- intensity[intensity > 0]
  if (length(positive) == 0L) {
    return(min_intensity_threshold)
  }

  bottom <- positive[positive <= stats::quantile(positive, 0.25)]
  max(min_intensity_threshold, stats::sd(bottom), na.rm = TRUE)
}

asari_peak_parameters <- function(parameters) {
  if (is.null(parameters$min_prominence_threshold)) {
    parameters$min_prominence_threshold <- 0.33 * parameters$min_peak_height
  }
  parameters$min_fwhm <- round(0.5 * parameters$min_timepoints)
  parameters$min_prominence_ratio <- 0.02
  parameters
}

audit_mass_track <- function(
  intensity,
  min_fwhm,
  min_intensity_threshold = 1000,
  min_peak_height = 100000,
  min_peak_ratio = 0.001
) {
  intensity <- as.numeric(intensity)
  low <- min_intensity_threshold
  high <- 1e8
  baseline <- low
  noise <- low
  scaling_factor <- 1
  max_intensity <- max(intensity)

  if (max_intensity > high) {
    scaling_factor <- max_intensity / high
    intensity <- intensity / scaling_factor
  }

  median_intensity <- stats::median(intensity)
  if (median_intensity > low) {
    active <- intensity[intensity > low]
    if (length(active) > length(intensity) * 0.5 && median_intensity > 10 * min_peak_height) {
      intensity <- as.numeric(stats::resid(stats::lm(intensity ~ seq_along(intensity))))
    }

    bottom <- intensity[intensity < low + stats::quantile(intensity, 0.25)]
    if (length(bottom) > 0L) {
      baseline <- mean(bottom)
      noise <- stats::sd(bottom)
    }
  }

  baseline <- max(baseline, low, na.rm = TRUE)
  noise <- max(noise, low, na.rm = TRUE)

  if (scaling_factor > 1 || max_intensity > 1000 * min_peak_height) {
    intensity <- smooth_moving_average(intensity, size = min_fwhm * 3L)
  } else if (100 * noise > max_intensity || max_intensity < 10 * min_peak_height) {
    intensity <- smooth_moving_average(intensity, size = min_fwhm + 2L)
  }

  cleaned <- intensity - baseline
  adjusted_min_peak_height <- max(min_peak_height, max_intensity * min_peak_ratio) - baseline

  list(
    baseline = baseline,
    noise = noise,
    scaling_factor = scaling_factor,
    min_peak_height = adjusted_min_peak_height,
    intensity = cleaned
  )
}

find_rois <- function(intensity, threshold, min_fwhm, gap_allowed = 2) {
  selected <- which(intensity > threshold)
  if (length(selected) == 0L) {
    return(list())
  }

  rois <- list()
  current <- selected[[1]]
  for (scan in selected[-1]) {
    if ((scan - current[length(current)]) <= gap_allowed + 1L) {
      current <- c(current, seq.int(current[length(current)] + 1L, scan))
    } else {
      if (length(current) >= min_fwhm + 2L) {
        rois[[length(rois) + 1L]] <- unique(current)
      }
      current <- scan
    }
  }

  if (length(current) >= min_fwhm + 2L) {
    rois[[length(rois) + 1L]] <- unique(current)
  }

  rois
}

extend_roi <- function(roi, number_of_scans) {
  left <- integer()
  right <- integer()
  if (roi[[1]] > 1L) {
    left <- seq.int(max(1L, roi[[1]] - 3L), roi[[1]] - 1L)
  }
  if (roi[[length(roi)]] < number_of_scans) {
    right <- seq.int(roi[[length(roi)]] + 1L, min(number_of_scans, roi[[length(roi)]] + 3L))
  }
  unique(c(left, roi, right))
}

local_peak_indices <- function(values, min_height, min_distance = 3) {
  if (length(values) < 3L) {
    return(integer())
  }

  candidates <- which(
    values >= min_height &
      c(-Inf, values[-length(values)]) < values &
      values >= c(values[-1], -Inf)
  )

  if (length(candidates) <= 1L) {
    return(candidates)
  }

  candidates <- candidates[order(values[candidates], decreasing = TRUE)]
  kept <- integer()
  for (candidate in candidates) {
    if (length(kept) == 0L || all(abs(candidate - kept) >= min_distance)) {
      kept <- c(kept, candidate)
    }
  }

  sort(kept)
}

local_maxima_asari <- function(values) {
  n <- length(values)
  if (n < 3L) {
    return(integer())
  }

  peaks <- integer()
  i <- 2L
  while (i < n) {
    if (values[[i - 1L]] < values[[i]]) {
      plateau_start <- i
      while (i < n && values[[i]] == values[[i + 1L]]) {
        i <- i + 1L
      }
      plateau_end <- i
      if (i < n && values[[i]] > values[[i + 1L]]) {
        peaks <- c(peaks, as.integer(floor((plateau_start + plateau_end) / 2)))
      }
    }
    i <- i + 1L
  }

  peaks
}

select_by_peak_distance <- function(peaks, values, distance) {
  if (length(peaks) <= 1L) {
    return(peaks)
  }

  distance <- max(1L, as.integer(distance))
  ordered <- peaks[order(values[peaks], decreasing = TRUE)]
  kept <- integer()
  for (peak in ordered) {
    if (length(kept) == 0L || all(abs(peak - kept) >= distance)) {
      kept <- c(kept, peak)
    }
  }

  sort(kept)
}

peak_prominence_asari <- function(values, peak, wlen = NULL) {
  peak_height <- values[[peak]]
  if (is.null(wlen)) {
    left_edge <- 1L
    right_edge <- length(values)
  } else {
    half_window <- max(1L, as.integer(floor(wlen / 2)))
    left_edge <- max(1L, peak - half_window)
    right_edge <- min(length(values), peak + half_window)
  }

  left_bound <- left_edge
  if (peak > left_edge) {
    for (i in seq.int(peak - 1L, left_edge)) {
      if (values[[i]] > peak_height) {
        left_bound <- i + 1L
        break
      }
      left_bound <- i
    }
  }

  right_bound <- right_edge
  if (peak < right_edge) {
    for (i in seq.int(peak + 1L, right_edge)) {
      if (values[[i]] > peak_height) {
        right_bound <- i - 1L
        break
      }
      right_bound <- i
    }
  }

  left_range <- seq.int(left_bound, peak)
  right_range <- seq.int(peak, right_bound)
  left_base <- left_range[[which.min(values[left_range])]]
  right_base <- right_range[[which.min(values[right_range])]]
  left_min <- values[[left_base]]
  right_min <- values[[right_base]]

  list(
    prominence = peak_height - max(left_min, right_min),
    left_base = left_base,
    right_base = right_base
  )
}

interpolated_crossing <- function(values, from, to, height) {
  if (from == to) {
    return(from)
  }

  step <- if (to > from) 1L else -1L
  current <- from
  while (current != to) {
    next_index <- current + step
    if (values[[next_index]] <= height) {
      y1 <- values[[current]]
      y2 <- values[[next_index]]
      if (y1 == y2) {
        return(next_index)
      }
      fraction <- (height - y1) / (y2 - y1)
      return(current + fraction * step)
    }
    current <- next_index
  }

  to
}

peak_width_asari <- function(values, peak, prominence, left_base, right_base, rel_height = 0.5) {
  width_height <- values[[peak]] - prominence * rel_height
  left_ip <- interpolated_crossing(values, peak, left_base, width_height)
  right_ip <- interpolated_crossing(values, peak, right_base, width_height)
  list(width = right_ip - left_ip, width_height = width_height, left_ip = left_ip, right_ip = right_ip)
}

find_peaks_asari <- function(values, height, distance, prominence, width, wlen = NULL) {
  peaks <- local_maxima_asari(values)
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(),
      prominences = numeric(),
      left_bases = integer(),
      right_bases = integer(),
      widths = numeric(),
      width_heights = numeric(),
      left_ips = numeric(),
      right_ips = numeric()
    )))
  }

  peaks <- peaks[values[peaks] >= height]
  peaks <- select_by_peak_distance(peaks, values, distance)
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(),
      prominences = numeric(),
      left_bases = integer(),
      right_bases = integer(),
      widths = numeric(),
      width_heights = numeric(),
      left_ips = numeric(),
      right_ips = numeric()
    )))
  }

  properties <- list(
    peak_heights = values[peaks],
    prominences = numeric(length(peaks)),
    left_bases = integer(length(peaks)),
    right_bases = integer(length(peaks)),
    widths = numeric(length(peaks)),
    width_heights = numeric(length(peaks)),
    left_ips = numeric(length(peaks)),
    right_ips = numeric(length(peaks))
  )

  keep <- logical(length(peaks))
  for (i in seq_along(peaks)) {
    prom <- peak_prominence_asari(values, peaks[[i]], wlen = wlen)
    width_info <- peak_width_asari(
      values,
      peak = peaks[[i]],
      prominence = prom$prominence,
      left_base = prom$left_base,
      right_base = prom$right_base
    )

    properties$prominences[[i]] <- prom$prominence
    properties$left_bases[[i]] <- prom$left_base
    properties$right_bases[[i]] <- prom$right_base
    properties$widths[[i]] <- width_info$width
    properties$width_heights[[i]] <- width_info$width_height
    properties$left_ips[[i]] <- width_info$left_ip
    properties$right_ips[[i]] <- width_info$right_ip
    keep[[i]] <- prom$prominence >= prominence && width_info$width >= width
  }

  for (name in names(properties)) {
    properties[[name]] <- properties[[name]][keep]
  }
  peaks <- peaks[keep]

  list(peaks = peaks, properties = properties)
}

peak_bounds <- function(values, apex_index, floor_value) {
  left <- apex_index
  while (left > 1L && values[[left - 1L]] > floor_value && values[[left - 1L]] <= values[[left]]) {
    left <- left - 1L
  }

  right <- apex_index
  while (right < length(values) && values[[right + 1L]] > floor_value && values[[right + 1L]] <= values[[right]]) {
    right <- right + 1L
  }

  c(left = left, right = right)
}

peak_cselectivity <- function(cleaned_intensity, left, right, height) {
  threshold <- 0.5 * height
  peak_points <- cleaned_intensity[left:right]
  peak_count <- sum(peak_points > threshold)
  background_count <- sum(cleaned_intensity > threshold)
  if (background_count >= peak_count && peak_count > 0L) {
    return(peak_count / background_count)
  }

  0
}

peaks_cselectivity_stats <- function(intensity, peaks) {
  if (length(peaks) == 0L) {
    return(list(cselectivity = numeric(), peak_datapoints = integer()))
  }

  peak_datapoints <- unique(unlist(lapply(peaks, function(peak) {
    seq.int(peak$left_base, peak$right_base)
  })))
  peak_values <- intensity[peak_datapoints]

  cselectivity <- vapply(peaks, function(peak) {
    threshold <- 0.5 * peak$height
    peak_count <- sum(peak_values > threshold)
    background_count <- sum(intensity > threshold)
    if (background_count >= peak_count && peak_count > 0L) {
      peak_count / background_count
    } else {
      0
    }
  }, numeric(1))

  list(cselectivity = cselectivity, peak_datapoints = peak_datapoints)
}

compute_noise_by_flanks <- function(
  peak,
  intensity,
  noise_data_points,
  min_intensity_threshold,
  old_noise_level
) {
  padding <- 30L
  n_points <- 100L
  noise_data_points <- sort(c(noise_data_points, peak$apex))
  center <- match(peak$apex, noise_data_points)
  if (is.na(center)) {
    return(old_noise_level)
  }

  left_positions <- seq.int(
    max(1L, center - padding - n_points),
    max(0L, center - padding - 1L)
  )
  right_positions <- seq.int(
    min(length(noise_data_points) + 1L, center + padding + 1L),
    min(length(noise_data_points), center + n_points + padding)
  )
  left_positions <- left_positions[left_positions >= 1L & left_positions <= length(noise_data_points)]
  right_positions <- right_positions[right_positions >= 1L & right_positions <= length(noise_data_points)]
  use_index <- noise_data_points[c(left_positions, right_positions)]

  if (length(use_index) == 0L) {
    return(old_noise_level)
  }

  max(min_intensity_threshold, mean(intensity[use_index] + 1))
}

gaussian_function <- function(x, a, mu, sigma) {
  a * exp(-((x - mu)^2) / (2 * sigma^2))
}

goodness_fitting_score <- function(observed, fitted) {
  denominator <- sum((observed - mean(observed))^2)
  if (denominator <= 0) {
    return(0)
  }

  1 - (sum((fitted - observed)^2) / denominator)
}

evaluate_gaussian_peak <- function(values, apex_index, left_index, right_index, min_fwhm) {
  if ((right_index - left_index + 1L) < 3L) {
    return(list(goodness_fitting = 0, sigma = NA_real_, left = left_index, right = right_index))
  }

  x <- seq.int(left_index, right_index)
  y <- values[x]
  initial_a <- max(y)
  initial_mu <- apex_index
  initial_sigma <- max(stats::sd(x), min_fwhm, 1)

  objective <- function(par) {
    a <- par[[1]]
    mu <- par[[2]]
    sigma <- exp(par[[3]])
    fitted <- gaussian_function(x, a, mu, sigma)
    sum((fitted - y)^2)
  }

  fit <- try(
    stats::optim(
      par = c(initial_a, initial_mu, log(initial_sigma)),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 200)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error") || fit$convergence > 1L) {
    return(list(goodness_fitting = 0, sigma = NA_real_, left = left_index, right = right_index))
  }

  a <- fit$par[[1]]
  mu <- fit$par[[2]]
  sigma <- exp(fit$par[[3]])
  fitted <- gaussian_function(x, a, mu, sigma)
  goodness <- goodness_fitting_score(y, fitted)

  if (!is.finite(sigma)) {
    return(list(goodness_fitting = 0, sigma = NA_real_, left = left_index, right = right_index))
  }

  halfwidth <- floor(max(abs(sigma), min_fwhm) * 6)
  halfwidth <- min(halfwidth, length(values))
  adjusted_left <- max(left_index, apex_index - halfwidth)
  adjusted_right <- min(right_index, apex_index + halfwidth)

  list(
    goodness_fitting = goodness,
    sigma = sigma,
    left = adjusted_left,
    right = adjusted_right
  )
}

evaluate_roi_peak_json <- function(
  ii,
  intensity_roi,
  rt_numbers_roi,
  peaks,
  properties,
  peakshape,
  min_fwhm
) {
  left_index <- properties$left_bases[[ii]]
  right_index <- properties$right_bases[[ii]]
  gaussian_fit <- evaluate_gaussian_peak(
    intensity_roi,
    apex_index = peaks[[ii]],
    left_index = left_index,
    right_index = right_index,
    min_fwhm = min_fwhm
  )

  if (gaussian_fit$goodness_fitting <= peakshape) {
    return(NULL)
  }

  left_base <- rt_numbers_roi[[gaussian_fit$left]]
  right_base <- rt_numbers_roi[[gaussian_fit$right]]
  peak_area <- sum(intensity_roi[gaussian_fit$left:gaussian_fit$right])

  list(
    apex = rt_numbers_roi[[peaks[[ii]]]],
    peak_area = as.integer(peak_area),
    height = as.integer(properties$peak_heights[[ii]]),
    left_base = left_base,
    right_base = right_base,
    goodness_fitting = gaussian_fit$goodness_fitting
  )
}

detect_evaluate_peaks_on_roi <- function(
  intensity_roi,
  rt_numbers_roi,
  min_peak_height,
  min_peak_ratio,
  min_fwhm,
  min_prominence_threshold,
  wlen,
  snr,
  peakshape,
  min_prominence_ratio,
  noise_level
) {
  max_intensity <- max(intensity_roi)
  min_prominence_threshold <- max(
    min_prominence_threshold,
    max_intensity * min_prominence_ratio,
    0.33 * min_peak_height,
    2 * noise_level
  )

  peak_result <- find_peaks_asari(
    intensity_roi,
    height = min_peak_height,
    distance = min_fwhm,
    prominence = min_prominence_threshold,
    width = min_fwhm,
    wlen = wlen
  )

  if (length(peak_result$peaks) == 0L && length(intensity_roi) > wlen) {
    peak_result <- find_peaks_asari(
      intensity_roi,
      height = min_peak_height,
      distance = min_fwhm,
      prominence = min_prominence_threshold,
      width = min_fwhm,
      wlen = NULL
    )
  }

  if (length(peak_result$peaks) == 0L) {
    return(list())
  }

  peaks <- list()
  for (ii in seq_along(peak_result$peaks)) {
    if ((peak_result$properties$right_bases[[ii]] - peak_result$properties$left_bases[[ii]]) >= min_fwhm + 2L) {
      peak <- evaluate_roi_peak_json(
        ii,
        intensity_roi,
        rt_numbers_roi,
        peak_result$peaks,
        peak_result$properties,
        peakshape,
        min_fwhm
      )
      if (!is.null(peak)) {
        peaks[[length(peaks) + 1L]] <- peak
      }
    }
  }

  check_overlap_peaks(peaks)
}

check_overlap <- function(peak1, peak2) {
  tuple1 <- c(peak1$left_base, peak1$right_base)
  tuple2 <- c(peak2$left_base, peak2$right_base)
  left_peak <- tuple1
  right_peak <- tuple2
  if (tuple1[[2]] > tuple2[[2]]) {
    left_peak <- tuple2
    right_peak <- tuple1
  }
  overlap <- max(0, left_peak[[2]] - right_peak[[1]])
  overlap > 3
}

merge_peak_cluster <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) {
    return(cluster_peaks[[1]])
  }

  peak_sizes <- vapply(cluster_peaks, function(peak) peak$right_base - peak$left_base, numeric(1))
  largest <- cluster_peaks[[which.max(peak_sizes)]]
  largest$left_base <- min(vapply(cluster_peaks, function(peak) peak$left_base, numeric(1)))
  largest$right_base <- max(vapply(cluster_peaks, function(peak) peak$right_base, numeric(1)))
  largest
}

cleanup_peak_cluster <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) {
    return(cluster_peaks)
  }

  if (length(cluster_peaks) == 2L) {
    peak1 <- cluster_peaks[[1]]
    peak2 <- cluster_peaks[[2]]
    bases <- sort(c(peak1$left_base, peak1$right_base, peak2$left_base, peak2$right_base))
    peak1$left_base <- bases[[1]]
    peak1$right_base <- bases[[2]]
    peak2$left_base <- bases[[3]]
    peak2$right_base <- bases[[4]]
    if (!(peak1$left_base < peak1$apex && peak1$apex < peak1$right_base)) {
      peak1$apex <- as.integer(0.5 * (peak1$left_base + peak1$right_base))
    }
    if (!(peak2$left_base < peak2$apex && peak2$apex < peak2$right_base)) {
      peak2$apex <- as.integer(0.5 * (peak2$left_base + peak2$right_base))
    }
    return(list(peak1, peak2))
  }

  list(merge_peak_cluster(cluster_peaks))
}

check_overlap_peaks <- function(peaks) {
  if (length(peaks) < 2L) {
    return(peaks)
  }

  clusters <- list()
  current <- list(peaks[[1]])
  for (peak in peaks[-1]) {
    if (check_overlap(peak, current[[length(current)]])) {
      current[[length(current) + 1L]] <- peak
    } else {
      clusters[[length(clusters) + 1L]] <- current
      current <- list(peak)
    }
  }
  clusters[[length(clusters) + 1L]] <- current

  unlist(lapply(clusters, cleanup_peak_cluster), recursive = FALSE)
}

iter_peak_detection_parameters <- function(list_mass_tracks, number_of_scans, parameters) {
  parameters <- asari_peak_parameters(parameters)
  Filter(
    Negate(is.null),
    lapply(list_mass_tracks, function(track) {
      if (any(track$intensity > parameters$min_intensity_threshold)) {
        list(track, number_of_scans, parameters)
      } else {
        NULL
      }
    })
  )
}

stats_detect_elution_peaks <- function(job) {
  mass_track <- job[[1]]
  number_of_scans <- job[[2]]
  parameters <- asari_peak_parameters(job[[3]])

  min_peak_height <- parameters$min_peak_height
  min_peak_ratio <- parameters$min_peak_ratio
  min_fwhm <- parameters$min_fwhm
  min_intensity_threshold <- parameters$min_intensity_threshold
  wlen <- parameters$wlen
  snr <- parameters$signal_noise_ratio
  peakshape <- parameters$gaussian_shape
  min_prominence_ratio <- parameters$min_prominence_ratio
  min_prominence_threshold <- parameters$min_prominence_threshold

  intensity <- mass_track$intensity
  audited <- audit_mass_track(
    mass_track$intensity,
    min_fwhm = min_fwhm,
    min_intensity_threshold = min_intensity_threshold,
    min_peak_height = min_peak_height,
    min_peak_ratio = min_peak_ratio
  )

  cleaned_intensity <- audited$intensity
  noise <- audited$noise
  min_peak_height <- audited$min_peak_height
  list_scans <- seq_len(number_of_scans)
  rois <- find_rois(cleaned_intensity, noise, min_fwhm = min_fwhm)
  if (length(rois) == 0L) {
    return(list())
  }

  list_json_peaks <- list()
  for (roi in rois) {
    if (length(roi) < 3L * min_fwhm) {
      roi <- extend_roi(roi, length(cleaned_intensity))
    }

    peaks <- detect_evaluate_peaks_on_roi(
      cleaned_intensity[roi],
      roi,
      min_peak_height,
      min_peak_ratio,
      min_fwhm,
      min_prominence_threshold,
      wlen,
      snr,
      peakshape,
      min_prominence_ratio,
      noise
    )
    list_json_peaks <- c(list_json_peaks, peaks)
  }

  if (length(list_json_peaks) == 0L) {
    return(list())
  }

  selectivity <- peaks_cselectivity_stats(cleaned_intensity, list_json_peaks)
  noise_data_points <- setdiff(list_scans, selectivity$peak_datapoints)

  list_peaks <- list()
  for (ii in seq_along(list_json_peaks)) {
    peak <- list_json_peaks[[ii]]
    peak$cSelectivity <- selectivity$cselectivity[[ii]]
    peak$parent_masstrack_id <- mass_track$id_number
    peak$mz <- mass_track$mz
    peak$height <- audited$scaling_factor * peak$height + audited$baseline
    height_for_snr <- min(peak$height, 1e100)
    peak$snr <- as.integer(height_for_snr / compute_noise_by_flanks(
      peak,
      intensity,
      noise_data_points,
      min_intensity_threshold,
      noise
    ))
    if (peak$snr >= snr) {
      peak$height <- as.integer(peak$height)
      list_peaks[[length(list_peaks) + 1L]] <- peak
    }
  }

  list_peaks
}

batch_deep_detect_elution_peaks <- function(list_mass_tracks, number_of_scans, parameters) {
  jobs <- iter_peak_detection_parameters(list_mass_tracks, number_of_scans, parameters)
  unlist(lapply(jobs, stats_detect_elution_peaks), recursive = FALSE)
}

feature_list_to_table <- function(feature_list, rt_table) {
  if (length(feature_list) == 0L) {
    return(data.frame())
  }

  rows <- lapply(feature_list, function(peak) {
    data.frame(
      mz = peak$mz,
      rtime = rt_table$rtime[[peak$apex]],
      rtime_left_base = rt_table$rtime[[peak$left_base]],
      rtime_right_base = rt_table$rtime[[peak$right_base]],
      parent_masstrack_id = peak$parent_masstrack_id,
      peak_area = peak$peak_area,
      cSelectivity = peak$cSelectivity,
      goodness_fitting = peak$goodness_fitting,
      snr = peak$snr,
      detection_counts = 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

detect_peaks <- function(
  mass_track_result,
  parameters
) {
  tracks <- mass_track_result$tracks
  if (length(tracks) == 0L) {
    return(empty_detected_feature_table())
  }

  list_mass_tracks <- vector("list", length(tracks))
  for (i in seq_along(tracks)) {
    list_mass_tracks[[i]] <- tracks[[i]]
    list_mass_tracks[[i]]$id_number <- i - 1L
  }

  feature_list <- batch_deep_detect_elution_peaks(
    list_mass_tracks,
    number_of_scans = nrow(mass_track_result$rt_table),
    parameters = parameters
  )

  if (length(feature_list) == 0L) {
    return(empty_detected_feature_table())
  }

  features <- feature_list_to_table(feature_list, mass_track_result$rt_table)
  features <- features[order(features$mz, features$rtime), , drop = FALSE]
  features$id_number <- paste0("F", seq_len(nrow(features)) - 1L)
  features <- features[, c(
    "id_number",
    "mz",
    "rtime",
    "rtime_left_base",
    "rtime_right_base",
    "parent_masstrack_id",
    "peak_area",
    "cSelectivity",
    "goodness_fitting",
    "snr",
    "detection_counts"
  )]

  rownames(features) <- NULL
  features
}

empty_detected_feature_table <- function() {
  data.frame(
    id_number = character(),
    mz = numeric(),
    rtime = numeric(),
    rtime_left_base = numeric(),
    rtime_right_base = numeric(),
    parent_masstrack_id = integer(),
    peak_area = numeric(),
    cSelectivity = numeric(),
    goodness_fitting = numeric(),
    snr = numeric(),
    detection_counts = integer(),
    stringsAsFactors = FALSE
  )
}
