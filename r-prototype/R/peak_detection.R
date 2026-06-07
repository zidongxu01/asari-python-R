source("r-prototype/R/mass_tracks.R")

estimate_track_noise <- function(intensity, min_intensity_threshold = 1000) {
  positive <- intensity[intensity > 0]
  if (length(positive) == 0L) {
    return(min_intensity_threshold)
  }

  bottom <- positive[positive <= stats::quantile(positive, 0.25)]
  max(min_intensity_threshold, stats::sd(bottom), na.rm = TRUE)
}

find_rois <- function(intensity, threshold, min_timepoints = 6, gap_allowed = 2) {
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
      if (length(current) >= min_timepoints) {
        rois[[length(rois) + 1L]] <- unique(current)
      }
      current <- scan
    }
  }

  if (length(current) >= min_timepoints) {
    rois[[length(rois) + 1L]] <- unique(current)
  }

  rois
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

detect_peaks_on_track <- function(
  track,
  track_id,
  rt_table,
  min_peak_height = 100000,
  min_timepoints = 6,
  min_intensity_threshold = 1000,
  signal_noise_ratio = 2,
  min_peak_ratio = 0.001
) {
  intensity <- track$intensity
  if (max(intensity) < min_peak_height) {
    return(data.frame())
  }

  noise <- estimate_track_noise(intensity, min_intensity_threshold)
  threshold <- max(noise, min_intensity_threshold)
  rois <- find_rois(intensity, threshold, min_timepoints = min_timepoints)
  if (length(rois) == 0L) {
    return(data.frame())
  }

  rows <- list()
  min_height <- max(min_peak_height, max(intensity) * min_peak_ratio)
  min_distance <- max(1L, round(0.5 * min_timepoints))

  for (roi in rois) {
    roi_values <- intensity[roi]
    local_peaks <- local_peak_indices(roi_values, min_height = min_height, min_distance = min_distance)

    for (local_apex in local_peaks) {
      bounds <- peak_bounds(roi_values, local_apex, floor_value = threshold)
      local_left <- bounds[["left"]]
      local_right <- bounds[["right"]]
      if ((local_right - local_left + 1L) < min_distance + 1L) {
        next
      }

      global_apex <- roi[[local_apex]]
      global_left <- roi[[local_left]]
      global_right <- roi[[local_right]]
      peak_height <- intensity[[global_apex]]
      snr <- floor(peak_height / max(noise, min_intensity_threshold))
      if (snr < signal_noise_ratio) {
        next
      }

      peak_points <- intensity[global_left:global_right]
      rows[[length(rows) + 1L]] <- data.frame(
        mz = track$mz,
        rtime = rt_table$rtime[[global_apex]],
        rtime_left_base = rt_table$rtime[[global_left]],
        rtime_right_base = rt_table$rtime[[global_right]],
        parent_masstrack_id = track_id,
        peak_area = sum(peak_points),
        cSelectivity = 1,
        goodness_fitting = NA_real_,
        snr = snr,
        detection_counts = 1L,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }

  do.call(rbind, rows)
}

detect_peaks <- function(
  mass_track_result,
  min_peak_height = 100000,
  min_timepoints = 6,
  min_intensity_threshold = 1000,
  signal_noise_ratio = 2,
  min_peak_ratio = 0.001
) {
  tracks <- mass_track_result$tracks
  if (length(tracks) == 0L) {
    return(empty_detected_feature_table())
  }

  rows <- list()
  for (i in seq_along(tracks)) {
    track_features <- detect_peaks_on_track(
      tracks[[i]],
      track_id = i - 1L,
      rt_table = mass_track_result$rt_table,
      min_peak_height = min_peak_height,
      min_timepoints = min_timepoints,
      min_intensity_threshold = min_intensity_threshold,
      signal_noise_ratio = signal_noise_ratio,
      min_peak_ratio = min_peak_ratio
    )
    if (nrow(track_features) > 0L) {
      rows[[length(rows) + 1L]] <- track_features
    }
  }

  if (length(rows) == 0L) {
    return(empty_detected_feature_table())
  }

  features <- do.call(rbind, rows)
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
