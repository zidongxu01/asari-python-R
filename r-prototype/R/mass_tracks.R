source("r-prototype/R/mzml_reader.R")

ppm_tolerance <- function(mz, ppm) {
  mz * ppm * 1e-6
}

has_consecutive_scans <- function(scan_ids, min_timepoints = 5, gap_allowed = 2) {
  scan_ids <- sort(unique(scan_ids))
  if (length(scan_ids) < min_timepoints) {
    return(FALSE)
  }

  spans <- scan_ids[min_timepoints:length(scan_ids)] - scan_ids[1:(length(scan_ids) - min_timepoints + 1)]
  min(spans) <= gap_allowed + min_timepoints - 1
}

asari_consensus_mz <- function(mzs, intensities) {
  0.5 * (stats::median(mzs) + mzs[which.max(intensities)])
}

build_mass_track <- function(points, rt_table) {
  rt_length <- nrow(rt_table)
  intensity_track <- numeric(rt_length)

  for (i in seq_len(nrow(points))) {
    scan_index <- points$scan_index[[i]]
    intensity_track[[scan_index]] <- max(intensity_track[[scan_index]], points$intensity[[i]])
  }

  list(
    mz = asari_consensus_mz(points$mz, points$intensity),
    intensity = intensity_track,
    points = points
  )
}

split_points_by_ppm <- function(points, mz_tolerance_ppm = 5) {
  points <- points[order(points$mz), , drop = FALSE]
  if (nrow(points) == 0L) {
    return(list())
  }

  groups <- list()
  start <- 1L

  if (nrow(points) > 1L) {
    for (i in seq.int(2, nrow(points))) {
      center_mz <- stats::median(points$mz[start:(i - 1L)])
      if (abs(points$mz[[i]] - center_mz) >= ppm_tolerance(points$mz[[i]], mz_tolerance_ppm)) {
        groups[[length(groups) + 1L]] <- points[start:(i - 1L), , drop = FALSE]
        start <- i
      }
    }
  }

  groups[[length(groups) + 1L]] <- points[start:nrow(points), , drop = FALSE]
  groups
}

points_for_rough_bin <- function(points_by_mz_bin, bin) {
  chunks <- points_by_mz_bin[as.character(bin)]
  chunks <- chunks[!vapply(chunks, is.null, logical(1))]
  if (length(chunks) == 0L) {
    return(data.frame())
  }
  if (length(chunks) == 1L) {
    return(chunks[[1]])
  }

  do.call(rbind, chunks)
}

rough_consecutive_scan_check <- function(points, min_timepoints = 5, gap_allowed = 2) {
  if (nrow(points) >= 4L * min_timepoints) {
    return(TRUE)
  }

  scan_ids <- sort(points$scan_index)
  if (length(scan_ids) < min_timepoints) {
    return(FALSE)
  }

  spans <- scan_ids[min_timepoints:length(scan_ids)] - scan_ids[1:(length(scan_ids) - min_timepoints + 1L)]
  min(spans) <= gap_allowed + min_timepoints - 1L
}

smooth_counts_nearest <- function(values, size) {
  size <- max(2L, as.integer(size))
  half_left <- floor((size - 1L) / 2L)
  half_right <- size - half_left - 1L
  smoothed <- numeric(length(values))

  for (i in seq_along(values)) {
    left <- max(1L, i - half_left)
    right <- min(length(values), i + half_right)
    window <- values[left:right]
    if (left == 1L && length(window) < size) {
      window <- c(rep(values[[1]], size - length(window)), window)
    }
    if (right == length(values) && length(window) < size) {
      window <- c(window, rep(values[[length(values)]], size - length(window)))
    }
    smoothed[[i]] <- mean(window)
  }

  smoothed
}

local_peak_indices_for_mz <- function(values, min_height = 0, min_distance = 1) {
  if (length(values) < 3L) {
    return(integer())
  }

  candidates <- which(values[-c(1, length(values))] >= values[-c(length(values) - 1L, length(values))] &
    values[-c(1, length(values))] >= values[-c(1, 2)] &
    values[-c(1, length(values))] > min_height) + 1L
  if (length(candidates) <= 1L) {
    return(candidates)
  }

  ordered <- candidates[order(values[candidates], decreasing = TRUE)]
  kept <- integer()
  for (candidate in ordered) {
    if (length(kept) == 0L || all(abs(candidate - kept) >= min_distance)) {
      kept <- c(kept, candidate)
    }
  }

  sort(kept)
}

identify_mass_peaks <- function(points, mz_tolerance, presorted = TRUE) {
  mz4 <- as.integer(points$mz * 10000)
  if (!presorted) {
    mz4 <- sort(mz4)
  }
  if (length(mz4) == 0L) {
    return(numeric())
  }

  tol4 <- as.integer(mz_tolerance * 10000)
  size <- max(2L, as.integer(0.5 * tol4))
  positions <- seq.int(min(mz4), max(mz4))
  counts <- tabulate(match(mz4, positions), nbins = length(positions))
  smoothed <- smooth_counts_nearest(counts, size = size)
  peak_positions <- local_peak_indices_for_mz(
    smoothed,
    min_height = 0,
    min_distance = max(1L, tol4)
  )

  if (length(peak_positions) == 0L) {
    return(numeric())
  }

  0.0001 * positions[peak_positions]
}

gap_divide_mz_cluster <- function(points) {
  points <- points[order(points$mz), , drop = FALSE]
  if (nrow(points) < 2L) {
    return(list(points))
  }

  gaps <- diff(points$mz)
  split_at <- which.max(gaps) + 1L
  list(
    points[seq_len(split_at - 1L), , drop = FALSE],
    points[split_at:nrow(points), , drop = FALSE]
  )
}

nn_cluster_by_mz_seeds <- function(points, mz_tolerance, presorted = TRUE) {
  seeds <- identify_mass_peaks(points, mz_tolerance, presorted = presorted)
  if (length(seeds) == 0L) {
    return(gap_divide_mz_cluster(points))
  }

  assignments <- vapply(points$mz, function(mz) which.min(abs(mz - seeds)), integer(1))
  unname(split(points, assignments))
}

bin_to_mass_tracks <- function(points, rt_table, mz_tolerance_ppm = 5) {
  points <- points[order(points$mz), , drop = FALSE]
  mz_range <- points$mz[[nrow(points)]] - points$mz[[1]]
  mz_tolerance <- points$mz[[1]] * mz_tolerance_ppm * 1e-6

  if (mz_range < mz_tolerance * 2) {
    return(list(build_mass_track(points, rt_table)))
  }

  lapply(
    nn_cluster_by_mz_seeds(points, mz_tolerance, presorted = FALSE),
    build_mass_track,
    rt_table = rt_table
  )
}

merge_two_mass_tracks <- function(track1, track2) {
  list(
    mz = 0.5 * (track1$mz + track2$mz),
    intensity = track1$intensity + track2$intensity,
    points = rbind(track1$points, track2$points)
  )
}

merge_close_mass_tracks <- function(tracks, mz_tolerance_ppm = 5) {
  if (length(tracks) < 2L) {
    return(tracks)
  }

  tracks <- tracks[order(vapply(tracks, function(track) track$mz, numeric(1)))]
  merged <- list()
  i <- 1L
  while (i <= length(tracks)) {
    if (i < length(tracks)) {
      mz_delta <- tracks[[i + 1L]]$mz - tracks[[i]]$mz
      mz_tolerance <- tracks[[i + 1L]]$mz * mz_tolerance_ppm * 1e-6
      if (mz_delta < mz_tolerance) {
        merged[[length(merged) + 1L]] <- merge_two_mass_tracks(tracks[[i + 1L]], tracks[[i]])
        i <- i + 2L
        next
      }
    }

    merged[[length(merged) + 1L]] <- tracks[[i]]
    i <- i + 1L
  }

  merged
}

extract_mass_tracks_from_points <- function(
  points,
  mz_tolerance_ppm = 5,
  min_timepoints = 5,
  min_peak_height = 1000
) {
  if (nrow(points) == 0L) {
    return(list(rt_table = data.frame(), tracks = list(), track_table = data.frame()))
  }

  points$scan_index <- match(points$scan_id, sort(unique(points$scan_id)))
  rt_table <- unique(points[, c("scan_id", "scan_index", "rtime")])
  rt_table <- rt_table[order(rt_table$scan_index), , drop = FALSE]

  points$mz_bin <- as.integer(floor(points$mz * 1000))
  point_counts <- table(points$mz_bin)
  bin_ids <- sort(as.integer(as.numeric(names(point_counts[point_counts >= min_timepoints]))))
  if (length(bin_ids) == 0L) {
    return(list(rt_table = rt_table, tracks = list(), track_table = data.frame()))
  }

  rough_bins <- list()
  current_ids <- bin_ids[[1]]

  if (length(bin_ids) > 1L) {
    for (i in seq.int(2, length(bin_ids))) {
      previous <- bin_ids[[i - 1L]]
      current <- bin_ids[[i]]
      if ((current - previous) == 1L || (current - previous) < current * mz_tolerance_ppm * 1e-6) {
        current_ids <- c(current_ids, current)
      } else {
        rough_bins[[length(rough_bins) + 1L]] <- current_ids
        current_ids <- current
      }
    }
  }
  rough_bins[[length(rough_bins) + 1L]] <- current_ids

  points <- points[points$mz_bin %in% bin_ids, , drop = FALSE]
  rough_bin_lookup <- integer(length(bin_ids))
  names(rough_bin_lookup) <- as.character(bin_ids)
  for (i in seq_along(rough_bins)) {
    rough_bin_lookup[as.character(rough_bins[[i]])] <- i
  }
  points$rough_bin <- unname(rough_bin_lookup[as.character(points$mz_bin)])
  points_by_rough_bin <- split(points, points$rough_bin)

  tracks <- list()
  for (i in seq_along(rough_bins)) {
    bin_points <- points_by_rough_bin[[as.character(i)]]
    if (max(bin_points$intensity) < min_peak_height) {
      next
    }
    if (!rough_consecutive_scan_check(bin_points, min_timepoints = min_timepoints)) {
      next
    }

    bin_tracks <- bin_to_mass_tracks(bin_points, rt_table, mz_tolerance_ppm = mz_tolerance_ppm)
    for (track in bin_tracks) {
      tracks[[length(tracks) + 1L]] <- track
    }
  }

  tracks <- merge_close_mass_tracks(tracks, mz_tolerance_ppm = mz_tolerance_ppm)

  if (length(tracks) == 0L) {
    track_table <- data.frame(
      track_id = integer(),
      mz = numeric(),
      scan_id = integer(),
      rtime = numeric(),
      intensity = numeric()
    )
  } else {
    track_rows <- vector("list", length(tracks))
    for (i in seq_along(tracks)) {
      track <- tracks[[i]]
      nonzero <- which(track$intensity > 0)
      track_rows[[i]] <- data.frame(
        track_id = i - 1L,
        mz = track$mz,
        scan_id = rt_table$scan_id[nonzero],
        rtime = rt_table$rtime[nonzero],
        intensity = track$intensity[nonzero],
        stringsAsFactors = FALSE
      )
    }
    track_table <- do.call(rbind, track_rows)
  }

  list(
    rt_table = rt_table,
    tracks = tracks,
    track_table = track_table
  )
}

extract_mass_tracks <- function(
  file,
  mz_tolerance_ppm = 5,
  min_intensity = 1000,
  min_timepoints = 5,
  min_peak_height = 1000,
  max_scans = NULL
) {
  points <- read_mzml_points(
    file,
    max_scans = max_scans,
    min_intensity = min_intensity
  )

  extract_mass_tracks_from_points(
    points,
    mz_tolerance_ppm = mz_tolerance_ppm,
    min_timepoints = min_timepoints,
    min_peak_height = min_peak_height
  )
}
