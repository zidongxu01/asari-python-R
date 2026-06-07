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
  current <- points[1, , drop = FALSE]

  for (i in seq.int(2, nrow(points))) {
    row <- points[i, , drop = FALSE]
    center_mz <- stats::median(current$mz)
    if (abs(row$mz - center_mz) < ppm_tolerance(row$mz, mz_tolerance_ppm)) {
      current <- rbind(current, row)
    } else {
      groups[[length(groups) + 1L]] <- current
      current <- row
    }
  }

  groups[[length(groups) + 1L]] <- current
  groups
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

  points$mz_bin <- floor(points$mz * 1000)
  bin_ids <- sort(unique(points$mz_bin))
  rough_bins <- list()
  current_ids <- bin_ids[[1]]

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
  rough_bins[[length(rough_bins) + 1L]] <- current_ids

  tracks <- list()
  for (bin in rough_bins) {
    bin_points <- points[points$mz_bin %in% bin, , drop = FALSE]
    if (max(bin_points$intensity) < min_peak_height) {
      next
    }
    if (!has_consecutive_scans(bin_points$scan_index, min_timepoints = min_timepoints)) {
      next
    }

    split_bins <- split_points_by_ppm(bin_points, mz_tolerance_ppm = mz_tolerance_ppm)
    for (split_bin in split_bins) {
      if (max(split_bin$intensity) < min_peak_height) {
        next
      }
      if (!has_consecutive_scans(split_bin$scan_index, min_timepoints = min_timepoints)) {
        next
      }
      tracks[[length(tracks) + 1L]] <- build_mass_track(split_bin, rt_table)
    }
  }

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

