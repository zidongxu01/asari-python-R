# Functions corresponding to asari/chromatograms.py.
#
# This file should contain mass track / EIC construction and RT alignment logic.
# The first implementation target is extract_mass_tracks().

extract_mass_tracks <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop(
      "The mzR package is required to read mzML files.",
      call. = FALSE
    )
  }
  if (!file.exists(infile)) {
    stop("mzML file does not exist: ", infile, call. = FALSE)
  }

  ms_expt <- mzR::openMSfile(infile)
  on.exit(mzR::close(ms_expt), add = TRUE)

  scan_header <- mzR::header(ms_expt)

  # pymzML reports scan time in minutes and the Python implementation converts
  # it to seconds. mzR's retentionTime column is already expressed in seconds.
  rt_times <- numeric()
  ms2_spectra <- list()
  ms1_rows <- list()
  ms1_scan_number <- 0L

  # Match Python's acquisition_time value: Unix seconds, or NULL if the mzML
  # start timestamp is unavailable or cannot be parsed.
  acquisition_time <- tryCatch({
    start_time <- mzR::runInfo(ms_expt)$startTimeStamp
    parsed_time <- as.POSIXct(
      start_time,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )
    if (length(parsed_time) != 1L || is.na(parsed_time)) {
      NULL
    } else {
      as.numeric(parsed_time)
    }
  }, error = function(e) NULL)

  for (scan_index in seq_len(nrow(scan_header))) {
    ms_level <- scan_header$msLevel[[scan_index]]

    if (identical(ms_level, 1L)) {
      rt_times <- c(rt_times, scan_header$retentionTime[[scan_index]])
      peak_matrix <- mzR::peaks(ms_expt, scan_index)

      if (nrow(peak_matrix) > 0L) {
        # Python first casts MS1 intensities to int and then applies a strict
        # greater-than filter. trunc() reproduces Python's conversion for the
        # non-negative intensity values found in mass spectra.
        intensities <- trunc(peak_matrix[, "intensity"])
        good_positions <- intensities > min_intensity

        if (any(good_positions)) {
          ms1_rows[[length(ms1_rows) + 1L]] <- data.frame(
            mz = peak_matrix[good_positions, "mz"],
            scan_num = rep.int(ms1_scan_number, sum(good_positions)),
            intensity = intensities[good_positions],
            stringsAsFactors = FALSE
          )
        }
      }

      # scan_num follows Python and is zero-based among MS1 scans only.
      ms1_scan_number <- ms1_scan_number + 1L
    } else if (identical(ms_level, 2L)) {
      precursor_mz <- scan_header$precursorMZ[[scan_index]]
      if (is.na(precursor_mz)) {
        next
      }

      peak_matrix <- mzR::peaks(ms_expt, scan_index)
      if (nrow(peak_matrix) == 0L) {
        next
      }

      keep <- peak_matrix[, "intensity"] >= min_intensity &
        peak_matrix[, "mz"] < precursor_mz - 1
      if (any(keep)) {
        ms2_spectra[[length(ms2_spectra) + 1L]] <- list(
          precursor_mz = precursor_mz,
          rtime = scan_header$retentionTime[[scan_index]],
          peaks = data.frame(
            mz = peak_matrix[keep, "mz"],
            intensity = peak_matrix[keep, "intensity"],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (length(ms1_rows) == 0L) {
    all_data <- data.frame(
      mz = numeric(),
      scan_num = integer(),
      intensity = numeric()
    )
  } else {
    all_data <- do.call(rbind, ms1_rows)
    rownames(all_data) <- NULL
  }

  # Python's mzTree is keyed by int(mz * 1000). For non-negative m/z values,
  # floor() gives the same thousandth-m/z bin.
  if (nrow(all_data) == 0L) {
    mz_tree <- list()
  } else {
    thousandth_mz <- floor(all_data$mz * 1000)
    mz_tree <- split(all_data, thousandth_mz, drop = TRUE)
  }

  rt_numbers <- if (length(rt_times) == 0L) {
    integer()
  } else {
    seq.int(0L, length(rt_times) - 1L)
  }
  rt_length <- length(rt_numbers)

  good_bins <- if (length(mz_tree) == 0L) {
    list()
  } else {
    get_thousandth_bins(
      mz_tree,
      mz_tolerance_ppm,
      min_timepoints,
      min_peak_height
    )
  }

  tracks <- list()
  for (bin in good_bins) {
    tracks <- c(
      tracks,
      bin_to_mass_tracks(bin, rt_length, mz_tolerance_ppm)
    )
  }

  # Track objects use list(mz = ..., intensity = ...). The helper
  # check_close_mzs() will return one-based index pairs for R.
  if (length(tracks) > 1L) {
    track_mzs <- vapply(tracks, function(track) track$mz, numeric(1))
    tracks <- tracks[order(track_mzs)]
    track_mzs <- vapply(tracks, function(track) track$mz, numeric(1))

    tracks_to_merge <- check_close_mzs(track_mzs, mz_tolerance_ppm)
    merged <- list()
    to_remove <- integer()

    for (pair in tracks_to_merge) {
      a <- pair[[1L]]
      b <- pair[[2L]]
      merged[[length(merged) + 1L]] <- merge_two_mass_tracks(
        tracks[[a]],
        tracks[[b]]
      )
      to_remove <- c(to_remove, a, b)
    }

    tracks <- c(tracks[setdiff(seq_along(tracks), to_remove)], merged)
  }

  list(
    rt_numbers = rt_numbers,
    rt_times = rt_times,
    tracks = tracks,
    ms2_spectra = ms2_spectra,
    acquisition_time = acquisition_time
  )
}

extract_single_track_full_rt_length <- function(bin, rt_length) {
  stop("Not implemented yet: extract_single_track_full_rt_length")
}

bin_to_mass_tracks <- function(bin_data_tuples, rt_length, mz_tolerance_ppm = 5) {
  stop("Not implemented yet: bin_to_mass_tracks")
}

get_thousandth_bins <- function(mz_tree,
                                mz_tolerance_ppm = 5,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  stop("Not implemented yet: get_thousandth_bins")
}

rt_lowess_calibration <- function(good_landmark_peaks,
                                  selected_reference_landmark_peaks,
                                  sample_rt_numbers,
                                  reference_rt_numbers,
                                  num_iterations,
                                  sample_name,
                                  outdir) {
  stop("Not implemented yet: rt_lowess_calibration")
}

remap_intensity_track <- function(intensity_track, new, rt_cal_dict) {
  stop("Not implemented yet: remap_intensity_track")
}

smooth_moving_average <- function(list_intensity, size = 9) {
  stats::filter(list_intensity, rep(1 / size, size), sides = 2)
}
