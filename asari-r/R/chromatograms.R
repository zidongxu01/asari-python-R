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
