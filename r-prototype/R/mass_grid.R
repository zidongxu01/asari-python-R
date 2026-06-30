source("r-prototype/R/mass_tracks.R")

normalize_sample_mass_tracks <- function(sample_mass_tracks) {
  if (is.null(names(sample_mass_tracks)) || any(names(sample_mass_tracks) == "")) {
    stop("sample_mass_tracks must be a named list.", call. = FALSE)
  }

  lapply(sample_mass_tracks, function(sample) {
    tracks <- if (!is.null(sample$tracks)) sample$tracks else sample
    lapply(seq_along(tracks), function(i) {
      track <- tracks[[i]]
      if (is.null(track$id_number)) {
        track$id_number <- i - 1L
      }
      track
    })
  })
}

track_mz_tuples <- function(sample_mass_tracks) {
  sample_names <- names(sample_mass_tracks)
  rows <- list()

  for (sample_index in seq_along(sample_mass_tracks)) {
    tracks <- sample_mass_tracks[[sample_index]]
    for (track in tracks) {
      rows[[length(rows) + 1L]] <- data.frame(
        mz = track$mz,
        track_id = track$id_number,
        sample_index = sample_index,
        sample_name = sample_names[[sample_index]],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      mz = numeric(),
      track_id = integer(),
      sample_index = integer(),
      sample_name = character(),
      stringsAsFactors = FALSE
    ))
  }

  tuples <- do.call(rbind, rows)
  tuples[order(tuples$mz), , drop = FALSE]
}

bin_mz_range_ok <- function(bin, mz_tolerance_ppm) {
  if (nrow(bin) <= 1L) {
    return(TRUE)
  }

  mz_tolerance <- bin$mz[[1]] * mz_tolerance_ppm * 1e-6
  (max(bin$mz) - min(bin$mz)) < 2 * mz_tolerance
}

bin_has_unique_samples <- function(bin) {
  length(unique(bin$sample_index)) == nrow(bin)
}

split_bin_by_sample_conflicts <- function(bin, mz_tolerance_ppm) {
  clusters <- list()
  current <- bin[1, , drop = FALSE]

  if (nrow(bin) > 1L) {
    for (i in seq.int(2L, nrow(bin))) {
      row <- bin[i, , drop = FALSE]
      candidate <- rbind(current, row)
      mz_tolerance <- stats::median(candidate$mz) * mz_tolerance_ppm * 1e-6
      if (
        bin_has_unique_samples(candidate) &&
          (max(candidate$mz) - min(candidate$mz)) < 2 * mz_tolerance
      ) {
        current <- candidate
      } else {
        clusters[[length(clusters) + 1L]] <- current
        current <- row
      }
    }
  }

  clusters[[length(clusters) + 1L]] <- current
  clusters
}

bin_track_mzs <- function(track_tuples, mz_tolerance_ppm = 5) {
  if (nrow(track_tuples) == 0L) {
    return(list())
  }

  tuples <- track_tuples[order(track_tuples$mz), , drop = FALSE]
  tol <- mz_tolerance_ppm * 1e-6
  bins_of_bins <- list()
  current <- tuples[1, , drop = FALSE]

  if (nrow(tuples) > 1L) {
    for (i in seq.int(2L, nrow(tuples))) {
      delta <- tuples$mz[[i]] - tuples$mz[[i - 1L]]
      if (delta < tol * tuples$mz[[i - 1L]]) {
        current <- rbind(current, tuples[i, , drop = FALSE])
      } else {
        bins_of_bins[[length(bins_of_bins) + 1L]] <- current
        current <- tuples[i, , drop = FALSE]
      }
    }
  }
  bins_of_bins[[length(bins_of_bins) + 1L]] <- current

  good_bins <- list()
  for (bin in bins_of_bins) {
    split_bins <- if (bin_mz_range_ok(bin, mz_tolerance_ppm) && bin_has_unique_samples(bin)) {
      list(bin)
    } else {
      split_bin_by_sample_conflicts(bin, mz_tolerance_ppm)
    }

    for (split_bin in split_bins) {
      good_bins[[length(good_bins) + 1L]] <- list(
        mz = stats::median(split_bin$mz),
        members = split_bin
      )
    }
  }

  good_bins
}

mass_grid_from_bins <- function(bins, sample_names) {
  grid <- data.frame(mz = vapply(bins, function(bin) bin$mz, numeric(1)))
  for (sample_name in sample_names) {
    grid[[sample_name]] <- NA_integer_
  }

  for (i in seq_along(bins)) {
    members <- bins[[i]]$members
    for (j in seq_len(nrow(members))) {
      grid[[members$sample_name[[j]]]][[i]] <- members$track_id[[j]]
    }
  }

  grid[order(grid$mz), , drop = FALSE]
}

build_grid_by_centroiding <- function(sample_mass_tracks, parameters) {
  sample_mass_tracks <- normalize_sample_mass_tracks(sample_mass_tracks)
  bins <- bin_track_mzs(
    track_mz_tuples(sample_mass_tracks),
    mz_tolerance_ppm = parameters$mz_tolerance_ppm
  )

  if (length(bins) == 0L) {
    grid <- data.frame(mz = numeric())
    for (sample_name in names(sample_mass_tracks)) {
      grid[[sample_name]] <- integer()
    }
    return(grid)
  }

  mass_grid_from_bins(bins, names(sample_mass_tracks))
}

nearest_ppm_match <- function(mz, reference_mz, mz_tolerance_ppm) {
  deltas <- abs(reference_mz - mz)
  nearest <- which.min(deltas)
  if (length(nearest) == 0L) {
    return(NA_integer_)
  }

  if (deltas[[nearest]] < reference_mz[[nearest]] * mz_tolerance_ppm * 1e-6) {
    nearest
  } else {
    NA_integer_
  }
}

add_sample_to_mass_grid <- function(grid, sample_name, tracks, mz_tolerance_ppm) {
  grid[[sample_name]] <- NA_integer_
  reference_mz <- grid$mz
  new_rows <- list()

  for (track in tracks) {
    match_index <- nearest_ppm_match(track$mz, reference_mz, mz_tolerance_ppm)
    if (is.na(match_index)) {
      row <- as.list(rep(NA_integer_, ncol(grid)))
      names(row) <- names(grid)
      row$mz <- track$mz
      row[[sample_name]] <- track$id_number
      new_rows[[length(new_rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    } else if (is.na(grid[[sample_name]][[match_index]])) {
      grid$mz[[match_index]] <- mean(c(grid$mz[[match_index]], track$mz))
      grid[[sample_name]][[match_index]] <- track$id_number
    } else {
      row <- as.list(rep(NA_integer_, ncol(grid)))
      names(row) <- names(grid)
      row$mz <- track$mz
      row[[sample_name]] <- track$id_number
      new_rows[[length(new_rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }
  }

  if (length(new_rows) > 0L) {
    grid <- rbind(grid, do.call(rbind, new_rows))
  }

  grid[order(grid$mz), , drop = FALSE]
}

build_grid_sample_wise <- function(sample_mass_tracks, parameters, reference_sample = NULL) {
  sample_mass_tracks <- normalize_sample_mass_tracks(sample_mass_tracks)
  sample_names <- names(sample_mass_tracks)
  if (is.null(reference_sample)) {
    reference_sample <- sample_names[[1]]
  }

  reference_tracks <- sample_mass_tracks[[reference_sample]]
  grid <- data.frame(
    mz = vapply(reference_tracks, function(track) track$mz, numeric(1)),
    stringsAsFactors = FALSE
  )
  grid[[reference_sample]] <- vapply(reference_tracks, function(track) track$id_number, integer(1))

  for (sample_name in sample_names[sample_names != reference_sample]) {
    grid <- add_sample_to_mass_grid(
      grid,
      sample_name = sample_name,
      tracks = sample_mass_tracks[[sample_name]],
      mz_tolerance_ppm = parameters$mz_tolerance_ppm
    )
  }

  missing_samples <- setdiff(sample_names, names(grid))
  for (sample_name in missing_samples) {
    grid[[sample_name]] <- NA_integer_
  }

  grid[, c("mz", sample_names), drop = FALSE]
}

construct_mass_grid <- function(sample_mass_tracks, parameters, reference_sample = NULL) {
  sample_mass_tracks <- normalize_sample_mass_tracks(sample_mass_tracks)
  number_of_samples <- length(sample_mass_tracks)

  if (number_of_samples <= parameters$project_sample_number_small) {
    build_grid_sample_wise(
      sample_mass_tracks,
      parameters = parameters,
      reference_sample = reference_sample
    )
  } else {
    build_grid_by_centroiding(
      sample_mass_tracks,
      parameters = parameters
    )
  }
}
