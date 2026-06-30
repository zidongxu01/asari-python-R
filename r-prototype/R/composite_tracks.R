source("r-prototype/R/mass_grid.R")

track_by_id <- function(tracks, track_id) {
  for (track in tracks) {
    if (identical(as.integer(track$id_number), as.integer(track_id))) {
      return(track)
    }
  }

  NULL
}

remap_intensity_track <- function(intensity, base_track, rt_cal_dict = NULL) {
  remapped <- base_track
  if (is.null(rt_cal_dict) || length(rt_cal_dict) == 0L) {
    limit <- min(length(intensity), length(remapped))
    remapped[seq_len(limit)] <- intensity[seq_len(limit)]
    return(remapped)
  }

  for (source_index in seq_along(intensity)) {
    target_index <- rt_cal_dict[[as.character(source_index)]]
    if (is.null(target_index)) {
      target_index <- source_index
    }
    if (!is.na(target_index) && target_index >= 1L && target_index <= length(remapped)) {
      remapped[[target_index]] <- remapped[[target_index]] + intensity[[source_index]]
    }
  }

  remapped
}

composite_rt_length <- function(sample_mass_tracks, rt_length = NULL) {
  if (!is.null(rt_length)) {
    return(rt_length)
  }

  lengths <- unlist(lapply(sample_mass_tracks, function(tracks) {
    vapply(tracks, function(track) length(track$intensity), integer(1))
  }))

  if (length(lengths) == 0L) {
    return(0L)
  }

  max(lengths)
}

build_composite_tracks <- function(
  mass_grid,
  sample_mass_tracks,
  rt_length = NULL,
  rt_cal_dicts = NULL,
  drop_unaligned_samples = FALSE,
  aligned_samples = NULL
) {
  sample_mass_tracks <- normalize_sample_mass_tracks(sample_mass_tracks)
  sample_names <- setdiff(names(mass_grid), "mz")
  rt_length <- composite_rt_length(sample_mass_tracks, rt_length = rt_length)
  base_track <- numeric(rt_length)

  if (is.null(rt_cal_dicts)) {
    rt_cal_dicts <- list()
  }
  if (is.null(aligned_samples)) {
    aligned_samples <- sample_names
  }

  composite_tracks <- vector("list", nrow(mass_grid))
  for (row_index in seq_len(nrow(mass_grid))) {
    composite_intensity <- base_track

    for (sample_name in sample_names) {
      if (drop_unaligned_samples && !sample_name %in% aligned_samples) {
        next
      }

      track_id <- mass_grid[[sample_name]][[row_index]]
      if (is.na(track_id)) {
        next
      }

      track <- track_by_id(sample_mass_tracks[[sample_name]], track_id)
      if (is.null(track)) {
        next
      }

      sample_track <- remap_intensity_track(
        track$intensity,
        base_track = numeric(rt_length),
        rt_cal_dict = rt_cal_dicts[[sample_name]]
      )
      composite_intensity <- composite_intensity + sample_track
    }

    composite_tracks[[row_index]] <- list(
      id_number = row_index - 1L,
      mz = mass_grid$mz[[row_index]],
      intensity = composite_intensity
    )
  }

  composite_tracks
}
