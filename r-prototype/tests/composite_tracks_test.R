source("r-prototype/R/composite_tracks.R")

make_track <- function(mz, id_number, intensity) {
  list(
    id_number = id_number,
    mz = mz,
    intensity = intensity
  )
}

sample_mass_tracks <- list(
  sample_a = list(
    make_track(100.0000, 0L, c(0, 10, 50, 10, 0)),
    make_track(200.0000, 1L, c(0, 5, 10, 5, 0))
  ),
  sample_b = list(
    make_track(100.0003, 0L, c(0, 20, 70, 20, 0)),
    make_track(300.0000, 1L, c(0, 1, 2, 1, 0))
  )
)

mass_grid <- data.frame(
  mz = c(100.00015, 200.0000, 300.0000),
  sample_a = c(0L, 1L, NA_integer_),
  sample_b = c(0L, NA_integer_, 1L)
)

composite_tracks <- build_composite_tracks(
  mass_grid,
  sample_mass_tracks,
  rt_length = 5
)

if (length(composite_tracks) != 3L) {
  stop("Expected one composite track per MassGrid row.")
}

if (!identical(composite_tracks[[1]]$intensity, c(0, 30, 120, 30, 0))) {
  stop("Expected shared m/z row to sum sample tracks.")
}

if (!identical(composite_tracks[[2]]$intensity, c(0, 5, 10, 5, 0))) {
  stop("Expected row missing sample_b to keep sample_a track only.")
}

if (!identical(composite_tracks[[3]]$intensity, c(0, 1, 2, 1, 0))) {
  stop("Expected row missing sample_a to keep sample_b track only.")
}

if (composite_tracks[[1]]$id_number != 0L || composite_tracks[[2]]$id_number != 1L) {
  stop("Expected composite track ids to follow MassGrid row order.")
}

shifted <- remap_intensity_track(
  c(10, 20, 30),
  base_track = numeric(4),
  rt_cal_dict = list("1" = 2L, "2" = 3L, "3" = 4L)
)

if (!identical(shifted, c(0, 10, 20, 30))) {
  stop("Expected remap_intensity_track to use RT mapping dictionary.")
}

cat("composite tracks test passed.\n")
