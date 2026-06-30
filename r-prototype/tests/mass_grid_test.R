source("r-prototype/R/mass_grid.R")
source("r-prototype/R/parameters.R")

make_track <- function(mz, id_number) {
  list(
    id_number = id_number,
    mz = mz,
    intensity = c(0, 10, 100, 10, 0)
  )
}

parameters <- default_parameters()
parameters$mz_tolerance_ppm <- 5
parameters$project_sample_number_small <- 10

sample_mass_tracks <- list(
  sample_a = list(
    make_track(100.00000, 0L),
    make_track(200.00000, 1L)
  ),
  sample_b = list(
    make_track(100.00030, 0L),
    make_track(300.00000, 1L)
  )
)

grid <- construct_mass_grid(sample_mass_tracks, parameters)

if (!identical(names(grid), c("mz", "sample_a", "sample_b"))) {
  stop("Expected MassGrid columns to match asari layout.")
}

if (nrow(grid) != 3L) {
  stop("Expected three MassGrid rows.")
}

matched_row <- which(!is.na(grid$sample_a) & !is.na(grid$sample_b))
if (length(matched_row) != 1L) {
  stop("Expected one matched m/z row.")
}

if (grid$sample_a[[matched_row]] != 0L || grid$sample_b[[matched_row]] != 0L) {
  stop("Expected matched row to preserve track ids.")
}

unmatched_a <- grid[!is.na(grid$sample_a) & is.na(grid$sample_b), , drop = FALSE]
unmatched_b <- grid[is.na(grid$sample_a) & !is.na(grid$sample_b), , drop = FALSE]
if (nrow(unmatched_a) != 1L || nrow(unmatched_b) != 1L) {
  stop("Expected one unmatched row per sample.")
}

centroid_parameters <- parameters
centroid_parameters$project_sample_number_small <- 1
centroid_grid <- construct_mass_grid(sample_mass_tracks, centroid_parameters)

if (nrow(centroid_grid) != 3L) {
  stop("Expected centroid MassGrid to preserve same row count.")
}

if (length(which(!is.na(centroid_grid$sample_a) & !is.na(centroid_grid$sample_b))) != 1L) {
  stop("Expected centroid MassGrid to align one shared m/z row.")
}

conflict_tracks <- list(
  sample_a = list(
    make_track(100.00000, 0L),
    make_track(100.00010, 1L)
  ),
  sample_b = list(
    make_track(100.00020, 0L)
  )
)

conflict_grid <- build_grid_by_centroiding(conflict_tracks, parameters)
if (anyDuplicated(stats::na.omit(conflict_grid$sample_a)) > 0L) {
  stop("Expected one track per sample per MassGrid row.")
}

cat("mass grid test passed.\n")
