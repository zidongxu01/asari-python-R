# Functions corresponding to asari/mass_functions.py.
#
# This file should contain m/z distance checks, clustering, and seed-based
# splitting used by chromatograms.R.

# Flatten a list of two-element tuples and retain each value once.
#
# This follows the Python helper's traversal order: all first elements are
# collected before all second elements. Unlike Python's set(), unique() keeps
# the first occurrence, so the R result has deterministic ordering.
flatten_tuplelist <- function(tuple_list) {
  if (!is.list(tuple_list)) {
    stop("tuple_list must be a list.", call. = FALSE)
  }
  if (length(tuple_list) == 0L) {
    return(NULL)
  }
  if (any(vapply(tuple_list, length, integer(1)) < 2L)) {
    stop("Each tuple must contain at least two elements.", call. = FALSE)
  }

  first_elements <- lapply(tuple_list, `[[`, 1L)
  second_elements <- lapply(tuple_list, `[[`, 2L)
  unique(unlist(c(first_elements, second_elements), use.names = FALSE))
}

# Find adjacent sorted m/z values separated by less than a ppm tolerance.
#
# The Python function returns zero-based pairs (ii, ii - 1). This R port
# returns the corresponding one-based positions so callers can use each pair
# directly to index R vectors and lists.
check_close_mzs <- function(mzs, mz_tolerance_ppm = 5) {
  if (!is.numeric(mzs) || any(!is.finite(mzs))) {
    stop("mzs must be a finite numeric vector.", call. = FALSE)
  }
  if (length(mz_tolerance_ppm) != 1L ||
      !is.finite(mz_tolerance_ppm) ||
      mz_tolerance_ppm < 0) {
    stop(
      "mz_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }
  if (is.unsorted(mzs, strictly = FALSE)) {
    stop("mzs must be sorted in ascending order.", call. = FALSE)
  }
  if (length(mzs) < 2L) {
    return(list())
  }

  close_pairs <- list()
  for (ii in 2:length(mzs)) {
    tolerance <- mzs[[ii]] * mz_tolerance_ppm * 1e-6
    distance <- mzs[[ii]] - mzs[[ii - 1L]]
    if (distance < tolerance) {
      close_pairs[[length(close_pairs) + 1L]] <- c(ii, ii - 1L)
    }
  }

  close_pairs
}

nn_cluster_by_mz_seeds <- function(datatuples, mz_tolerance, presorted = FALSE) {
  stop("Not implemented yet: nn_cluster_by_mz_seeds")
}

identify_mass_peaks <- function(datatuples, mz_tolerance) {
  stop("Not implemented yet: identify_mass_peaks")
}
