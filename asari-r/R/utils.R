# Shared helper functions for asariR.

ppm_to_mz_tolerance <- function(mz, ppm) {
  mz * ppm * 1e-6
}
