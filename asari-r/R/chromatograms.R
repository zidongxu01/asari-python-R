# Functions corresponding to asari/chromatograms.py.
#
# This file should contain mass track / EIC construction and RT alignment logic.
# The first implementation target is extract_mass_tracks().

extract_mass_tracks <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  stop("Not implemented yet: extract_mass_tracks")
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
