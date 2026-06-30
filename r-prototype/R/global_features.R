source("r-prototype/R/composite_tracks.R")
source("r-prototype/R/peak_detection.R")
source("r-prototype/R/feature_tables.R")

get_peak_area_sum <- function(track_intensity, left_base, right_base) {
  left_base <- as.integer(left_base)
  right_base <- as.integer(right_base)
  if (is.na(left_base) || is.na(right_base)) {
    return(0)
  }

  left_base <- max(1L, left_base)
  right_base <- min(length(track_intensity), right_base)
  if (left_base > right_base) {
    return(0)
  }

  sum(track_intensity[left_base:right_base])
}

get_peak_area_auc <- function(track_intensity, left_base, right_base) {
  left_base <- as.integer(left_base)
  right_base <- as.integer(right_base)
  if (is.na(left_base) || is.na(right_base)) {
    return(0)
  }

  left_base <- max(1L, left_base)
  right_base <- min(length(track_intensity), right_base)
  if (left_base > right_base) {
    return(0)
  }

  values <- track_intensity[left_base:right_base]
  if (length(values) <= 1L) {
    return(sum(values))
  }

  sum(pmax(values, c(values[-1], 0)))
}

get_peak_area_gaussian <- function(track_intensity, left_base, right_base) {
  left_base <- as.integer(left_base)
  right_base <- as.integer(right_base)
  if (is.na(left_base) || is.na(right_base)) {
    return(0)
  }

  left_base <- max(1L, left_base)
  right_base <- min(length(track_intensity), right_base)
  if ((right_base - left_base + 1L) < 3L) {
    return(0)
  }

  values <- track_intensity[left_base:right_base]
  apex_index <- which.max(values)
  fit <- evaluate_gaussian_peak(
    values,
    apex_index = apex_index,
    left_index = 1L,
    right_index = length(values),
    min_fwhm = 1L
  )
  if (is.na(fit$sigma) || fit$goodness_fitting <= 0) {
    return(0)
  }

  height <- max(values)
  as.integer(height * sqrt(2 * pi * fit$sigma^2))
}

peak_area_function <- function(parameters) {
  methods <- list(
    sum = get_peak_area_sum,
    auc = get_peak_area_auc,
    gauss = get_peak_area_gaussian
  )

  method <- parameters$peak_area
  if (is.null(method) || !method %in% names(methods)) {
    method <- "sum"
  }
  methods[[method]]
}

reverse_rt_base <- function(rt_cal_dict, base) {
  if (is.null(rt_cal_dict) || length(rt_cal_dict) == 0L) {
    return(base)
  }

  mapped <- rt_cal_dict[[as.character(base)]]
  if (is.null(mapped)) {
    base
  } else {
    mapped
  }
}

extract_features_per_sample <- function(
  sample_name,
  sample_tracks,
  mass_grid,
  feature_list,
  peak_area_fn,
  reverse_rt_cal_dict = NULL
) {
  f_list <- numeric(length(feature_list))

  for (i in seq_along(feature_list)) {
    peak <- feature_list[[i]]
    track_number <- mass_grid[[sample_name]][[peak$parent_masstrack_id + 1L]]
    if (is.na(track_number)) {
      f_list[[i]] <- 0
      next
    }

    mass_track <- track_by_id(sample_tracks, track_number)
    if (is.null(mass_track)) {
      f_list[[i]] <- 0
      next
    }

    left_base <- reverse_rt_base(reverse_rt_cal_dict, peak$left_base)
    right_base <- reverse_rt_base(reverse_rt_cal_dict, peak$right_base)
    f_list[[i]] <- peak_area_fn(mass_track$intensity, left_base, right_base)
  }

  f_list
}

generate_feature_table <- function(
  feature_list,
  mass_grid,
  sample_mass_tracks,
  parameters,
  rt_table = NULL,
  reverse_rt_cal_dicts = NULL,
  drop_unaligned_samples = FALSE,
  aligned_samples = NULL
) {
  if (is.null(rt_table)) {
    rt_table <- default_rt_table(if (length(feature_list) == 0L) 0L else max(vapply(feature_list, function(peak) peak$right_base, integer(1))))
  }

  feature_table <- feature_list_to_global_table(feature_list, rt_table)
  if (nrow(feature_table) == 0L) {
    for (sample_name in names(sample_mass_tracks)) {
      feature_table[[sample_name]] <- numeric()
    }
    return(feature_table)
  }

  if (is.null(reverse_rt_cal_dicts)) {
    reverse_rt_cal_dicts <- list()
  }
  if (is.null(aligned_samples)) {
    aligned_samples <- names(sample_mass_tracks)
  }

  peak_area_fn <- peak_area_function(parameters)
  for (sample_name in names(sample_mass_tracks)) {
    if (drop_unaligned_samples && !sample_name %in% aligned_samples) {
      next
    }

    feature_table[[sample_name]] <- extract_features_per_sample(
      sample_name = sample_name,
      sample_tracks = sample_mass_tracks[[sample_name]],
      mass_grid = mass_grid,
      feature_list = feature_list,
      peak_area_fn = peak_area_fn,
      reverse_rt_cal_dict = reverse_rt_cal_dicts[[sample_name]]
    )
  }

  feature_table
}

default_rt_table <- function(number_of_scans) {
  data.frame(
    scan_id = seq_len(number_of_scans),
    scan_index = seq_len(number_of_scans),
    rtime = seq_len(number_of_scans),
    stringsAsFactors = FALSE
  )
}

feature_list_to_global_table <- function(feature_list, rt_table) {
  if (length(feature_list) == 0L) {
    return(empty_detected_feature_table())
  }

  features <- feature_list_to_table(feature_list, rt_table)
  features <- features[order(features$mz, features$rtime), , drop = FALSE]
  features$id_number <- paste0("F", seq_len(nrow(features)) - 1L)
  features <- features[, c(
    "id_number",
    "mz",
    "rtime",
    "rtime_left_base",
    "rtime_right_base",
    "parent_masstrack_id",
    "peak_area",
    "cSelectivity",
    "goodness_fitting",
    "snr",
    "detection_counts"
  )]
  rownames(features) <- NULL
  features
}

global_peak_detection <- function(composite_tracks, number_of_scans, parameters, rt_table = NULL) {
  if (is.null(rt_table)) {
    rt_table <- default_rt_table(number_of_scans)
  }

  feature_list <- batch_deep_detect_elution_peaks(
    composite_tracks,
    number_of_scans = number_of_scans,
    parameters = parameters
  )

  list(
    feature_list = feature_list,
    feature_table = feature_list_to_global_table(feature_list, rt_table)
  )
}

process_global_features <- function(
  sample_mass_tracks,
  parameters,
  rt_table = NULL,
  reference_sample = NULL,
  rt_cal_dicts = NULL,
  drop_unaligned_samples = FALSE,
  aligned_samples = NULL
) {
  sample_mass_tracks <- normalize_sample_mass_tracks(sample_mass_tracks)
  mass_grid <- construct_mass_grid(
    sample_mass_tracks,
    parameters = parameters,
    reference_sample = reference_sample
  )

  rt_length <- composite_rt_length(sample_mass_tracks)
  composite_tracks <- build_composite_tracks(
    mass_grid,
    sample_mass_tracks,
    rt_length = rt_length,
    rt_cal_dicts = rt_cal_dicts,
    drop_unaligned_samples = drop_unaligned_samples,
    aligned_samples = aligned_samples
  )

  peak_result <- global_peak_detection(
    composite_tracks,
    number_of_scans = rt_length,
    parameters = parameters,
    rt_table = rt_table
  )

  feature_table <- generate_feature_table(
    peak_result$feature_list,
    mass_grid = mass_grid,
    sample_mass_tracks = sample_mass_tracks,
    parameters = parameters,
    rt_table = rt_table,
    reverse_rt_cal_dicts = rt_cal_dicts,
    drop_unaligned_samples = drop_unaligned_samples,
    aligned_samples = aligned_samples
  )
  sample_names <- if (drop_unaligned_samples) {
    intersect(names(sample_mass_tracks), aligned_samples)
  } else {
    names(sample_mass_tracks)
  }
  feature_tables <- prepare_feature_tables(
    feature_table,
    sample_names = sample_names
  )

  list(
    mass_grid = mass_grid,
    composite_tracks = composite_tracks,
    feature_list = peak_result$feature_list,
    feature_table = feature_table,
    full_feature_table = feature_tables$full,
    preferred_feature_table = feature_tables$preferred,
    global_feature_table = peak_result$feature_table
  )
}
