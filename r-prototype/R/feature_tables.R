asari_feature_columns <- function(sample_names) {
  c(
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
    "detection_counts",
    sample_names
  )
}

compute_detection_counts <- function(feature_table, sample_names) {
  if (length(sample_names) == 0L || nrow(feature_table) == 0L) {
    return(integer(nrow(feature_table)))
  }

  rowSums(feature_table[, sample_names, drop = FALSE] > 1)
}

round_asari_feature_columns <- function(feature_table) {
  if (nrow(feature_table) == 0L) {
    return(feature_table)
  }

  round_columns <- list(
    mz = 4,
    rtime = 2,
    rtime_left_base = 2,
    rtime_right_base = 2,
    cSelectivity = 2,
    goodness_fitting = 2
  )

  for (column in names(round_columns)) {
    if (column %in% names(feature_table)) {
      feature_table[[column]] <- round(feature_table[[column]], round_columns[[column]])
    }
  }

  feature_table
}

prepare_full_feature_table <- function(feature_table, sample_names) {
  feature_table$detection_counts <- compute_detection_counts(feature_table, sample_names)
  feature_table <- feature_table[, asari_feature_columns(sample_names), drop = FALSE]
  round_asari_feature_columns(feature_table)
}

prepare_preferred_feature_table <- function(
  full_feature_table,
  snr = 2,
  peak_shape = 0.7,
  cselectivity = 0.7
) {
  if (nrow(full_feature_table) == 0L) {
    return(full_feature_table)
  }

  preferred <- full_feature_table[full_feature_table$detection_counts > 0, , drop = FALSE]
  preferred <- preferred[
    preferred$snr > snr &
      preferred$goodness_fitting > peak_shape &
      preferred$cSelectivity > cselectivity,
    ,
    drop = FALSE
  ]
  rownames(preferred) <- NULL
  preferred
}

prepare_feature_tables <- function(
  feature_table,
  sample_names,
  snr = 2,
  peak_shape = 0.7,
  cselectivity = 0.7
) {
  full <- prepare_full_feature_table(feature_table, sample_names)
  preferred <- prepare_preferred_feature_table(
    full,
    snr = snr,
    peak_shape = peak_shape,
    cselectivity = cselectivity
  )

  list(
    full = full,
    preferred = preferred
  )
}
