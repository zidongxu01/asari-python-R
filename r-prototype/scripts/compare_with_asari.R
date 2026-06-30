source("r-prototype/R/parameters.R")
source("r-prototype/R/global_features.R")
source("r-prototype/R/preview_utils.R")

usage <- function() {
  paste(
    "Usage:",
    "  Rscript r-prototype/scripts/compare_with_asari.R path/to/file.mzML",
    "  Rscript r-prototype/scripts/compare_with_asari.R path/to/file.mzML --python-output path/to/asari_output",
    "",
    "Options:",
    "  --python-output DIR  Directory containing Python asari output tables.",
    "  --max-scans N        Only read the first N scans for a quick check.",
    "  --mz-ppm N           Matching tolerance for mz, default 10.",
    "  --rt-sec N           Matching tolerance for retention time, default 2.",
    sep = "\n"
  )
}

parse_args <- function(args) {
  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat(usage(), "\n")
    quit(status = if (length(args) == 0L) 1L else 0L)
  }

  options <- list(
    mzml_file = NULL,
    python_output = NULL,
    max_scans = NULL,
    mz_ppm = 10,
    rt_sec = 2
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--python-output") {
      i <- i + 1L
      if (i > length(args)) stop("--python-output requires a directory.", call. = FALSE)
      options$python_output <- args[[i]]
    } else if (arg == "--max-scans") {
      i <- i + 1L
      if (i > length(args)) stop("--max-scans requires a number.", call. = FALSE)
      options$max_scans <- as.integer(args[[i]])
    } else if (arg == "--mz-ppm") {
      i <- i + 1L
      if (i > length(args)) stop("--mz-ppm requires a number.", call. = FALSE)
      options$mz_ppm <- as.numeric(args[[i]])
    } else if (arg == "--rt-sec") {
      i <- i + 1L
      if (i > length(args)) stop("--rt-sec requires a number.", call. = FALSE)
      options$rt_sec <- as.numeric(args[[i]])
    } else if (is.null(options$mzml_file)) {
      options$mzml_file <- arg
    } else {
      stop("Unexpected argument: ", arg, "\n\n", usage(), call. = FALSE)
    }
    i <- i + 1L
  }

  if (is.null(options$mzml_file)) {
    stop("Missing mzML file path.\n\n", usage(), call. = FALSE)
  }
  if (!file.exists(options$mzml_file)) {
    stop("Input mzML file does not exist: ", options$mzml_file, call. = FALSE)
  }
  if (!is.null(options$python_output) && !dir.exists(options$python_output)) {
    stop("Python output directory does not exist: ", options$python_output, call. = FALSE)
  }
  if (!is.null(options$max_scans) && (is.na(options$max_scans) || options$max_scans <= 0L)) {
    stop("--max-scans must be a positive integer.", call. = FALSE)
  }

  options$mzml_file <- normalizePath(options$mzml_file, mustWork = TRUE)
  if (!is.null(options$python_output)) {
    options$python_output <- normalizePath(options$python_output, mustWork = TRUE)
  }
  options
}

table_path <- function(root, filename) {
  direct <- file.path(root, filename)
  export <- file.path(root, "export", filename)
  if (file.exists(direct)) {
    return(direct)
  }
  if (file.exists(export)) {
    return(export)
  }

  matches <- list.files(root, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
  if (length(matches) == 0L) {
    return(NULL)
  }
  matches[[1]]
}

read_feature_table <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  utils::read.delim(path, check.names = FALSE)
}

count_feature_filters <- function(full_table, snr = 2, peak_shape = 0.7, cselectivity = 0.7) {
  if (is.null(full_table) || nrow(full_table) == 0L) {
    return(c(
      detection_counts = 0L,
      snr = 0L,
      goodness_fitting = 0L,
      cSelectivity = 0L,
      all_preferred = 0L
    ))
  }

  pass_detection <- full_table$detection_counts > 0
  pass_snr <- full_table$snr > snr
  pass_shape <- full_table$goodness_fitting > peak_shape
  pass_selectivity <- full_table$cSelectivity > cselectivity

  c(
    detection_counts = sum(pass_detection),
    snr = sum(pass_snr),
    goodness_fitting = sum(pass_shape),
    cSelectivity = sum(pass_selectivity),
    all_preferred = sum(pass_detection & pass_snr & pass_shape & pass_selectivity)
  )
}

feature_match_summary <- function(python_table, r_table, mz_ppm = 10, rt_sec = 2) {
  if (is.null(python_table) || is.null(r_table)) {
    return(NULL)
  }
  if (nrow(python_table) == 0L || nrow(r_table) == 0L) {
    return(list(
      python_rows = nrow(python_table),
      r_rows = nrow(r_table),
      matched = 0L,
      unmatched_python = nrow(python_table),
      unmatched_r = nrow(r_table),
      mz_ppm = numeric(),
      rt_delta = numeric(),
      area_ratio = numeric(),
      missing_python = python_table,
      extra_r = r_table
    ))
  }

  used_r <- rep(FALSE, nrow(r_table))
  matched_python <- rep(FALSE, nrow(python_table))
  matched_r <- rep(FALSE, nrow(r_table))
  mz_errors <- numeric()
  rt_errors <- numeric()
  area_ratios <- numeric()

  for (i in seq_len(nrow(python_table))) {
    mz_tolerance <- python_table$mz[[i]] * mz_ppm * 1e-6
    candidates <- which(
      !used_r &
        abs(r_table$mz - python_table$mz[[i]]) <= mz_tolerance &
        abs(r_table$rtime - python_table$rtime[[i]]) <= rt_sec
    )
    if (length(candidates) == 0L) {
      next
    }

    score <- abs(r_table$mz[candidates] - python_table$mz[[i]]) / mz_tolerance +
      abs(r_table$rtime[candidates] - python_table$rtime[[i]]) / rt_sec
    best <- candidates[[which.min(score)]]
    used_r[[best]] <- TRUE
    matched_python[[i]] <- TRUE
    matched_r[[best]] <- TRUE

    mz_errors <- c(mz_errors, abs(r_table$mz[[best]] - python_table$mz[[i]]) / python_table$mz[[i]] * 1e6)
    rt_errors <- c(rt_errors, abs(r_table$rtime[[best]] - python_table$rtime[[i]]))
    if ("peak_area" %in% names(python_table) && "peak_area" %in% names(r_table) && python_table$peak_area[[i]] > 0) {
      area_ratios <- c(area_ratios, r_table$peak_area[[best]] / python_table$peak_area[[i]])
    }
  }

  matched <- length(mz_errors)
  list(
    python_rows = nrow(python_table),
    r_rows = nrow(r_table),
    matched = matched,
    unmatched_python = nrow(python_table) - matched,
    unmatched_r = nrow(r_table) - matched,
    mz_ppm = mz_errors,
    rt_delta = rt_errors,
    area_ratio = area_ratios,
    missing_python = python_table[!matched_python, , drop = FALSE],
    extra_r = r_table[!matched_r, , drop = FALSE]
  )
}

print_quantiles <- function(label, values) {
  if (length(values) == 0L) {
    cat(label, "none\n")
    return(invisible(NULL))
  }

  qs <- stats::quantile(values, probs = c(0.5, 0.95, 1), na.rm = TRUE, names = FALSE)
  cat(label, "median", round(qs[[1]], 4), "p95", round(qs[[2]], 4), "max", round(qs[[3]], 4), "\n")
}

print_match_summary <- function(label, summary) {
  if (is.null(summary)) {
    cat(label, "Python table missing, skipped.\n")
    return(invisible(NULL))
  }

  cat(label, "\n")
  cat("  Python rows:", summary$python_rows, "\n")
  cat("  R rows:", summary$r_rows, "\n")
  cat("  matched:", summary$matched, "\n")
  cat("  unmatched Python:", summary$unmatched_python, "\n")
  cat("  unmatched R:", summary$unmatched_r, "\n")
  cat("  ")
  print_quantiles("mz ppm error:", summary$mz_ppm)
  cat("  ")
  print_quantiles("rt seconds error:", summary$rt_delta)
  cat("  ")
  print_quantiles("R/Python peak_area:", summary$area_ratio)
  print_unmatched_profile("  Python rows missing in R", summary$missing_python)
  print_unmatched_profile("  R rows not matched to Python", summary$extra_r)
}

print_unmatched_profile <- function(label, table) {
  if (is.null(table) || nrow(table) == 0L) {
    cat(label, ": none\n")
    return(invisible(NULL))
  }

  cat(label, ":", nrow(table), "\n")
  for (column in intersect(c("mz", "rtime", "peak_area", "cSelectivity", "goodness_fitting", "snr", "detection_counts"), names(table))) {
    values <- table[[column]]
    qs <- stats::quantile(values, probs = c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE)
    cat(
      "    ", column,
      " p10 ", round(qs[[1]], 3),
      " median ", round(qs[[2]], 3),
      " p90 ", round(qs[[3]], 3),
      "\n",
      sep = ""
    )
  }
}

run_r_checkpoints <- function(options) {
  parameters <- default_parameters()
  sample_name <- tools::file_path_sans_ext(basename(options$mzml_file))

  timed_tracks <- system.time({
    mass_tracks <- extract_mass_tracks(
      options$mzml_file,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm,
      min_intensity = parameters$min_intensity_threshold,
      min_timepoints = parameters$min_timepoints,
      min_peak_height = parameters$min_peak_height,
      max_scans = options$max_scans
    )
  })

  sample_mass_tracks <- list()
  sample_mass_tracks[[sample_name]] <- mass_tracks

  timed_features <- system.time({
    feature_result <- process_global_features(
      sample_mass_tracks,
      parameters = parameters,
      rt_table = mass_tracks$rt_table,
      drop_unaligned_samples = parameters$drop_unaligned_samples
    )
  })

  list(
    parameters = parameters,
    sample_name = sample_name,
    mass_tracks = mass_tracks,
    feature_result = feature_result,
    track_seconds = unname(timed_tracks[["elapsed"]]),
    feature_seconds = unname(timed_features[["elapsed"]])
  )
}

main <- function() {
  options <- parse_args(commandArgs(trailingOnly = TRUE))
  print_input_file_summary(options$mzml_file)

  if (!is.null(options$max_scans)) {
    cat("Max scans:", options$max_scans, "\n\n")
  }

  r_result <- run_r_checkpoints(options)
  mass_tracks <- r_result$mass_tracks
  feature_result <- r_result$feature_result
  full_table <- feature_result$full_feature_table
  preferred_table <- feature_result$preferred_feature_table

  cat("R checkpoints:\n")
  cat("  scans:", nrow(mass_tracks$rt_table), "\n")
  cat("  mass tracks:", length(mass_tracks$tracks), "\n")
  cat("  mass grid rows:", nrow(feature_result$mass_grid), "\n")
  cat("  composite tracks:", length(feature_result$composite_tracks), "\n")
  cat("  global peak candidates:", length(feature_result$feature_list), "\n")
  cat("  full Feature table rows:", nrow(full_table), "\n")
  cat("  preferred Feature table rows:", nrow(preferred_table), "\n")
  cat("  mass track seconds:", round(r_result$track_seconds, 3), "\n")
  cat("  feature seconds:", round(r_result$feature_seconds, 3), "\n\n")

  filter_counts <- count_feature_filters(full_table)
  cat("R preferred filters, counted independently:\n")
  for (name in names(filter_counts)) {
    cat("  pass", name, ":", filter_counts[[name]], "\n")
  }
  cat("\n")

  if (is.null(options$python_output)) {
    cat("Python asari comparison skipped: no --python-output directory was provided.\n")
    cat("Run again with --python-output path/to/asari_output when you have Python asari output.\n")
    return(invisible(NULL))
  }

  python_full_path <- table_path(options$python_output, "full_Feature_table.tsv")
  python_preferred_path <- table_path(options$python_output, "preferred_Feature_table.tsv")
  python_full <- read_feature_table(python_full_path)
  python_preferred <- read_feature_table(python_preferred_path)

  cat("Python asari output:\n")
  cat("  full table:", ifelse(is.null(python_full_path), "missing", python_full_path), "\n")
  cat("  preferred table:", ifelse(is.null(python_preferred_path), "missing", python_preferred_path), "\n\n")

  full_summary <- feature_match_summary(python_full, full_table, mz_ppm = options$mz_ppm, rt_sec = options$rt_sec)
  preferred_summary <- feature_match_summary(python_preferred, preferred_table, mz_ppm = options$mz_ppm, rt_sec = options$rt_sec)

  cat("Match tolerance:", options$mz_ppm, "ppm,", options$rt_sec, "seconds\n\n")
  print_match_summary("Full table match:", full_summary)
  cat("\n")
  print_match_summary("Preferred table match:", preferred_summary)
}

main()
