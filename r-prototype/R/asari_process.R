source("r-prototype/R/parameters.R")
source("r-prototype/R/global_features.R")

read_project_dir <- function(directory, file_pattern = "\\.mzML$") {
  if (!dir.exists(directory)) {
    stop("Input directory does not exist: ", directory, call. = FALSE)
  }

  files <- list.files(
    path = directory,
    pattern = file_pattern,
    full.names = TRUE,
    ignore.case = FALSE
  )

  sort(normalizePath(files, mustWork = TRUE))
}

register_samples <- function(input_files) {
  data.frame(
    sample_id = seq_along(input_files) - 1L,
    input_file = input_files,
    sample_name = tools::file_path_sans_ext(basename(input_files)),
    status = "registered",
    stringsAsFactors = FALSE
  )
}

create_export_folders <- function(parameters, timestamp = format(Sys.time(), "%Y%m%d%H%M%S")) {
  outdir <- paste(parameters$outdir, parameters$project_name, timestamp, sep = "_")
  export_outdir <- file.path(outdir, "export")

  dir.create(export_outdir, recursive = TRUE, showWarnings = FALSE)

  parameters$outdir <- normalizePath(outdir, mustWork = TRUE)
  parameters$export_outdir <- normalizePath(export_outdir, mustWork = TRUE)
  parameters
}

select_single_input_file <- function(input_files, selected_file = NULL) {
  if (is.null(selected_file)) {
    if (length(input_files) == 1L) {
      return(input_files[[1]])
    }
    stop(
      "selected_file is required when process_features = TRUE and input_dir contains multiple mzML files.",
      call. = FALSE
    )
  }

  selected_file <- normalizePath(selected_file, mustWork = TRUE)
  if (!selected_file %in% input_files) {
    stop("selected_file must be one of the mzML files under input_dir.", call. = FALSE)
  }

  selected_file
}

extract_project_mass_tracks <- function(samples, parameters, max_scans = NULL) {
  sample_mass_tracks <- list()

  for (i in seq_len(nrow(samples))) {
    sample_name <- samples$sample_name[[i]]
    input_file <- samples$input_file[[i]]
    mass_tracks <- extract_mass_tracks(
      input_file,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm,
      min_intensity = parameters$min_intensity_threshold,
      min_timepoints = parameters$min_timepoints,
      min_peak_height = parameters$min_peak_height,
      max_scans = max_scans
    )
    sample_mass_tracks[[sample_name]] <- mass_tracks
  }

  sample_mass_tracks
}

process_project_features <- function(samples, parameters, max_scans = NULL) {
  sample_mass_tracks <- extract_project_mass_tracks(
    samples,
    parameters = parameters,
    max_scans = max_scans
  )

  rt_table <- NULL
  if (length(sample_mass_tracks) > 0L && !is.null(sample_mass_tracks[[1]]$rt_table)) {
    rt_table <- sample_mass_tracks[[1]]$rt_table
  }

  process_global_features(
    sample_mass_tracks,
    parameters = parameters,
    rt_table = rt_table,
    drop_unaligned_samples = parameters$drop_unaligned_samples
  )
}

empty_feature_table <- function(samples) {
  feature_table <- data.frame(
    id_number = character(),
    mz = numeric(),
    rtime = numeric(),
    rtime_left_base = numeric(),
    rtime_right_base = numeric(),
    parent_masstrack_id = integer(),
    peak_area = numeric(),
    cSelectivity = numeric(),
    goodness_fitting = numeric(),
    snr = numeric(),
    detection_counts = integer(),
    stringsAsFactors = FALSE
  )
  for (sample_name in samples$sample_name) {
    feature_table[[sample_name]] <- numeric()
  }
  feature_table
}

preview_single_sample_features <- function(input_file, parameters, max_scans = NULL) {
  mass_tracks <- extract_mass_tracks(
    file = input_file,
    mz_tolerance_ppm = parameters$mz_tolerance_ppm,
    min_intensity = parameters$min_intensity_threshold,
    min_timepoints = parameters$min_timepoints,
    min_peak_height = parameters$min_peak_height,
    max_scans = max_scans
  )

  list(
    mass_tracks = mass_tracks
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

write_project_json <- function(parameters, samples) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to write project.json.", call. = FALSE)
  }

  payload <- parameters
  payload$number_of_samples <- nrow(samples)
  payload$samples <- lapply(seq_len(nrow(samples)), function(i) as.list(samples[i, ]))

  jsonlite::write_json(
    payload,
    path = file.path(parameters$outdir, "project.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

asari_process <- function(
  input_dir,
  output_dir = "output",
  project_name = "asari_project",
  process_features = FALSE,
  selected_file = NULL,
  max_scans = NULL
) {
  parameters <- default_parameters()
  parameters$input <- normalizePath(input_dir, mustWork = TRUE)
  parameters$outdir <- normalizePath(output_dir, mustWork = FALSE)
  parameters$project_name <- project_name
  parameters$timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  parameters$process_features <- process_features
  parameters$max_scans <- max_scans

  input_files <- read_project_dir(parameters$input)
  if (length(input_files) == 0L) {
    stop("No valid mzML files are found in the input directory.", call. = FALSE)
  }

  selected_input_file <- NULL
  sample_files <- input_files
  if (process_features) {
    selected_input_file <- select_single_input_file(input_files, selected_file = selected_file)
    sample_files <- selected_input_file
  }

  samples <- register_samples(sample_files)
  parameters <- create_export_folders(parameters)

  full_feature_table <- empty_feature_table(samples)
  preferred_feature_table <- full_feature_table
  feature_result <- NULL
  if (process_features) {
    parameters$selected_file <- selected_input_file

    feature_result <- process_project_features(
      samples,
      parameters = parameters,
      max_scans = max_scans
    )

    full_feature_table <- feature_result$full_feature_table
    preferred_feature_table <- feature_result$preferred_feature_table
  }

  write_tsv(preferred_feature_table, file.path(parameters$outdir, "preferred_Feature_table.tsv"))
  write_tsv(full_feature_table, file.path(parameters$export_outdir, "full_Feature_table.tsv"))
  write_project_json(parameters, samples)

  list(
    parameters = parameters,
    samples = samples,
    feature_table = preferred_feature_table,
    full_feature_table = full_feature_table,
    preferred_feature_table = preferred_feature_table,
    feature_result = feature_result
  )
}
