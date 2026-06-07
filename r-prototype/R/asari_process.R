default_parameters <- function() {
  list(
    project_name = "asari_project",
    outdir = "output",
    mode = "pos",
    mz_tolerance_ppm = 5,
    min_timepoints = 6,
    min_intensity_threshold = 1000,
    min_peak_height = 100000,
    workflow = "LC",
    anno = FALSE,
    asari_r_prototype = TRUE
  )
}

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

empty_feature_table <- function(samples) {
  fixed_columns <- data.frame(
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
    fixed_columns[[sample_name]] <- numeric()
  }

  fixed_columns
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

asari_process <- function(input_dir, output_dir = "output", project_name = "asari_project") {
  parameters <- default_parameters()
  parameters$input <- normalizePath(input_dir, mustWork = TRUE)
  parameters$outdir <- normalizePath(output_dir, mustWork = FALSE)
  parameters$project_name <- project_name
  parameters$timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  input_files <- read_project_dir(parameters$input)
  if (length(input_files) == 0L) {
    stop("No valid mzML files are found in the input directory.", call. = FALSE)
  }

  samples <- register_samples(input_files)
  parameters <- create_export_folders(parameters)

  feature_table <- empty_feature_table(samples)
  write_tsv(feature_table, file.path(parameters$outdir, "preferred_Feature_table.tsv"))
  write_tsv(feature_table, file.path(parameters$export_outdir, "full_Feature_table.tsv"))
  write_project_json(parameters, samples)

  list(
    parameters = parameters,
    samples = samples,
    feature_table = feature_table
  )
}

