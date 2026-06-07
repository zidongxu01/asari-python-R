preview_input_file <- function(default_file) {
  args <- commandArgs(trailingOnly = TRUE)
  file <- if (length(args) > 0L) args[[1]] else default_file

  if (!file.exists(file)) {
    stop("Input mzML file does not exist: ", file, call. = FALSE)
  }

  normalizePath(file, mustWork = TRUE)
}

print_input_file_summary <- function(file) {
  info <- file.info(file)
  size_mb <- round(info$size / 1024^2, 2)

  cat("Input mzML:\n")
  cat(file, "\n")
  cat("File size:", size_mb, "MB\n\n")
}

