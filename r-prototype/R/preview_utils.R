preview_input_file <- function(script_name = "preview script") {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    stop(
      "Missing mzML file path.\n",
      "Usage: Rscript ", script_name, " \"path/to/file.mzML\"",
      call. = FALSE
    )
  }

  file <- args[[1]]

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
