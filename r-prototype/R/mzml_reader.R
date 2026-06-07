require_mzR <- function() {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop(
      "The mzR package is required to read mzML files. ",
      "Install it with BiocManager::install('mzR').",
      call. = FALSE
    )
  }
}

read_mzml_header <- function(file) {
  require_mzR()

  if (!file.exists(file)) {
    stop("mzML file does not exist: ", file, call. = FALSE)
  }

  handle <- mzR::openMSfile(file)
  on.exit(mzR::close(handle), add = TRUE)

  mzR::header(handle)
}

read_mzml_points <- function(file, max_scans = NULL, min_intensity = 0) {
  require_mzR()

  if (!file.exists(file)) {
    stop("mzML file does not exist: ", file, call. = FALSE)
  }
  if (!is.null(max_scans) && max_scans <= 0) {
    stop("max_scans must be positive when provided.", call. = FALSE)
  }

  handle <- mzR::openMSfile(file)
  on.exit(mzR::close(handle), add = TRUE)

  scan_header <- mzR::header(handle)
  scan_ids <- which(scan_header$msLevel == 1L)
  if (!is.null(max_scans)) {
    scan_ids <- head(scan_ids, max_scans)
  }

  rows <- vector("list", length(scan_ids))
  for (i in seq_along(scan_ids)) {
    scan_id <- scan_ids[[i]]
    peak_matrix <- mzR::peaks(handle, scan_id)

    if (nrow(peak_matrix) == 0L) {
      rows[[i]] <- data.frame()
      next
    }

    keep <- peak_matrix[, "intensity"] > min_intensity
    peak_matrix <- peak_matrix[keep, , drop = FALSE]

    rows[[i]] <- data.frame(
      scan_id = scan_id,
      rtime = scan_header$retentionTime[[scan_id]],
      mz = peak_matrix[, "mz"],
      intensity = peak_matrix[, "intensity"],
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(
      scan_id = integer(),
      rtime = numeric(),
      mz = numeric(),
      intensity = numeric()
    ))
  }

  do.call(rbind, rows)
}

