# 对应 Python asari/tools/extract_masstracks.py：旧版独立质量轨迹导出工具。

# 对应 read_project_dir：不排序，保留目录枚举顺序并按子串过滤。
read_project_dir_masstracks <- function(directory, file_pattern = ".mzML") {
  entries <- list.files(directory, full.names = FALSE, all.files = FALSE, no.. = TRUE)
  file.path(directory, entries[grepl(file_pattern, entries, fixed = TRUE)])
}

# 提供原函数名的显式别名；workflow.R加载后可能覆盖同名主流程函数。
read_project_dir_extract_masstracks <- read_project_dir_masstracks

# 从环境或list读取字段。
.extract_masstracks_get <- function(object, name) {
  if (is.environment(object)) return(object[[name]])
  object[[name]]
}

# 对应 process_single_file：提取轨迹并导出每条轨迹的摘要统计。
process_single_file <- function(infile, outdir = "") {
  # 原Python使用过时的SimpleSample(input_file=...)接口；R版允许注入同功能工厂。
  factory <- getOption("asariR.simple_sample_from_file")
  if (!is.function(factory)) {
    factory <- get0("SimpleSample_from_file", mode = "function", inherits = TRUE)
  }
  if (!is.function(factory)) {
    stop(
      "process_single_file requires options(asariR.simple_sample_from_file = function(path) ...).",
      call. = FALSE
    )
  }
  sample <- factory(infile)
  extractor <- .extract_masstracks_get(sample, "get_mass_tracks_")
  if (!is.function(extractor)) stop("Sample object lacks get_mass_tracks_.")
  extractor(mz_tolerance_ppm = 5, min_intensity = 100, min_timepoints = 5)
  tracks <- .extract_masstracks_get(sample, "list_mass_tracks")
  rows <- lapply(tracks, function(track) {
    intensity <- as.numeric(track$intensity)
    positive <- intensity[intensity > 0]
    median_positive <- if (length(positive)) as.integer(stats::median(positive)) else NA_integer_
    c(
      as.character(track$id_number), as.character(track$mz),
      as.character(length(intensity)), as.character(max(intensity)),
      as.character(median_positive)
    )
  })
  lines <- c(
    "id_number\tmz\tnumber_scans\tmax_intensity\tmedian_intensity",
    vapply(rows, paste, collapse = "\t", FUN.VALUE = "")
  )
  outfile <- file.path(outdir, sub(".mzML", "_massTrack.tsv", basename(infile), fixed = TRUE))
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}
