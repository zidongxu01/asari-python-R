# 对应 Python asari/utils.py：并行任务、文件下载、校验和及 mzML 工具。

# 将 R matrix 按 Python ndarray.tolist() 的行顺序转换为嵌套 list。
.utils_matrix_to_list <- function(value) {
  lapply(seq_len(nrow(value)), function(row) {
    as.list(unname(value[row, , drop = TRUE]))
  })
}

# 下载 URL 到临时文件；file:// 分支仅用于无网络的本地验证。
.utils_download <- function(url, destination) {
  if (startsWith(url, "file://")) {
    source <- sub("^file://", "", url)
    if (!file.copy(source, destination, overwrite = TRUE)) {
      stop("Failed to copy local archive: ", source, call. = FALSE)
    }
    return(invisible(destination))
  }
  status <- utils::download.file(
    url, destination, mode = "wb", quiet = TRUE
  )
  if (!identical(status, 0L)) {
    stop("Failed to download archive: ", url, call. = FALSE)
  }
  invisible(destination)
}

# 定位 R package 资源目录；目录路径输入用于本地开发和测试。
.utils_package_path <- function(package) {
  if (dir.exists(package)) return(normalizePath(package))
  path <- system.file(package = package)
  if (!nzchar(path)) stop("R package is not installed: ", package, call. = FALSE)
  path
}

# 从路径或测试提供的 header 中获取 mzML 扫描头信息。
.utils_mzml_header <- function(mzml_file) {
  if (is.list(mzml_file) && !is.null(mzml_file$header)) {
    return(mzml_file$header)
  }
  if (is.environment(mzml_file) && exists(
      "header", envir = mzml_file, inherits = FALSE)) {
    return(mzml_file$header)
  }
  path <- if (is.character(mzml_file)) mzml_file else {
    if (is.environment(mzml_file)) mzml_file$path else mzml_file$path
  }
  handle <- mzR::openMSfile(path)
  on.exit(mzR::close(handle), add = TRUE)
  mzR::header(handle)
}

# Windows 没有 fork，单独使用临时 PSOCK cluster 执行并行任务。
.utils_psock_apply <- function(arguments, command, workers) {
  worker_cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(worker_cluster), add = TRUE)
  parallel::parLapply(worker_cluster, arguments, command)
}

# 创建与 Python json.JSONEncoder 子类相对应的轻量可变对象。
NpEncoder <- function() {
  self <- new.env(parent = emptyenv())
  class(self) <- c("NpEncoder", "environment")
  self$default <- function(obj) NpEncoder_default(self, obj)
  self
}

# 对应 NpEncoder.default：把 numpy 对应类型转换成 JSON 可编码的基础类型。
NpEncoder_default <- function(self, obj) {
  # self 与 Python JSONEncoder 实例对应，转换本身不依赖实例状态。
  invisible(self)
  if (is.integer(obj) && length(obj) == 1L) return(as.integer(obj))
  if (is.numeric(obj) && length(obj) == 1L) return(as.numeric(obj))
  if (is.matrix(obj)) return(.utils_matrix_to_list(obj))
  if (is.array(obj) || (is.atomic(obj) && length(obj) > 1L)) {
    return(as.list(unname(obj)))
  }
  stop(
    "Object of class ", paste(class(obj), collapse = "/"),
    " is not JSON serializable.",
    call. = FALSE
  )
}

# 对应 bulk_process：把 command 逐个应用到 arguments，并按输入顺序返回结果。
bulk_process <- function(command, arguments, jobs_per_worker = FALSE) {
  if (!is.function(command)) stop("command must be a function.", call. = FALSE)
  if (!is.list(arguments)) arguments <- as.list(arguments)
  workers <- if (
      isFALSE(jobs_per_worker) || identical(jobs_per_worker, "auto")) {
    parallel::detectCores()
  } else {
    as.integer(jobs_per_worker)
  }
  if (is.na(workers) || workers < 1L) stop("Invalid worker count.")
  if (length(arguments) == 0L) return(list())
  if (workers == 1L) return(lapply(arguments, command))

  # Unix 使用 fork；Windows 使用临时 PSOCK cluster，均保留 starmap 的输入顺序。
  if (.Platform$OS.type != "windows") {
    parallel::mclapply(arguments, command, mc.cores = workers)
  } else {
    .utils_psock_apply(arguments, command, workers)
  }
}

# 对应 download_and_unzip_to_pkg_resources：下载 ZIP 到 package 父目录的子目录。
download_and_unzip_to_pkg_resources <- function(url, package, subdir = "data") {
  message("HERE")
  package_dir <- dirname(.utils_package_path(package))
  extract_to <- file.path(package_dir, subdir)
  dir.create(extract_to, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".zip")
  on.exit(unlink(archive), add = TRUE)
  .utils_download(url, archive)
  utils::unzip(archive, exdir = extract_to)
  message("Extracted to: ", extract_to)
  invisible(NULL)
}

# 对应 download_and_unzip：下载 ZIP 并解压到明确目标目录。
download_and_unzip <- function(url, extract_to) {
  dir.create(extract_to, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".zip")
  on.exit(unlink(archive), add = TRUE)
  .utils_download(url, archive)
  utils::unzip(archive, exdir = extract_to)
  message("Extracted to: ", extract_to)
  invisible(NULL)
}

# 对应 validate_mzml_file：完整遍历 mzML 扫描，任何读取异常均返回 FALSE。
validate_mzml_file <- function(file) {
  tryCatch({
    handle <- mzR::openMSfile(file)
    on.exit(mzR::close(handle), add = TRUE)
    scan_header <- mzR::header(handle)
    if (nrow(scan_header) > 0L) {
      for (scan in seq_len(nrow(scan_header))) mzR::peaks(handle, scan)
    }
    TRUE
  }, error = function(error) FALSE)
}

# 对应 build_boolean_dict：保留 Python 接受的 8 个布尔输入键。
build_boolean_dict <- function() {
  list(
    T = TRUE,
    F = FALSE,
    `1` = TRUE,
    `0` = FALSE,
    True = TRUE,
    False = FALSE,
    `TRUE` = TRUE,
    `FALSE` = FALSE,
    true = TRUE,
    false = FALSE
  )
}

# 对应 sizeof_fmt：将字节数格式化为 Python 使用的 IEC 单位文本。
sizeof_fmt <- function(num, suffix = "B") {
  units <- c("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi")
  value <- as.numeric(num)
  for (unit in units) {
    if (abs(value) < 1024) return(sprintf("%3.1f%s%s", value, unit, suffix))
    value <- value / 1024
  }
  sprintf("%.1fYi%s", value, suffix)
}

# 对应 checksum_file：计算文件 MD5；chunksize 保留 Python 接口。
checksum_file <- function(file, chunksize = 16384) {
  if (!file.exists(file) || dir.exists(file)) stop("file must be a regular file.")
  if (length(chunksize) != 1L || chunksize <= 0) stop("Invalid chunksize.")
  unname(tools::md5sum(file))
}

# 对应 wait_with_pbar：每秒更新一次等待进度条。
wait_with_pbar <- function(wait = 5) {
  if (length(wait) != 1L || is.na(wait) || wait < 0 || wait != as.integer(wait)) {
    stop("wait must be a non-negative integer.")
  }
  wait <- as.integer(wait)
  if (wait == 0L) return(invisible(NULL))
  progress <- utils::txtProgressBar(min = 0, max = wait, style = 3)
  on.exit(close(progress), add = TRUE)
  for (ii in seq_len(wait)) {
    Sys.sleep(1)
    utils::setTxtProgressBar(progress, ii)
  }
  invisible(NULL)
}

# 对应 get_ionization_mode_mzml：检查前 limit+1 个扫描的正负离子模式。
get_ionization_mode_mzml <- function(mzml_file, limit = 50) {
  header <- .utils_mzml_header(mzml_file)
  if (!("polarity" %in% names(header))) stop("mzML polarity is unavailable.")
  modes <- character()
  scan_limit <- min(nrow(header), as.integer(limit) + 1L)
  if (scan_limit > 0L) {
    for (ii in seq_len(scan_limit)) {
      polarity <- header$polarity[[ii]]
      if (is.na(polarity) || !(polarity %in% c(0, 1))) {
        stop("Scan has neither positive nor negative polarity.")
      }
      modes <- unique(c(modes, if (polarity == 1) "pos" else "neg"))
      if (length(modes) > 1L) return("mixed")
    }
  }
  if (length(modes) == 0L) stop("mzML contains no spectra.")
  modes[[1L]]
}

# 保留早期 R 模块使用的 ppm 换算辅助函数；它不属于 Python 的 10 个 def。
ppm_to_mz_tolerance <- function(mz, ppm) {
  mz * ppm * 1e-6
}
