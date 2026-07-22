# 面向普通用户的附加入口：把已经完成的纯R模块接到稳定、可校验的公开函数。

# 检查一个参数是否为单个正数，并返回数值形式。
.asari_public_positive_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
    stop(name, " must be one positive number.", call. = FALSE)
  }
  as.numeric(value)
}

# 从显式路径或data.frame读取特征表，不寻找任何默认文件。
.asari_public_feature_table <- function(value, name = "feature_table") {
  if (is.data.frame(value)) return(value)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop(name, " must be a data.frame or one explicit TSV file path.", call. = FALSE)
  }
  path <- path.expand(value)
  if (!file.exists(path) || dir.exists(path)) {
    stop(name, " does not exist or is not a file: ", path, call. = FALSE)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# 检查特征表必须具备的列，错误中直接列出缺少的名称。
.asari_public_require_columns <- function(table, columns, name = "feature_table") {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0L) {
    stop(
      name, " is missing required column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

# 准备一个明确的输出文件，并创建它的父目录。
.asari_public_output_file <- function(output, name = "output") {
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(trimws(output))) {
    stop(name, " must be one explicit file path.", call. = FALSE)
  }
  path <- path.expand(output)
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", parent, call. = FALSE)
  }
  normalizePath(path, mustWork = FALSE)
}

# 检查并返回一个明确的mzML文件。
.asari_public_one_mzml <- function(input) {
  files <- .asari_collect_input_files(input, recursive = FALSE)
  if (length(files) != 1L) {
    stop("This function requires exactly one explicit .mzML file.", call. = FALSE)
  }
  files[[1L]]
}

#' 使用DIMS流程处理直接进样质谱数据
#'
#' 该函数复用[asari_process()]的输入检查、质量轨迹提取和结果导出，
#' 但固定使用`DIMS`流程并关闭保留时间对齐。
#'
#' @inheritParams asari_process
#' @return 一个`asari_result`对象。DIMS的主要结果位于`full_feature_table`。
#' @export
asari_process_dims <- function(
    input,
    output = "asari-output",
    project_name = "asari_dims",
    mode = NULL,
    ppm = NULL,
    multicores = NULL,
    database_mode = NULL,
    recursive = FALSE,
    parameters = list()) {
  # DIMS不使用色谱峰边界或跨样本保留时间对齐。
  parameters$workflow <- "DIMS"
  parameters$rt_align_on <- FALSE
  asari_process(
    input = input,
    output = output,
    project_name = project_name,
    mode = mode,
    ppm = ppm,
    multicores = multicores,
    rt_align = FALSE,
    database_mode = database_mode,
    recursive = recursive,
    parameters = parameters
  )
}

#' 分析一个mzML文件的扫描和质量轨迹
#'
#' 该函数只执行本地mzML读取、质量轨迹提取和同位素锚点统计，
#' 不连接已知化合物数据库，也不执行外部注释。
#'
#' @param input 一个明确的`.mzML`文件。
#' @param ppm m/z容差，单位为ppm。
#' @param min_intensity 最低数据点强度。
#' @param min_timepoints 一条质量轨迹至少包含的数据点数。
#' @param min_peak_height 最低峰高。
#' @param parameters 其他高级参数的具名覆盖列表。
#' @return 包含扫描摘要、锚点统计、质量轨迹表和原始轨迹的具名列表。
#' @export
asari_analyze <- function(
    input,
    ppm = 5,
    min_intensity = 1000,
    min_timepoints = 6L,
    min_peak_height = 100000,
    parameters = list()) {
  # 所有分析都从用户给出的唯一文件开始，不读取工作目录中的其他mzML。
  infile <- .asari_public_one_mzml(input)
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Analyzing mzML requires the Bioconductor package mzR.", call. = FALSE)
  }
  effective <- .asari_merge_parameters(parameters)
  effective$mz_tolerance_ppm <- .asari_public_positive_number(ppm, "ppm")
  effective$min_intensity_threshold <- .asari_public_positive_number(
    min_intensity, "min_intensity"
  )
  effective$min_timepoints <- as.integer(
    .asari_public_positive_number(min_timepoints, "min_timepoints")
  )
  effective$min_peak_height <- .asari_public_positive_number(
    min_peak_height, "min_peak_height"
  )

  # 只读取一次扫描头，并把它注入现有分析模块以避免重复打开大文件。
  header <- .analyze_scan_header(infile, effective)
  effective$scan_header <- header
  sample <- get_file_masstrack_stats(infile, effective, return_sample = TRUE)
  tracks <- sample$list_mass_tracks
  track_table <- if (length(tracks) == 0L) {
    data.frame(
      id_number = integer(), mz = numeric(), number_scans = integer(),
      max_intensity = numeric(), median_positive_intensity = numeric()
    )
  } else {
    do.call(rbind, lapply(tracks, function(track) {
      intensity <- as.numeric(track$intensity)
      positive <- intensity[intensity > 0]
      data.frame(
        id_number = as.integer(track$id_number),
        mz = as.numeric(track$mz),
        number_scans = length(intensity),
        max_intensity = max(intensity),
        median_positive_intensity = if (length(positive)) {
          stats::median(positive)
        } else {
          NA_real_
        }
      )
    }))
  }

  # 返回既适合快速查看、又足以继续下游计算的结构化结果。
  levels <- as.integer(header$msLevel)
  result <- list(
    input_file = infile,
    ion_mode = sample$ion_mode,
    ms1_scans = sum(levels == 1L, na.rm = TRUE),
    ms2_scans = sum(levels == 2L, na.rm = TRUE),
    max_retention_time = if (length(sample$list_retention_time)) {
      max(sample$list_retention_time)
    } else {
      NA_real_
    },
    mass_track_count = length(tracks),
    isotope_pair_count = length(sample$anchor_mz_pairs),
    landmark_min_height = sample$landmark_min_height,
    landmark_max_height = sample$landmark_max_height,
    mass_track_table = track_table,
    mass_tracks = tracks,
    parameters = effective[names(effective) != "scan_header"]
  )
  class(result) <- "asari_analysis"
  result
}

#' 导出一个mzML文件的质量轨迹摘要
#'
#' @inheritParams asari_analyze
#' @param output 可选的明确TSV输出路径；为`NULL`时只返回data.frame。
#' @return 质量轨迹摘要data.frame。
#' @export
asari_extract_mass_tracks <- function(
    input,
    output = NULL,
    ppm = 5,
    min_intensity = 1000,
    min_timepoints = 6L,
    min_peak_height = 100000,
    parameters = list()) {
  # 复用单文件分析，保证轨迹提取参数和统计口径只有一套实现。
  analysis <- asari_analyze(
    input = input,
    ppm = ppm,
    min_intensity = min_intensity,
    min_timepoints = min_timepoints,
    min_peak_height = min_peak_height,
    parameters = parameters
  )
  table <- analysis$mass_track_table
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(
      table, outfile, sep = "\t", row.names = FALSE,
      col.names = TRUE, quote = FALSE, na = ""
    )
    attr(table, "output_file") <- outfile
  }
  attr(table, "analysis") <- analysis
  table
}

#' 读取asariR已经生成的项目结果
#'
#' @param project 一个项目目录，或[asari_process()]返回的`asari_result`。
#' @param table 要读取`"preferred"`、`"full"`或`"both"`表。
#' @param max_samples 可选的最大样本列数；`NULL`表示读取全部样本列。
#' @return 包含项目描述和所选特征表的`asari_project_results`对象。
#' @export
asari_read_results <- function(
    project,
    table = c("both", "preferred", "full"),
    max_samples = NULL) {
  table <- match.arg(table)
  project_dir <- if (inherits(project, "asari_result")) {
    project$project_dir
  } else {
    project
  }
  if (!is.character(project_dir) || length(project_dir) != 1L ||
      is.na(project_dir) || !dir.exists(path.expand(project_dir))) {
    stop("project must be an existing project directory or asari_result.", call. = FALSE)
  }
  project_dir <- normalizePath(path.expand(project_dir), mustWork = TRUE)
  description_file <- file.path(project_dir, "project.json")
  if (!file.exists(description_file)) {
    stop("Project is missing project.json: ", project_dir, call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Reading project.json requires jsonlite.", call. = FALSE)
  }
  if (!is.null(max_samples)) {
    if (!is.numeric(max_samples) || length(max_samples) != 1L ||
        is.na(max_samples) || max_samples < 0 || max_samples != as.integer(max_samples)) {
      stop("max_samples must be NULL or one non-negative integer.", call. = FALSE)
    }
    max_samples <- as.integer(max_samples)
  }

  # asari表的前11列是固定特征元数据，后面的列才是样本丰度。
  read_one <- function(path) {
    if (!file.exists(path)) stop("Project is missing feature table: ", path, call. = FALSE)
    value <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(max_samples) && ncol(value) > 11L + max_samples) {
      value <- value[, seq_len(11L + max_samples), drop = FALSE]
    }
    value
  }
  project_description <- jsonlite::fromJSON(
    description_file, simplifyVector = FALSE
  )
  # 使用项目实际记录的表名，兼容output_feature_table高级覆盖值。
  feature_table_name <- if (!is.null(project_description$output_feature_table)) {
    as.character(project_description$output_feature_table)
  } else {
    "Feature_table.tsv"
  }
  result <- list(
    project_dir = project_dir,
    project = project_description,
    preferred = NULL,
    full = NULL
  )
  if (table %in% c("both", "preferred")) {
    result$preferred <- read_one(file.path(
      project_dir, paste0("preferred_", feature_table_name)
    ))
  }
  if (table %in% c("both", "full")) {
    result$full <- read_one(file.path(
      project_dir, "export", paste0("full_", feature_table_name)
    ))
  }
  class(result) <- "asari_project_results"
  result
}

#' 比较两张LC-MS特征表
#'
#' 两个方向都选择ppm和RT窗口内m/z偏差最小的候选，只保留相互最佳匹配。
#'
#' @param left,right TSV路径或data.frame。
#' @param ppm 正数m/z容差。
#' @param rt_tolerance 正数保留时间容差，单位为秒，使用严格小于关系。
#' @param left_rt_unit,right_rt_unit 各表保留时间单位。
#' @param verbose 是否打印底层双向匹配统计。
#' @return 每行一个相互最佳匹配的data.frame。
#' @export
asari_compare_features <- function(
    left,
    right,
    ppm = 5,
    rt_tolerance = 5,
    left_rt_unit = c("seconds", "minutes"),
    right_rt_unit = c("seconds", "minutes"),
    verbose = FALSE) {
  ppm <- .asari_public_positive_number(ppm, "ppm")
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  left_rt_unit <- match.arg(left_rt_unit)
  right_rt_unit <- match.arg(right_rt_unit)
  left_table <- .asari_public_feature_table(left, "left")
  right_table <- .asari_public_feature_table(right, "right")
  .asari_public_require_columns(left_table, c("mz", "rtime"), "left")
  .asari_public_require_columns(right_table, c("mz", "rtime"), "right")

  # 没有id_number时使用稳定的行号ID，避免要求用户先修改自己的表。
  if (!("id_number" %in% names(left_table))) {
    left_table$id_number <- paste0("row", seq_len(nrow(left_table)))
  }
  if (!("id_number" %in% names(right_table))) {
    right_table$id_number <- paste0("row", seq_len(nrow(right_table)))
  }
  left_table$id_number <- as.character(left_table$id_number)
  right_table$id_number <- as.character(right_table$id_number)
  if (anyDuplicated(left_table$id_number) || anyDuplicated(right_table$id_number)) {
    stop("id_number values must be unique within each feature table.", call. = FALSE)
  }
  left_seconds <- as.numeric(left_table$rtime) * if (left_rt_unit == "minutes") 60 else 1
  right_seconds <- as.numeric(right_table$rtime) * if (right_rt_unit == "minutes") 60 else 1
  if (anyNA(c(as.numeric(left_table$mz), left_seconds)) ||
      anyNA(c(as.numeric(right_table$mz), right_seconds))) {
    stop("mz and rtime columns must contain numeric values without NA.", call. = FALSE)
  }
  left_features <- lapply(seq_len(nrow(left_table)), function(ii) list(
    id = left_table$id_number[[ii]], mz = as.numeric(left_table$mz[[ii]]),
    rtime = left_seconds[[ii]]
  ))
  right_features <- lapply(seq_len(nrow(right_table)), function(ii) list(
    id = right_table$id_number[[ii]], mz = as.numeric(right_table$mz[[ii]]),
    rtime = right_seconds[[ii]]
  ))

  # 直接调用已经逐def翻译并测试的双向最佳匹配算法。
  matched <- if (isTRUE(verbose)) {
    bidirectional_best_match(left_features, right_features, ppm, rt_tolerance)
  } else {
    invisible(capture.output(
      value <- bidirectional_best_match(
        left_features, right_features, ppm, rt_tolerance
      )
    ))
    value
  }
  pairs <- matched[[1L]]
  if (length(pairs) == 0L) {
    result <- data.frame(
      left_id = character(), right_id = character(),
      left_mz = numeric(), right_mz = numeric(), ppm_error = numeric(),
      left_rtime = numeric(), right_rtime = numeric(), rt_delta = numeric()
    )
  } else {
    result <- do.call(rbind, lapply(pairs, function(pair) {
      li <- match(pair[[1L]], left_table$id_number)
      ri <- match(pair[[2L]], right_table$id_number)
      data.frame(
        left_id = pair[[1L]], right_id = pair[[2L]],
        left_mz = as.numeric(left_table$mz[[li]]),
        right_mz = as.numeric(right_table$mz[[ri]]),
        ppm_error = (as.numeric(right_table$mz[[ri]]) -
          as.numeric(left_table$mz[[li]])) /
          as.numeric(left_table$mz[[li]]) * 1e6,
        left_rtime = left_seconds[[li]], right_rtime = right_seconds[[ri]],
        rt_delta = right_seconds[[ri]] - left_seconds[[li]],
        stringsAsFactors = FALSE
      )
    }))
  }
  attr(result, "summary") <- c(
    left_features = nrow(left_table), right_features = nrow(right_table),
    mutual_matches = nrow(result), left_unmatched = nrow(left_table) - nrow(result),
    right_unmatched = nrow(right_table) - nrow(result)
  )
  result
}

#' 按目标m/z从特征表中提取候选峰
#'
#' @param feature_table asari特征表路径或data.frame。
#' @param targets 一个或多个目标m/z。
#' @param ppm 正数m/z容差。
#' @param output 可选的明确TSV输出路径。
#' @return 命中的特征行；首列为查询质量，并附带ppm误差。
#' @export
asari_extract_targets <- function(feature_table, targets, ppm = 5, output = NULL) {
  table <- .asari_public_feature_table(feature_table)
  .asari_public_require_columns(table, "mz")
  ppm <- .asari_public_positive_number(ppm, "ppm")
  if (!is.numeric(targets) || length(targets) == 0L || anyNA(targets) ||
      any(targets <= 0)) {
    stop("targets must contain one or more positive numeric m/z values.", call. = FALSE)
  }
  feature_mz <- as.numeric(table$mz)
  if (anyNA(feature_mz)) stop("feature_table$mz must be numeric without NA.", call. = FALSE)
  mapping <- all_mass_paired_mapping(feature_mz, as.numeric(targets), ppm)
  pairs <- mapping$mapped
  if (length(pairs) == 0L) {
    result <- table[FALSE, , drop = FALSE]
    result <- cbind(
      query_target = numeric(), ppm_error = numeric(), result,
      stringsAsFactors = FALSE
    )
  } else {
    feature_index <- vapply(pairs, `[[`, integer(1), 1L)
    target_index <- vapply(pairs, `[[`, integer(1), 2L)
    query <- as.numeric(targets)[target_index]
    result <- cbind(
      query_target = query,
      ppm_error = (feature_mz[feature_index] - query) / query * 1e6,
      table[feature_index, , drop = FALSE],
      stringsAsFactors = FALSE
    )
  }
  attr(result, "unmatched_targets") <- as.numeric(targets)[mapping$list2_unmapped]
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(
      result, outfile, sep = "\t", row.names = FALSE,
      col.names = TRUE, quote = FALSE, na = ""
    )
    attr(result, "output_file") <- outfile
  }
  result
}

#' 为特征表生成QC图
#'
#' @param feature_table asari完整特征表路径或data.frame。
#' @param output 明确的PDF输出路径。
#' @param height PDF高度。
#' @param aspect PDF宽高比。
#' @return 生成的PDF绝对路径。
#' @export
asari_feature_qc <- function(feature_table, output, height = 12, aspect = 0.7) {
  table <- .asari_public_feature_table(feature_table)
  .asari_public_require_columns(
    table,
    c("snr", "peak_area", "cSelectivity", "goodness_fitting")
  )
  outfile <- .asari_public_output_file(output)
  asari_qc_plot(table, outfile = outfile, height = height, aspect = aspect)
  normalizePath(outfile, mustWork = TRUE)
}

#' 为一个mzML文件生成HTML质量检查报告
#'
#' @param input 一个明确的`.mzML`文件。
#' @param output 明确的HTML输出路径。
#' @param spikeins `NULL`、目标list或目标JSON文件。
#' @return 生成的HTML绝对路径。
#' @export
asari_qc_report <- function(input, output, spikeins = NULL) {
  infile <- .asari_public_one_mzml(input)
  outfile <- .asari_public_output_file(output)
  generate_qc_report(list(infile, outfile, spikeins))
  normalizePath(outfile, mustWork = TRUE)
}

#' 为一个mzML文件生成MS1/MS2扫描摘要PDF
#'
#' @param input 一个明确的`.mzML`文件。
#' @param output 明确的PDF输出路径。
#' @param nspec_plot 要展示的示例扫描数。
#' @param offset_n 从第几个扫描附近开始展示。
#' @param width,height PDF尺寸。
#' @param title 图标题。
#' @return 生成的PDF绝对路径。
#' @export
asari_scan_summary <- function(
    input,
    output,
    nspec_plot = 50L,
    offset_n = 100L,
    width = 10,
    height = 6,
    title = "") {
  infile <- .asari_public_one_mzml(input)
  outfile <- .asari_public_output_file(output)
  summarize_ms_file(
    infile = infile, nspec_plot = nspec_plot, offset_n = offset_n,
    width = width, height = height, title = title, outfile = outfile
  )
  normalizePath(outfile, mustWork = TRUE)
}

#' 查看当前工作流的可用状态
#'
#' @return 一张说明本地可直接使用和外部依赖状态的data.frame。
#' @export
asari_available_workflows <- function() {
  # 状态描述安装后普通用户现在能直接调用的公开入口。
  data.frame(
    workflow = c("LC", "DIMS", "GC", "LCMSMS", "LC_annotation", "dashboard", "RAW_conversion"),
    direct_use = rep(TRUE, 7L),
    public_function = c(
      "asari_process", "asari_process_dims", "asari_annotate_gc",
      "asari_extract_ms2 / asari_match_ms2 / asari_search_ms2",
      "asari_annotate_lc", "asari_dashboard", "asari_convert_raw"
    ),
    note = c(
      "Local core workflow",
      "Local direct-infusion workflow; use full_feature_table",
      "Requires explicit Kovats standards and an MSP/JSON EI library",
      "Local extraction and search; entropy uses a documented local fallback",
      "Mass/RT candidate annotation; hmdb4 first load needs Python pickle decoding",
      "Self-contained offline HTML; no Shiny or dashboard pickle required",
      "Requires ThermoRawFileParser and, on macOS/Linux, a compatible runtime"
    ),
    stringsAsFactors = FALSE
  )
}

# 打印单文件分析摘要，不自动展开可能很大的轨迹list。
#' @export
print.asari_analysis <- function(x, ...) {
  invisible(list(...))
  cat("asariR single-file analysis.\n")
  cat("Input: ", x$input_file, "\n", sep = "")
  cat("Ion mode: ", x$ion_mode, "\n", sep = "")
  cat("MS1 scans: ", x$ms1_scans, "; MS2 scans: ", x$ms2_scans, "\n", sep = "")
  cat("Mass tracks: ", x$mass_track_count,
      "; isotope pairs: ", x$isotope_pair_count, "\n", sep = "")
  invisible(x)
}

# 打印已经读取的项目摘要，表格本身仍保存在对象中。
#' @export
print.asari_project_results <- function(x, ...) {
  invisible(list(...))
  cat("asariR project results.\n")
  cat("Project directory: ", x$project_dir, "\n", sep = "")
  if (!is.null(x$preferred)) cat("Preferred features: ", nrow(x$preferred), "\n", sep = "")
  if (!is.null(x$full)) cat("Full features: ", nrow(x$full), "\n", sep = "")
  invisible(x)
}
