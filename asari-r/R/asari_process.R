# 面向普通用户的R包入口：校验显式输入，组织参数并调用已经完成的asari主流程。

# 检查实际处理mzML和写出项目记录所必需的R依赖。
.asari_check_runtime_dependencies <- function() {
  required <- c("mzR", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Install mzR with BiocManager and jsonlite with install.packages().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# 从用户明确给出的文件或目录收集mzML；不搜索隐藏的默认数据位置。
.asari_collect_input_files <- function(input, recursive = FALSE) {
  if (!is.character(input) || length(input) == 0L || anyNA(input) ||
      any(!nzchar(trimws(input)))) {
    stop("input must contain at least one explicit mzML file or directory.", call. = FALSE)
  }
  if (!is.logical(recursive) || length(recursive) != 1L || is.na(recursive)) {
    stop("recursive must be TRUE or FALSE.", call. = FALSE)
  }

  # 逐项展开路径，目录只读取标准的.mzML文件，单文件则严格检查扩展名。
  collected <- character()
  for (item in path.expand(input)) {
    if (dir.exists(item)) {
      collected <- c(
        collected,
        list.files(
          item,
          pattern = "[.]mzML$",
          full.names = TRUE,
          recursive = recursive,
          ignore.case = FALSE
        )
      )
    } else if (file.exists(item)) {
      if (!grepl("[.]mzML$", basename(item))) {
        stop("Input file is not a .mzML file: ", item, call. = FALSE)
      }
      collected <- c(collected, item)
    } else {
      stop("Input path does not exist: ", item, call. = FALSE)
    }
  }

  # 规范化并排序，保证相同输入在不同运行中的样本顺序稳定。
  collected <- sort(unique(normalizePath(collected, mustWork = TRUE)))
  if (length(collected) == 0L) {
    stop("No .mzML files were found in the supplied input path(s).", call. = FALSE)
  }

  # asari使用文件名登记样本，重复文件名会造成输出列难以区分，因此提前报错。
  duplicated_names <- unique(basename(collected)[duplicated(basename(collected))])
  if (length(duplicated_names) > 0L) {
    stop(
      "Duplicate mzML file name(s) are not supported: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  collected
}

# 把用户覆盖项合并进一份新的默认参数，避免修改模块级PARAMETERS对象。
.asari_merge_parameters <- function(overrides = list()) {
  if (!is.list(overrides)) {
    stop("parameters must be a named list.", call. = FALSE)
  }
  if (length(overrides) > 0L &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("Every parameters entry must have a name.", call. = FALSE)
  }
  parameters <- default_parameters()
  if (length(overrides) > 0L) parameters[names(overrides)] <- overrides
  parameters
}

# 验证公开入口中最常用的参数，并写回内部参数列表。
.asari_prepare_parameters <- function(
    parameters, output, project_name, mode = NULL, ppm = NULL,
    multicores = NULL, rt_align = NULL, database_mode = NULL) {
  if (!is.character(project_name) || length(project_name) != 1L ||
      is.na(project_name) || !nzchar(trimws(project_name)) ||
      grepl("[/\\\\]", project_name)) {
    stop("project_name must be one non-empty name without path separators.", call. = FALSE)
  }
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(trimws(output))) {
    stop("output must be one explicit directory path.", call. = FALSE)
  }

  # 创建用户指定的父目录；每次运行仍由workflow生成独立的时间戳项目目录。
  output <- path.expand(output)
  if (file.exists(output) && !dir.exists(output)) {
    stop("output exists but is not a directory: ", output, call. = FALSE)
  }
  if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", output, call. = FALSE)
  }
  output <- normalizePath(output, mustWork = TRUE)

  # 显式参数优先于parameters列表；未显式给出时保留用户列表或asari默认值。
  if (!is.null(mode)) parameters$mode <- mode
  if (!parameters$mode %in% c("pos", "neg")) {
    stop("mode must be either 'pos' or 'neg'.", call. = FALSE)
  }
  if (!is.null(ppm)) parameters$mz_tolerance_ppm <- ppm
  if (!is.numeric(parameters$mz_tolerance_ppm) ||
      length(parameters$mz_tolerance_ppm) != 1L ||
      is.na(parameters$mz_tolerance_ppm) || parameters$mz_tolerance_ppm <= 0) {
    stop("ppm must be one positive number.", call. = FALSE)
  }
  if (!is.null(multicores)) parameters$multicores <- multicores
  if (!is.numeric(parameters$multicores) || length(parameters$multicores) != 1L ||
      is.na(parameters$multicores) || parameters$multicores < 1 ||
      parameters$multicores != as.integer(parameters$multicores)) {
    stop("multicores must be one positive integer.", call. = FALSE)
  }
  parameters$multicores <- as.integer(parameters$multicores)
  if (!is.null(rt_align)) parameters$rt_align_on <- rt_align
  if (!is.logical(parameters$rt_align_on) || length(parameters$rt_align_on) != 1L ||
      is.na(parameters$rt_align_on)) {
    stop("rt_align must be TRUE or FALSE.", call. = FALSE)
  }

  # 新的便捷入口默认使用内存模式，避免普通LC处理依赖Python pickle转换。
  if (is.null(database_mode)) {
    database_mode <- if ("database_mode" %in% names(parameters)) {
      parameters$database_mode
    } else {
      "memory"
    }
  }
  if (!database_mode %in% c("memory", "ondisk", "auto")) {
    stop("database_mode must be 'memory', 'ondisk', or 'auto'.", call. = FALSE)
  }
  parameters$database_mode <- database_mode
  parameters$project_name <- project_name
  parameters$outdir <- file.path(output, "run")
  parameters
}

# 按Python命令行入口的顺序补齐峰检测派生阈值，同时保留用户明确给出的高级覆盖值。
.asari_finalize_peak_parameters <- function(parameters, input_files, override_names = character()) {
  if (!is.logical(parameters$autoheight) || length(parameters$autoheight) != 1L ||
      is.na(parameters$autoheight)) {
    stop("parameters$autoheight must be TRUE or FALSE.", call. = FALSE)
  }

  # 自动峰高直接使用已经收集的显式文件，避免再次扫描任何默认目录。
  if (isTRUE(parameters$autoheight)) {
    if (!exists("estimate_min_peak_height", mode = "function", inherits = TRUE)) {
      stop(
        "autoheight requires the optional analyze module; set autoheight = FALSE ",
        "or install the complete asariR source tree.",
        call. = FALSE
      )
    }
    estimator <- get("estimate_min_peak_height", mode = "function", inherits = TRUE)
    parameters$min_peak_height <- estimator(input_files, parameters)
    parameters$min_intensity_threshold <- parameters$min_peak_height / 10
  }

  # Python入口会根据最终峰高计算这两个值；高级列表中的明确值具有最高优先级。
  if (!("min_prominence_threshold" %in% override_names)) {
    parameters$min_prominence_threshold <- as.integer(0.33 * parameters$min_peak_height)
  }
  if (!("cal_min_peak_height" %in% override_names)) {
    parameters$cal_min_peak_height <- 10 * parameters$min_peak_height
  }

  # 在进入耗时处理前一次性检查所有关键阈值，给用户明确的参数名称。
  positive_fields <- c(
    "min_peak_height", "min_prominence_threshold",
    "cal_min_peak_height", "min_intensity_threshold"
  )
  for (name in positive_fields) {
    value <- parameters[[name]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
      stop("parameters$", name, " must be one positive number.", call. = FALSE)
    }
  }
  parameters
}

#' 使用asariR处理中心化mzML数据
#'
#' 这是面向普通用户的主入口。函数接收用户明确指定的mzML文件或目录，
#' 运行已经翻译的asari流程，并返回主要结果表的位置。
#'
#' @param input 一个或多个明确的`.mzML`文件，或者包含`.mzML`文件的目录。
#' @param output 存放带时间戳项目目录的父目录。
#' @param project_name 用于输出目录名称的简短项目名。
#' @param mode 可选离子模式，只能是`"pos"`或`"neg"`。
#' @param ppm 可选的正数m/z容差，单位为ppm。
#' @param multicores 可选的正整数工作进程数。
#' @param rt_align 可选逻辑值，控制是否进行保留时间对齐。
#' @param database_mode 数据存储模式。便捷入口默认使用`"memory"`；
#'   `"ondisk"`和部分`"auto"`运行需要Python解释器写兼容pickle。
#' @param recursive 是否递归查找输入目录。
#' @param parameters 来自[asari_default_parameters()]的高级参数覆盖列表；
#'   显式参数的优先级更高。
#'
#' @return 一个`asari_result`对象，包含项目目录、特征表路径、输入文件和有效参数。
#' @export
asari_process <- function(
    input,
    output = "asari-output",
    project_name = "asari_project",
    mode = NULL,
    ppm = NULL,
    multicores = NULL,
    rt_align = NULL,
    database_mode = NULL,
    recursive = FALSE,
    parameters = list()) {
  # 先完成轻量输入校验，再检查处理真实mzML所需的软件依赖。
  input_files <- .asari_collect_input_files(input, recursive = recursive)
  .asari_check_runtime_dependencies()

  # 合并高级参数并应用公开入口中最常用、最容易理解的显式参数。
  effective <- .asari_merge_parameters(parameters)
  if (is.null(database_mode) && !("database_mode" %in% names(parameters))) {
    database_mode <- "memory"
  }
  effective <- .asari_prepare_parameters(
    effective,
    output = output,
    project_name = project_name,
    mode = mode,
    ppm = ppm,
    multicores = multicores,
    rt_align = rt_align,
    database_mode = database_mode
  )
  effective <- .asari_finalize_peak_parameters(
    effective,
    input_files = input_files,
    override_names = names(parameters)
  )

  # 磁盘中间模式需要Python标准库写出与原asari兼容的pickle，提前给出明确错误。
  if (!identical(effective$database_mode, "memory")) {
    .samples_find_python()
  }

  # 记录运行前已有目录，用差集准确定位本次生成的时间戳项目。
  project_pattern <- paste0(
    effective$outdir, "_", effective$project_name, "_*"
  )
  existing_projects <- Sys.glob(project_pattern)
  process_project(input_files, effective)
  generated_projects <- setdiff(Sys.glob(project_pattern), existing_projects)
  if (length(generated_projects) != 1L) {
    stop(
      "Processing finished but the new project directory could not be identified.",
      call. = FALSE
    )
  }

  # 返回最常用结果路径，同时保留输入和有效参数便于复查与复现。
  project_dir <- normalizePath(generated_projects[[1L]], mustWork = TRUE)
  result <- list(
    project_dir = project_dir,
    preferred_feature_table = file.path(project_dir, "preferred_Feature_table.tsv"),
    full_feature_table = file.path(project_dir, "export", "full_Feature_table.tsv"),
    input_files = input_files,
    parameters = effective
  )
  class(result) <- "asari_result"
  print(result)
  invisible(result)
}

#' 返回一份独立的asariR默认参数
#'
#' @return 一份可以安全修改的具名列表。
#' @export
asari_default_parameters <- function() {
  default_parameters()
}

# 打印简洁结果摘要，避免用户自己寻找时间戳输出目录。
#' @export
print.asari_result <- function(x, ...) {
  # 接收并忽略额外打印参数，以符合R的S3打印方法约定。
  dots <- list(...)
  invisible(dots)
  cat("asariR processing completed.\n")
  cat("Project directory: ", x$project_dir, "\n", sep = "")
  cat("Preferred feature table: ", x$preferred_feature_table, "\n", sep = "")
  cat("Full feature table: ", x$full_feature_table, "\n", sep = "")
  invisible(x)
}
