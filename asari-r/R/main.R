# 对应 Python asari/main.py：命令行参数、子命令分派和入口流程。

# 包版本与 Python asari.__version__ 保持一致。
ASARI_VERSION <- "1.17.0"
SUBCOMMANDS <- c("analyze", "process", "annotate", "viz", "join", "list_workflows")
# 在入口模块内保留独立映射，使main.R不受R包文件加载顺序影响。
booleandict <- list(
  T = TRUE, F = FALSE, `1` = TRUE, `0` = FALSE,
  True = TRUE, False = FALSE, `TRUE` = TRUE, `FALSE` = FALSE,
  true = TRUE, false = FALSE
)

# 判断命令行值是否具有 Python 风格的真值。
.main_truthy <- function(value) {
  !is.null(value) && length(value) > 0L && !identical(value, FALSE) &&
    !(is.numeric(value) && length(value) == 1L && value == 0) &&
    !(is.character(value) && length(value) == 1L && !nzchar(value))
}

# 对应 process：发现输入目录中的mzML文件并启动项目处理。
process <- function(parameters) {
  list_input_files <- read_project_dir(parameters$input)
  if (length(list_input_files) == 0L) {
    cat("No valid mzML files are found in the input directory :(\n")
    return(invisible(NULL))
  }
  process_project(list_input_files, parameters)
}

# 对应 analyze：分析一个显式指定的mzML文件。
analyze <- function(parameters, args) {
  analyze_single_sample(args$input, parameters = parameters)
}

# 对应 annotate：根据参数注释用户输入。
annotate <- function(parameters, args) {
  annotate_project(args$input, parameters = parameters)
}

# 对应 join：Python原版仅保留未实现占位行为。
join <- function(parameters, args) {
  invisible(parameters); invisible(args)
  cat("NOT IMPLEMENTED\n")
  invisible(NULL)
}

# 对应 viz：读取项目并选择完整或推荐特征表。
viz <- function(parameters, args) {
  invisible(parameters)
  tryCatch({
    project <- read_project(args$input)
    project_desc <- project[[1L]]; cmap <- project[[2L]]; epd <- project[[3L]]
    Ftable <- project[[4L]]; Ptable <- project[[5L]]
    if (identical(args$table_for_viz, "full")) {
      dashboard(project_desc, cmap, epd, Ftable)
    } else if (identical(args$table_for_viz, "preferred")) {
      dashboard(project_desc, cmap, epd, Ptable)
    } else {
      cat("Table for viz not recognized! Using Full table.\n")
      dashboard(project_desc, cmap, epd, Ftable)
    }
  }, error = function(error) {
    cat("Error loading Dashboard, likely missing libraries.\n")
    invisible(NULL)
  })
}

# 对应 update_peak_detection_params：更新自动峰高及三个派生阈值。
update_peak_detection_params <- function(parameters, args = NULL) {
  if (isTRUE(parameters$autoheight)) {
    tryCatch({
      parameters$min_peak_height <- estimate_min_peak_height(
        read_project_dir(args$input), parameters
      )
      parameters$min_intensity_threshold <- parameters$min_peak_height / 10
    }, error = function(error) {
      cat("Problems with input files: ", conditionMessage(error),
          ". Back to default min_peak_height.\n", sep = "")
    })
  } else if (!is.null(args) && .main_truthy(args$min_peak_height)) {
    value <- suppressWarnings(as.numeric(args$min_peak_height))
    if (is.na(value)) {
      cat("Problems with specified min_height. Back to default min_peak_height.\n")
    } else {
      parameters$min_peak_height <- value
    }
  }

  # Python int() 对正数向零取整，在R中用as.integer对应。
  parameters$min_prominence_threshold <- as.integer(0.33 * parameters$min_peak_height)
  parameters$cal_min_peak_height <- 10 * parameters$min_peak_height

  for (name in c("min_prominence_threshold", "cal_min_peak_height",
                 "min_intensity_threshold")) {
    if (!is.null(args) && .main_truthy(args[[name]])) {
      value <- suppressWarnings(as.numeric(args[[name]]))
      if (is.na(value)) {
        cat("Problems with specified ", name, ". Back to default ", name, ".\n", sep = "")
      } else {
        parameters[[name]] <- value
      }
    }
  }
  parameters
}

# 对应 update_params_from_CLI 中的内部 __debug_print。
`__debug_print` <- function(debug_print, to_print) {
  if (isTRUE(debug_print)) cat(to_print, "\n", sep = "")
  invisible(NULL)
}

# 读取JSON或YAML参数文件。
.main_read_parameter_file <- function(path) {
  if (!file.exists(path)) stop("Parameter file does not exist: ", path)
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "json")) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading JSON requires jsonlite.")
    return(jsonlite::fromJSON(path, simplifyVector = TRUE))
  }
  if (extension %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) stop("Reading YAML requires the yaml package.")
    return(yaml::read_yaml(path))
  }
  # 未知扩展名先按YAML尝试，失败再按JSON，对应CLI声明的两种格式。
  if (requireNamespace("yaml", quietly = TRUE)) {
    parsed <- tryCatch(yaml::read_yaml(path), error = function(error) NULL)
    if (is.list(parsed)) return(parsed)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading JSON requires jsonlite.")
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

# 要求参数为正数，对应Python的assert。
.main_positive <- function(value, message) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) stop(message)
  value
}

# 对应 update_params_from_CLI：CLI > 参数文件 > 默认参数。
update_params_from_CLI <- function(parameters, args, debug_print = FALSE) {
  if (is.null(args)) stop("No arguments provided.")
  if (is.null(parameters)) stop("No parameters provided.")
  debug <- function(text) `__debug_print`(debug_print, text)

  if (.main_truthy(args$parameters)) {
    supplied <- tryCatch(
      .main_read_parameter_file(args$parameters),
      error = function(error) stop("Failure parsing provided YAML/JSON parameters file: ", conditionMessage(error))
    )
    if (!is.list(supplied)) stop("Parameters file must contain an object/mapping.")
    parameters[names(supplied)] <- supplied
    debug(paste("Updating default parameters from", args$parameters))
  } else debug("Using default parameters")

  if (.main_truthy(args$mode)) {
    if (!args$mode %in% c("pos", "neg")) stop("Mode must be either pos or neg.")
    parameters$mode <- args$mode
  }
  if (.main_truthy(args$ppm)) {
    parameters$ppm <- .main_positive(args$ppm, "PPM must be greater than 0.")
  }

  # 单文件输入在Python中会转成所在目录；文件本身仍保留在args$input供analyze使用。
  if (.main_truthy(args$input)) {
    if (dir.exists(args$input)) parameters$input <- args$input
    else if (file.exists(args$input)) parameters$input <- dirname(args$input)
  }
  if (.main_truthy(args$output)) parameters$outdir <- normalizePath(args$output, mustWork = FALSE)
  else parameters$outdir <- normalizePath(file.path(".", parameters$outdir), mustWork = FALSE)
  if (.main_truthy(args$project)) parameters$project_name <- args$project

  if (!is.null(args$multicores)) {
    available <- parallel::detectCores(logical = TRUE)
    if (is.na(available)) available <- 1L
    parameters$multicores <- if (args$multicores == 0) available else min(available, args$multicores)
  }
  if (.main_truthy(args$reference)) parameters$reference <- args$reference
  if (.main_truthy(args$database_mode)) {
    if (!args$database_mode %in% c("auto", "ondisk", "memory")) {
      stop("Database mode must be either auto, ondisk, or memory.")
    }
    parameters$database_mode <- args$database_mode
  }

  positive_fields <- c(
    wlen = "Wlen must be greater than 0.",
    max_retention_shift = "Max retention shift must be greater than 0.",
    num_lowess_iterations = "Num lowess iterations must be greater than 0.",
    min_peak_height = "Min peak height must be greater than 0.",
    min_prominence_threshold = "Min prominence threshold must be greater than 0.",
    cal_min_peak_height = "Cal min peak height must be greater than 0.",
    min_intensity_threshold = "Min intensity threshold must be greater than 0.",
    vizualization_max_samples = "Visualization max samples must be greater than 0."
  )
  for (name in names(positive_fields)) {
    if (.main_truthy(args[[name]])) {
      parameters[[name]] <- .main_positive(args[[name]], positive_fields[[name]])
    }
  }

  boolean_fields <- c("autoheight", "anno", "debug_rtime_align", "drop_unaligned_samples")
  for (name in boolean_fields) {
    if (.main_truthy(args[[name]])) {
      value <- booleandict[[as.character(args[[name]])]]
      if (is.null(value)) stop(name, " must be either True or False.")
      parameters[[name]] <- value
    }
  }

  if (.main_truthy(args$peak_area)) {
    if (!args$peak_area %in% c("sum", "auc", "gauss")) stop("Peak area must be either sum, auc, or gauss.")
    parameters$peak_area <- args$peak_area
  }
  if (.main_truthy(args$kovats)) parameters$kovats <- args$kovats
  if (.main_truthy(args$db)) parameters$db <- args$db
  if (.main_truthy(args$denovo)) {
    value <- booleandict[[as.character(args$denovo)]]
    if (is.null(value)) stop("denovo must be either True or False.")
    parameters$denovo <- value
  }
  if (.main_truthy(args$denovo_limit)) parameters$max_core_features <- args$denovo_limit
  if (.main_truthy(args$table_for_viz)) {
    if (!args$table_for_viz %in% c("preferred", "full")) stop("Table for viz must be either preferred or full.")
    parameters$table_for_viz <- args$table_for_viz
  }
  if (.main_truthy(args$workflow)) {
    if (!args$workflow %in% c("LC", "GC", "DIMS")) stop("Workflow must be either LC, GC or DIMS.")
    parameters$workflow <- args$workflow
  }
  if (.main_truthy(args$run)) parameters$run <- sub("[[:space:]]+$", "", args$run)

  debug(paste("Using run", parameters$run))
  parameters
}

# 对应 initialize_parameters：写入版本号与当前时间戳。
initialize_parameters <- function(parameters, args) {
  invisible(args)
  parameters$asari_version <- ASARI_VERSION
  parameters$timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  parameters
}

# CLI选项元数据，所有名称与Python argparse定义对应。
.main_cli_spec <- function() {
  list(
    mode = "character", ppm = "integer", input = "character", output = "character",
    project = "character", parameters = "character", multicores = "integer",
    reference = "character", database_mode = "character", wlen = "integer",
    max_retention_shift = "numeric", num_lowess_iterations = "integer",
    autoheight = "character", min_peak_height = "integer",
    min_prominence_threshold = "integer", cal_min_peak_height = "integer",
    min_intensity_threshold = "integer", peak_area = "character",
    keep_intermediates = "character", anno = "character",
    debug_rtime_align = "character", drop_unaligned_samples = "character",
    kovats = "character", db = "character", denovo = "character",
    denovo_limit = "integer", table_for_viz = "character",
    vizualization_max_samples = "integer", workflow = "character"
  )
}

# 转换一个CLI字符串，并在类型无效时给出明确错误。
.main_convert_cli <- function(value, type, name) {
  converted <- switch(type,
    integer = suppressWarnings(as.integer(value)),
    numeric = suppressWarnings(as.numeric(value)),
    character = value
  )
  if ((type %in% c("integer", "numeric")) && is.na(converted)) {
    stop("Invalid value for --", name, ": ", value)
  }
  converted
}

# 对应 build_parser：在不引入额外R包的情况下解析同名CLI参数。
build_parser <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (any(argv %in% c("-v", "--version"))) {
    cat(ASARI_VERSION, "\n", sep = "")
    return(structure(list(version = ASARI_VERSION), class = "asari_version_result"))
  }
  if (length(argv) == 0L) stop("A subcommand is required: ", paste(SUBCOMMANDS, collapse = ", "))

  short_names <- c(m = "mode", i = "input", o = "output", j = "project",
                   p = "parameters", c = "multicores", f = "reference")
  spec <- .main_cli_spec()
  args <- stats::setNames(vector("list", length(spec)), names(spec))
  positional <- character()
  index <- 1L
  while (index <= length(argv)) {
    token <- argv[[index]]
    if (!startsWith(token, "-")) {
      positional <- c(positional, token)
      index <- index + 1L
      next
    }
    inline <- grepl("=", token, fixed = TRUE)
    pieces <- strsplit(sub("^--?", "", token), "=", fixed = TRUE)[[1L]]
    raw_name <- pieces[[1L]]
    name <- if (nchar(raw_name) == 1L && raw_name %in% names(short_names)) short_names[[raw_name]] else raw_name
    if (!name %in% names(spec)) stop("Unknown option: ", token)
    if (inline) value <- paste(pieces[-1L], collapse = "=") else {
      if (index == length(argv)) stop("Option requires a value: ", token)
      index <- index + 1L
      value <- argv[[index]]
    }
    args[[name]] <- .main_convert_cli(value, spec[[name]], name)
    index <- index + 1L
  }
  if (length(positional) != 1L) stop("Exactly one subcommand is required.")
  args$run <- positional[[1L]]
  args
}

# 对应 run_asari：按run字段调用对应子命令。
run_asari <- function(parameters, args = NULL) {
  run <- parameters$run
  if (identical(run, "process")) process(parameters)
  else if (identical(run, "analyze")) analyze(parameters, args)
  else if (identical(run, "annotate")) annotate(parameters, args)
  else if (identical(run, "join")) join(parameters, args)
  else if (identical(run, "viz")) viz(parameters, args)
  else if (identical(run, "list_workflows")) {
    print(asari_available_workflows(), row.names = FALSE)
    invisible(asari_available_workflows())
  }
  else cat("Expecting one of the subcommands: ", paste(SUBCOMMANDS, collapse = ", "), ".\n", sep = "")
}

# 对应 main：创建独立默认参数、应用CLI覆盖并运行子命令。
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  cat("\n\n~~~~~~~ Hello from Asari (", ASARI_VERSION, ") ~~~~~~~~~\n\n", sep = "")
  parameters <- default_parameters()
  args <- build_parser(argv)
  if (inherits(args, "asari_version_result")) return(invisible(args))
  parameters <- initialize_parameters(parameters, args)
  parameters <- update_params_from_CLI(parameters, args)
  parameters <- update_peak_detection_params(parameters, args)
  run_asari(parameters, args)
  invisible(parameters)
}
