# 面向普通用户的项目页面、交互式向导、命令行和RAW转换入口。

# 转义将要写入HTML的文本，避免样本名或路径破坏页面。
.asari_html_escape <- function(value) {
  value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  gsub("'", "&#39;", value, fixed = TRUE)
}

# 把有限行特征表转为可搜索HTML表格。
.asari_dashboard_table <- function(table, max_rows) {
  shown <- utils::head(table, max_rows)
  headers <- paste0("<th>", .asari_html_escape(names(shown)), "</th>", collapse = "")
  rows <- if (nrow(shown) == 0L) "" else vapply(seq_len(nrow(shown)), function(ii) {
    values <- vapply(shown[ii, , drop = FALSE], function(value) {
      text <- if (length(value) == 0L || is.na(value[[1L]])) "" else as.character(value[[1L]])
      paste0("<td>", .asari_html_escape(text), "</td>")
    }, "")
    paste0("<tr>", paste0(values, collapse = ""), "</tr>")
  }, "")
  paste0(
    "<div class='table-wrap'><table id='feature-table'><thead><tr>", headers,
    "</tr></thead><tbody>", paste0(rows, collapse = ""), "</tbody></table></div>"
  )
}

# 生成一张不依赖外部JavaScript库的m/z-RT散点SVG。
.asari_dashboard_scatter <- function(table, max_points = 3000L) {
  if (!all(c("mz", "rtime") %in% names(table)) || nrow(table) == 0L) {
    return("<p class='empty'>No m/z and retention-time data are available.</p>")
  }
  points <- table[seq_len(min(nrow(table), as.integer(max_points))), , drop = FALSE]
  mz <- suppressWarnings(as.numeric(points$mz))
  rt <- suppressWarnings(as.numeric(points$rtime))
  valid <- is.finite(mz) & is.finite(rt)
  mz <- mz[valid]
  rt <- rt[valid]
  if (length(mz) == 0L) return("<p class='empty'>No numeric m/z and RT values are available.</p>")
  scale_axis <- function(value, low, high) {
    limits <- range(value)
    if (diff(limits) == 0) return(rep((low + high) / 2, length(value)))
    low + (value - limits[[1L]]) / diff(limits) * (high - low)
  }
  x <- scale_axis(mz, 50, 930)
  y <- 330 - scale_axis(rt, 0, 280)
  circles <- paste0(
    "<circle cx='", round(x, 2), "' cy='", round(y, 2),
    "' r='2.2'><title>m/z ", round(mz, 5), "; RT ", round(rt, 2), "</title></circle>",
    collapse = ""
  )
  paste0(
    "<svg viewBox='0 0 980 370' role='img' aria-label='m/z retention time scatter'>",
    "<line x1='50' y1='330' x2='930' y2='330'/><line x1='50' y1='50' x2='50' y2='330'/>",
    circles,
    "<text x='490' y='362'>m/z</text><text x='14' y='190' transform='rotate(-90 14 190)'>Retention time</text>",
    "</svg>"
  )
}

#' 为asariR项目生成自包含的可视化页面
#'
#' 该页面直接读取`project.json`和特征表，不要求Python pickle、
#' Shiny或外部Web服务。
#'
#' @param project 项目目录或[asari_process()]返回的结果。
#' @param table `"preferred"`或`"full"`。
#' @param output 可选HTML输出路径；默认写入项目目录。
#' @param open 是否生成后使用默认浏览器打开。
#' @param max_rows HTML表格最多写入的特征行数。
#' @return 生成的HTML绝对路径。
#' @export
asari_dashboard <- function(
    project,
    table = c("preferred", "full"),
    output = NULL,
    open = FALSE,
    max_rows = 5000L) {
  table <- match.arg(table)
  max_rows <- as.integer(.asari_public_positive_number(max_rows, "max_rows"))
  result <- asari_read_results(project, table = table)
  features <- result[[table]]
  if (is.null(output)) output <- file.path(result$project_dir, paste0("asari_dashboard_", table, ".html"))
  outfile <- .asari_public_output_file(output)
  project_name <- result$project$project_name %||% basename(result$project_dir)
  samples <- max(0L, ncol(features) - 11L)
  mz_range <- if ("mz" %in% names(features) && nrow(features)) {
    paste(round(range(as.numeric(features$mz), na.rm = TRUE), 4), collapse = " - ")
  } else "not available"
  rt_range <- if ("rtime" %in% names(features) && nrow(features)) {
    paste(round(range(as.numeric(features$rtime), na.rm = TRUE), 2), collapse = " - ")
  } else "not available"
  cards <- c(
    paste0("<article><strong>", nrow(features), "</strong><span>Features</span></article>"),
    paste0("<article><strong>", samples, "</strong><span>Sample columns</span></article>"),
    paste0("<article><strong>", .asari_html_escape(mz_range), "</strong><span>m/z range</span></article>"),
    paste0("<article><strong>", .asari_html_escape(rt_range), "</strong><span>RT range</span></article>")
  )
  # 所有CSS和JavaScript直接写入文件，离线打开也能搜索。
  html <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    paste0("<title>", .asari_html_escape(project_name), " - asariR</title>"),
    "<style>body{margin:0;background:#f4f7f8;color:#183033;font:15px system-ui,sans-serif}main{max-width:1200px;margin:auto;padding:32px}h1{margin-bottom:4px}.muted{color:#66797c}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:24px 0}.cards article{background:white;border-radius:12px;padding:18px;box-shadow:0 3px 16px #18303312}.cards strong{display:block;font-size:22px}.cards span{color:#66797c}section{background:white;border-radius:12px;padding:20px;margin:16px 0;box-shadow:0 3px 16px #18303312}svg{width:100%;height:auto}svg line{stroke:#9eb0b3}svg circle{fill:#0b8390;opacity:.55}svg text{fill:#52686b;text-anchor:middle}.table-wrap{overflow:auto;max-height:620px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 9px;border-bottom:1px solid #e4ebec;white-space:nowrap;text-align:right}th{position:sticky;top:0;background:#eaf2f3}th:first-child,td:first-child{text-align:left}input{box-sizing:border-box;width:100%;padding:11px;border:1px solid #b9c8ca;border-radius:8px;margin-bottom:12px}</style></head><body><main>",
    paste0("<h1>", .asari_html_escape(project_name), "</h1>"),
    paste0("<p class='muted'>", .asari_html_escape(result$project_dir), " - ", table, " table</p>"),
    paste0("<div class='cards'>", paste0(cards, collapse = ""), "</div>"),
    paste0("<section><h2>Feature map</h2>", .asari_dashboard_scatter(features), "</section>"),
    "<section><h2>Feature table</h2><input id='filter' placeholder='Search any visible value'>",
    .asari_dashboard_table(features, max_rows),
    if (nrow(features) > max_rows) paste0("<p class='muted'>Showing first ", max_rows, " rows.</p>") else "",
    "</section><script>const box=document.getElementById('filter');box.addEventListener('input',()=>{const q=box.value.toLowerCase();document.querySelectorAll('#feature-table tbody tr').forEach(r=>r.style.display=r.innerText.toLowerCase().includes(q)?'':'none')});</script>",
    "</main></body></html>"
  )
  writeLines(html, outfile, useBytes = TRUE)
  outfile <- normalizePath(outfile, mustWork = TRUE)
  if (isTRUE(open)) utils::browseURL(outfile)
  outfile
}

#' 使用交互式向导准备或运行asariR
#'
#' @param input 明硤的mzML文件或目录；交互会话中为`NULL`时会询问。
#' @param output 输出根目录。
#' @param project_name 项目名。
#' @param mode `"pos"`或`"neg"`。
#' @param workflow `"LC"`、`"GC"`或`"DIMS"`。
#' @param ppm m/z容差。
#' @param multicores 使用的CPU核数。
#' @param run `TRUE`直接运行；`FALSE`只返回将要使用的配置。
#' @param parameters 其他高级参数。
#' @return `run = TRUE`时返回处理结果，否则返回配置list。
#' @export
asari_gui <- function(
    input = NULL,
    output = "asari-output",
    project_name = "asari_project",
    mode = c("pos", "neg"),
    workflow = c("LC", "GC", "DIMS"),
    ppm = 5,
    multicores = 1L,
    run = TRUE,
    parameters = list()) {
  # 基础R的交互式向导可在RStudio和Terminal工作，不强制安装Shiny。
  if (is.null(input)) {
    if (!interactive()) stop("input is required outside an interactive R session.", call. = FALSE)
    input <- readline("mzML file or directory: ")
  }
  mode <- match.arg(mode)
  workflow <- match.arg(workflow)
  configuration <- list(
    input = path.expand(input), output = path.expand(output), project_name = project_name,
    mode = mode, workflow = workflow, ppm = ppm, multicores = as.integer(multicores),
    parameters = parameters
  )
  if (!isTRUE(run)) return(configuration)
  parameters$workflow <- workflow
  if (workflow == "DIMS") {
    return(asari_process_dims(
      input = configuration$input, output = configuration$output,
      project_name = project_name, mode = mode, ppm = ppm,
      multicores = multicores, parameters = parameters
    ))
  }
  asari_process(
    input = configuration$input, output = configuration$output,
    project_name = project_name, mode = mode, ppm = ppm,
    multicores = multicores, parameters = parameters
  )
}

# 打印一份简短命令行帮助，避免用户必须查找长文档。
.asari_cli_help <- function() {
  cat(paste0(
    "asariR commands:\n",
    "  process         Process one mzML file or a directory\n",
    "  analyze         Inspect one mzML file\n",
    "  annotate        Annotate an LC or GC feature table\n",
    "  viz             Create an offline HTML dashboard\n",
    "  join            Merge comma-separated feature tables\n",
    "  list_workflows  Show available workflows\n\n",
    "Common options: --input PATH --output PATH --project NAME --mode pos|neg ",
    "--ppm 5 --multicores 1 --workflow LC|GC|DIMS\n",
    "Annotation options: --db PATH|hmdb4; GC also requires --kovats PATH\n"
  ))
  invisible(NULL)
}

#' 运行asariR命令行入口
#'
#' @param argv 命令行参数字符向量。
#' @return 子命令的结果，不可见返回。
#' @export
asari_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
    .asari_cli_help()
    return(invisible(NULL))
  }
  args <- build_parser(argv)
  if (inherits(args, "asari_version_result")) return(invisible(args))
  if (identical(args$run, "list_workflows")) {
    print(asari_available_workflows(), row.names = FALSE)
    return(invisible(asari_available_workflows()))
  }
  if (is.null(args$input) || !nzchar(args$input)) stop("--input is required for this command.", call. = FALSE)
  parameters <- initialize_parameters(default_parameters(), args)
  parameters <- update_params_from_CLI(parameters, args)
  parameters <- update_peak_detection_params(parameters, args)
  ppm <- if (!is.null(args$ppm)) as.numeric(args$ppm) else parameters$mz_tolerance_ppm
  output <- if (!is.null(args$output)) args$output else parameters$outdir
  project_name <- if (!is.null(args$project)) args$project else parameters$project_name
  result <- if (identical(args$run, "join")) {
    inputs <- trimws(strsplit(args$input, ",", fixed = TRUE)[[1L]])
    if (length(inputs) < 2L) stop("join requires comma-separated paths in --input.")
    asari_merge_feature_tables(inputs, output = output, ppm = ppm)
  } else if (identical(args$run, "process")) {
    if (identical(parameters$workflow, "DIMS")) {
      asari_process_dims(
        args$input, output, project_name, parameters$mode, ppm,
        parameters$multicores, parameters$database_mode, parameters = parameters
      )
    } else {
      asari_process(
        args$input, output, project_name, parameters$mode, ppm,
        parameters$multicores, database_mode = parameters$database_mode,
        parameters = parameters
      )
    }
  } else if (identical(args$run, "analyze")) {
    asari_analyze(
      args$input, ppm = ppm,
      min_intensity = parameters$min_intensity_threshold,
      min_timepoints = parameters$min_timepoints,
      min_peak_height = parameters$min_peak_height,
      parameters = parameters
    )
  } else if (identical(args$run, "viz")) {
    asari_dashboard(
      args$input, table = parameters$table_for_viz,
      output = if (is.null(args$output)) NULL else file.path(output, "asari_dashboard.html")
    )
  } else if (identical(args$run, "annotate")) {
    if (identical(parameters$workflow, "GC")) {
      if (is.null(args$kovats) || is.null(args$db)) stop("GC annotation requires --kovats and --db.")
      asari_annotate_gc(args$input, args$kovats, args$db, output, project_name, parameters$denovo)
    } else {
      asari_annotate_lc(
        args$input, database = args$db %||% "hmdb4",
        output = file.path(output, paste0(project_name, "_LC_annotation.tsv")),
        mode = parameters$mode, ppm = ppm,
        rt_tolerance = parameters$khipu_rtime_tolerance
      )
    }
  } else {
    stop("Unknown or unavailable command: ", args$run, call. = FALSE)
  }
  invisible(result)
}

#' 安装Thermo RAW到mzML转换器
#'
#' @param version ThermoRawFileParser版本。
#' @return 转换器安装目录。
#' @export
asari_install_raw_converter <- function(version = "1.4.5") {
  mzMLconverter_install_converter(version)
}

#' 把Thermo RAW文件批量转换为mzML
#'
#' @param input 一个或多个明硤RAW文件路径。
#' @param output_dir 输出目录；`NULL`时各自写入RAW所在目录。
#' @param version ThermoRawFileParser版本。
#' @param multicores 并行任务数。
#' @param dry_run `TRUE`只返回将执行的命令，不运行转换。
#' @return 转换命令和预期输出路径。
#' @export
asari_convert_raw <- function(
    input,
    output_dir = NULL,
    version = "1.4.5",
    multicores = 1L,
    dry_run = FALSE) {
  if (!is.character(input) || length(input) == 0L) stop("input must contain explicit RAW file paths.")
  files <- vapply(input, function(path) {
    path <- path.expand(path)
    if (!file.exists(path) || dir.exists(path)) stop("RAW file does not exist: ", path)
    normalizePath(path, mustWork = TRUE)
  }, "")
  if (!all(tolower(tools::file_ext(files)) == "raw")) stop("Every input file must use the .raw extension.")
  if (!is.null(output_dir)) {
    output_dir <- path.expand(output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Unable to create output directory: ", output_dir)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
  }
  # dry_run不触发网络下载；它只能使用已经安装的转换器。
  if (isTRUE(dry_run)) {
    root <- getOption("asariR.converter_dir", tools::R_user_dir("asariR", "data"))
    installed <- file.path(root, paste0("ThermoRawFileParser", version))
    executable <- file.path(installed, "ThermoRawFileParser.exe")
    if (!file.exists(executable)) {
      stop(
        "Converter is not installed. Run asari_install_raw_converter() before dry_run.",
        call. = FALSE
      )
    }
    engine <- mzMLconverter__determine_executable(executable)
    command_template <- c(engine, executable, "-i", "INPUT", "-o", "OUTPUT", "-f", "2")
  } else {
    converter <- mzMLconverter(version = version, multicores = multicores)
    command_template <- converter$command_template
  }
  commands <- lapply(files, function(path) {
    command <- gsub("INPUT", path, command_template, fixed = TRUE)
    gsub("OUTPUT", output_dir %||% dirname(path), command, fixed = TRUE)
  })
  outputs <- file.path(
    if (is.null(output_dir)) dirname(files) else output_dir,
    paste0(tools::file_path_sans_ext(basename(files)), ".mzML")
  )
  if (!isTRUE(dry_run)) {
    runner <- getOption("asariR.converter_run")
    if (!is.function(runner)) runner <- function(command) system2(command[[1L]], command[-1L])
    status <- bulk_process(runner, commands, jobs_per_worker = as.integer(multicores))
    failed <- which(vapply(status, function(value) !identical(as.integer(value), 0L), FALSE))
    if (length(failed)) stop("RAW conversion failed for: ", paste(basename(files[failed]), collapse = ", "))
  }
  list(commands = commands, output_files = outputs, dry_run = isTRUE(dry_run))
}

#' 卸载指定版本的Thermo RAW转换器
#'
#' @inheritParams asari_install_raw_converter
#' @return 不可见的`NULL`。
#' @export
asari_uninstall_raw_converter <- function(version = "1.4.5") {
  mzMLconverter_uninstall_converter(version)
}
