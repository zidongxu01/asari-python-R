# 对应 Python asari/workflow.py：项目调度、样本登记和质量轨迹批量提取。
# 本模块只负责流程编排；实际 mzML 读取与质量轨迹算法位于 chromatograms.R。

# 读取 list 或 environment 参数，兼容 Python 可变 dict 与 R 具名 list。
.workflow_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && name %in% names(object)) return(object[[name]])
  default
}

# 修改参数；environment 原地更新，list 返回更新后的副本供 R 调用者继续使用。
.workflow_set <- function(object, name, value) {
  if (is.environment(object)) {
    assign(name, value, envir = object)
    return(object)
  }
  object[[name]] <- value
  object
}

# 按 sample ID 读取或写入 registry，样本 ID 始终保留 Python 0-based 值。
.workflow_registry_get <- function(registry, sample_id) {
  key <- as.character(sample_id)
  if (!is.null(names(registry)) && key %in% names(registry)) {
    return(registry[[key]])
  }
  registry[[as.integer(sample_id) + 1L]]
}

.workflow_registry_set <- function(registry, sample_id, value) {
  key <- as.character(sample_id)
  if (!is.null(names(registry)) && key %in% names(registry)) {
    registry[[key]] <- value
  } else {
    registry[[as.integer(sample_id) + 1L]] <- value
  }
  registry
}

# 获取与 Python time.localtime()[1:6] 一致的月、日、时、分、秒字符串。
.workflow_time_parts <- function(now = Sys.time()) {
  local <- as.POSIXlt(now)
  as.character(c(
    local$mon + 1L, local$mday, local$hour, local$min,
    as.integer(local$sec)
  ))
}

# 从 parameters 注入项或 chromatograms.R 查找质量轨迹提取器。
.workflow_mass_track_extractor <- function(parameters) {
  callback <- .workflow_get(parameters, "extract_mass_tracks")
  if (is.function(callback)) return(callback)
  if (exists("extract_mass_tracks", mode = "function")) return(extract_mass_tracks)
  stop("extract_mass_tracks is unavailable; load chromatograms.R first.")
}

# 调用 mass2chem 等价的 m/z 差值锚点查找器。
.workflow_find_anchor_pairs <- function(list_mass_tracks, parameters) {
  callback <- .workflow_get(parameters, "find_mzdiff_pairs")
  if (is.function(callback)) {
    return(callback(
      list_mass_tracks,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm
    ))
  }
  if (exists(".constructors_find_mzdiff_pairs", mode = "function")) {
    return(.constructors_find_mzdiff_pairs(
      list_mass_tracks,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm
    ))
  }
  stop("m/z anchor-pair finder is unavailable; load constructors.R first.")
}

# 将单样本中间对象写成 Python pickle，保证 samples.R 和 Python 都能读取。
.workflow_write_intermediate <- function(value, path) {
  if (!exists(".experiment_write_pickle", mode = "function")) {
    stop("Python pickle writer is unavailable; load experiment.R first.")
  }
  .experiment_write_pickle(value, path)
}

# 对应 workflow_setup：建立输出目录、批量提取 EIC 并构造 ext_Experiment。
workflow_setup <- function(list_input_files, parameters) {
  sample_registry <- register_samples(list_input_files)

  # auto 模式按项目大小选择 memory 或 ondisk，阈值与 Python 完全相同。
  if (identical(parameters$database_mode, "auto")) {
    database_mode <- if (
        length(list_input_files) <= parameters$project_sample_number_small) {
      "memory"
    } else {
      "ondisk"
    }
    parameters <- .workflow_set(parameters, "database_mode", database_mode)
  }

  time_parts <- .workflow_time_parts()
  parameters <- .workflow_set(
    parameters, "time_stamp", paste(time_parts, collapse = ":")
  )
  parameters <- create_export_folders(
    parameters, paste(time_parts, collapse = "")
  )
  shared_dict <- batch_EIC_from_samples_(sample_registry, parameters)
  if (length(shared_dict) == 0L) {
    stop("No data was processed, check the input files.", call. = FALSE)
  }

  for (sample_id in seq_along(sample_registry) - 1L) {
    sample <- .workflow_registry_get(sample_registry, sample_id)
    result <- .workflow_registry_get(shared_dict, sample_id)
    if (is.null(result)) {
      stop("Missing extraction result for sample ID ", sample_id, call. = FALSE)
    }

    # Python 当前源码把最后三个 tuple 字段错位赋值；这里按 single_sample_EICs_
    # 的实际返回顺序恢复 acquisition_time、sample_data、sparsified，保证流程可运行。
    field_names <- c(
      "status:mzml_parsing", "status:eic", "data_location",
      "max_scan_number", "list_scan_numbers", "list_retention_time",
      "track_mzs", "number_anchor_mz_pairs", "anchor_mz_pairs",
      "acquisition_time", "sample_data", "sparsified"
    )
    for (ii in seq_along(field_names)) sample[[field_names[[ii]]]] <- result[[ii]]
    sample$name <- sub("\\.mzML", "", basename(sample$input_file), fixed = FALSE)
    sample_registry <- .workflow_registry_set(sample_registry, sample_id, sample)
  }

  ext_Experiment__init__(sample_registry, parameters)
}

# 对应 process_project：按 workflow 执行实验处理、导出并清理磁盘中间文件。
process_project <- function(list_input_files, parameters = default_parameters()) {
  experiment <- workflow_setup(list_input_files, parameters)
  workflow <- experiment$parameters$workflow
  message("Processing Experiment Using ", workflow, " Workflow...")

  if (workflow %in% c("LC", "GC")) {
    experiment$process_all_LC()
  } else if (identical(workflow, "DIMS")) {
    experiment$process_all_DIMS()
  } else if (identical(workflow, "LCMSMS")) {
    experiment$process_all_LCMSMS()
  } else {
    message("Workflow not supported. Error 100.")
  }
  experiment$export_all(anno = isTRUE(experiment$parameters$anno))

  if (!identical(experiment$parameters$database_mode, "memory")) {
    remove_intermediate_pickles(experiment$parameters)
  }
  invisible(NULL)
}

# 对应 read_project_dir：返回目录中名称包含 file_pattern 的文件路径。
read_project_dir <- function(directory, file_pattern = ".mzML") {
  message("Working on ~~ ", directory, " ~~ \n")
  entries <- list.files(directory, all.files = FALSE, no.. = TRUE)
  entries <- entries[grepl(file_pattern, entries, fixed = TRUE)]
  file.path(directory, entries)
}

# 对应 read_project_file：从文本文件读取包含 file_pattern 的绝对路径。
read_project_file <- function(project_file, file_pattern = ".mzML") {
  message("Working on ~~ ", project_file, " ~~ \n")
  lines <- readLines(project_file, warn = FALSE)
  selected <- trimws(lines[grepl(file_pattern, lines, fixed = TRUE)])
  vapply(selected, function(path) {
    expanded <- path.expand(path)
    is_absolute <- grepl("^/", expanded) ||
      grepl("^[A-Za-z]:[/\\\\]", expanded)
    candidate <- if (is_absolute) expanded else file.path(getwd(), expanded)
    normalizePath(candidate, mustWork = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

# 对应 register_samples：为输入文件依次分配 Python 0-based sample_id。
register_samples <- function(list_input_files) {
  registry <- lapply(seq_along(list_input_files), function(ii) {
    list(sample_id = ii - 1L, input_file = list_input_files[[ii]])
  })
  names(registry) <- as.character(seq_along(list_input_files) - 1L)
  registry
}

# 对应 create_export_folders：创建项目、export 和 pickle 目录并更新参数。
create_export_folders <- function(parameters, time_stamp = NULL) {
  if (grepl(parameters$project_name, parameters$outdir, fixed = TRUE)) {
    message("Export folders already exist, will not overwrite.")
    return(parameters)
  }

  suffix <- if (is.null(time_stamp)) {
    parameters$time_stamp_for_dir
  } else {
    time_stamp
  }
  outdir <- paste(parameters$outdir, parameters$project_name, suffix, sep = "_")
  parameters <- .workflow_set(parameters, "outdir", outdir)

  # Python os.makedirs 在目标已存在时抛错；R 同样不静默覆盖已有项目。
  if (!dir.create(outdir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Failed to create output directory: ", outdir, call. = FALSE)
  }
  export_dir <- file.path(outdir, "export")
  if (!dir.create(export_dir, showWarnings = FALSE) && !dir.exists(export_dir)) {
    stop("Failed to create export directory: ", export_dir, call. = FALSE)
  }
  pickle_dir <- file.path(outdir, "pickle")
  if (!dir.create(pickle_dir, showWarnings = FALSE) && !dir.exists(pickle_dir)) {
    stop("Failed to create pickle directory: ", pickle_dir, call. = FALSE)
  }
  parameters <- .workflow_set(parameters, "export_outdir", export_dir)
  parameters <- .workflow_set(parameters, "tmp_pickle_dir", pickle_dir)
  parameters
}

# 对应 remove_intermediate_pickles：删除项目 pickle 目录下的临时文件及空目录。
remove_intermediate_pickles <- function(parameters) {
  if (!is.null(parameters$reuse_intermediates)) {
    stop("Cannot remove when reuse_intermediates is set.", call. = FALSE)
  }
  message("Removing temporary pickle files...")
  pickle_dir <- parameters$tmp_pickle_dir
  entries <- list.files(pickle_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  for (entry in entries) {
    if (dir.exists(entry)) {
      stop("Intermediate directory contains a subdirectory: ", entry)
    }
    if (!file.remove(entry)) stop("Failed to remove intermediate: ", entry)
  }
  # 前面已逐项确认并删除文件，此处只移除已经为空的精确 pickle_dir。
  removal_status <- unlink(pickle_dir, recursive = TRUE, force = FALSE)
  if (removal_status != 0L) {
    message("Failed to remove directory ", pickle_dir, ".")
  }
  invisible(NULL)
}

# 对应 make_iter_parameters：生成每个样本的提取任务 tuple。
make_iter_parameters <- function(sample_registry, parameters) {
  lapply(sample_registry, function(sample) {
    sample_name <- sub("\\.mzML", "", basename(sample$input_file))
    outfile <- file.path(
      parameters$outdir, "pickle", paste0(sample_name, ".pickle")
    )
    list(sample$sample_id, sample$input_file, outfile, parameters)
  })
}

# 对应 batch_EIC_from_samples_：批量执行单样本质量轨迹提取并合并结果。
batch_EIC_from_samples_ <- function(sample_registry, parameters) {
  jobs <- make_iter_parameters(sample_registry, parameters)
  workers <- parameters$multicores
  use_parallel <- .Platform$OS.type != "windows" && is.numeric(workers) &&
    length(workers) == 1L && workers > 1L
  results <- if (use_parallel) {
    parallel::mclapply(jobs, single_sample_EICs_, mc.cores = as.integer(workers))
  } else {
    lapply(jobs, single_sample_EICs_)
  }
  merged <- list()
  for (result in results) {
    for (name in names(result)) merged[[name]] <- result[[name]]
  }
  merged
}

# 对应 single_sample_EICs_：提取一个 mzML 的质量轨迹、锚点和扫描信息。
single_sample_EICs_ <- function(job) {
  sample_id <- job[[1L]]
  infile <- job[[2L]]
  outfile <- job[[3L]]
  parameters <- job[[4L]]

  tryCatch({
    extractor <- .workflow_mass_track_extractor(parameters)
    xdict <- extractor(
      infile,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm,
      min_intensity = parameters$min_intensity_threshold,
      min_timepoints = parameters$min_timepoints,
      min_peak_height = parameters$min_peak_height
    )
    tracks <- .workflow_get(xdict, "tracks", list())
    list_mass_tracks <- lapply(seq_along(tracks), function(ii) {
      list(id_number = ii - 1L, mz = tracks[[ii]][[1L]],
           intensity = tracks[[ii]][[2L]])
    })
    message(
      "Extracted ", basename(infile), " to ", length(tracks), " mass tracks."
    )
    anchor_pairs <- .workflow_find_anchor_pairs(list_mass_tracks, parameters)
    rt_numbers <- .workflow_get(xdict, "rt_numbers")
    if (is.null(rt_numbers) || length(rt_numbers) == 0L) {
      stop("No MS1 scan number was extracted.")
    }
    new <- list(
      sample_id = sample_id,
      input_file = infile,
      ion_mode = parameters$mode,
      list_mass_tracks = list_mass_tracks,
      anchor_mz_pairs = anchor_pairs,
      number_anchor_mz_pairs = length(anchor_pairs),
      track_mzs = lapply(seq_along(tracks), function(ii) {
        c(tracks[[ii]][[1L]], ii - 1L)
      }),
      ms2_spectra = .workflow_get(xdict, "ms2_spectra", list()),
      max_scan_number = max(rt_numbers),
      acquisition_time = .workflow_get(xdict, "acquisition_time")
    )

    data_filepath <- outfile
    if (identical(parameters$database_mode, "ondisk")) {
      .workflow_write_intermediate(new, outfile)
    }
    tuple <- list(
      "passed", "passed", data_filepath, new$max_scan_number,
      rt_numbers, .workflow_get(xdict, "rt_times"), new$track_mzs,
      new$number_anchor_mz_pairs, anchor_pairs, new$acquisition_time,
      if (identical(parameters$database_mode, "memory")) new else list(),
      parameters$compress
    )
    result <- list(tuple)
    names(result) <- as.character(sample_id)
    result
  }, error = function(error) {
    message("Failed to extract ", basename(infile), ": ", conditionMessage(error))
    tuple <- list(
      "failed", "failed", NULL, NULL, NULL, NULL, NULL, NULL,
      NULL, NULL, NULL, parameters$compress
    )
    result <- list(tuple)
    names(result) <- as.character(sample_id)
    result
  })
}

# 对应 process_xics：仅提取并保存 XIC，不运行 Experiment/FeatureTable。
process_xics <- function(list_input_files, parameters) {
  parameters <- .workflow_set(parameters, "database_mode", "ondisk")
  time_stamp <- paste(.workflow_time_parts(), collapse = "")
  parameters <- create_export_folders(parameters, time_stamp)

  # Python 当前源码误写为 register_samples(register_samples(files))；
  # R 使用一次登记，避免把 sample ID 当作 mzML 路径。
  batch_EIC_from_samples_(register_samples(list_input_files), parameters)
  message(
    "XICs were stored as pickle objects under ",
    file.path(parameters$outdir, "pickle")
  )
  invisible(NULL)
}

# 对应 get_mz_list：跳过表头，读取首列 tab/comma 分隔的目标 m/z。
get_mz_list <- function(infile) {
  lines <- readLines(infile, warn = FALSE)
  if (length(lines) <= 1L) return(numeric())
  values <- vapply(lines[-1L], function(line) {
    first_tab_field <- strsplit(line, "\t", fixed = TRUE)[[1L]][[1L]]
    first_comma_field <- strsplit(first_tab_field, ",", fixed = TRUE)[[1L]][[1L]]
    value <- suppressWarnings(as.numeric(first_comma_field))
    if (is.na(value)) stop("Invalid m/z value: ", first_comma_field)
    value
  }, numeric(1))
  unname(values)
}
