# 对应 Python asari/samples.py：轻量级样本对象与中间数据读取。
#
# 数据结构约定：
# 1. Python 的可变 SimpleSample 对象在 R 中使用 environment，保证字段修改会保留。
# 2. Python dict 对应具名 list 或 environment，None 对应 NULL，空 list 对应 list()。
# 3. scan、track 等业务编号保持 Python 原始值；本模块不擅自进行 0/1-based 转换。
# 4. Python pickle 不是 R 原生格式，因此通过可配置的 Python 解释器转换为 JSON 后读取。

# 读取具名 list 或 environment 的字段，缺失时返回给定默认值。
.samples_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && name %in% names(object)) return(object[[name]])
  default
}

# 按 Python registry['field'] 语义读取必需字段；缺失时立即报错。
.samples_registry_get_required <- function(registry, name) {
  sentinel <- new.env(parent = emptyenv())
  value <- .samples_get(registry, name, sentinel)
  if (identical(value, sentinel)) {
    stop("Missing registry field: ", name, call. = FALSE)
  }
  value
}

# 判断质量轨迹容器是否等价于 Python 的非空 list。
.samples_has_items <- function(value) {
  !is.null(value) && length(value) > 0L
}

# 修改实验的 parameters 字段；实验对象必须可变才能保留 Python 的副作用。
.samples_set_experiment_parameter <- function(experiment, name, value) {
  if (!is.environment(experiment)) {
    stop("experiment must be an environment to update parameters.", call. = FALSE)
  }
  parameters <- .samples_get(experiment, "parameters", list())
  parameters[[name]] <- value
  experiment$parameters <- parameters
  invisible(value)
}

# 递归还原 json_tricks 对 numpy ndarray 和 set 使用的标记结构。
.samples_decode_json_tricks <- function(value) {
  if (!is.list(value)) return(value)

  if ("__ndarray__" %in% names(value)) {
    data <- unlist(
      lapply(value$`__ndarray__`, .samples_decode_json_tricks),
      use.names = FALSE
    )
    dtype <- .samples_get(value, "dtype", "")
    if (grepl("^(u?int)", dtype)) data <- as.integer(data)
    if (grepl("^bool", dtype)) data <- as.logical(data)

    # numpy 按行保存多维数组；先反转维度填充，再转置回 Python 轴顺序。
    shape <- as.integer(unlist(.samples_get(value, "shape", length(data))))
    if (length(shape) <= 1L) return(data)
    reversed <- array(data, dim = rev(shape))
    return(aperm(reversed, rev(seq_along(shape))))
  }
  if ("__set__" %in% names(value)) {
    return(lapply(value$`__set__`, .samples_decode_json_tricks))
  }

  lapply(value, .samples_decode_json_tricks)
}

# 将 JSON 文本解析成 R list，并还原 Python json_tricks 的常见特殊类型。
.samples_parse_json <- function(text) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Reading JSON intermediate files requires the optional jsonlite package.",
      call. = FALSE
    )
  }
  parsed <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  .samples_decode_json_tricks(parsed)
}

# 检查文件是否具有 ZIP 的 PK 文件头，避免对普通文件调用 unzip 产生警告。
.samples_is_zipfile <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- readBin(connection, what = "raw", n = 4L)
  length(signature) >= 2L &&
    identical(as.integer(signature[1:2]), c(0x50L, 0x4bL))
}

# 读取 ZIP 第一项的原始字节；顺序与 Python z.namelist()[0] 一致。
.samples_read_first_zip_entry <- function(path) {
  entries <- utils::unzip(path, list = TRUE)
  if (nrow(entries) == 0L) {
    return(list(name = NULL, raw = raw()))
  }

  entry_name <- entries$Name[[1L]]
  entry_size <- as.integer(entries$Length[[1L]])
  connection <- unz(path, entry_name, open = "rb")
  on.exit(close(connection), add = TRUE)
  list(
    name = entry_name,
    raw = readBin(connection, what = "raw", n = entry_size)
  )
}

# 定位 pickle 转换使用的 Python；用户可通过 option 或环境变量明确指定。
.samples_find_python <- function() {
  configured <- getOption("asariR.python", "")
  if (nzchar(configured)) return(configured)

  configured <- Sys.getenv("ASARIR_PYTHON", unset = "")
  if (nzchar(configured)) return(configured)

  candidates <- c(
    file.path(".venv", "bin", "python"),
    file.path("..", ".venv", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Reading Python pickle files requires a Python interpreter. ",
      "Set options(asariR.python = '/path/to/python').",
      call. = FALSE
    )
  }
  existing[[1L]]
}

# 通过 Python 标准库解码 pickle，再转换成 R 能读取的 JSON。
.samples_read_pickle <- function(data_location, zip_entry = NULL) {
  python <- .samples_find_python()
  script_path <- tempfile(fileext = ".py")
  output_path <- tempfile(fileext = ".json")
  on.exit(unlink(c(script_path, output_path)), add = TRUE)

  # 转换器处理 asari 中常见的 numpy 数组、numpy 标量和带 __dict__ 的谱对象。
  script <- c(
    "import json, pickle, sys, zipfile",
    "source, entry, output = sys.argv[1:4]",
    "if entry:",
    "    with zipfile.ZipFile(source, 'r') as archive:",
    "        obj = pickle.loads(archive.read(entry))",
    "else:",
    "    with open(source, 'rb') as handle:",
    "        obj = pickle.load(handle)",
    "def convert(value):",
    "    try:",
    "        import numpy as np",
    "        if isinstance(value, np.ndarray): return value.tolist()",
    "        if isinstance(value, np.generic): return value.item()",
    "    except ImportError:",
    "        pass",
    "    if isinstance(value, (set, tuple)): return list(value)",
    "    if isinstance(value, bytes): return {'__bytes__': list(value)}",
    "    if hasattr(value, '__dict__'):",
    "        payload = dict(value.__dict__)",
    "        scan_time = getattr(value, 'scan_time_in_minutes', None)",
    "        if callable(scan_time): payload['scan_time_in_minutes'] = scan_time()",
    "        for name in ('mz', 'intensity', 'precursor_mz'):",
    "            try: payload[name] = getattr(value, name)",
    "            except Exception: pass",
    "        return payload",
    "    raise TypeError(f'Unsupported pickle value: {type(value)!r}')",
    "with open(output, 'w', encoding='utf-8') as handle:",
    "    json.dump(obj, handle, default=convert, allow_nan=True)"
  )
  writeLines(script, script_path, useBytes = TRUE)

  status <- system2(
    python,
    # system2 会交给 shell 解析参数，因此显式引用含空格的文件路径和空字符串。
    args = shQuote(c(
      script_path,
      data_location,
      if (is.null(zip_entry)) "" else zip_entry,
      output_path
    )),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop(
      "Failed to decode Python pickle: ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }
  .samples_parse_json(paste(readLines(output_path, warn = FALSE), collapse = "\n"))
}

# 从实验对象或全局环境查找 Python 版 extract_ms2 所依赖的可选回调。
.samples_resolve_callback <- function(self, name) {
  callback <- .samples_get(self$experiment, name)
  if (is.function(callback)) return(callback)
  get0(name, envir = .GlobalEnv, mode = "function", inherits = TRUE)
}

# 对应 SimpleSample.__init__：根据 registry 创建可变样本对象。
SimpleSample__init__ <- function(registry = list(),
                                 experiment = NULL,
                                 database_mode = "ondisk",
                                 mode = "pos",
                                 is_reference = FALSE) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("SimpleSample", "environment")

  # 保存构造参数和私有 registry；字段名逐一对应 Python 对象。
  self$experiment <- experiment
  self$mode <- mode
  self$database_mode <- database_mode
  self$is_reference <- is_reference
  self$`.__registry` <- registry

  self$input_file <- .samples_registry_get_required(registry, "input_file")
  self$name <- .samples_registry_get_required(registry, "name")
  self$sample_id <- .samples_registry_get_required(registry, "sample_id")
  self$data_location <- .samples_registry_get_required(registry, "data_location")
  self$track_mzs <- .samples_registry_get_required(registry, "track_mzs")
  self$max_scan_number <- .samples_registry_get_required(
    registry, "max_scan_number"
  )
  self$anchor_mz_pairs <- .samples_registry_get_required(
    registry, "anchor_mz_pairs"
  )
  self$rt_numbers <- .samples_registry_get_required(
    registry, "list_scan_numbers"
  )
  self$list_retention_time <- .samples_registry_get_required(
    registry, "list_retention_time"
  )

  # memory 模式直接持有质量轨迹；ondisk 模式只在需要时读取。
  if (identical(database_mode, "memory")) {
    sample_data <- .samples_registry_get_required(registry, "sample_data")
    self$list_mass_tracks <- .samples_registry_get_required(
      sample_data, "list_mass_tracks"
    )
  } else {
    self$list_mass_tracks <- list()
  }

  # m/z landmark 展平逻辑复用 mass_functions.R，并初始化 RT 校准状态。
  self$`_mz_landmarks_` <- flatten_tuplelist(self$anchor_mz_pairs)
  self$rt_landmarks <- list()
  self$rt_cal_dict <- NULL
  self$reverse_rt_cal_dict <- NULL
  self$is_rt_aligned <- is_reference
  self$mz_calibration_function <- NULL

  # active binding 对应 Python @property：每次访问都重新读取私有 registry。
  makeActiveBinding(
    "list_scan_numbers",
    local({
      sample <- self
      function(value) {
        if (!missing(value)) {
          stop("list_scan_numbers is a read-only property.", call. = FALSE)
        }
        SimpleSample_list_scan_numbers(sample)
      }
    }),
    self
  )

  # 为 environment 安装方法闭包，使 constructors.R 可按 Python 对象方式调用。
  self$get_mass_tracks_for_sample <- function(sample) {
    SimpleSample_get_mass_tracks_for_sample(sample)
  }
  self$get_masstracks_and_anchors <- function() {
    SimpleSample_get_masstracks_and_anchors(self)
  }
  self$get_rt_calibration_records <- function() {
    SimpleSample_get_rt_calibration_records(self)
  }
  self$extract_ms2 <- function(export_format = "msp") {
    SimpleSample_extract_ms2(self, export_format)
  }
  self$`_get_sample_data` <- function() SimpleSample__get_sample_data(self)
  self$`_retrieve_from_disk` <- function() SimpleSample__retrieve_from_disk(self)

  self
}

# 对应 @property SimpleSample.list_scan_numbers：始终从原始 registry 读取。
SimpleSample_list_scan_numbers <- function(self) {
  .samples_registry_get_required(self$`.__registry`, "list_scan_numbers")
}

# 对应 @staticmethod SimpleSample.get_mass_tracks_for_sample。
SimpleSample_get_mass_tracks_for_sample <- function(sample) {
  SimpleSample_get_masstracks_and_anchors(sample)
}

# 对应 SimpleSample.get_masstracks_and_anchors：优先使用内存，否则按需读盘。
SimpleSample_get_masstracks_and_anchors <- function(self) {
  if (.samples_has_items(self$list_mass_tracks)) {
    return(self$list_mass_tracks)
  }
  sample_data <- SimpleSample__get_sample_data(self)
  .samples_registry_get_required(sample_data, "list_mass_tracks")
}

# 对应 SimpleSample.get_rt_calibration_records：导出 RT 校准所需的四个字段。
SimpleSample_get_rt_calibration_records <- function(self) {
  list(
    sample_id = self$sample_id,
    name = self$name,
    rt_landmarks = self$rt_landmarks,
    reverse_rt_cal_dict = self$reverse_rt_cal_dict
  )
}

# 对应 SimpleSample.extract_ms2：转换 MS2 谱并调用可选保存器。
SimpleSample_extract_ms2 <- function(self, export_format = "msp") {
  tryCatch({
    # Python 版在 memory 模式从 registry 读取，其他模式从中间文件读取。
    if (identical(self$database_mode, "memory")) {
      sample_data <- .samples_registry_get_required(
        self$`.__registry`, "sample_data"
      )
    } else {
      sample_data <- SimpleSample__get_sample_data(self)
    }
    ms2_data <- .samples_registry_get_required(sample_data, "ms2_spectra")

    Spectrum <- .samples_resolve_callback(self, "Spectrum")
    save_spectra <- .samples_resolve_callback(self, "save_spectra")

    spectra <- lapply(ms2_data, function(spec) {
      # Python 只在实际遍历谱图时解析 Spectrum；空谱列表不需要该构造器。
      if (!is.function(Spectrum)) stop("Spectrum callback is unavailable.")
      scan_time <- .samples_get(spec, "scan_time_in_minutes")
      if (is.null(scan_time)) stop("MS2 spectrum has no scan time.")
      if (is.function(scan_time)) scan_time <- scan_time()
      rtime <- scan_time * 60
      precursor_mz <- tryCatch(
        .samples_get(spec, "precursor_mz", NULL),
        error = function(error) NULL
      )

      # 与 Python zip(spec.mz, spec.intensity) 一样，只保留两者共同长度。
      mzs <- .samples_get(spec, "mz", numeric())
      intensities <- .samples_get(spec, "intensity", numeric())
      common_length <- min(length(mzs), length(intensities))
      if (common_length == 0L) {
        mzs <- intensities <- numeric()
      } else {
        mzs <- as.numeric(mzs[seq_len(common_length)])
        intensities <- as.numeric(intensities[seq_len(common_length)])
      }

      Spectrum(
        mz = mzs,
        intensities = intensities,
        metadata = list(
          scan_time = rtime,
          origin = self$name,
          precursor_mz = precursor_mz
        )
      )
    })

    # Python 仅移除开头的一个点，并把最终格式写回实验参数。
    if (!is.character(export_format) || length(export_format) != 1L ||
        !nzchar(export_format)) {
      stop("export_format must contain at least one character.")
    }
    if (startsWith(export_format, ".")) {
      export_format <- substring(export_format, 2L)
    }
    .samples_set_experiment_parameter(
      self$experiment, "ms2_export_format", export_format
    )
    output_dir <- .samples_registry_get_required(
      self$experiment$parameters, "ms2_spectra_outdir"
    )
    path <- file.path(
      output_dir,
      sprintf("ms2_%s.%s", self$name, export_format)
    )
    if (!is.function(save_spectra)) {
      stop("save_spectra callback is unavailable.")
    }
    save_spectra(spectra, path, export_style = "matchms")
    invisible(NULL)
  }, error = function(error) {
    # Python 捕获所有异常，只报告当前样本，不阻断主流程。
    message("Error Extracting MS2 for: ", self$name)
    invisible(NULL)
  })
}

# 对应 SimpleSample._get_sample_data：保留未来数据库模式扩展所需的包装层。
SimpleSample__get_sample_data <- function(self) {
  SimpleSample_load_intermediate(self$data_location)
}

# 对应 SimpleSample._retrieve_from_disk：从当前样本的数据位置读取中间文件。
SimpleSample__retrieve_from_disk <- function(self) {
  SimpleSample_load_intermediate(self$data_location)
}

# 对应 @staticmethod SimpleSample.load_intermediate：读取 JSON、pickle 或其 ZIP。
SimpleSample_load_intermediate <- function(data_location) {
  if (length(data_location) != 1L || !is.character(data_location) ||
      !file.exists(data_location)) {
    stop("Failed to load sample data from: ", data_location, call. = FALSE)
  }

  sample_data <- NULL
  if (.samples_is_zipfile(data_location)) {
    first_entry <- .samples_read_first_zip_entry(data_location)
    entry_name <- first_entry$name

    # Python 根据 ZIP 第一项的扩展名决定解码器。
    if (!is.null(entry_name) && endsWith(entry_name, ".pickle")) {
      sample_data <- .samples_read_pickle(data_location, entry_name)
    } else if (!is.null(entry_name) && endsWith(entry_name, ".json")) {
      json_text <- rawToChar(first_entry$raw)
      sample_data <- .samples_parse_json(json_text)
    }
  } else if (endsWith(data_location, ".pickle")) {
    sample_data <- .samples_read_pickle(data_location)
  } else if (endsWith(data_location, ".json")) {
    json_text <- paste(readLines(data_location, warn = FALSE), collapse = "\n")
    sample_data <- .samples_parse_json(json_text)
  } else {
    stop("Unknown file format: ", data_location, call. = FALSE)
  }

  if (is.null(sample_data)) {
    stop("Failed to load sample data from: ", data_location, call. = FALSE)
  }
  sample_data
}

# 使用接近 Python 类名的构造器别名，供 constructors.R 的 sample_factory 直接使用。
SimpleSample <- SimpleSample__init__

# 保留早期 smoke_check.R 使用的兼容构造器；它不属于 Python 的 9 个 def。
new_sample <- function(name, infile) {
  list(
    name = name,
    infile = infile,
    mass_tracks = NULL,
    peaks = NULL
  )
}
