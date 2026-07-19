# samples.py 的 9 个 def 必须逐一存在明确的 R 对应函数。
test_that("all 9 Python samples defs have explicit R counterparts", {
  expected_functions <- c(
    "SimpleSample__init__",
    "SimpleSample_list_scan_numbers",
    "SimpleSample_get_mass_tracks_for_sample",
    "SimpleSample_get_masstracks_and_anchors",
    "SimpleSample_get_rt_calibration_records",
    "SimpleSample_extract_ms2",
    "SimpleSample__get_sample_data",
    "SimpleSample__retrieve_from_disk",
    "SimpleSample_load_intermediate"
  )

  # 数量必须恰好为 9，并且每个名称都解析为函数。
  expect_length(expected_functions, 9L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# 创建 tests 使用的最小 registry；字段与 workflow.register_samples 输出一致。
make_sample_registry <- function(data_location = tempfile(fileext = ".json"),
                                 sample_data = list(
                                   list_mass_tracks = list(
                                     list(id_number = 0L, mz = 100)
                                   ),
                                   ms2_spectra = list()
                                 )) {
  list(
    input_file = "sample_1.mzML",
    name = "sample_1",
    sample_id = 0L,
    data_location = data_location,
    track_mzs = list(c(100, 0L)),
    max_scan_number = 3L,
    anchor_mz_pairs = list(c(0L, 1L), c(1L, 2L)),
    list_scan_numbers = 0:2,
    list_retention_time = c(0, 5, 10),
    sample_data = sample_data
  )
}

# 创建可变实验对象，以保留 extract_ms2 对 parameters 的更新。
make_sample_experiment <- function(output_dir = tempdir()) {
  experiment <- new.env(parent = emptyenv())
  experiment$parameters <- list(ms2_spectra_outdir = output_dir)
  experiment
}

# 验证构造器逐字段复制 registry，并正确初始化 Python 对象状态。
test_that("SimpleSample constructor preserves fields and mutable state", {
  registry <- make_sample_registry()
  experiment <- make_sample_experiment()
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory",
    mode = "neg",
    is_reference = TRUE
  )

  expect_s3_class(sample, "SimpleSample")
  expect_identical(sample$experiment, experiment)
  expect_equal(sample$mode, "neg")
  expect_true(sample$is_reference)
  expect_equal(sample$input_file, registry$input_file)
  expect_equal(sample$name, registry$name)
  expect_equal(sample$sample_id, 0L)
  expect_equal(sample$rt_numbers, 0:2)
  expect_equal(sample$list_scan_numbers, 0:2)
  expect_equal(sample$list_mass_tracks, registry$sample_data$list_mass_tracks)
  expect_equal(sample$`_mz_landmarks_`, c(0L, 1L, 2L))
  expect_equal(sample$rt_landmarks, list())
  expect_null(sample$rt_cal_dict)
  expect_null(sample$reverse_rt_cal_dict)
  expect_true(sample$is_rt_aligned)
  expect_null(sample$mz_calibration_function)
})

# Python property 每次从私有 registry 读取，而不是使用可能过期的公开副本。
test_that("list_scan_numbers property reads the current registry value", {
  registry <- list2env(make_sample_registry(), parent = emptyenv())
  sample <- SimpleSample__init__(registry, database_mode = "memory")

  registry$list_scan_numbers <- c(4L, 5L, 6L)

  expect_equal(SimpleSample_list_scan_numbers(sample), c(4L, 5L, 6L))
  expect_equal(sample$list_scan_numbers, c(4L, 5L, 6L))
  expect_error(
    sample$list_scan_numbers <- 7:9,
    "read-only property"
  )
})

# memory 模式应直接返回同一个质量轨迹，静态入口和实例入口结果一致。
test_that("mass tracks are returned directly when already in memory", {
  registry <- make_sample_registry()
  sample <- SimpleSample__init__(registry, database_mode = "memory")

  expect_identical(
    SimpleSample_get_masstracks_and_anchors(sample),
    registry$sample_data$list_mass_tracks
  )
  expect_identical(
    SimpleSample_get_mass_tracks_for_sample(sample),
    registry$sample_data$list_mass_tracks
  )
  expect_identical(
    sample$get_masstracks_and_anchors(),
    registry$sample_data$list_mass_tracks
  )
})

# ondisk 模式不缓存质量轨迹，而是通过 _get_sample_data 按需加载 JSON。
test_that("mass tracks are loaded lazily in ondisk mode", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    '{"list_mass_tracks":[{"id_number":0,"mz":101.25}]}',
    json_path
  )
  registry <- make_sample_registry(data_location = json_path)
  sample <- SimpleSample__init__(registry, database_mode = "ondisk")

  expect_equal(sample$list_mass_tracks, list())
  tracks <- SimpleSample_get_masstracks_and_anchors(sample)
  expect_equal(tracks[[1L]]$id_number, 0L)
  expect_equal(tracks[[1L]]$mz, 101.25)
  expect_equal(SimpleSample__get_sample_data(sample),
               SimpleSample__retrieve_from_disk(sample))
})

# 校准记录必须包含 Python 原函数返回的四个键，并反映对象的最新状态。
test_that("RT calibration records contain current sample state", {
  sample <- SimpleSample__init__(make_sample_registry(), database_mode = "memory")
  sample$rt_landmarks <- c(10L, 20L)
  sample$reverse_rt_cal_dict <- c(`10` = 11L, `20` = 19L)

  records <- SimpleSample_get_rt_calibration_records(sample)

  expect_named(
    records,
    c("sample_id", "name", "rt_landmarks", "reverse_rt_cal_dict")
  )
  expect_equal(records$sample_id, 0L)
  expect_equal(records$name, "sample_1")
  expect_equal(records$rt_landmarks, c(10L, 20L))
  expect_equal(records$reverse_rt_cal_dict, c(`10` = 11L, `20` = 19L))
})

# 普通 JSON 和 ZIP 中第一项 JSON 都应恢复相同的 Python dict/list 结构。
test_that("load_intermediate reads JSON and zipped JSON", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    '{"list_mass_tracks":[{"id_number":0,"mz":100.5}],"flag":true}',
    json_path
  )
  direct <- SimpleSample_load_intermediate(json_path)

  zip_dir <- tempfile("sample_zip_")
  dir.create(zip_dir)
  zip_json <- file.path(zip_dir, "sample.json")
  file.copy(json_path, zip_json)
  zip_path <- tempfile(fileext = ".zip")
  old_directory <- setwd(zip_dir)
  on.exit(setwd(old_directory), add = TRUE)
  suppressWarnings(utils::zip(zip_path, "sample.json"))
  zipped <- SimpleSample_load_intermediate(zip_path)

  expect_equal(direct, zipped)
  expect_equal(direct$list_mass_tracks[[1L]]$mz, 100.5)
  expect_true(direct$flag)
})

# json_tricks 保存的 numpy 标记必须恢复成 R 向量，而不是暴露编码器内部字典。
test_that("load_intermediate decodes common json_tricks numpy markers", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    paste0(
      '{"intensity":{"__ndarray__":[1,2,3],',
      '"dtype":"int64","shape":[3]}}'
    ),
    json_path
  )

  decoded <- SimpleSample_load_intermediate(json_path)

  expect_equal(decoded$intensity, c(1L, 2L, 3L))
})

# Python pickle 和 ZIP 中 pickle 都通过转换桥接恢复常见 asari 数据类型。
test_that("load_intermediate reads Python pickle and zipped pickle", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle parity testing.")

  work_dir <- tempfile("pickle_test_")
  dir.create(work_dir)
  pickle_path <- file.path(work_dir, "sample.pickle")
  zip_path <- file.path(work_dir, "sample.zip")
  script_path <- file.path(work_dir, "make_pickle.py")
  writeLines(c(
    "import pickle, sys, zipfile",
    "path, zip_path = sys.argv[1:3]",
    "data = {'list_mass_tracks': [{'id_number': 0, 'mz': 100.25}], 'values': (1, 2)}",
    "with open(path, 'wb') as handle: pickle.dump(data, handle)",
    "with zipfile.ZipFile(zip_path, 'w') as archive: archive.write(path, 'sample.pickle')"
  ), script_path)
  output <- system2(
    python,
    shQuote(c(script_path, pickle_path, zip_path)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"))

  direct <- SimpleSample_load_intermediate(pickle_path)
  zipped <- SimpleSample_load_intermediate(zip_path)

  expect_equal(direct, zipped)
  expect_equal(direct$list_mass_tracks[[1L]]$mz, 100.25)
  expect_equal(unlist(direct$values), c(1L, 2L))
})

# 使用可注入 Spectrum/save_spectra 回调验证 MS2 转换、路径和参数副作用。
test_that("extract_ms2 converts spectra and calls the configured exporter", {
  captured <- new.env(parent = emptyenv())
  experiment <- make_sample_experiment(tempdir())
  experiment$Spectrum <- function(mz, intensities, metadata) {
    list(mz = mz, intensities = intensities, metadata = metadata)
  }
  experiment$save_spectra <- function(spectra, path, export_style) {
    captured$spectra <- spectra
    captured$path <- path
    captured$export_style <- export_style
    invisible(NULL)
  }
  ms2_spectrum <- list(
    scan_time_in_minutes = function() 1.5,
    mz = c(50, 75),
    intensity = c(1000, 500),
    precursor_mz = 150
  )
  registry <- make_sample_registry(sample_data = list(
    list_mass_tracks = list(),
    ms2_spectra = list(ms2_spectrum)
  ))
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory"
  )

  expect_null(SimpleSample_extract_ms2(sample, ".msp"))
  expect_equal(experiment$parameters$ms2_export_format, "msp")
  expect_equal(captured$path, file.path(tempdir(), "ms2_sample_1.msp"))
  expect_equal(captured$export_style, "matchms")
  expect_equal(captured$spectra[[1L]]$mz, c(50, 75))
  expect_equal(captured$spectra[[1L]]$intensities, c(1000, 500))
  expect_equal(captured$spectra[[1L]]$metadata$scan_time, 90)
  expect_equal(captured$spectra[[1L]]$metadata$origin, "sample_1")
  expect_equal(captured$spectra[[1L]]$metadata$precursor_mz, 150)
})

# 与 Python 的 try/except 一致：缺少 MS2 导出依赖只提示，不抛出错误。
test_that("extract_ms2 reports optional dependency failure without stopping", {
  experiment <- make_sample_experiment()
  registry <- make_sample_registry(sample_data = list(
    list_mass_tracks = list(),
    ms2_spectra = list()
  ))
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory"
  )

  expect_message(
    expect_null(SimpleSample_extract_ms2(sample)),
    "Error Extracting MS2 for: sample_1"
  )
  # 空谱列表会先更新格式，再在 Python 的 save_spectra 名称解析处失败。
  expect_equal(experiment$parameters$ms2_export_format, "msp")
})

# 未知扩展名、ZIP 第一项未知以及缺失 registry 字段都必须明确失败。
test_that("invalid intermediate files and incomplete registries fail clearly", {
  unknown_path <- tempfile(fileext = ".txt")
  writeLines("not an intermediate", unknown_path)
  expect_error(
    SimpleSample_load_intermediate(unknown_path),
    "Unknown file format"
  )
  expect_error(
    SimpleSample_load_intermediate(tempfile(fileext = ".json")),
    "Failed to load sample data"
  )
  expect_error(
    SimpleSample__init__(list()),
    "Missing registry field: input_file"
  )
})
