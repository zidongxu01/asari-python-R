# experiment.py 的 20 个 def 必须逐一存在明确的 R 对应函数。
test_that("all 20 Python experiment defs have explicit R counterparts", {
  expected_functions <- c(
    "ext_Experiment__init__",
    "ext_Experiment_get_reference_sample_id",
    "ext_Experiment_get_valid_sample_ids",
    "ext_Experiment_determine_acquisition_order",
    "ext_Experiment_get_max_scan_number",
    "ext_Experiment_process_all_LC",
    "ext_Experiment_process_all_DIMS",
    "ext_Experiment_process_all_LCMSMS",
    "ext_Experiment_export_all",
    "ext_Experiment_annotate",
    "ext_Experiment_generate_qc_plot_pdf",
    "ext_Experiment_export_CMAP_pickle",
    "ext_Experiment_load_annotation_db",
    "ext_Experiment_db_mass_calibrate",
    "ext_Experiment_append_orphans_to_epmCpds",
    "ext_Experiment_export_peak_annotation",
    "ext_Experiment_select_unique_compound_features",
    "ext_Experiment_export_feature_tables",
    "ext_Experiment_export_log",
    "ext_Experiment_export_readme"
  )

  expect_length(expected_functions, 20L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# 生成 experiment 测试使用的完整参数，避免依赖尚未扩充的默认参数占位文件。
make_experiment_parameters <- function(outdir = tempfile("experiment_out_"),
                                       reference = NULL,
                                       workflow = "LC") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "export"), showWarnings = FALSE)
  list(
    outdir = outdir,
    mode = "pos",
    mz_tolerance_ppm = 5,
    check_isotope_ratio = FALSE,
    database_mode = "memory",
    reference = reference,
    workflow = workflow,
    rt_align_on = FALSE,
    drop_unaligned_samples = FALSE,
    output_feature_table = "Feature_table.tsv",
    mass_grid_mapping = "_mass_grid_mapping.csv",
    anno = FALSE,
    project_sample_number_small = 10,
    correction_tolerance_ppm = 1,
    cal_min_peak_height = 1000,
    peak_number_rt_calibration = 2,
    num_lowess_iterations = 3,
    max_retention_shift = NULL,
    debug_rtime_align = FALSE,
    peak_area = "sum",
    min_peak_height = 1000,
    min_peak_ratio = 0.001,
    min_timepoints = 3,
    min_intensity_threshold = 100,
    wlen = 5,
    signal_noise_ratio = 2,
    gaussian_shape = 0.5,
    min_prominence_threshold = 100
  )
}

# 创建 Python workflow.register_samples 风格的两个样本登记项。
make_experiment_registry <- function() {
  make_one <- function(id, name, anchors, status = "passed", time = NULL) {
    mz_offset <- id * 0.0001
    mzs <- 100:104 + mz_offset
    intensities <- lapply(seq_along(mzs), function(ii) {
      list(
        id_number = ii - 1L,
        mz = mzs[[ii]],
        intensity = as.integer(c(0L, 10L, 100L, 10L, 0L) / ii)
      )
    })
    list(
      input_file = paste0(name, ".mzML"),
      name = name,
      sample_id = id,
      data_location = tempfile(fileext = ".pickle"),
      track_mzs = lapply(seq_along(mzs), function(ii) {
        c(mzs[[ii]], ii - 1L)
      }),
      max_scan_number = 4L,
      anchor_mz_pairs = list(c(0L, 1L), c(2L, 3L)),
      number_anchor_mz_pairs = anchors,
      list_scan_numbers = 0:4,
      list_retention_time = seq(0, 20, 5),
      `status:eic` = status,
      acquisition_time = time,
      sample_data = list(list_mass_tracks = intensities, ms2_spectra = list())
    )
  }
  list(
    `0` = make_one(0L, "sample_a", 2L, time = "2024-01-02"),
    `1` = make_one(1L, "sample_b", 5L, time = "2024-01-01"),
    # 提取失败的 Python registry 通常没有 anchor 数，因此不参与参考样本候选。
    `2` = make_one(2L, "failed", NULL, status = "failed", time = "2024-01-03")
  )
}

# 构造器应建立正式类身份、有效样本清单、共享样本别名和自动参考样本。
test_that("experiment constructor initializes Python fields and aliases", {
  registry <- make_experiment_registry()
  parameters <- make_experiment_parameters()
  experiment <- ext_Experiment__init__(registry, parameters)

  expect_s3_class(experiment, "ext_Experiment")
  expect_equal(unlist(experiment$valid_sample_ids), c(0L, 1L))
  expect_equal(experiment$number_of_samples, 2L)
  expect_equal(experiment$number_scans, 5)
  expect_equal(experiment$reference_sample_id, 1L)
  expect_equal(experiment$parameters$reference, "sample_b.mzML")
  expect_equal(experiment$output_dir, parameters$outdir)
  expect_equal(experiment$mode, "pos")
  expect_equal(experiment$mz_tolerance_ppm, 5)
  expect_equal(experiment$database_mode, "memory")

  marker <- new.env(parent = emptyenv())
  experiment$all_samples <- list(marker)
  expect_identical(experiment$all_sample_instances[[1L]], marker)
  experiment$all_sample_instances <- list()
  expect_length(experiment$all_samples, 0L)
})

# 用户指定参考文件时支持完整文件名及省略 .mzML 的名称。
test_that("reference selection honors user file names", {
  registry <- make_experiment_registry()
  exact <- ext_Experiment__init__(
    registry,
    make_experiment_parameters(reference = "sample_a.mzML")
  )
  without_extension <- ext_Experiment__init__(
    registry,
    make_experiment_parameters(reference = "sample_b")
  )

  expect_equal(exact$reference_sample_id, 0L)
  expect_equal(without_extension$reference_sample_id, 1L)
})

# 有效样本、采集顺序和最大 scan 必须与 Python 的筛选和排序规则一致。
test_that("valid IDs acquisition order and max scan follow registry", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  order_by_time <- ext_Experiment_determine_acquisition_order(experiment)

  expect_equal(
    vapply(order_by_time, function(pair) as.character(pair[[1L]]), character(1)),
    c("1", "0", "2")
  )
  expect_equal(ext_Experiment_get_max_scan_number(
    experiment, experiment$sample_registry
  ), 5)

  experiment$sample_registry[["0"]]$acquisition_time <- NULL
  fallback <- ext_Experiment_determine_acquisition_order(experiment)
  expect_equal(
    vapply(fallback, function(pair) as.character(pair[[1L]]), character(1)),
    c("0", "1", "2")
  )
})

# mock CompositeMap 记录 LC 和 DIMS 主流程调用顺序，同时验证 LCMSMS 空占位。
test_that("processing methods call CompositeMap stages in Python order", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  calls <- character()
  experiment$CompositeMap_factory <- function(experiment) {
    cmap <- new.env(parent = emptyenv())
    cmap$construct_mass_grid <- function() calls <<- c(calls, "grid")
    cmap$mock_rentention_alignment <- function() calls <<- c(calls, "mock_rt")
    cmap$build_composite_tracks <- function() calls <<- c(calls, "tracks")
    cmap$global_peak_detection <- function() calls <<- c(calls, "peaks")
    cmap$get_DIMS_feature_table <- function() calls <<- c(calls, "dims_table")
    cmap
  }

  expect_null(ext_Experiment_process_all_LC(experiment))
  expect_equal(calls, c("grid", "mock_rt", "tracks", "peaks"))
  calls <- character()
  expect_null(ext_Experiment_process_all_DIMS(experiment))
  expect_equal(calls, c("grid", "dims_table"))
  expect_null(ext_Experiment_process_all_LCMSMS(experiment))
})

# 不使用 mock，验证 experiment、samples 和 constructors 能完成真实 DIMS 串联。
test_that("real DIMS integration connects samples constructors and experiment", {
  registry <- make_experiment_registry()[c("0", "1")]
  experiment <- ext_Experiment__init__(
    registry,
    make_experiment_parameters(workflow = "DIMS")
  )

  expect_null(ext_Experiment_process_all_DIMS(experiment))

  expect_s3_class(experiment$CMAP, "CompositeMap")
  expect_equal(nrow(experiment$CMAP$MassGrid), 5L)
  expect_equal(nrow(experiment$CMAP$FeatureTable), 5L)
  expect_length(experiment$all_samples, 2L)
  expect_true(all(c("sample_a", "sample_b") %in%
                  names(experiment$CMAP$FeatureTable)))
})

# 使用满足 selectivity 前提的真实轨迹完成 LC 质量网格、合成轨迹和峰检测。
test_that("real LC integration reaches global feature table", {
  registry <- make_experiment_registry()[c("0", "1")]
  experiment <- ext_Experiment__init__(
    registry,
    make_experiment_parameters(workflow = "LC")
  )

  expect_null(ext_Experiment_process_all_LC(experiment))

  expect_equal(nrow(experiment$CMAP$MassGrid), 5L)
  expect_length(experiment$CMAP$composite_mass_tracks, 5L)
  expect_true(is.data.frame(experiment$CMAP$FeatureTable))
})

# 创建 FeatureTable 和样本对象，供表格、日志及 export_all 测试共用。
attach_export_fixture <- function(experiment) {
  sample_a <- new.env(parent = emptyenv())
  sample_a$name <- "sample_a"
  sample_a$is_rt_aligned <- TRUE
  sample_b <- new.env(parent = emptyenv())
  sample_b$name <- "sample_b"
  sample_b$is_rt_aligned <- FALSE
  experiment$all_samples <- list(sample_a, sample_b)

  cmap <- new.env(parent = emptyenv())
  cmap$MassGrid <- data.frame(mz = c(100, 101), sample_a = c(0, 1),
                              sample_b = c(0, 1), check.names = FALSE)
  cmap$FeatureList <- list(
    list(id_number = "F0", mz = 100.123456, rtime = 5.678,
         rtime_left_base = 4.444, rtime_right_base = 6.666,
         parent_masstrack_id = 0L, peak_area = 1000,
         cSelectivity = 0.95, goodness_fitting = 0.9, snr = 5,
         apex = 1L),
    list(id_number = "F1", mz = 101.987654, rtime = 10.111,
         rtime_left_base = 9.111, rtime_right_base = 11.111,
         parent_masstrack_id = 1L, peak_area = 100,
         cSelectivity = 0.5, goodness_fitting = 0.6, snr = 1,
         apex = 2L)
  )
  cmap$FeatureTable <- data.frame(
    id_number = c("F0", "F1"), mz = c(100.123456, 101.987654),
    rtime = c(5.678, 10.111), rtime_left_base = c(4.444, 9.111),
    rtime_right_base = c(6.666, 11.111), parent_masstrack_id = 0:1,
    peak_area = c(1000, 100), cSelectivity = c(0.95, 0.5),
    goodness_fitting = c(0.9, 0.6), snr = c(5, 1),
    sample_a = c(10, 0), sample_b = c(20, 0), check.names = FALSE
  )
  experiment$CMAP <- cmap
  invisible(experiment)
}

# full/preferred 表应正确计算 detection count、舍入并按严格阈值过滤。
test_that("feature table export reproduces filtering and sample dropping", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$parameters$drop_unaligned_samples <- TRUE

  expect_null(ext_Experiment_export_feature_tables(experiment))

  full_path <- file.path(
    experiment$parameters$outdir, "export", "full_Feature_table.tsv"
  )
  preferred_path <- file.path(
    experiment$parameters$outdir, "preferred_Feature_table.tsv"
  )
  full <- utils::read.delim(full_path, check.names = FALSE)
  preferred <- utils::read.delim(preferred_path, check.names = FALSE)
  expect_equal(nrow(full), 2L)
  expect_equal(nrow(preferred), 1L)
  expect_equal(full$detection_counts, c(1L, 0L))
  expect_equal(full$mz[[1L]], 100.1235)
  expect_equal(names(full)[[length(names(full))]], "sample_a")
  expect_equal(experiment$dropped_sample_names, "sample_b")
  expect_equal(experiment$number_of_samples, 1L)
})

# unique 表在注释模式下插入 empCpd、formula 和 ion relation 三列。
test_that("unique compound table uses selected representative features", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$parameters$anno <- TRUE
  experiment$selected_unique_features <- list(
    F0 = list(100001, "C6H12O6", "M+H")
  )

  ext_Experiment_export_feature_tables(experiment)
  path <- file.path(
    experiment$parameters$outdir, "export",
    "unique_compound__Feature_table.tsv"
  )
  table <- utils::read.delim(path, check.names = FALSE)

  expect_equal(nrow(table), 1L)
  expect_equal(table$empCpd, 100001)
  expect_equal(table$neutral_formula, "C6H12O6")
  expect_equal(table$ion_relation, "M+H")
})

# export_log 和 export_readme 写出项目统计与 Python 默认说明。
test_that("log and README exports contain project metadata", {
  skip_if_not_installed("jsonlite")
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$dropped_sample_names <- "sample_b"

  expect_null(ext_Experiment_export_log(experiment))
  expect_null(ext_Experiment_export_readme(experiment))
  project <- jsonlite::fromJSON(
    file.path(experiment$parameters$outdir, "project.json")
  )
  readme <- paste(readLines(
    file.path(experiment$parameters$outdir, "README.txt"), warn = FALSE
  ), collapse = "\n")

  expect_equal(project$number_of_samples, 2L)
  expect_equal(project$number_scans, 5L)
  expect_equal(project$dropped_samples, "sample_b")
  expect_match(readme, "preferred_Feature_table.tsv", fixed = TRUE)
})

# 非注释 LC 和 DIMS 的 export_all 应建立质量网格、特征表和日志文件。
test_that("export_all dispatches LC and DIMS outputs", {
  skip_if_not_installed("jsonlite")
  lc <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters(workflow = "LC")
  )
  attach_export_fixture(lc)
  ext_Experiment_export_all(lc, anno = FALSE)

  expect_true(file.exists(file.path(
    lc$parameters$outdir, "export", "_mass_grid_mapping.csv"
  )))
  expect_true(file.exists(file.path(lc$parameters$outdir, "README.txt")))

  dims <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters(workflow = "DIMS")
  )
  attach_export_fixture(dims)
  ext_Experiment_export_all(dims, anno = FALSE)
  expect_true(file.exists(file.path(
    dims$parameters$outdir, "export", "_mass_grid_mapping.csv"
  )))
  expect_false(file.exists(file.path(dims$parameters$outdir, "README.txt")))
})

# 数据库估计偏差超过 2 ppm 时，全部 feature m/z 应使用相同比例除法校正。
test_that("database mass calibration updates feature mz values", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$CMAP$good_reference_landmark_peaks <- list(
    list(ref_id_num = 0L), list(ref_id_num = 1L)
  )
  KCD <- new.env(parent = emptyenv())
  KCD$evaluate_mass_accuracy_ratio <- function(mzs, mode, mz_tolerance_ppm) {
    expect_equal(mzs, c(100, 101))
    expect_equal(mode, "pos")
    expect_equal(mz_tolerance_ppm, 10)
    5e-6
  }
  experiment$KCD <- KCD

  before <- vapply(experiment$CMAP$FeatureList, `[[`, numeric(1), "mz")
  ext_Experiment_db_mass_calibrate(experiment)
  after <- vapply(experiment$CMAP$FeatureList, `[[`, numeric(1), "mz")

  expect_equal(after, before / 1.000005)
  expect_true(all(vapply(experiment$CMAP$FeatureList, function(feature) {
    identical(feature$mz_corrected_by_division, 1.000005)
  }, logical(1))))
})

# 没有数据库匹配时不校正 m/z，且函数只报告跳过。
test_that("database mass calibration skips missing matches", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$CMAP$good_reference_landmark_peaks <- list()
  KCD <- new.env(parent = emptyenv())
  KCD$evaluate_mass_accuracy_ratio <- function(...) NULL
  experiment$KCD <- KCD

  expect_message(
    expect_null(ext_Experiment_db_mass_calibrate(experiment)),
    "Mass accuracy check is skipped"
  )
})

# 第三阶段注释把 FeatureList 中尚未出现的 peak 包装成孤儿 empCpd。
test_that("orphan features are appended as empirical compounds", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  compounds <- list(`100000` = list(
    interim_id = 100000L,
    neutral_formula = "C1",
    neutral_formula_mass = 12,
    MS1_pseudo_Spectra = list(experiment$CMAP$FeatureList[[1L]])
  ))

  result <- ext_Experiment_append_orphans_to_epmCpds(experiment, compounds)

  expect_length(result, 2L)
  expect_equal(result[["100001"]]$interim_id, 100001L)
  expect_equal(result[["100001"]]$MS1_pseudo_Spectra[[1L]]$id_number, "F1")
  expect_null(result[["100001"]]$ion_relation)
})

# 注释 TSV 应包含数据库首选名称、短报告、feature 与 empCpd 字段。
test_that("peak annotation export writes matched database details", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  KCD <- new.env(parent = emptyenv())
  KCD$mass_indexed_compounds <- list(db1 = list(
    compounds = list(list(name = "Glucose"))
  ))
  KCD$short_report_emp_cpd <- function(id) paste0("short-", id)
  compounds <- list(`7` = list(
    interim_id = 7L,
    neutral_formula = "C6H12O6",
    neutral_formula_mass = 180.063,
    list_matches = list(list("db1", "M+H", 2L)),
    MS1_pseudo_Spectra = list(experiment$CMAP$FeatureList[[1L]])
  ))

  ext_Experiment_export_peak_annotation(
    experiment, compounds, KCD, "Feature_annotation"
  )
  text <- paste(readLines(file.path(
    experiment$parameters$outdir, "Feature_annotation.tsv"
  )), collapse = "\n")

  expect_match(text, "Glucose", fixed = TRUE)
  expect_match(text, "short-db1", fixed = TRUE)
  expect_match(text, "C6H12O6", fixed = TRUE)
})

# 每个 empCpd 选择 peak_area、goodness_fitting 最大的峰；singleton 使用固定关系。
test_that("unique feature selection follows peak quality ordering", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  compounds <- list(
    `10` = list(
      interim_id = 10L, neutral_formula = "C1",
      MS1_pseudo_Spectra = list(
        list(id_number = "F0", peak_area = 10, goodness_fitting = 0.9),
        list(id_number = "F1", peak_area = 20, goodness_fitting = 0.8,
             ion_relation = "M+H")
      )
    ),
    `11` = list(
      interim_id = 11L, neutral_formula = "C2",
      MS1_pseudo_Spectra = list(list(id_number = "F2"))
    )
  )

  ext_Experiment_select_unique_compound_features(experiment, compounds)

  expect_equal(experiment$selected_unique_features$F1,
               list(10L, "C1", "M+H"))
  expect_equal(experiment$selected_unique_features$F2,
               list(11L, "C2", "singleton"))
})

# QC 绘图器存在时收到 FeatureTable 和 export 路径；缺失时只提示跳过。
test_that("QC plot callback is optional and receives Python arguments", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  captured <- new.env(parent = emptyenv())
  experiment$asari_qc_plot <- function(table, outfile) {
    captured$table <- table
    captured$outfile <- outfile
  }

  ext_Experiment_generate_qc_plot_pdf(experiment, "custom.pdf")
  expect_equal(captured$table, experiment$CMAP$FeatureTable)
  expect_equal(captured$outfile, file.path(
    experiment$parameters$outdir, "export", "custom.pdf"
  ))
  rm("asari_qc_plot", envir = experiment)
  expect_message(
    expect_null(ext_Experiment_generate_qc_plot_pdf(experiment)),
    "cannot import matplotlib"
  )
})

# 外部 JMS 数据库通过明确 factory 接入；未提供时应报告缺少适配器。
test_that("annotation database loader uses configured JMS adapter", {
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  expect_error(
    ext_Experiment_load_annotation_db(experiment),
    "knownCompoundDatabase_factory"
  )
  marker <- new.env(parent = emptyenv())
  experiment$knownCompoundDatabase_factory <- function(src) {
    expect_equal(src, "custom")
    marker
  }
  expect_null(ext_Experiment_load_annotation_db(experiment, "custom"))
  expect_identical(experiment$KCD, marker)
})

# 完整 annotate 顺序使用 mock JMS 对象，并输出 JSON、pickle 和 annotation TSV。
test_that("annotation workflow connects KCD and EED adapters", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle export testing.")
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$CMAP$good_reference_landmark_peaks <- list()
  KCD <- new.env(parent = emptyenv())
  KCD$mass_indexed_compounds <- list()
  KCD$evaluate_mass_accuracy_ratio <- function(...) NULL
  experiment$knownCompoundDatabase_factory <- function(src) KCD

  calls <- character()
  experiment$ExperimentalEcpdDatabase_factory <- function(
      mode, mz_tolerance_ppm, rt_tolerance) {
    expect_equal(c(mode, mz_tolerance_ppm, rt_tolerance), c("pos", 5, 2))
    EED <- new.env(parent = emptyenv())
    EED$dict_empCpds <- list()
    EED$build_from_list_peaks <- function(peaks) calls <<- c(calls, "build")
    EED$extend_empCpd_annotation <- function(KCD) calls <<- c(calls, "extend")
    EED$annotate_singletons <- function(KCD) calls <<- c(calls, "singletons")
    EED
  }

  expect_message(
    expect_null(ext_Experiment_annotate(experiment)),
    "Mass accuracy check is skipped"
  )
  expect_equal(calls, c("build", "extend", "singletons"))
  expect_true(file.exists(file.path(
    experiment$parameters$outdir, "Annotated_empiricalCompounds.json"
  )))
  expect_true(file.exists(file.path(
    experiment$parameters$outdir, "export", "epd.pickle"
  )))
  expect_length(experiment$selected_unique_features, 2L)
})

# CMAP pickle 应能由 samples.R 的 Python pickle 读取器恢复核心字段。
test_that("CMAP export writes a readable Python pickle", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle export testing.")
  experiment <- ext_Experiment__init__(
    make_experiment_registry(), make_experiment_parameters()
  )
  attach_export_fixture(experiment)
  experiment$CMAP$`_number_of_samples_` <- 2L
  experiment$CMAP$rt_length <- 5L
  experiment$CMAP$good_reference_landmark_peaks <- list(list(apex = 2L))
  experiment$CMAP$dict_scan_rtime <- c(`0` = 0, `1` = 5)
  experiment$CMAP$composite_mass_tracks <- list(
    list(id_number = 0L, mz = 100, intensity = 1:5)
  )
  for (sample in experiment$all_samples) {
    sample$sample_id <- 0L
    sample$rt_landmarks <- list()
    sample$reverse_rt_cal_dict <- NULL
    sample$get_rt_calibration_records <- function() list(
      sample_id = 0L, name = sample$name, rt_landmarks = list(),
      reverse_rt_cal_dict = NULL
    )
  }

  ext_Experiment_export_CMAP_pickle(experiment)
  restored <- SimpleSample_load_intermediate(file.path(
    experiment$parameters$outdir, "export", "cmap.pickle"
  ))

  expect_equal(restored$`_number_of_samples_`, 2L)
  expect_equal(restored$rt_length, 5L)
  expect_equal(unlist(restored$rt_reference_landmarks), 2L)
  expect_length(restored$rt_records, 2L)
})
