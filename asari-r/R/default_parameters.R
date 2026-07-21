# 对应 Python asari/default_parameters.py：默认参数和注释搜索规则。

# 质子与电子的精确质量常数；名称和值与 Python 原版保持一致。
PROTON <- 1.00727646677
electron <- 0.000549

# 正离子模式的初始加合物搜索模式，质量差均相对于 M+H+。
adduct_search_patterns_pos <- list(
  list(21.98194, "Na/H"),
  list(41.026549, "ACN"),
  list(67.987424, "NaCOOH"),
  list(37.955882, "K/H")
)

# Python 中 adduct_search_patterns 是正离子模式列表的同一对象。
adduct_search_patterns <- adduct_search_patterns_pos

# 负离子模式的初始加合物搜索模式。
adduct_search_patterns_neg <- list(
  list(21.98194, "Na/H"),
  list(67.987424, "NaCOOH"),
  list(82.0030, "C2HF3"),
  list(1.99566, "F <-> OH"),
  list(41.026549, "ACN"),
  list(37.955882, "K/H")
)

# 同位素搜索模式：质量差、标签和允许的相对强度范围。
isotope_search_patterns <- list(
  list(1.003355, "13C/12C", c(0, 0.8)),
  list(2.00671, "13C/12C*2", c(0, 0.8)),
  list(3.010065, "13C/12C*3", c(0, 0.8)),
  list(1.9970, "37Cl/35Cl", c(0.1, 0.8))
)

# 扩展加合物及中性丢失规则，用于经验化合物注释。
extended_adducts <- list(
  list(1.0078, "H"),
  list(-1.0078, "-H"),
  list(10.991, "Na/H, double charged"),
  list(0.5017, "13C/12C, double charged"),
  list(117.02655, "-NH3"),
  list(17.02655, "NH3"),
  list(-18.0106, "-H2O"),
  list(18.0106, "H2O"),
  list(18.033823, "NH4"),
  list(27.01089904, "HCN"),
  list(37.94694, "Ca/H2"),
  list(32.026215, "MeOH"),
  list(43.96389, "Na2/H2"),
  list(67.987424, "NaCOOH"),
  list(83.961361, "KCOOH"),
  list(97.96737927, "H2SO4"),
  list(97.97689507, "H3PO4")
)

# Python 原版写入输出目录的说明文本。
readme_doc_str <- paste0(
  "\nThe recommended feature table is `preferred_Feature_table.tsv`. \n\n",
  "All peaks are kept in `export/full_Feature_table.tsv` \n",
  "if they meet signal (snr) and shape standards \n",
  "(part of input parameters but default values are fine for most people). \n",
  "The filtering decisions are left to end users.\n\n",
  "Annotation is in JSON (`Annotated_empricalCompounds.json`) \n",
  "and in tab delimited text (`Feature_annotation.tsv`).\n\n",
  "The processing parameters and history are in `project.json`.\n\n",
  "Please refer to https://github.com/shuzhao-li-lab/asari for details, \n",
  "report bugs or request features.\n"
)

# 构造一份新的默认参数列表，避免调用者之间共享可变状态。
default_parameters <- function() {
  parameters <- list(
    project_name = "asari_project",
    outdir = "output",
    keep_intermediates = FALSE,
    reuse_intermediates = NULL,
    database_mode = "ondisk",
    multicores = 4L,
    mode = "pos",
    mass_range = c(50, 2000),
    mz_tolerance_ppm = 5,
    correction_tolerance_ppm = 1,
    min_timepoints = 6L,
    signal_noise_ratio = 2,
    min_intensity_threshold = 1000,
    min_peak_height = 100000,
    min_peak_ratio = 0.001,
    wlen = 25L,
    autoheight = FALSE,
    gaussian_shape = 0.5,
    peak_area = "sum",
    sample_metadata = NULL,
    workflow = "LC",
    min_min_peak_height = 10,
    dynamic_range = 1000,
    num_files_to_check = 20L,
    reference = NULL,
    rt_align_method = "lowess",
    rt_align_on = TRUE,
    debug_rtime_align = FALSE,
    drop_unaligned_samples = FALSE,
    rtime_tolerance = 50,
    cal_min_peak_height = 100000,
    peak_number_rt_calibration = 20L,
    max_retention_shift = NULL,
    num_lowess_iterations = 3L,
    project_sample_number_small = 10L,
    ms2_tolerance_in_ppm = 10,
    ms2_tolerance_in_da = 0.005,
    ri_tolerance = 30,
    score_cutoff_cosine = 0.4,
    score_cutoff_entropy = 0.4,
    corr_cutoff = 0.6,
    max_ri_delta = 100,
    do_mirror_plot = TRUE,
    max_core_features = 20000L,
    denovo = FALSE,
    anno = FALSE,
    check_isotope_ratio = FALSE,
    khipu_rtime_tolerance = 2,
    compress = FALSE,
    output_feature_table = "Feature_table.tsv",
    mass_grid_mapping = "_mass_grid_mapping.csv",
    annotation_filename = "Annotation_table.tsv",
    json_empricalCompounds = "_empCpd_json.json",
    table_for_viz = "preferred",
    vizualization_max_samples = 20L
  )

  # 兼容早期R原型曾使用的字段；真实提取仍以 min_intensity_threshold 为准。
  parameters$min_intensity <- parameters$min_intensity_threshold
  parameters
}

# 对应 Python 模块级 PARAMETERS；R调用者如需修改应先复制该列表。
PARAMETERS <- default_parameters()
