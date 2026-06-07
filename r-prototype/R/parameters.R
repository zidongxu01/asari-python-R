default_parameters <- function() {
  list(
    project_name = "asari_project",
    outdir = "output",
    keep_intermediates = FALSE,
    reuse_intermediates = NULL,
    database_mode = "ondisk",
    multicores = 4,

    mode = "pos",
    mass_range = c(50, 2000),
    mz_tolerance_ppm = 5,
    correction_tolerance_ppm = 1,

    min_timepoints = 6,
    signal_noise_ratio = 2,
    min_intensity_threshold = 1000,
    min_peak_height = 100000,
    min_peak_ratio = 0.001,
    wlen = 25,
    autoheight = FALSE,
    gaussian_shape = 0.5,
    peak_area = "sum",

    sample_metadata = NULL,
    workflow = "LC",

    min_min_peak_height = 10,
    dynamic_range = 1000,
    num_files_to_check = 20,

    reference = NULL,
    rt_align_method = "lowess",
    rt_align_on = TRUE,
    debug_rtime_align = FALSE,
    drop_unaligned_samples = FALSE,
    rtime_tolerance = 50,
    cal_min_peak_height = 100000,
    peak_number_rt_calibration = 20,
    max_retention_shift = NULL,
    num_lowess_iterations = 3,
    project_sample_number_small = 10,

    ms2_tolerance_in_ppm = 10,
    ms2_tolerance_in_da = 0.005,
    ri_tolerance = 30,
    score_cutoff_cosine = 0.4,
    score_cutoff_entropy = 0.4,
    corr_cutoff = 0.6,
    max_ri_delta = 100,
    do_mirror_plot = TRUE,
    max_core_features = 20000,
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
    vizualization_max_samples = 20,

    asari_r_prototype = TRUE
  )
}

asari_constants <- function() {
  list(
    proton = 1.00727646677,
    electron = 0.000549,

    adduct_search_patterns_pos = data.frame(
      mass_shift = c(21.98194, 41.026549, 67.987424, 37.955882),
      name = c("Na/H", "ACN", "NaCOOH", "K/H"),
      stringsAsFactors = FALSE
    ),

    adduct_search_patterns_neg = data.frame(
      mass_shift = c(21.98194, 67.987424, 82.0030, 1.99566, 41.026549, 37.955882),
      name = c("Na/H", "NaCOOH", "C2HF3", "F <-> OH", "ACN", "K/H"),
      stringsAsFactors = FALSE
    ),

    isotope_search_patterns = data.frame(
      mass_shift = c(1.003355, 2.00671, 3.010065, 1.9970),
      name = c("13C/12C", "13C/12C*2", "13C/12C*3", "37Cl/35Cl"),
      ratio_min = c(0, 0, 0, 0.1),
      ratio_max = c(0.8, 0.8, 0.8, 0.8),
      stringsAsFactors = FALSE
    ),

    extended_adducts = data.frame(
      mass_shift = c(
        1.0078, -1.0078, 10.991, 0.5017, 117.02655, 17.02655,
        -18.0106, 18.0106, 18.033823, 27.01089904, 37.94694,
        32.026215, 43.96389, 67.987424, 83.961361, 97.96737927,
        97.97689507
      ),
      name = c(
        "H", "-H", "Na/H, double charged", "13C/12C, double charged",
        "-NH3", "NH3", "-H2O", "H2O", "NH4", "HCN", "Ca/H2",
        "MeOH", "Na2/H2", "NaCOOH", "KCOOH", "H2SO4", "H3PO4"
      ),
      stringsAsFactors = FALSE
    )
  )
}

asari_readme_text <- function() {
  paste(
    "The recommended feature table is `preferred_Feature_table.tsv`.",
    "",
    "All peaks are kept in `export/full_Feature_table.tsv` if they meet",
    "signal (snr) and shape standards. The filtering decisions are left to",
    "end users.",
    "",
    "Annotation is in JSON (`Annotated_empricalCompounds.json`) and in tab",
    "delimited text (`Feature_annotation.tsv`).",
    "",
    "The processing parameters and history are in `project.json`.",
    "",
    "Please refer to https://github.com/shuzhao-li-lab/asari for details,",
    "report bugs or request features.",
    sep = "\n"
  )
}

