# 对应 Python asari/tools/entropy_search.py：MS2谱图清理和谱图熵数据库检索。

# Python模块级默认参数。
entropy_search_param <- list(
  ion_mode = "positive",
  mz_tol_ms1 = 0.01,
  mz_tol_ms2 = 0.02,
  ms2_sim_tol = 0.6,
  precursor_mz_offset = 1.6
)
param <- entropy_search_param

# 调用可注入的熵检索后端方法；环境和具名list均支持。
.entropy_backend_method <- function(entropy_search, name) {
  method <- if (is.environment(entropy_search)) {
    get0(name, envir = entropy_search, mode = "function", inherits = FALSE)
  } else {
    entropy_search[[name]]
  }
  if (!is.function(method)) stop("Entropy backend is missing method: ", name)
  method
}

# 读取最佳数据库条目；后端可提供entries、get_entry或直接索引函数。
.entropy_backend_entry <- function(entropy_search, index0) {
  if (is.environment(entropy_search) && exists(
      "get_entry", envir = entropy_search, inherits = FALSE)) {
    return(entropy_search$get_entry(index0))
  }
  entries <- if (is.environment(entropy_search)) entropy_search$entries else entropy_search$entries
  if (!is.null(entries)) return(entries[[index0 + 1L]])
  if (is.function(entropy_search)) return(entropy_search(index0))
  stop("Entropy backend cannot return database entries.")
}

# 对应 clean_list_ms2spectra：清理普通MS2谱图列表。
clean_list_ms2spectra <- function(list_ms2_spectra, entropy_search, params = list()) {
  cleaner <- .entropy_backend_method(entropy_search, "clean_spectrum_for_search")
  cleaned <- list()
  for (spectrum in list_ms2_spectra) {
    peaks <- do.call(rbind, lapply(spectrum$peaks, as.numeric))
    cleaned_peaks <- cleaner(
      spectrum$precursor_mz,
      peaks,
      precursor_ions_removal_da = params$precursor_mz_offset
    )
    if (length(cleaned_peaks) > 0L && NROW(cleaned_peaks) > 0L) {
      cleaned[[length(cleaned) + 1L]] <- list(
        spectrum$precursor_mz, cleaned_peaks, spectrum$rtime, spectrum$id
      )
    }
  }
  cleaned
}

# 执行一次identity_search并规范化后端返回值。
.entropy_identity_search <- function(entropy_search, precursor_mz, peaks, params) {
  search <- .entropy_backend_method(entropy_search, "identity_search")
  search(
    precursor_mz = precursor_mz,
    peaks = peaks,
    ms1_tolerance_in_da = params$mz_tol_ms1,
    ms2_tolerance_in_da = params$mz_tol_ms2,
    output_matched_peak_number = TRUE
  )
}

# 去掉数据库条目的peaks字段，使结果可JSON序列化。
.entropy_result_entry <- function(entropy_search, best_index0) {
  entry <- .entropy_backend_entry(entropy_search, best_index0)
  entry$peaks <- NULL
  entry
}

# 对应 search_ms2_spectra：返回超过阈值的最佳数据库条目。
search_ms2_spectra <- function(cleaned_spectra, entropy_search, params = list()) {
  results <- list()
  for (spectrum in cleaned_spectra) {
    searched <- .entropy_identity_search(entropy_search, spectrum[[1L]], spectrum[[2L]], params)
    similarity <- as.numeric(searched[[1L]])
    matched_num <- as.integer(searched[[2L]])
    if (length(similarity) > 0L && max(similarity) > params$ms2_sim_tol) {
      best <- which.max(similarity)
      results[[length(results) + 1L]] <- list(
        .entropy_result_entry(entropy_search, best - 1L),
        spectrum[[4L]], spectrum[[1L]], spectrum[[3L]],
        as.numeric(similarity[[best]]), as.integer(matched_num[[best]])
      )
    }
  }
  results
}

# 对应 get_cleaned_experimental_ms2：清理已经按MS1特征分组的最佳MS2谱图。
get_cleaned_experimental_ms2 <- function(dict_ms2_spectra, entropy_search, params = list()) {
  cleaner <- .entropy_backend_method(entropy_search, "clean_spectrum_for_search")
  cleaned <- list()
  for (feature_id in names(dict_ms2_spectra)) {
    spectrum <- dict_ms2_spectra[[feature_id]][[1L]]
    peaks <- do.call(rbind, lapply(spectrum$peaks, as.numeric))
    cleaned_peaks <- cleaner(
      spectrum$precursor_mz,
      peaks,
      precursor_ions_removal_da = params$precursor_mz_offset
    )
    if (length(cleaned_peaks) > 0L && NROW(cleaned_peaks) > 0L) {
      cleaned[[length(cleaned) + 1L]] <- list(feature_id, spectrum$precursor_mz, cleaned_peaks)
    }
  }
  cleaned
}

# 对应 search_ms1matchedms2_spectra：以feature id为键保存最佳匹配。
search_ms1matchedms2_spectra <- function(cleaned_spectra, entropy_search, params = list()) {
  results <- list()
  for (spectrum in cleaned_spectra) {
    searched <- .entropy_identity_search(entropy_search, spectrum[[2L]], spectrum[[3L]], params)
    similarity <- as.numeric(searched[[1L]])
    matched_num <- as.integer(searched[[2L]])
    if (length(similarity) > 0L && max(similarity) > params$ms2_sim_tol) {
      best <- which.max(similarity)
      results[[as.character(spectrum[[1L]])]] <- list(
        .entropy_result_entry(entropy_search, best - 1L),
        as.numeric(similarity[[best]]),
        as.integer(matched_num[[best]])
      )
    }
  }
  results
}
