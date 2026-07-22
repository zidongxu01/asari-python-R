# 对应 Python asari/tools/merge.py：合并同一研究中的多个asari特征表。

# Python模块文档中保留的旧接口；依赖调用方提供全局容差。
match2 <- function(F1, F2) {
  ppm <- get0("PPM_tolerance", ifnotfound = 4, inherits = TRUE)
  rt <- get0("RTime_tolerance", ifnotfound = 15, inherits = TRUE)
  abs(F1[[1L]] - F2[[1L]]) / F1[[1L]] < ppm &&
    abs(F1[[2L]] - F2[[2L]]) < rt
}

# Python模块文档中的旧RT分箱辅助函数。
`__bin_by_median_rt__` <- function(List_of_peaks, tolerance) {
  tuples <- lapply(List_of_peaks, function(peak) {
    value <- if (is.environment(peak)) peak$cal_rtime else peak$cal_rtime
    list(value, peak)
  })
  tuples <- tuples[order(vapply(tuples, `[[`, 0, 1L))]
  tools_merge_bin_by_median(tuples, function(value) max(tolerance, 0.1 * value))
}

# Python模块文档中的旧选择性统计接口。
`_formula_selectivity_` <- function(L) {
  result <- as.list(as.integer(table(L)))
  names(result) <- names(table(L))
  result
}

# 对应可执行的 bin_by_median：动态中位数决定是否开始新分箱。
tools_merge_bin_by_median <- function(List_of_tuples, func_tolerance) {
  if (length(List_of_tuples) == 0L) stop("List_of_tuples must be non-empty.")
  bins <- list(list(List_of_tuples[[1L]]))
  if (length(List_of_tuples) > 1L) {
    for (tuple in List_of_tuples[-1L]) {
      current <- bins[[length(bins)]]
      median_value <- stats::median(vapply(current, `[[`, 0, 1L))
      if (tuple[[1L]] - median_value < func_tolerance(tuple[[1L]])) {
        bins[[length(bins)]][[length(current) + 1L]] <- tuple
      } else {
        bins[[length(bins) + 1L]] <- list(tuple)
      }
    }
  }
  lapply(bins, function(bin) lapply(bin, `[[`, 2L))
}


# 保留Python对外公共名，并让模块归属别名始终可审计。
bin_by_median <- tools_merge_bin_by_median

# 对应 _read_asari_ftables_：读取特征键、质量和原始文本行。
`_read_asari_ftables_` <- function(infile, make_feature_id = FALSE) {
  lines <- readLines(infile, warn = FALSE)
  if (length(lines) == 0L) stop("Feature table is empty.")
  key_dict <- list()
  data_dict <- list(header = lines[[1L]])
  for (line in lines[-1L]) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (isTRUE(make_feature_id)) {
      formula_mass <- fields[[1L]]
      mz <- fields[[2L]]
      rtime <- fields[[3L]]
      feature_id <- paste0(mz, "@", rtime)
    } else {
      feature_id <- fields[[1L]]
      formula_mass <- fields[[2L]]
      mz <- fields[[3L]]
      rtime <- fields[[4L]]
    }
    key_dict[[feature_id]] <- list(formula_mass, as.numeric(mz), as.numeric(rtime))
    data_dict[[feature_id]] <- sub("[\r\n]+$", "", line)
  }
  list(key_dict, data_dict)
}

# 内部函数，对应Python嵌套def find_min_delta。
.merge_find_min_delta <- function(values) {
  values <- sort(as.numeric(values))
  differences <- diff(c(0, values))
  min(differences[-1L])
}

# 对应 _masstrace2features_：在同一质量轨迹中按RT合并不同表的特征。
`_masstrace2features_` <- function(LL, RTime_tolerance) {
  number_tables <- length(LL)
  rt_tuples <- list()
  for (table_index in seq_along(LL)) {
    for (feature_id in LL[[table_index]]) {
      rt <- as.numeric(strsplit(feature_id, "@", fixed = TRUE)[[1L]][[2L]])
      rt_tuples[[length(rt_tuples) + 1L]] <- list(rt, list(feature_id, table_index))
    }
  }
  multiplet_gaps <- vapply(LL[lengths(LL) > 1L], function(ids) {
    .merge_find_min_delta(vapply(ids, function(id) {
      as.numeric(strsplit(id, "@", fixed = TRUE)[[1L]][[2L]])
    }, 0))
  }, 0)
  if (length(multiplet_gaps) > 0L) {
    RTime_tolerance <- min(RTime_tolerance, min(multiplet_gaps) - 1)
  }
  rt_tuples <- rt_tuples[order(vapply(rt_tuples, `[[`, 0, 1L))]
  grouped_features <- bin_by_median(rt_tuples, function(value) RTime_tolerance)
  lapply(grouped_features, function(group) {
    row <- rep(list(NULL), number_tables)
    for (item in group) row[[item[[2L]]]] <- item[[1L]]
    row
  })
}

# 对应 merge_feature_tables：统一质量轨迹、按RT对齐并写出合并表。
merge_feature_tables <- function(
    list_of_feature_tables,
    make_feature_id = FALSE,
    PPM_tolerance = 4,
    RTime_tolerance = 15,
    outfile = "merged_featuretable.tsv") {
  number_tables <- length(list_of_feature_tables)
  formula_mass_dict <- list()
  mass_pointer <- list()
  all_keys <- list()
  all_tables <- list()
  unassigned <- list()

  for (file in list_of_feature_tables) {
    parsed <- `_read_asari_ftables_`(file, make_feature_id)
    key_dict <- parsed[[1L]]
    data_dict <- parsed[[2L]]
    all_keys[[length(all_keys) + 1L]] <- key_dict
    all_tables[[length(all_tables) + 1L]] <- data_dict
    for (key in names(key_dict)) {
      value <- key_dict[[key]]
      if (startsWith(value[[1L]], "_M_")) {
        unassigned[[length(unassigned) + 1L]] <- list(value[[2L]], value[[1L]])
      } else {
        formula_mass_dict[[value[[1L]]]] <- rep(list(character()), number_tables)
      }
    }
  }

  if (length(unassigned) > 0L) {
    unassigned <- unassigned[order(vapply(unassigned, `[[`, 0, 1L))]
    mass_bins <- bin_by_median(
      unassigned,
      function(value) PPM_tolerance * 1e-6 * value
    )
    for (bin in mass_bins) {
      representative <- bin[[1L]]
      for (mass_id in bin) mass_pointer[[mass_id]] <- representative
    }
    for (representative in unique(unlist(mass_pointer, use.names = FALSE))) {
      formula_mass_dict[[representative]] <- rep(list(character()), number_tables)
    }
  }

  mass_to_features <- lapply(all_keys, function(key_dict) {
    result <- list()
    for (feature_id in names(key_dict)) {
      mass_id <- key_dict[[feature_id]][[1L]]
      result[[mass_id]] <- c(result[[mass_id]], feature_id)
    }
    result
  })
  for (table_index in seq_len(number_tables)) {
    for (mass_id in names(mass_to_features[[table_index]])) {
      unified <- if (startsWith(mass_id, "_M_")) mass_pointer[[mass_id]] else mass_id
      formula_mass_dict[[unified]][[table_index]] <- mass_to_features[[table_index]][[mass_id]]
    }
  }

  new_features <- lapply(names(formula_mass_dict), function(mass_id) {
    list(mass_id, `_masstrace2features_`(formula_mass_dict[[mass_id]], RTime_tolerance))
  })
  column_counts <- vapply(all_tables, function(table) {
    length(strsplit(table$header, "\t", fixed = TRUE)[[1L]])
  }, 0L)
  header <- paste0(
    "new_feature_id\tnumber_src_tables\told_ids\tformula_mass\t",
    paste(vapply(all_tables, `[[`, "", "header"), collapse = "\t")
  )
  output_lines <- header
  for (mass_group in new_features) {
    mass_id <- mass_group[[1L]]
    for (row in mass_group[[2L]]) {
      present <- !vapply(row, is.null, FALSE)
      old_ids <- vapply(row, function(id) if (is.null(id)) "_" else id, "")
      line <- paste(c(old_ids[present][[1L]], sum(present), paste(old_ids, collapse = ","), mass_id), collapse = "\t")
      for (table_index in seq_len(number_tables)) {
        line <- if (present[[table_index]]) {
          paste(line, all_tables[[table_index]][[row[[table_index]]]], sep = "\t")
        } else {
          paste0(line, paste(rep("\t", column_counts[[table_index]]), collapse = ""))
        }
      }
      output_lines <- c(output_lines, line)
    }
  }
  writeLines(output_lines, outfile, useBytes = TRUE)
  cat("Merged result was written to ", outfile, "\n")
  invisible(NULL)
}
