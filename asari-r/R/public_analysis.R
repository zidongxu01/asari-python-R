# 面向普通用户的合并、GC特征图和绘图入口。

# 从标准asari特征表中识别样本强度列。
.asari_sample_columns <- function(table) {
  metadata <- c(
    "id_number", "mz", "rtime", "rtime_left_base", "rtime_right_base",
    "parent_masstrack_id", "peak_area", "cSelectivity", "goodness_fitting",
    "snr", "detection_counts"
  )
  setdiff(names(table), metadata)
}

# 为一个新特征在当前共识组中找到最近且可用的组。
.asari_merge_candidate <- function(clusters, mz, rtime, table_index, ppm, rt_tolerance) {
  if (length(clusters) == 0L) return(NA_integer_)
  centers_mz <- vapply(clusters, function(cluster) stats::median(cluster$mz), 0)
  delta_da <- mz * ppm * 1e-6
  ordering <- order(centers_mz)
  sorted <- centers_mz[ordering]
  left <- findInterval(mz - delta_da, sorted) + 1L
  right <- findInterval(mz + delta_da, sorted)
  if (left > right || left > length(sorted) || right < 1L) return(NA_integer_)
  candidates <- ordering[seq.int(max(1L, left), min(length(sorted), right))]
  candidates <- candidates[vapply(candidates, function(index) {
    cluster <- clusters[[index]]
    !(table_index %in% cluster$table_index) &&
      abs(rtime - stats::median(cluster$rtime)) <= rt_tolerance
  }, FALSE)]
  if (length(candidates) == 0L) return(NA_integer_)
  distance <- vapply(candidates, function(index) {
    cluster <- clusters[[index]]
    mass_part <- abs(mz - stats::median(cluster$mz)) / max(delta_da, .Machine$double.eps)
    rt_part <- abs(rtime - stats::median(cluster$rtime)) / rt_tolerance
    mass_part + rt_part
  }, 0)
  candidates[[which.min(distance)]]
}

#' 合并多张asari特征表
#'
#' `method = "standard"`直接合并普通asari表；`method = "legacy"`
#' 调用Python `tools/merge.py`对应的旧式formula-mass表合并器。
#'
#' @param feature_tables 两个或更多明硤TSV路径。
#' @param output 明硤的合并TSV输出路径。
#' @param ppm m/z容差。
#' @param rt_tolerance RT容差，单位为秒。
#' @param method 标准asari表或旧式formula-mass表。
#' @param make_feature_id 仅传给旧式合并器。
#' @return `standard`返回合并data.frame；`legacy`返回输出路径。
#' @export
asari_merge_feature_tables <- function(
    feature_tables,
    output,
    ppm = 4,
    rt_tolerance = 15,
    method = c("standard", "legacy"),
    make_feature_id = FALSE) {
  method <- match.arg(method)
  ppm <- .asari_public_positive_number(ppm, "ppm")
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  if (!is.character(feature_tables) || length(feature_tables) < 2L) {
    stop("feature_tables must contain at least two explicit TSV paths.", call. = FALSE)
  }
  paths <- vapply(feature_tables, function(path) {
    path <- path.expand(path)
    if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
    normalizePath(path, mustWork = TRUE)
  }, "")
  outfile <- .asari_public_output_file(output)
  if (method == "legacy") {
    merge_feature_tables(paths, make_feature_id, ppm, rt_tolerance, outfile)
    return(normalizePath(outfile, mustWork = TRUE))
  }

  tables <- lapply(paths, .asari_public_feature_table)
  invisible(lapply(seq_along(tables), function(ii) {
    .asari_public_require_columns(tables[[ii]], c("mz", "rtime"), paste0("feature_tables[[", ii, "]]"))
  }))
  for (ii in seq_along(tables)) {
    if (!("id_number" %in% names(tables[[ii]]))) {
      tables[[ii]]$id_number <- paste0("row", seq_len(nrow(tables[[ii]])))
    }
  }

  # 按输入表顺序将特征加入共识组，每张表在一组中最多占一行。
  clusters <- list()
  for (table_index in seq_along(tables)) {
    table <- tables[[table_index]]
    row_order <- order(as.numeric(table$mz), as.numeric(table$rtime))
    for (row_index in row_order) {
      mz <- as.numeric(table$mz[[row_index]])
      rtime <- as.numeric(table$rtime[[row_index]])
      if (!is.finite(mz) || !is.finite(rtime)) next
      candidate <- .asari_merge_candidate(
        clusters, mz, rtime, table_index, ppm, rt_tolerance
      )
      member <- list(table_index = table_index, row_index = row_index)
      if (is.na(candidate)) {
        clusters[[length(clusters) + 1L]] <- list(
          mz = mz, rtime = rtime, table_index = table_index, members = list(member)
        )
      } else {
        clusters[[candidate]]$mz <- c(clusters[[candidate]]$mz, mz)
        clusters[[candidate]]$rtime <- c(clusters[[candidate]]$rtime, rtime)
        clusters[[candidate]]$table_index <- c(clusters[[candidate]]$table_index, table_index)
        clusters[[candidate]]$members <- c(clusters[[candidate]]$members, list(member))
      }
    }
  }
  clusters <- clusters[order(
    vapply(clusters, function(cluster) stats::median(cluster$mz), 0),
    vapply(clusters, function(cluster) stats::median(cluster$rtime), 0)
  )]

  # 所有样本列加上来源表前缀，避免不同项目使用相同样本名时覆盖。
  table_labels <- make.unique(tools::file_path_sans_ext(basename(paths)))
  sample_names <- unlist(lapply(seq_along(tables), function(ii) {
    paste0(table_labels[[ii]], "::", .asari_sample_columns(tables[[ii]]))
  }), use.names = FALSE)
  rows <- lapply(seq_along(clusters), function(ii) {
    cluster <- clusters[[ii]]
    intensities <- stats::setNames(as.list(rep(NA_real_, length(sample_names))), sample_names)
    source_ids <- character()
    for (member in cluster$members) {
      table_index <- member$table_index
      row_index <- member$row_index
      table <- tables[[table_index]]
      source_ids <- c(source_ids, paste0(
        table_labels[[table_index]], ":", as.character(table$id_number[[row_index]])
      ))
      for (sample in .asari_sample_columns(table)) {
        key <- paste0(table_labels[[table_index]], "::", sample)
        intensities[[key]] <- suppressWarnings(as.numeric(table[[sample]][[row_index]]))
      }
    }
    c(list(
      consensus_id = paste0("C", ii - 1L),
      mz = stats::median(cluster$mz), rtime = stats::median(cluster$rtime),
      detection_tables = length(unique(cluster$table_index)),
      source_ids = paste(source_ids, collapse = ";")
    ), intensities)
  })
  result <- if (length(rows)) {
    as.data.frame(do.call(rbind, lapply(rows, function(row) {
      as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
    })), check.names = FALSE)
  } else data.frame(
    consensus_id = character(), mz = numeric(), rtime = numeric(),
    detection_tables = integer(), source_ids = character()
  )
  # rbind会将混合list转为文本，这里恢复共识数值列。
  for (name in c("mz", "rtime", "detection_tables", sample_names)) {
    if (name %in% names(result)) result[[name]] <- suppressWarnings(as.numeric(result[[name]]))
  }
  utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  attr(result, "output_file") <- outfile
  result
}

#' 按共洗脱关系构建GC特征图
#'
#' @param feature_table 明硤的asari特征表TSV。
#' @param rt_tolerance 构建边的RT容差。
#' @param output 可选的特征-簇TSV输出路径。
#' @return 每行一个特征和其谱簇ID的data.frame。
#' @export
asari_feature_graph <- function(feature_table, rt_tolerance = 0.5, output = NULL) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  graph <- FeatureGraph_ftgraph_from_ft(normalizePath(path, mustWork = TRUE))$filter_graph(rt_tolerance)
  graph$find_spectral_clusters()
  rows <- unlist(lapply(names(graph$clusters), function(cluster_id) {
    lapply(graph$clusters[[cluster_id]], function(feature_id) {
      data.frame(feature_id = feature_id, cluster_id = as.integer(cluster_id))
    })
  }), recursive = FALSE)
  result <- if (length(rows)) do.call(rbind, rows) else {
    data.frame(feature_id = character(), cluster_id = integer())
  }
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
    attr(result, "output_file") <- outfile
  }
  result
}

#' 绘制两张谱图的镜像图
#'
#' @param query,reference 两列`m/z, intensity` matrix或data.frame。
#' @param output 明硤的PDF输出路径。
#' @param tolerance 可选的碎片m/z匹配容差Da。
#' @param title 图标题。
#' @return PDF绝对路径。
#' @export
asari_mirror_plot <- function(query, reference, output, tolerance = NULL, title = "MS Mirror Plot") {
  outfile <- .asari_public_output_file(output)
  mirror_plot(
    as.matrix(query), as.matrix(reference), match_tol = tolerance,
    title = title, outfile = outfile
  )
  normalizePath(outfile, mustWork = TRUE)
}

#' 绘制mzML指定扫描号和m/z区域
#'
#' @param input 一个明硤mzML文件。
#' @param output 明硤PDF输出路径。
#' @param scan_range 长度为2的扫描号开区间。
#' @param mz_range 长度为2的m/z开区间。
#' @param ms_level 质谱级别。
#' @param style `"scatter"`、`"double"`或`"line"`。
#' @return PDF绝对路径。
#' @export
asari_plot_mz_region <- function(
    input,
    output,
    scan_range,
    mz_range,
    ms_level = 1L,
    style = c("scatter", "double", "line")) {
  infile <- .asari_public_one_mzml(input)
  if (!is.numeric(scan_range) || length(scan_range) != 2L || diff(scan_range) <= 0) {
    stop("scan_range must be two increasing numbers.")
  }
  if (!is.numeric(mz_range) || length(mz_range) != 2L || diff(mz_range) <= 0) {
    stop("mz_range must be two increasing numbers.")
  }
  style <- match.arg(style)
  points <- get_plot_region_from_file(
    infile, scan_range[[1L]], scan_range[[2L]], mz_range[[1L]], mz_range[[2L]],
    as.integer(ms_level)
  )
  if (length(points) == 0L) stop("No points were found in the requested region.")
  outfile <- .asari_public_output_file(output)
  grDevices::pdf(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (style == "scatter") plot_scatter_map_region(points)
  else if (style == "double") double_scatter_map_region(points)
  else with_line_scatter_map_region(points)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  normalizePath(outfile, mustWork = TRUE)
}

#' 对特征表样本强度执行PCA并绘图
#'
#' @param feature_table 明硤的asari特征表TSV。
#' @param output 可选PDF输出路径；`NULL`使用当前绘图设备。
#' @return 样本的PC1和PC2坐标matrix。
#' @export
asari_pca <- function(feature_table, output = NULL) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  opened <- !is.null(output)
  if (opened) {
    outfile <- .asari_public_output_file(output)
    grDevices::pdf(outfile)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  scores <- pca_ftable(path)
  if (opened) {
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
    attr(scores, "output_file") <- normalizePath(outfile, mustWork = TRUE)
  }
  scores
}

#' 绘制样本汇总强度与peak_area的相关图
#'
#' @param feature_table 明硤的asari特征表TSV。
#' @param output 明硤PDF输出路径。
#' @return PDF绝对路径。
#' @export
asari_plot_correlations <- function(feature_table, output) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  outfile <- .asari_public_output_file(output)
  grDevices::pdf(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_correlations(path)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  normalizePath(outfile, mustWork = TRUE)
}
