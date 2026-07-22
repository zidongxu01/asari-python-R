# 对应 Python asari/tools/feature_graph.py：把GC特征表转换成共洗脱关系图。

# 读取data.frame一行并转换成标量具名list。
.feature_graph_row <- function(data, row) {
  lapply(data[row, , drop = FALSE], function(value) value[[1L]])
}

# 创建Python FeatureGraph对应的可变R环境对象。
FeatureGraph <- function(ft_path, graph = NULL) FeatureGraph__init__(ft_path, graph)

# 对应 FeatureGraph.__init__。
FeatureGraph__init__ <- function(ft_path, graph = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("FeatureGraph", "environment")
  self$ft_path <- ft_path
  self$df <- if (!is.null(ft_path) && nzchar(ft_path)) {
    utils::read.delim(ft_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else NULL
  self$graph <- if (!is.null(graph)) graph else FeatureGraph_ft_to_graph(self)
  self$clusters <- NULL
  self$reverse_clusters <- NULL
  if (is.null(self$graph)) stop("Either a graph or a feature table path must be provided.")
  if (is.null(self$df)) stop("Either a graph or a feature table path must be provided.")

  # 安装与Python实例方法同名的闭包。
  self$ft_to_graph <- function() FeatureGraph_ft_to_graph(self)
  self$graph_to_ft <- function() FeatureGraph_graph_to_ft(self)
  self$filter_graph <- function(drt = 0.5) FeatureGraph_filter_graph(self, drt)
  self$find_spectral_clusters <- function() FeatureGraph_find_spectral_clusters(self)
  self$extract_fragmentation_spectrum <- function(find_clusters = FALSE, MIN_PEAKS_EXTRACTION = 3L) {
    FeatureGraph_extract_fragmentation_spectrum(self, find_clusters, MIN_PEAKS_EXTRACTION)
  }
  self$map_annotations <- function(matches, to_extract = c("compound_name", "inchikey", "formula")) {
    FeatureGraph_map_annotations(self, matches, to_extract)
  }
  self
}

# 对应静态方法 ftgraph_from_ft。
FeatureGraph_ftgraph_from_ft <- function(ft_path) FeatureGraph(ft_path)

# 对应静态方法 ftgraph_from_graph。
FeatureGraph_ftgraph_from_graph <- function(ft_path, graph) FeatureGraph(ft_path, graph)

# 对应静态方法 metric。
FeatureGraph_metric <- function(x, y) {
  c(abs(as.numeric(x$mz) - as.numeric(y$mz)), abs(as.numeric(x$rtime) - as.numeric(y$rtime)))
}

# 对应 ft_to_graph：建立包含全部特征对的无向完全图。
FeatureGraph_ft_to_graph <- function(self) {
  if (is.null(self$df)) return(NULL)
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  node_ids <- vapply(records, function(record) as.character(record$id_number), "")
  nodes <- stats::setNames(records, node_ids)
  edges <- data.frame(u = character(), v = character(), dmz = numeric(), drt = numeric())
  if (length(records) > 1L) {
    pairs <- utils::combn(seq_along(records), 2L)
    edges <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(ii) {
      left <- pairs[1L, ii]
      right <- pairs[2L, ii]
      metric <- FeatureGraph_metric(records[[left]], records[[right]])
      data.frame(
        u = node_ids[[left]], v = node_ids[[right]],
        dmz = metric[[1L]], drt = metric[[2L]], stringsAsFactors = FALSE
      )
    }))
  }
  list(nodes = nodes, edges = edges)
}

# 对应 graph_to_ft：把图节点属性恢复成data.frame。
FeatureGraph_graph_to_ft <- function(self) {
  records <- unname(self$graph$nodes)
  if (length(records) == 0L || all(lengths(records) == 0L)) return(data.frame())
  keys <- unique(unlist(lapply(records, names), use.names = FALSE))
  rows <- lapply(records, function(record) {
    stats::setNames(lapply(keys, function(key) record[[key]] %||% NA), keys)
  })
  as.data.frame(do.call(rbind, lapply(rows, as.data.frame)), stringsAsFactors = FALSE)
}

# 本文件私有的空值替代运算符。
`%||%` <- function(value, fallback) if (is.null(value)) fallback else value

# 对应 filter_graph：仅保留drt严格小于阈值的边及其端点。
FeatureGraph_filter_graph <- function(self, drt = 0.5) {
  selected <- self$graph$edges[self$graph$edges$drt < drt, , drop = FALSE]
  node_ids <- unique(c(selected$u, selected$v))
  # NetworkX add_edges_from只带(u,v)，因此过滤图不保留原边属性。
  filtered_edges <- if (nrow(selected)) {
    data.frame(u = selected$u, v = selected$v, dmz = NA_real_, drt = NA_real_)
  } else data.frame(u = character(), v = character(), dmz = numeric(), drt = numeric())
  filtered_nodes <- stats::setNames(rep(list(list()), length(node_ids)), node_ids)
  FeatureGraph_ftgraph_from_graph(
    self$ft_path,
    list(nodes = filtered_nodes, edges = filtered_edges)
  )
}

# 计算无向图连接分量，返回节点id字符向量列表。
.feature_graph_components <- function(graph) {
  nodes <- names(graph$nodes)
  components <- list()
  unseen <- nodes
  while (length(unseen) > 0L) {
    queue <- unseen[[1L]]
    component <- character()
    while (length(queue) > 0L) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      if (current %in% component) next
      component <- c(component, current)
      edges <- graph$edges
      neighbors <- unique(c(edges$v[edges$u == current], edges$u[edges$v == current]))
      queue <- unique(c(queue, setdiff(neighbors, component)))
    }
    components[[length(components) + 1L]] <- component
    unseen <- setdiff(unseen, component)
  }
  components
}

# 对应 find_spectral_clusters：连接分量编号从0开始。
FeatureGraph_find_spectral_clusters <- function(self) {
  if (is.null(self$graph)) stop("Must have a graph to find spectral clusters.")
  components <- .feature_graph_components(self$graph)
  cluster_to_feature <- list()
  feature_to_cluster <- list()
  for (ii in seq_along(components)) {
    cluster_id <- ii - 1L
    cluster_to_feature[[as.character(cluster_id)]] <- components[[ii]]
    for (feature in components[[ii]]) {
      self$graph$nodes[[feature]]$clique_id <- cluster_id
      feature_to_cluster[[feature]] <- cluster_id
    }
  }
  cat(sprintf("Clustering finds: %d spectral clusters!\n", length(cluster_to_feature)))
  self$clusters <- cluster_to_feature
  self$reverse_clusters <- feature_to_cluster
  invisible(NULL)
}

# 创建轻量Spectrum兼容对象。
.feature_graph_spectrum <- function(mzs, intensities, metadata) {
  list(
    mz = as.numeric(mzs), intensities = as.numeric(intensities),
    metadata = metadata,
    peaks = cbind(mz = as.numeric(mzs), intensity = as.numeric(intensities))
  )
}

# 默认谱图过滤/归一化；用户可通过options注入matchms等价后端。
.feature_graph_process_spectrum <- function(spectrum) {
  filter <- getOption("asariR.spectrum_default_filters")
  normalize <- getOption("asariR.spectrum_normalize")
  if (is.function(filter)) spectrum <- filter(spectrum)
  if (is.function(normalize)) return(normalize(spectrum))
  maximum <- max(spectrum$intensities)
  if (is.finite(maximum) && maximum > 0) {
    spectrum$intensities <- spectrum$intensities / maximum
    spectrum$peaks[, 2L] <- spectrum$intensities
  }
  spectrum
}

# 对应 extract_fragmentation_spectrum：每个cluster选择总强度最高且峰数足够的样本。
FeatureGraph_extract_fragmentation_spectrum <- function(
    self, find_clusters = FALSE, MIN_PEAKS_EXTRACTION = 3L) {
  if (isTRUE(find_clusters)) self$find_spectral_clusters()
  else if (is.null(self$clusters)) stop("Must find clusters before extracting fragmentation spectra.")
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  feature_table <- stats::setNames(records, vapply(records, function(x) as.character(x$id_number), ""))
  samples <- if (ncol(self$df) > 11L) names(self$df)[12:ncol(self$df)] else character()
  processed <- list()
  for (cluster_id in names(self$clusters)) {
    feature_ids <- self$clusters[[cluster_id]]
    spectra <- list()
    for (sample in samples) {
      values <- lapply(feature_ids, function(feature) {
        c(as.numeric(feature_table[[feature]]$mz), as.numeric(feature_table[[feature]][[sample]]))
      })
      values <- values[order(vapply(values, `[[`, 0, 1L), vapply(values, `[[`, 0, 2L))]
      mzs <- vapply(values, `[[`, 0, 1L)
      intensities <- vapply(values, `[[`, 0, 2L)
      if (sum(intensities) > 0) {
        spectra[[length(spectra) + 1L]] <- .feature_graph_spectrum(
          mzs, intensities, list(cluster_id = as.integer(cluster_id), sample = sample)
        )
      }
    }
    if (length(spectra) == 1L) {
      # Python原版只赋值但没有append，保留这一实际行为。
      selected_spectrum <- spectra[[1L]]
      invisible(selected_spectrum)
    } else if (length(spectra) > 1L) {
      totals <- vapply(spectra, function(spectrum) sum(spectrum$intensities), 0)
      for (index in order(-totals, seq_along(totals))) {
        selected <- .feature_graph_process_spectrum(spectra[[index]])
        if (nrow(selected$peaks) >= MIN_PEAKS_EXTRACTION) {
          processed[[length(processed) + 1L]] <- selected
          break
        }
      }
    }
  }
  processed
}

# 对应 map_annotations：把cluster匹配传播到其中每个feature并输出注释表。
FeatureGraph_map_annotations <- function(
    self, matches, to_extract = c("compound_name", "inchikey", "formula")) {
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  feature_table <- stats::setNames(records, vapply(records, function(x) as.character(x$id_number), ""))
  for (match in matches) {
    cluster_id <- as.character(match$extract$metadata$cluster_id)
    annotation <- unlist(lapply(to_extract, function(field) {
      c(match$library$metadata[[field]] %||% "Not Found", match$similarity, match$match_peaks)
    }), use.names = FALSE)
    for (feature in self$clusters[[cluster_id]]) {
      existing <- feature_table[[feature]]$annotations %||% list()
      # Python缩进导致同一最终annotation每个字段追加一次。
      feature_table[[feature]]$annotations <- c(
        existing, rep(list(annotation), length(to_extract))
      )
    }
  }
  # 原始列保持data.frame类型，注释单独作为list-column，避免rbind将嵌list错转为字符。
  ids <- as.character(self$df$id_number)
  self$df$annotations <- I(lapply(ids, function(id) {
    feature_table[[id]]$annotations %||% list()
  }))
  output <- self$df
  if ("annotations" %in% names(output)) {
    output$annotations <- vapply(output$annotations, function(value) paste(capture.output(dput(value)), collapse = ""), "")
  }
  outfile <- gsub(".tsv", "_annotated_gc_beta.tsv", self$ft_path, fixed = TRUE)
  utils::write.table(output, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(NULL)
}
