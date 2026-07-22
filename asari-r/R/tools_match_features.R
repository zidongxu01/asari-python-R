# 对应 Python asari/tools/match_features.py：按m/z和保留时间匹配LC-MS特征。

# 在R中直接执行ppm窗口搜索，返回Python centurion tree检索得到的记录顺序。
.match_features_mz_candidates <- function(query_mz, features, mz_ppm) {
  Filter(function(feature) {
    abs(feature$mz - query_mz) / query_mz * 1e6 <= mz_ppm
  }, features)
}

# 对应 get_featureList：读取表格并使用Python行号生成row id。
get_featureList <- function(infile, start_row = 1L, mz_col = 0L, rt_col = 1L, sep = "\t") {
  lines <- readLines(infile, warn = FALSE)
  # start_row和列号均沿用Python的0基语义。
  selected <- if (length(lines) > start_row) lines[seq.int(start_row + 1L, length(lines))] else character()
  lapply(seq_along(selected), function(ii) {
    fields <- strsplit(selected[[ii]], sep, fixed = TRUE)[[1L]]
    list(
      id = paste0("row", ii + start_row),
      mz = as.numeric(fields[[as.integer(mz_col) + 1L]]),
      rtime = as.numeric(fields[[as.integer(rt_col) + 1L]])
    )
  })
}

# 对应 list_match_lcms_features：保存所有同时满足ppm和RT窗口的匹配。
list_match_lcms_features <- function(list1, list2, mz_ppm = 5, rt_tolerance = 5) {
  mapped <- list()
  for (p1 in list1) {
    matched <- character()
    candidates <- .match_features_mz_candidates(p1$mz, list2, mz_ppm)
    for (p2 in candidates) {
      if (abs(p1$rtime - p2$rtime) < rt_tolerance) matched <- c(matched, p2$id)
    }
    if (length(matched) > 0L) mapped[[as.character(p1$id)]] <- matched
  }
  cat(sprintf(
    "Of %d list1 features, number of uni-direction matched features is %d.\n",
    length(list1), length(mapped)
  ))
  mapped
}

# 对应 bidirectional_match：分别返回两个方向的一对多映射。
bidirectional_match <- function(list1, list2, mz_ppm = 5, rt_tolerance = 5) {
  dict1 <- list_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
  dict2 <- list_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)
  cat("    ~~~ match_numbers ~~~     \n\n")
  unique1 <- names(dict1)[lengths(dict1) == 1L]
  unique2 <- names(dict2)[lengths(dict2) == 1L]
  cat("Unique Number of matched features in table 1: ", length(unique1), "\n")
  cat("Unique Number of matched features in table 2: ", length(unique2), "\n")
  bi_unique <- Filter(function(id) dict1[[id]][[1L]] %in% unique2, unique1)
  cat("Biodirectional, unique Number of matched feature pairs: ", length(bi_unique), "\n")
  list(dict1, dict2)
}

# 选择m/z偏差最小的候选；偏差相同按id排序以复现Python tuple排序。
best_mz_match_lcms_features <- function(list1, list2, mz_ppm = 5, rt_tolerance = 5) {
  mapped <- list()
  for (p1 in list1) {
    candidates <- Filter(function(p2) {
      abs(p1$rtime - p2$rtime) < rt_tolerance
    }, .match_features_mz_candidates(p1$mz, list2, mz_ppm))
    if (length(candidates) > 0L) {
      deltas <- vapply(candidates, function(p2) abs(p1$mz - p2$mz), 0)
      ids <- vapply(candidates, `[[`, "", "id")
      mapped[[as.character(p1$id)]] <- ids[order(deltas, ids)][[1L]]
    }
  }
  cat(sprintf("Of %d list1 features, number of uni-direction matched features is %d.\n", length(list1), length(mapped)))
  mapped
}

# 选择保留时间偏差最小的候选；偏差相同按id排序。
best_rt_match_lcms_features <- function(list1, list2, mz_ppm = 5, rt_tolerance = 5) {
  mapped <- list()
  for (p1 in list1) {
    candidates <- .match_features_mz_candidates(p1$mz, list2, mz_ppm)
    candidates <- Filter(function(p2) abs(p1$rtime - p2$rtime) < rt_tolerance, candidates)
    if (length(candidates) > 0L) {
      deltas <- vapply(candidates, function(p2) abs(p1$rtime - p2$rtime), 0)
      ids <- vapply(candidates, `[[`, "", "id")
      mapped[[as.character(p1$id)]] <- ids[order(deltas, ids)][[1L]]
    }
  }
  cat(sprintf("Of %d list1 features, number of uni-direction matched features is %d.\n", length(list1), length(mapped)))
  mapped
}

# 对应 bidirectional_best_match：返回双向最佳m/z匹配及两个方向字典。
bidirectional_best_match <- function(list1, list2, mz_ppm = 5, rt_tolerance = 5) {
  cat("\n    ~~~ By best rtime matches ~~~     \n\n")
  rt1 <- best_rt_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
  rt2 <- best_rt_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)
  rt_pairs <- c(
    Map(c, names(rt1), unname(rt1)),
    Map(c, unname(rt2), names(rt2))
  )
  pair_keys <- vapply(rt_pairs, paste, collapse = "\r", FUN.VALUE = "")
  cat("~~~ Biodirectional, unique Number of matched feature pairs: ~~~\n", length(pair_keys) - length(unique(pair_keys)), "\n")

  cat("\n\n########################################################################\n")
  cat("    ~~~ By best m/z matches ~~~     \n\n")
  dict1 <- best_mz_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
  dict2 <- best_mz_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)
  reverse_keys <- paste(unname(dict2), names(dict2), sep = "\r")
  # Map会沿用dict1的名称，Python的配对list没有这层额外名称。
  valid <- unname(Map(c, names(dict1), unname(dict1)))
  valid <- valid[vapply(valid, function(pair) paste(pair, collapse = "\r") %in% reverse_keys, FALSE)]
  cat("~~~ Biodirectional, unique Number of matched feature pairs: ~~~\n    ", length(valid), "\n")
  cat("########################################################################\n\n")
  list(valid, dict1, dict2)
}

# 对应 convert_min2secs：把每个特征的分钟保留时间转成秒。
convert_min2secs <- function(LL) {
  lapply(LL, function(point) {
    point$rtime <- point$rtime * 60
    point
  })
}

# 对应 convert_sec2mins：把每个特征的秒保留时间转成分钟。
convert_sec2mins <- function(LL) {
  lapply(LL, function(point) {
    point$rtime <- point$rtime / 60
    point
  })
}
