# 对应 Python asari/tools/cosine.py：质谱峰配对与余弦相似度。

# 对应 find_matches：返回与NumPy一致的0基峰索引。
find_matches <- function(ref_spec_mz, qry_spec_mz, tolerance, shift = 0.0) {
  ref_spec_mz <- as.numeric(ref_spec_mz)
  qry_spec_mz <- as.numeric(qry_spec_mz)
  matches1 <- integer()
  matches2 <- integer()
  lowest_idx <- 0L
  if (length(ref_spec_mz) == 0L || length(qry_spec_mz) == 0L) {
    return(list(matches1, matches2))
  }
  for (peak1_idx in 0:(length(ref_spec_mz) - 1L)) {
    mz <- ref_spec_mz[[peak1_idx + 1L]]
    low_bound <- mz - tolerance
    high_bound <- mz + tolerance
    if (lowest_idx > length(qry_spec_mz) - 1L) next
    for (peak2_idx in lowest_idx:(length(qry_spec_mz) - 1L)) {
      mz2 <- qry_spec_mz[[peak2_idx + 1L]] - shift
      if (mz2 > high_bound) break
      if (mz2 < low_bound) {
        # 原Python保留当前位置，而不是移动到下一位置。
        lowest_idx <- peak2_idx
      } else {
        matches1 <- c(matches1, peak1_idx)
        matches2 <- c(matches2, peak2_idx)
      }
    }
  }
  list(as.integer(matches1), as.integer(matches2))
}

# 对应 collect_peak_pairs：收集精确/偏移匹配并按贡献分数降序排列。
collect_peak_pairs <- function(
    ref_spec,
    qry_spec,
    min_matched_peak,
    sqrt_transform,
    tolerance,
    shift = 0.0) {
  ref_spec <- as.matrix(ref_spec)
  qry_spec <- as.matrix(qry_spec)
  if (nrow(ref_spec) == 0L || nrow(qry_spec) == 0L) {
    return(list(integer(), integer(), numeric()))
  }
  matches <- find_matches(ref_spec[, 1L], qry_spec[, 1L], tolerance, 0.0)
  idx1 <- matches[[1L]]
  idx2 <- matches[[2L]]
  if (abs(shift) > 1e-6) {
    shifted <- find_matches(ref_spec[, 1L], qry_spec[, 1L], tolerance, shift)
    idx1 <- c(idx1, shifted[[1L]])
    idx2 <- c(idx2, shifted[[2L]])
  }
  if (length(idx1) < as.integer(min_matched_peak)) {
    return(list(integer(), integer(), numeric()))
  }
  products <- ref_spec[idx1 + 1L, 2L] * qry_spec[idx2 + 1L, 2L]
  scores <- if (isTRUE(sqrt_transform)) sqrt(products) else products
  ordering <- order(-scores, seq_along(scores))
  list(as.integer(idx1[ordering]), as.integer(idx2[ordering]), as.numeric(scores[ordering]))
}

# 对应 score_matches：贪心选择不重叠峰对并计算带惩罚的最终分数。
score_matches <- function(
    matches_idx1,
    matches_idx2,
    scores,
    ref_spec,
    qry_spec,
    sqrt_transform,
    penalty) {
  ref_spec <- as.matrix(ref_spec)
  qry_spec <- as.matrix(qry_spec)
  used1 <- rep(FALSE, nrow(ref_spec))
  used2 <- rep(FALSE, nrow(qry_spec))
  total_score <- 0
  used_matches <- 0L
  if (length(matches_idx1) > 0L) {
    for (ii in seq_along(matches_idx1)) {
      idx1 <- matches_idx1[[ii]] + 1L
      idx2 <- matches_idx2[[ii]] + 1L
      if (!used1[[idx1]] && !used2[[idx2]]) {
        total_score <- total_score + scores[[ii]]
        used1[[idx1]] <- TRUE
        used2[[idx2]] <- TRUE
        used_matches <- used_matches + 1L
      }
    }
  }
  if (used_matches == 0L) return(list(0.0, 0L))

  ref_intensity <- ref_spec[, 2L]
  norm1 <- if (isTRUE(sqrt_transform)) {
    sqrt(sum(sqrt(ref_intensity * ref_intensity)))
  } else {
    sqrt(sum(ref_intensity * ref_intensity))
  }
  qry_intensity <- qry_spec[, 2L]
  new_qry_intensity <- ifelse(used2, qry_intensity, qry_intensity * (1 - penalty))
  norm2 <- if (isTRUE(sqrt_transform)) {
    sqrt(sum(sqrt(new_qry_intensity * new_qry_intensity)))
  } else {
    sqrt(sum(new_qry_intensity * new_qry_intensity))
  }
  if (norm1 == 0 || norm2 == 0) return(list(0.0, used_matches))
  list(min(as.numeric(total_score / (norm1 * norm2)), 1.0), used_matches)
}

# 对应 cosine_similarity：先分别按最大强度归一化，再进行峰匹配和评分。
cosine_similarity <- function(
    qry_spec,
    ref_spec,
    tolerance = 0.1,
    min_matched_peak = 1L,
    sqrt_transform = TRUE,
    penalty = 0.0,
    shift = 0.0) {
  qry_spec <- as.matrix(qry_spec)
  ref_spec <- as.matrix(ref_spec)
  if (length(qry_spec) == 0L || length(ref_spec) == 0L ||
      nrow(qry_spec) == 0L || nrow(ref_spec) == 0L) {
    return(list(0.0, 0L))
  }
  storage.mode(qry_spec) <- "double"
  storage.mode(ref_spec) <- "double"
  ref_spec[, 2L] <- ref_spec[, 2L] / max(ref_spec[, 2L])
  qry_spec[, 2L] <- qry_spec[, 2L] / max(qry_spec[, 2L])
  pairs <- collect_peak_pairs(
    ref_spec, qry_spec, min_matched_peak, sqrt_transform,
    as.numeric(tolerance), as.numeric(shift)
  )
  if (length(pairs[[1L]]) == 0L) return(list(0.0, 0L))
  score_matches(
    pairs[[1L]], pairs[[2L]], pairs[[3L]], ref_spec, qry_spec,
    sqrt_transform, as.numeric(penalty)
  )
}
