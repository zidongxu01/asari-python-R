# 对应 asari/mass_functions.py 的质量操作函数。
# 本文件放置 m/z 距离检查、选择性评分、聚类和基于 seed 的分组逻辑。
# 兼容约定：Python 0-based 位置转为 R 1-based，tuple 返回值转为命名 list，
# Python None 转为 NA_integer_；有效科学输入下保留 Python 的计算、排序和边界规则。

# 将二元 tuple 列表展平，并且每个值只保留一次。
#
# 遍历顺序与 Python 版一致：先收集所有第一个元素，再收集所有第二个元素。
# R 的 unique() 会保留首次出现的顺序，因此返回顺序比 Python set() 更稳定。
flatten_tuplelist <- function(tuple_list) {
  if (!is.list(tuple_list)) {
    stop("tuple_list must be a list.", call. = FALSE)
  }
  if (length(tuple_list) == 0L) {
    return(NULL)
  }
  if (any(vapply(tuple_list, length, integer(1)) < 2L)) {
    stop("Each tuple must contain at least two elements.", call. = FALSE)
  }

  first_elements <- lapply(tuple_list, `[[`, 1L)
  second_elements <- lapply(tuple_list, `[[`, 2L)
  unique(unlist(c(first_elements, second_elements), use.names = FALSE))
}

# 在已按升序排列的 m/z 中，找出距离小于 ppm 容差的相邻值。
#
# Python 版返回从 0 开始的 (ii, ii - 1)；这里返回对应的 R 1-based 位置，
# 使调用者可以直接用它索引 R 向量或列表。
check_close_mzs <- function(mzs, mz_tolerance_ppm = 5) {
  if (!is.numeric(mzs) || any(!is.finite(mzs))) {
    stop("mzs must be a finite numeric vector.", call. = FALSE)
  }
  if (length(mz_tolerance_ppm) != 1L ||
      !is.finite(mz_tolerance_ppm) ||
      mz_tolerance_ppm < 0) {
    stop(
      "mz_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }
  if (is.unsorted(mzs, strictly = FALSE)) {
    stop("mzs must be sorted in ascending order.", call. = FALSE)
  }
  if (length(mzs) < 2L) {
    return(list())
  }

  close_pairs <- list()
  for (ii in 2:length(mzs)) {
    tolerance <- mzs[[ii]] * mz_tolerance_ppm * 1e-6
    distance <- mzs[[ii]] - mzs[[ii - 1L]]
    if (distance < tolerance) {
      close_pairs[[length(close_pairs) + 1L]] <- c(ii, ii - 1L)
    }
  }

  close_pairs
}

# 计算排序 m/z 列表中每个 m/z 的质量选择性（mSelectivity）。
#
# 每个分数位于 0 到 1 之间：越接近 1，表示该 m/z 与附近 m/z 越容易区分。
# 实现保留 Python 版的计算方式：先计算相邻 ppm 距离，再将当前 m/z
# 与最多两个较低、两个较高邻居的选择性分数相乘。
calculate_selectivity <- function(sorted_mz_list, std_ppm = 5) {
  if (!is.numeric(sorted_mz_list) || any(!is.finite(sorted_mz_list))) {
    stop("sorted_mz_list must be a finite numeric vector.", call. = FALSE)
  }
  if (length(sorted_mz_list) <= 3L) {
    stop("sorted_mz_list must contain more than three m/z values.", call. = FALSE)
  }
  if (is.unsorted(sorted_mz_list, strictly = FALSE)) {
    stop("sorted_mz_list must be sorted in ascending order.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm <= 0) {
    stop("std_ppm must be one finite, positive number.", call. = FALSE)
  }

  # 将两个 m/z 的 ppm 距离转成 0 到 1 的选择性分数。
  selectivity_from_ppm <- function(ppm_distance) {
    if (ppm_distance > 100) {
      1
    } else if (ppm_distance < 0.1) {
      0
    } else {
      1 - exp(-ppm_distance / std_ppm)
    }
  }

  mz_count <- length(sorted_mz_list)
  ppm_distances <- 1e6 *
    diff(sorted_mz_list) /
    sorted_mz_list[-mz_count]

  # 前两个 m/z 只有不完整的左侧邻居，按 Python 版单独计算。
  selectivities <- c(
    selectivity_from_ppm(ppm_distances[[1L]]) *
      selectivity_from_ppm(ppm_distances[[1L]] + ppm_distances[[2L]]),
    selectivity_from_ppm(ppm_distances[[1L]]) *
      selectivity_from_ppm(ppm_distances[[2L]]) *
      selectivity_from_ppm(ppm_distances[[2L]] + ppm_distances[[3L]])
  )

  # 中间的 m/z 同时考虑两个较低和两个较高的邻居。
  if (mz_count > 4L) {
    for (ii in 3:(mz_count - 2L)) {
      selectivities <- c(
        selectivities,
        selectivity_from_ppm(
          ppm_distances[[ii - 2L]] + ppm_distances[[ii - 1L]]
        ) *
          selectivity_from_ppm(ppm_distances[[ii - 1L]]) *
          selectivity_from_ppm(ppm_distances[[ii]]) *
          selectivity_from_ppm(
            ppm_distances[[ii]] + ppm_distances[[ii + 1L]]
          )
      )
    }
  }

  # 最后两个 m/z 只有不完整的右侧邻居，按 Python 版单独计算。
  distance_count <- length(ppm_distances)
  selectivities <- c(
    selectivities,
    selectivity_from_ppm(
      ppm_distances[[distance_count - 2L]] +
        ppm_distances[[distance_count - 1L]]
    ) *
      selectivity_from_ppm(ppm_distances[[distance_count - 1L]]) *
      selectivity_from_ppm(ppm_distances[[distance_count]]),
    selectivity_from_ppm(
      ppm_distances[[distance_count - 1L]] +
        ppm_distances[[distance_count]]
    ) *
      selectivity_from_ppm(ppm_distances[[distance_count]])
  )

  selectivities
}

# 在两个 m/z 列表之间寻找没有邻近候选冲突的一对一匹配。
#
# 函数先保留每个 m/z 的列表来源和原始位置，再合并排序两个列表。
# 只有来自不同列表、距离小于 ppm 容差，且下一个相邻值不会造成冲突时，
# 当前两个 m/z 才会被记录为匹配。R 版返回 1-based 原始位置。
mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  list1_count <- length(list1)
  list2_count <- length(list2)
  combined_count <- list1_count + list2_count

  if (combined_count < 2L) {
    return(list(mapped = list(), ratio_deltas = numeric()))
  }

  combined_mz <- c(list1, list2)
  list_origin <- c(rep(1L, list1_count), rep(2L, list2_count))
  original_index <- c(seq_along(list1), seq_along(list2))

  # Python tuple 的排序依次使用 m/z、列表来源和原始索引。
  sorted_order <- order(combined_mz, list_origin, original_index)
  combined_mz <- combined_mz[sorted_order]
  list_origin <- list_origin[sorted_order]
  original_index <- original_index[sorted_order]

  # 追加与 Python 版一致的哨兵值，使最后一个真实元素也能检查下一个距离。
  combined_mz <- c(combined_mz, 999999)
  list_origin <- c(list_origin, 2L)
  original_index <- c(original_index, NA_integer_)

  mapped <- list()
  ratio_deltas <- numeric()

  for (ii in 2:combined_count) {
    if (list_origin[[ii]] != list_origin[[ii - 1L]]) {
      tolerance <- combined_mz[[ii]] * std_ppm * 1e-6
      distance <- combined_mz[[ii]] - combined_mz[[ii - 1L]]

      if (distance < tolerance &&
          combined_mz[[ii + 1L]] - combined_mz[[ii]] > tolerance) {
        if (list_origin[[ii]] > list_origin[[ii - 1L]]) {
          mapped[[length(mapped) + 1L]] <- c(
            original_index[[ii - 1L]],
            original_index[[ii]]
          )
          ratio_deltas <- c(ratio_deltas, distance / combined_mz[[ii]])
        } else {
          mapped[[length(mapped) + 1L]] <- c(
            original_index[[ii]],
            original_index[[ii - 1L]]
          )
          ratio_deltas <- c(ratio_deltas, -distance / combined_mz[[ii]])
        }
      }
    }
  }

  list(mapped = mapped, ratio_deltas = ratio_deltas)
}

# 在两个 m/z 列表之间尽可能完成一对一匹配。
#
# 函数先收集来自不同列表且距离小于 ppm 容差的相邻候选对。
# 如果同一个 m/z 出现多个候选，保留绝对 m/z 距离最小的一对。
# 返回成功匹配的 1-based 原始位置，以及两个列表各自未匹配的位置。
complete_mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  list1_count <- length(list1)
  list2_count <- length(list2)
  combined_count <- list1_count + list2_count
  candidates <- list()

  if (combined_count >= 2L) {
    combined_mz <- c(list1, list2)
    list_origin <- c(rep(1L, list1_count), rep(2L, list2_count))
    original_index <- c(seq_along(list1), seq_along(list2))

    # Python tuple 的排序依次使用 m/z、列表来源和原始索引。
    sorted_order <- order(combined_mz, list_origin, original_index)
    combined_mz <- combined_mz[sorted_order]
    list_origin <- list_origin[sorted_order]
    original_index <- original_index[sorted_order]

    for (ii in 2:combined_count) {
      if (list_origin[[ii]] != list_origin[[ii - 1L]]) {
        tolerance <- combined_mz[[ii]] * std_ppm * 1e-6
        distance <- combined_mz[[ii]] - combined_mz[[ii - 1L]]

        if (distance < tolerance) {
          if (list_origin[[ii]] > list_origin[[ii - 1L]]) {
            candidates[[length(candidates) + 1L]] <- list(
              list1_index = original_index[[ii - 1L]],
              distance = distance,
              list2_index = original_index[[ii]]
            )
          } else {
            candidates[[length(candidates) + 1L]] <- list(
              list1_index = original_index[[ii]],
              distance = distance,
              list2_index = original_index[[ii - 1L]]
            )
          }
        }
      }
    }
  }

  # 按 Python 版的顺序化冲突处理：共享任一原始索引时，距离更小者获胜。
  selected <- list()
  if (length(candidates) > 0L) {
    staged <- candidates[[1L]]

    if (length(candidates) > 1L) {
      for (ii in 2:length(candidates)) {
        current <- candidates[[ii]]
        has_conflict <-
          current$list1_index == staged$list1_index ||
          current$list2_index == staged$list2_index

        if (has_conflict) {
          if (current$distance < staged$distance) {
            staged <- current
          }
        } else {
          selected[[length(selected) + 1L]] <- staged
          staged <- current
        }
      }
    }

    selected[[length(selected) + 1L]] <- staged
  }

  mapped <- lapply(
    selected,
    function(candidate) {
      c(candidate$list1_index, candidate$list2_index)
    }
  )

  mapped_list1 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_list2 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  list(
    mapped = mapped,
    list1_unmapped = setdiff(seq_along(list1), mapped_list1),
    list2_unmapped = setdiff(seq_along(list2), mapped_list2)
  )
}

# 返回查询 m/z 在 Centurion 索引中 ppm 容差内的全部候选。
#
# 这是 all_mass_paired_mapping() 的内部搜索辅助函数。索引以
# floor(m/z * 100) 作为桶，查询时检查当前桶及左右相邻桶。
.find_all_mzmatches_centurion_indexed_list <- function(
    query_mz,
    mz_centurion_tree,
    limit_ppm = 5) {
  query_key <- as.integer(query_mz * 100)
  mz_tolerance <- query_mz * limit_ppm * 1e-6
  results <- list()

  for (bucket_key in (query_key - 1L):(query_key + 1L)) {
    bucket <- mz_centurion_tree[[as.character(bucket_key)]]
    if (is.null(bucket)) {
      next
    }

    for (entry in bucket) {
      if (abs(entry[[1L]] - query_mz) < mz_tolerance) {
        results[[length(results) + 1L]] <- entry
      }
    }
  }

  results
}

# 返回两个 m/z 列表之间 ppm 容差内的全部候选配对。
#
# 与前两个配对函数不同，这里不强制一对一，也不在多个候选中只选最近值。
# 一个 list1 m/z 可以对应多个 list2 m/z，反之亦然。函数返回全部
# 1-based 配对位置，以及两个列表中完全没有参与匹配的位置。
all_mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  # 构建与 mass2chem build_centurion_tree_mzlist() 相同概念的百分之一 m/z 索引。
  mz_centurion_tree <- list()
  for (ii in seq_along(list1)) {
    bucket_key <- as.character(as.integer(list1[[ii]] * 100))
    mz_centurion_tree[[bucket_key]] <- c(
      mz_centurion_tree[[bucket_key]],
      list(c(list1[[ii]], ii))
    )
  }

  mapped <- list()
  for (ii in seq_along(list2)) {
    matches <- .find_all_mzmatches_centurion_indexed_list(
      list2[[ii]],
      mz_centurion_tree,
      limit_ppm = std_ppm
    )

    for (match in matches) {
      mapped[[length(mapped) + 1L]] <- c(as.integer(match[[2L]]), ii)
    }
  }

  mapped_list1 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_list2 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  list(
    mapped = mapped,
    list1_unmapped = setdiff(seq_along(list1), mapped_list1),
    list2_unmapped = setdiff(seq_along(list2), mapped_list2)
  )
}

# 在可靠 m/z 配对的基础上估计 list2 的整体质量偏移，必要时校正后重新配对。
#
# 函数先调用 mass_paired_mapping() 获得高可信配对及其相对偏移，
# 再以平均相对偏移作为 list2 的校正比例。为与 Python 原版完全对齐，
# 只有正方向偏移严格大于 correction_tolerance_ppm 时才执行校正；
# 负方向偏移不触发校正。返回 R 1-based 配对和校正前估计的偏移比例。
mass_paired_mapping_with_correction <- function(
    list1,
    list2,
    std_ppm = 5,
    correction_tolerance_ppm = 1) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }
  if (length(correction_tolerance_ppm) != 1L ||
      !is.finite(correction_tolerance_ppm) ||
      correction_tolerance_ppm < 0) {
    stop(
      "correction_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }

  initial_mapping <- mass_paired_mapping(list1, list2, std_ppm)
  mapped <- initial_mapping$mapped
  correction_ratio <- mean(initial_mapping$ratio_deltas)

  if (isTRUE(correction_ratio > correction_tolerance_ppm * 1e-6)) {
    corrected_list2 <- list2 / (1 + correction_ratio)
    mapped <- mass_paired_mapping(
      list1,
      corrected_list2,
      std_ppm
    )$mapped
  }

  list(mapped = mapped, correction_ratio = correction_ratio)
}

# 以高可信 landmark m/z 为锚点，将一个新样本的 m/z 列表加入当前参考列表。
#
# 函数先对齐参考和新样本的 landmarks，在高可信配对数量足够时估计整体
# m/z 偏移，并在正或负偏移绝对值超过阈值时校正新样本。随后对齐剩余 m/z，
# 用配对两端的平均值更新参考 m/z，并把新样本中未配对的 m/z 追加为新参考行。
# R 版使用 1-based landmark 和映射位置，用 NA_integer_ 表示 Python 中的 None。
landmark_guided_mapping <- function(
    REF_reference_mzlist,
    REF_landmarks,
    SM_mzlist,
    SM_landmarks,
    std_ppm = 5,
    correction_tolerance_ppm = 1) {
  if (!is.numeric(REF_reference_mzlist) ||
      any(!is.finite(REF_reference_mzlist)) ||
      !is.numeric(SM_mzlist) ||
      any(!is.finite(SM_mzlist))) {
    stop(
      "REF_reference_mzlist and SM_mzlist must be finite numeric vectors.",
      call. = FALSE
    )
  }
  if (!is.numeric(REF_landmarks) || any(!is.finite(REF_landmarks)) ||
      any(REF_landmarks != as.integer(REF_landmarks)) ||
      any(REF_landmarks < 1L) ||
      any(REF_landmarks > length(REF_reference_mzlist))) {
    stop(
      "REF_landmarks must contain valid 1-based reference positions.",
      call. = FALSE
    )
  }
  if (!is.numeric(SM_landmarks) || any(!is.finite(SM_landmarks)) ||
      any(SM_landmarks != as.integer(SM_landmarks)) ||
      any(SM_landmarks < 1L) ||
      any(SM_landmarks > length(SM_mzlist))) {
    stop(
      "SM_landmarks must contain valid 1-based sample positions.",
      call. = FALSE
    )
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }
  if (length(correction_tolerance_ppm) != 1L ||
      !is.finite(correction_tolerance_ppm) ||
      correction_tolerance_ppm < 0) {
    stop(
      "correction_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }

  REF_landmarks <- as.integer(REF_landmarks)
  SM_landmarks <- as.integer(SM_landmarks)
  original_reference_count <- length(REF_reference_mzlist)
  new_reference_map2 <- rep(NA_integer_, original_reference_count)
  correction_ratio <- NULL

  # 第一步：只使用高可信 landmarks 估计样本间 m/z 偏移。
  reference_anchors <- REF_reference_mzlist[REF_landmarks]
  sample_anchors <- SM_mzlist[SM_landmarks]
  anchor_result <- mass_paired_mapping(
    reference_anchors,
    sample_anchors,
    std_ppm
  )
  anchor_mapping <- anchor_result$mapped

  if (length(anchor_mapping) > 0.2 * original_reference_count) {
    correction_ratio <- mean(anchor_result$ratio_deltas)

    if (abs(correction_ratio) > correction_tolerance_ppm * 1e-6) {
      SM_mzlist <- SM_mzlist / (1 + correction_ratio)
      sample_anchors <- SM_mzlist[SM_landmarks]
      anchor_result <- mass_paired_mapping(
        reference_anchors,
        sample_anchors,
        std_ppm
      )
      anchor_mapping <- anchor_result$mapped
    }
  }

  # 将 landmark 子列表中的位置还原为完整 m/z 列表位置。
  mapped <- lapply(
    anchor_mapping,
    function(pair) {
      c(REF_landmarks[[pair[[1L]]]], SM_landmarks[[pair[[2L]]]])
    }
  )
  mapped_reference <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_sample <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  # 第二步：在排除 landmark 配对后，尽可能配对剩余 m/z。
  remaining_reference <- setdiff(
    seq_along(REF_reference_mzlist),
    mapped_reference
  )
  remaining_sample <- setdiff(seq_along(SM_mzlist), mapped_sample)
  remaining_result <- complete_mass_paired_mapping(
    REF_reference_mzlist[remaining_reference],
    SM_mzlist[remaining_sample],
    std_ppm
  )
  sample_unmatched <- remaining_sample[remaining_result$list2_unmapped]
  remaining_mapped <- lapply(
    remaining_result$mapped,
    function(pair) {
      c(
        remaining_reference[[pair[[1L]]]],
        remaining_sample[[pair[[2L]]]]
      )
    }
  )
  mapped_pairs <- c(mapped, remaining_mapped)

  message(
    sprintf(
      "    mapped pairs = %d / %d ",
      length(mapped_pairs),
      length(SM_mzlist)
    )
  )

  # 成功配对时，记录新样本位置，并用配对两端的平均值更新参考 m/z。
  for (pair in mapped_pairs) {
    new_reference_map2[[pair[[1L]]]] <- pair[[2L]]
    REF_reference_mzlist[[pair[[1L]]]] <- 0.5 * (
      REF_reference_mzlist[[pair[[1L]]]] + SM_mzlist[[pair[[2L]]]]
    )
  }

  # 新样本中没有参考配对的 m/z 追加为 MassGrid 的新参考行。
  for (ii in seq_along(sample_unmatched)) {
    sample_index <- sample_unmatched[[ii]]
    new_reference_map2 <- c(new_reference_map2, sample_index)

    if (sample_index %in% SM_landmarks) {
      REF_landmarks <- c(REF_landmarks, original_reference_count + ii)
    }
  }

  new_reference_mzlist <- c(
    REF_reference_mzlist,
    SM_mzlist[sample_unmatched]
  )

  list(
    new_reference_mzlist = new_reference_mzlist,
    new_reference_map2 = new_reference_map2,
    REF_landmarks = REF_landmarks,
    correction_ratio = correction_ratio
  )
}

# 按当前分箱中位数和用户提供的容差函数，将已排序的 (value, object) 分成连续箱。
#
# 每个新 value 都与当前最后一箱的 value 中位数比较。距离严格小于
# func_tolerance(value) 时放入当前箱，否则开始新箱。返回时只保留每个 tuple 的 object。
bin_by_median <- function(List_of_tuples, func_tolerance) {
  if (!is.list(List_of_tuples) || length(List_of_tuples) == 0L) {
    stop("List_of_tuples must be a non-empty list.", call. = FALSE)
  }
  if (any(vapply(List_of_tuples, length, integer(1)) < 2L)) {
    stop("Each tuple must contain a value and an object.", call. = FALSE)
  }
  if (!is.function(func_tolerance)) {
    stop("func_tolerance must be a function.", call. = FALSE)
  }

  bins <- list(list(List_of_tuples[[1L]]))
  if (length(List_of_tuples) > 1L) {
    for (ii in 2:length(List_of_tuples)) {
      current <- List_of_tuples[[ii]]
      current_bin <- bins[[length(bins)]]
      current_median <- stats::median(
        vapply(current_bin, function(tuple) tuple[[1L]], numeric(1))
      )

      if (current[[1L]] - current_median < func_tolerance(current[[1L]])) {
        bins[[length(bins)]][[length(current_bin) + 1L]] <- current
      } else {
        bins[[length(bins) + 1L]] <- list(current)
      }
    }
  }

  lapply(
    bins,
    function(bin) lapply(bin, `[[`, 2L)
  )
}

# 在已排序 m/z tuple 中找到最大相邻间隔，并从该位置将数据分成左右两组。
#
# mz_tolerance 参数与 Python 原版一样保留在函数接口中，但当前算法并不使用它。
# 如果有多个相同的最大间隔，与 numpy.argmax() 一样选择第一个。
gap_divide_mz_cluster <- function(bin_data_tuples, mz_tolerance) {
  if (!is.list(bin_data_tuples) || length(bin_data_tuples) < 2L) {
    stop("bin_data_tuples must contain at least two tuples.", call. = FALSE)
  }
  if (any(vapply(bin_data_tuples, length, integer(1)) < 1L)) {
    stop("Each tuple must contain an m/z value.", call. = FALSE)
  }

  # 按第一个最大 m/z 间隔分割列表，对应 Python 内嵌函数 __divide_by_largest_gap__()。
  divide_by_largest_gap <- function(tuple_list) {
    mzs <- vapply(tuple_list, function(tuple) tuple[[1L]], numeric(1))
    split_position <- which.max(diff(mzs)) + 1L
    list(
      tuple_list[seq_len(split_position - 1L)],
      tuple_list[split_position:length(tuple_list)]
    )
  }

  divide_by_largest_gap(bin_data_tuples)
}

# 用最近边界扩展复制 SciPy uniform_filter1d(size, mode="nearest") 的一维均值滤波。
#
# asari 向 SciPy 传入整数计数，SciPy 会保留整数 dtype 并截断小数部分；
# 这里显式使用 as.integer() 保留该行为。偶数窗口向左多取一个位置。
.uniform_filter1d_nearest_integer <- function(values, size) {
  if (length(values) == 0L || length(size) != 1L || size < 1L) {
    stop("values must be non-empty and size must be positive.", call. = FALSE)
  }

  size <- as.integer(size)
  left_width <- size %/% 2L
  right_width <- size - left_width - 1L
  padded <- c(
    rep(values[[1L]], left_width),
    values,
    rep(values[[length(values)]], right_width)
  )
  cumulative <- c(0, cumsum(padded))
  value_count <- length(values)
  window_sums <- cumulative[(size + 1L):(size + value_count)] -
    cumulative[seq_len(value_count)]

  as.integer(window_sums / size)
}

# 复制 scipy.signal.find_peaks(values, distance=...) 在 asari 所用参数下的峰位置选择。
#
# 端点不作为峰；平顶峰返回中点并向下取整。峰距离冲突时先保留较高峰，
# 等高时按 SciPy 的顺序先处理较右侧峰，最后将保留位置按升序返回。
.find_peaks_with_distance <- function(values, distance) {
  if (length(distance) != 1L || !is.finite(distance) || distance < 1) {
    stop("distance must be at least 1.", call. = FALSE)
  }
  value_count <- length(values)
  if (value_count < 3L) {
    return(integer())
  }

  candidates <- integer()
  ii <- 2L
  while (ii <= value_count - 1L) {
    if (values[[ii]] > values[[ii - 1L]]) {
      plateau_end <- ii
      while (plateau_end < value_count &&
             values[[plateau_end + 1L]] == values[[ii]]) {
        plateau_end <- plateau_end + 1L
      }

      if (plateau_end < value_count &&
          values[[plateau_end]] > values[[plateau_end + 1L]]) {
        candidates <- c(candidates, (ii + plateau_end) %/% 2L)
      }
      ii <- plateau_end + 1L
    } else {
      ii <- ii + 1L
    }
  }

  if (length(candidates) <= 1L) {
    return(candidates)
  }

  minimum_distance <- ceiling(distance)
  priority_order <- order(-values[candidates], -candidates)
  kept <- integer()
  for (candidate in candidates[priority_order]) {
    if (length(kept) == 0L ||
        all(abs(candidate - kept) >= minimum_distance)) {
      kept <- c(kept, candidate)
    }
  }

  sort(kept)
}

# 将 m/z tuple 投影到万分之一 m/z 整数网格，并识别可作为聚类 seeds 的质量峰。
#
# 函数逐步复制 Python 原版：m/z 乘 10000 后截断为整数、计算每个网格的频数、
# 使用 nearest 模式整数均值滤波，再按 tol4 最小距离找峰。返回值按 0.0001 m/z 还原。
identify_mass_peaks <- function(
    bin_data_tuples,
    mz_tolerance,
    presorted = TRUE) {
  if (!is.list(bin_data_tuples) || length(bin_data_tuples) == 0L) {
    stop("bin_data_tuples must be a non-empty list.", call. = FALSE)
  }
  if (any(vapply(bin_data_tuples, length, integer(1)) < 1L)) {
    stop("Each tuple must contain an m/z value.", call. = FALSE)
  }
  if (length(mz_tolerance) != 1L ||
      !is.finite(mz_tolerance) ||
      mz_tolerance <= 0) {
    stop("mz_tolerance must be one finite, positive number.", call. = FALSE)
  }

  mz4 <- as.integer(
    vapply(bin_data_tuples, function(tuple) tuple[[1L]], numeric(1)) * 10000
  )
  if (!presorted) {
    mz4 <- sort(mz4)
  }

  tol4 <- as.integer(mz_tolerance * 10000)
  if (tol4 < 1L) {
    stop("mz_tolerance must correspond to at least 0.0001 m/z.", call. = FALSE)
  }
  filter_size <- max(2L, as.integer(0.5 * tol4))
  positioned <- seq.int(mz4[[1L]], mz4[[length(mz4)]])
  count_indices <- mz4 - positioned[[1L]] + 1L
  counts <- tabulate(count_indices, nbins = length(positioned))
  filtered_counts <- .uniform_filter1d_nearest_integer(
    counts,
    size = filter_size
  )
  peak_indices <- .find_peaks_with_distance(
    filtered_counts,
    distance = tol4
  )

  0.0001 * positioned[peak_indices]
}

# 根据 identify_mass_peaks() 识别的 m/z seeds，将每个 tuple 分配给距离最近的 seed。
#
# 距离相同时选择位置更靠前的 seed，各聚类内保留原 tuple 顺序。
# 如果没有识别到 seed，与 Python 原版一样改用 gap_divide_mz_cluster() 在最大 m/z 间隔处分成两组。
nn_cluster_by_mz_seeds <- function(
    bin_data_tuples,
    mz_tolerance,
    presorted = TRUE) {
  mz_seeds <- identify_mass_peaks(
    bin_data_tuples,
    mz_tolerance,
    presorted = presorted
  )

  if (length(mz_seeds) > 0L) {
    clusters <- lapply(mz_seeds, function(seed) list())
    for (tuple in bin_data_tuples) {
      distances <- abs(tuple[[1L]] - mz_seeds)
      cluster_index <- which.min(distances)
      clusters[[cluster_index]][[length(clusters[[cluster_index]]) + 1L]] <-
        tuple
    }
    clusters
  } else {
    gap_divide_mz_cluster(bin_data_tuples, mz_tolerance)
  }
}
