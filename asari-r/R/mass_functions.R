# 对应 asari/mass_functions.py 的质量操作函数。
# 本文件放置 m/z 距离检查、选择性评分、聚类和基于 seed 的分组逻辑。

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

# 根据 m/z seeds 将数据点分配到最近聚类；当前仍是待实现骨架。
nn_cluster_by_mz_seeds <- function(datatuples, mz_tolerance, presorted = FALSE) {
  stop("Not implemented yet: nn_cluster_by_mz_seeds")
}

# 从 m/z 数据中识别可用作聚类中心的质量峰；当前仍是待实现骨架。
identify_mass_peaks <- function(datatuples, mz_tolerance) {
  stop("Not implemented yet: identify_mass_peaks")
}
