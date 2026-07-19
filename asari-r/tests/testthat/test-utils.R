# utils.py 的 10 个 def 必须逐一存在明确的 R 对应函数。
test_that("all 10 Python utils defs have explicit R counterparts", {
  expected_functions <- c(
    "NpEncoder_default",
    "bulk_process",
    "download_and_unzip_to_pkg_resources",
    "download_and_unzip",
    "validate_mzml_file",
    "build_boolean_dict",
    "sizeof_fmt",
    "checksum_file",
    "wait_with_pbar",
    "get_ionization_mode_mzml"
  )

  expect_length(expected_functions, 10L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# NpEncoder 保留类身份，并把标量和 ndarray 对应结构转换为基础 R 类型。
test_that("NpEncoder default converts numeric and array values", {
  encoder <- NpEncoder()

  expect_s3_class(encoder, "NpEncoder")
  expect_identical(encoder$default(3L), 3L)
  expect_equal(encoder$default(3.5), 3.5)
  expect_equal(
    encoder$default(matrix(1:6, nrow = 2L, byrow = TRUE)),
    list(as.list(1:3), as.list(4:6))
  )
  expect_equal(encoder$default(c(1, 2)), list(1, 2))
  expect_error(encoder$default(list(a = 1)), "not JSON serializable")
})

# bulk_process 必须保持输入顺序，并把每个 argument 作为唯一参数传给 command。
test_that("bulk_process applies command in argument order", {
  result <- bulk_process(
    function(value) value * value,
    as.list(1:5),
    jobs_per_worker = 1
  )

  expect_equal(unlist(result), c(1, 4, 9, 16, 25))
  expect_equal(bulk_process(identity, list(), jobs_per_worker = 1), list())
  expect_error(
    bulk_process("not a function", list(1), jobs_per_worker = 1),
    "command must be a function"
  )
})

# 创建本地 ZIP，避免下载测试访问外网。
make_utils_zip <- function() {
  source_dir <- tempfile("utils_zip_source_")
  dir.create(source_dir)
  writeLines("payload", file.path(source_dir, "payload.txt"))
  archive <- tempfile(fileext = ".zip")
  old_directory <- setwd(source_dir)
  on.exit(setwd(old_directory), add = TRUE)
  suppressWarnings(utils::zip(archive, "payload.txt"))
  archive
}

# download_and_unzip 应下载并恢复 ZIP 中的文件内容。
test_that("download_and_unzip extracts a local archive", {
  archive <- make_utils_zip()
  target <- tempfile("utils_extract_")

  expect_message(
    expect_null(download_and_unzip(paste0("file://", archive), target)),
    "Extracted to"
  )
  expect_true(file.exists(file.path(target, "payload.txt")))
  expect_equal(readLines(file.path(target, "payload.txt")), "payload")
})

# package resource 版本按 Python 逻辑解压到 package 目录父级下的 subdir。
test_that("download_and_unzip_to_pkg_resources resolves package parent", {
  archive <- make_utils_zip()
  package_parent <- tempfile("fake_library_")
  package_path <- file.path(package_parent, "fakePackage")
  dir.create(package_path, recursive = TRUE)
  subdir <- "downloaded_data"
  expected_path <- file.path(package_parent, subdir, "payload.txt")

  expect_message(
    expect_null(download_and_unzip_to_pkg_resources(
      paste0("file://", archive), package_path, subdir
    )),
    "Extracted to"
  )
  expect_true(file.exists(expected_path))
})

# 无效或不存在的 mzML 必须返回 FALSE，而不是把解析异常传播到调用者。
test_that("validate_mzml_file rejects invalid files", {
  invalid <- tempfile(fileext = ".mzML")
  writeLines("not mzML", invalid)

  expect_false(validate_mzml_file(invalid))
  expect_false(validate_mzml_file(tempfile(fileext = ".mzML")))
})

# Python 支持的十种布尔输入形式必须全部映射到正确逻辑值。
test_that("build_boolean_dict reproduces all Python keys", {
  booleans <- build_boolean_dict()

  expect_named(
    booleans,
    c("T", "F", "1", "0", "True", "False", "TRUE", "FALSE", "true", "false")
  )
  expect_equal(
    unlist(booleans, use.names = FALSE),
    c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
})

# sizeof_fmt 的单位进位、负数和自定义后缀与 Python 格式一致。
test_that("sizeof_fmt uses Python IEC units and formatting", {
  expect_equal(sizeof_fmt(0), "0.0B")
  expect_equal(sizeof_fmt(1023), "1023.0B")
  expect_equal(sizeof_fmt(1024), "1.0KiB")
  expect_equal(sizeof_fmt(-1024^2), "-1.0MiB")
  expect_equal(sizeof_fmt(1536, suffix = "bytes"), "1.5Kibytes")
  expect_equal(sizeof_fmt(1024^8), "1.0YiB")
})

# checksum_file 返回标准小写 MD5，并拒绝目录或不存在的路径。
test_that("checksum_file computes MD5", {
  path <- tempfile()
  writeChar("hello", path, eos = NULL)

  expect_equal(checksum_file(path), "5d41402abc4b2a76b9719d911017c592")
  expect_error(checksum_file(tempdir()), "regular file")
  expect_error(checksum_file(tempfile()), "regular file")
  expect_error(checksum_file(path, chunksize = 0), "Invalid chunksize")
})

# wait=0 与 Python range(0) 一样立即完成，非法计数应明确失败。
test_that("wait_with_pbar handles zero and invalid waits", {
  expect_null(wait_with_pbar(0))
  expect_error(wait_with_pbar(-1), "non-negative integer")
  expect_error(wait_with_pbar(1.5), "non-negative integer")
})

# polarity header 覆盖正、负、混合、limit 和未设置模式。
test_that("get_ionization_mode_mzml identifies scan polarity", {
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(1, 1)))),
    "pos"
  )
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(0, 0)))),
    "neg"
  )
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(1, 0)))),
    "mixed"
  )
  expect_equal(
    get_ionization_mode_mzml(
      list(header = data.frame(polarity = c(1, 0))), limit = 0
    ),
    "pos"
  )
  expect_error(
    get_ionization_mode_mzml(list(header = data.frame(polarity = -1))),
    "neither positive nor negative"
  )
})

# 保留原 R 项目已使用的 ppm 换算兼容函数。
test_that("ppm_to_mz_tolerance converts ppm to absolute m/z tolerance", {
  expect_equal(ppm_to_mz_tolerance(100, 5), 0.0005)
  expect_equal(ppm_to_mz_tolerance(800, 5), 0.004)
})
