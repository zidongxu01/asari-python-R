# 从任意工作目录定位asari-r根目录，避免testthat切换目录后找不到R源码。
arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (length(script_argument) == 1L) {
  script_file <- sub("^--file=", "", script_argument)
  package_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
} else if (file.exists("DESCRIPTION")) {
  package_root <- normalizePath(".", mustWork = TRUE)
} else {
  stop("Unable to locate the asari-r package root.", call. = FALSE)
}
setwd(package_root)

library(testthat)

# 把测试中的隐式绘图设备定向到临时文件，避免在源码目录留下Rplots.pdf。
old_device <- getOption("device")
options(device = function(...) {
  grDevices::pdf(file = tempfile("asariR-test-", fileext = ".pdf"), ...)
})
on.exit(options(device = old_device), add = TRUE)

# 先按R包加载顺序读取全部模块；部分旧测试只包含断言，没有重复source目标模块。
r_files <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
invisible(lapply(r_files, source, local = .GlobalEnv))

# 手工登记每个测试文件，使旧的逐文件source测试和统一汇总可以同时工作。
test_files <- sort(list.files("tests/testthat", pattern = "^test-.*[.]R$", full.names = TRUE))
reporter <- ListReporter$new()
with_reporter(reporter, {
  for (test_file in test_files) {
    reporter$start_file(basename(test_file))
    source(test_file, local = .GlobalEnv)
    reporter$end_file()
  }
})
results <- reporter$get_results()
results_frame <- as.data.frame(results)

# 输出稳定的总数，失败时返回非零状态，便于Terminal和CI判断。
passed <- sum(results_frame$passed)
failed <- sum(results_frame$failed) + sum(results_frame$error)
warnings <- sum(results_frame$warning)
skipped <- sum(results_frame$skipped)
cat(sprintf(
  "\nPASS %d | FAIL %d | WARN %d | SKIP %d\n",
  passed, failed, warnings, skipped
))
if (failed > 0L || warnings > 0L) {
  print(results_frame[
    results_frame$failed > 0L | results_frame$error |
      results_frame$warning > 0L,
    , drop = FALSE
  ])
  stop("The asariR test suite did not pass cleanly.", call. = FALSE)
}
