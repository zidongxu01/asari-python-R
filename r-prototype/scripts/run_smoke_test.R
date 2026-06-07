source("r-prototype/R/asari_process.R")

result <- asari_process(
  input_dir = "实验数据",
  output_dir = file.path(tempdir(), "asari-r-prototype"),
  project_name = "smoke"
)

cat("Created output folder:\n")
cat(result$parameters$outdir, "\n")
cat("Registered samples:\n")
print(result$samples[, c("sample_id", "sample_name")])

