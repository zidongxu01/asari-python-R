# asariR

`asariR`是Python [asari](https://github.com/shuzhao-li-lab/asari) 1.17.0的R实现，
用于从中心化（centroided）mzML数据中构建质量轨迹、对齐样本、
检测色谱峰并导出特征表。

本项目已经提供LC-MS、DIMS、GC-MS注释、MS2处理、QC、结果比较、
离线Dashboard、交互式向导和Terminal命令入口。所有mzML示例都要求
用户明确写出自己的文件路径，不会暗中选择测试数据。

## 当前状态

- R包可安装，公开函数可直接通过`library(asariR)`使用。
- LC-MS核心流程已使用真实mzML运行并生成完整/推荐特征表。
- 同一文件与Python asari对照时，大部分峰可相互匹配，但数量和编号不是完全相同。
- GC注释需要用户提供Kovats校准表和EI谱库。
- LC的本地候选注释不等于JMS/Khipu经验化合物分组，详见本页的限制说明。

## 安装

先安装mzML读取依赖：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("mzR")
install.packages(c("jsonlite", "remotes"))
```

然后从GitHub安装`asari-r`子目录：

```r
remotes::install_github(
  "zidongxu01/asari-python-R",
  subdir = "asari-r",
  upgrade = "never"
)
```

加载并检查：

```r
library(asariR)
packageVersion("asariR")
asari_available_workflows()
```

## 快速开始：LC-MS

将下面两个路径替换为自己的路径：

```r
library(asariR)

input_file <- "/absolute/path/to/sample.mzML"
output_dir <- "/absolute/path/to/asari-results"

result <- asari_process(
  input = input_file,
  output = output_dir,
  project_name = "my_project",
  mode = "pos",
  ppm = 5,
  multicores = 1,
  rt_align = FALSE,
  database_mode = "memory"
)

result
```

多样本时，把`input`改为包含mzML的目录，并开启RT对齐：

```r
result <- asari_process(
  input = "/absolute/path/to/mzml-folder",
  output = "/absolute/path/to/asari-results",
  project_name = "batch_project",
  mode = "pos",
  ppm = 5,
  multicores = 4,
  rt_align = TRUE
)
```

## 查看结果

```r
tables <- asari_read_results(result, table = "both")

View(tables$preferred)
View(tables$full)

dashboard_file <- asari_dashboard(
  result,
  table = "preferred",
  open = TRUE
)
```

项目目录中最重要的文件是：

- `preferred_Feature_table.tsv`：默认质量标准筛选后的推荐特征。
- `export/full_Feature_table.tsv`：保留更多峰的完整特征表。
- `project.json`：项目参数、样本和运行信息。
- `asari_dashboard_preferred.html`：调用`asari_dashboard()`后生成的离线页面。

## 常用功能

| 目的 | 函数 |
|---|---|
| LC/GC数据预处理 | `asari_process()` |
| DIMS直接进样 | `asari_process_dims()` |
| 分析一个mzML | `asari_analyze()` |
| 导出质量轨迹摘要 | `asari_extract_mass_tracks()` |
| 读取项目表格 | `asari_read_results()` |
| 比较两张特征表 | `asari_compare_features()` |
| 合并多张特征表 | `asari_merge_feature_tables()` |
| 提取目标m/z | `asari_extract_targets()` |
| LC质量/RT候选注释 | `asari_annotate_lc()` |
| GC Kovats + EI谱库注释 | `asari_annotate_gc()` |
| 提取MS2 | `asari_extract_ms2()` |
| MS2匹配MS1特征 | `asari_match_ms2()` |
| MS2本地谱库检索 | `asari_search_ms2()` |
| 特征表QC PDF | `asari_feature_qc()` |
| mzML QC HTML | `asari_qc_report()` |
| MS1/MS2扫描摘要PDF | `asari_scan_summary()` |
| 离线项目页面 | `asari_dashboard()` |
| GC共洗脱特征图 | `asari_feature_graph()` |
| PCA和特征表绘图 | `asari_pca()` / `asari_plot_correlations()` |
| 谱图镜像图 | `asari_mirror_plot()` |
| 交互式操作向导 | `asari_gui()` |
| Terminal命令 | `asari_cli()` |
| Thermo RAW转换 | `asari_convert_raw()` |

查看某个函数的全部参数：

```r
?asari_process
?asari_annotate_lc
?asari_extract_ms2
?asari_dashboard
```

## DIMS

```r
dims_result <- asari_process_dims(
  input = "/absolute/path/to/sample.mzML",
  output = "/absolute/path/to/dims-results",
  project_name = "dims_project",
  mode = "pos",
  ppm = 5
)

dims_table <- asari_read_results(dims_result, table = "full")$full
```

DIMS没有色谱峰形和RT对齐，因此主要结果是`full`表，
不应把空的`preferred`表解释为没有信号。

## LC注释

用户自己的标准品表（至少包含`mz`，可选`rtime`、`id`和`name`）：

```r
annotation <- asari_annotate_lc(
  feature_table = result$preferred_feature_table,
  database = "/absolute/path/to/standards.tsv",
  output = "/absolute/path/to/LC_annotation.tsv",
  mode = "pos",
  ppm = 5,
  rt_tolerance = 5
)
```

使用Python asari安装中的HMDB4 pickle：

```r
options(asariR.python = "/absolute/path/to/python")
options(asariR.db_dir = "/absolute/path/to/asari/db")

annotation <- asari_annotate_lc(
  result$preferred_feature_table,
  database = "hmdb4",
  output = "/absolute/path/to/HMDB_candidates.tsv",
  mode = "pos"
)
```

第一次会将pickle转为RDS缓存；后续检索不再重复解码。

## GC-MS注释

```r
gc_result <- asari_annotate_gc(
  feature_table = "/absolute/path/to/full_Feature_table.tsv",
  kovats = "/absolute/path/to/kovats.tsv",
  library = "/absolute/path/to/library.msp",
  output_dir = "/absolute/path/to/gc-annotation",
  project_name = "gc_project",
  denovo = FALSE,
  mirror_plots = FALSE
)
```

EI谱库必须是包含保留指数和谱峰的MSP或JSON文件。

## MS2

提取一个mzML中的MS2：

```r
ms2 <- asari_extract_ms2(
  input = "/absolute/path/to/sample.mzML",
  output = "/absolute/path/to/extracted_ms2.json",
  min_intensity = 1000
)
```

将MS2匹配到已生成的MS1特征：

```r
matched <- asari_match_ms2(
  feature_table = result$full_feature_table,
  ms2_files = c(
    "/absolute/path/to/sample_01.mzML",
    "/absolute/path/to/sample_02.mzML"
  ),
  output = "/absolute/path/to/matched_ms2.json"
)
```

检索本地MSP、MGF或JSON谱库：

```r
hits <- asari_search_ms2(
  spectra = ms2,
  library = "/absolute/path/to/library.msp",
  output = "/absolute/path/to/ms2_hits.tsv",
  method = "cosine"
)
```

## QC与比较

```r
asari_feature_qc(
  result$full_feature_table,
  "/absolute/path/to/feature_qc.pdf"
)

asari_qc_report(
  "/absolute/path/to/sample.mzML",
  "/absolute/path/to/sample_qc.html"
)

comparison <- asari_compare_features(
  left = "/absolute/path/to/R_table.tsv",
  right = "/absolute/path/to/Python_table.tsv",
  ppm = 5,
  rt_tolerance = 5
)
```

特征编号只是各自运行的顺序标识。比较两个版本时，应使用m/z和RT匹配，
不要仅比较`F33`这样的编号。

## 交互式向导与Terminal

在RStudio中先查看配置：

```r
configuration <- asari_gui(
  input = "/absolute/path/to/sample.mzML",
  output = "/absolute/path/to/results",
  workflow = "LC",
  run = FALSE
)
configuration
```

确认后把`run`改为`TRUE`。

Terminal中可以这样运行：

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- \
  process \
  --input "/absolute/path/to/sample.mzML" \
  --output "/absolute/path/to/results" \
  --project "my_project" \
  --mode pos \
  --ppm 5 \
  --multicores 1
```

其他子命令是`analyze`、`annotate`、`viz`、`join`和`list_workflows`。

## Thermo RAW转换

第一次需要网络下载ThermoRawFileParser：

```r
asari_install_raw_converter()

converted <- asari_convert_raw(
  input = c(
    "/absolute/path/to/sample_01.raw",
    "/absolute/path/to/sample_02.raw"
  ),
  output_dir = "/absolute/path/to/mzml",
  multicores = 2
)
```

macOS/Linux通常还需要可执行`.exe`的Mono环境。已经mzML的用户不需要这个工具。

## 依赖和网络

| 功能 | 额外要求 |
|---|---|
| 处理mzML | `mzR`、`jsonlite`；运行时不需要网络 |
| 本地TSV/CSV/JSON注释 | 不需要网络 |
| 读取Python pickle | 第一次需要可用Python解释器 |
| 离线Dashboard | 无额外R包，不需要网络 |
| RAW转换器首次安装 | 需要网络 |
| 开发测试 | 可选`testthat` |

## 重要限制

1. R版目标是行为尽可能接近Python asari，但不应宣称所有输出逐位完全相同。
2. 峰编号可以不同；请根据m/z和RT判定是否是同一特征。
3. `asari_annotate_lc()`是可直接使用的质量/RT候选匹配，不完全复制外部JMS/Khipu分组。
4. `asari_search_ms2(method = "entropy")`在没有外部后端时使用本地JS相似度，与`ms_entropy`不保证逐值相同。
5. GC、MS2和注释结果取决于用户提供的校准表、谱库版本和容差。
6. 真实多样本项目应先用小批次确认参数和结果，再运行全部数据。

## 开发检查

`testthat`只用于开发测试，不参与mzML计算。

```bash
cd "/absolute/path/to/asari-r"
Rscript scripts/run_tests.R
R CMD check --no-manual .
```

自动测试通过说明已覆盖的函数和边界条件未发现错误，
不单独证明所有真实数据都与Python版完全相同。

## License

MIT License。
