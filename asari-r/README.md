# asariR

`asariR`是Python asari 1.17.0的R实现，用于高分辨率质谱代谢组学数据预处理。
当前R包已经接入LC-MS核心流程，并在真实单样本数据上完成主流程验证。

## 当前能力

- 从中心化mzML文件提取质量轨道/EIC。
- 构建单样本或多样本质量网格。
- 进行保留时间对齐和复合质量轨道构建。
- 检测色谱峰并计算峰面积、信噪比和峰质量指标。
- 导出完整特征表和推荐特征表。

GC-MS、MS2、注释、Dashboard和部分外部数据库功能已经有对应代码，
但仍需要相应数据库、可选依赖或后端适配器才能完整使用。

## 安装依赖

处理mzML需要Bioconductor的`mzR`和CRAN的`jsonlite`：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("mzR")
install.packages("jsonlite")
```

如果需要读取YAML参数文件，可以另外安装：

```r
install.packages("yaml")
```

## 安装asariR

从本项目根目录执行：

```bash
R CMD INSTALL asari-r
```

代码上传到GitHub后，也可以通过以下方式安装：

```r
install.packages("remotes")
remotes::install_github("zidongxu01/asari-python-R", subdir = "asari-r")
```

## 最简单的运行方式

输入路径必须由用户明确提供；程序不会偷偷搜索默认mzML文件：

```r
library(asariR)

result <- asari_process(
  input = "/absolute/path/to/mzML-files",
  output = "/absolute/path/to/results",
  project_name = "my_experiment",
  multicores = 1
)
```

`input`既可以是一个mzML文件，也可以是包含多个mzML文件的目录。
运行成功后会直接打印：

- 本次项目目录；
- `preferred_Feature_table.tsv`的位置；
- `export/full_Feature_table.tsv`的位置。

返回对象也保存了这些路径：

```r
result$project_dir
result$preferred_feature_table
result$full_feature_table
```

## 常用参数

```r
result <- asari_process(
  input = "/absolute/path/to/sample.mzML",
  output = "/absolute/path/to/results",
  mode = "pos",
  ppm = 5,
  multicores = 1,
  rt_align = FALSE,
  database_mode = "memory"
)
```

便捷入口默认使用`database_mode = "memory"`，普通LC处理因此不需要Python。
对于大型项目，可以选择`"ondisk"`；为了保留与Python asari兼容的pickle中间文件，
该模式需要本机存在Python解释器。

完整参数可以这样查看和修改：

```r
parameters <- asari_default_parameters()
parameters$min_peak_height <- 50000

result <- asari_process(
  input = "/absolute/path/to/mzML-files",
  output = "/absolute/path/to/results",
  parameters = parameters
)
```

## 输出目录

每次运行都会在指定的`output`下面创建独立的时间戳目录，例如：

```text
results/
└── run_my_experiment_719143012/
    ├── preferred_Feature_table.tsv
    ├── project.json
    ├── README.txt
    └── export/
        ├── full_Feature_table.tsv
        └── _mass_grid_mapping.csv
```

## 测试

安装`testthat`后，从`asari-r`目录运行：

```r
install.packages("testthat")
```

```bash
Rscript scripts/run_tests.R
Rscript scripts/smoke_check.R
```

正式发布前仍需要继续完成真实多样本项目、GC-MS、MS2和注释工作流的端到端验证。
