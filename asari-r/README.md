# asariR

`asariR`是[Python asari](https://github.com/shuzhao-li-lab/asari) 1.17.0的R语言实现。
它读取中心化（centroided）的mzML文件，构建质量轨迹、检测色谱峰，
并生成可以继续统计分析的特征表。

主要流程：

```text
mzML文件 → 质量轨迹 → 峰检测 → 特征表
```

## 当前状态

- R包可以安装，LC-MS核心流程可以处理真实mzML并输出结果。
- 在一份完成严格对照的真实文件中，Python的2234个推荐峰都被R找到，
  但R额外报告了4个峰；因此目前不能宣称两个版本完全相同。
- 示例文件在当前开发电脑上的R处理时间约为10秒；实际时间取决于数据和电脑。
- LC、GC、DIMS、MS2、QC和离线结果页面均有公开R函数入口。

## 输入要求

输入必须是中心化（centroided）的`.mzML`文件。

```r
input_file <- "/absolute/path/to/sample.mzML"
file.exists(input_file)
```

结果必须是`TRUE`。不要使用不存在的路径、文件夹路径或0字节文件。

## 安装

先安装mzML读取依赖：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("mzR")
install.packages(c("jsonlite", "remotes"))
```

再从GitHub安装asariR：

```r
remotes::install_github(
  "zidongxu01/asari-python-R",
  subdir = "asari-r",
  upgrade = "never"
)
```

加载软件：

```r
library(asariR)
packageVersion("asariR")
```

## 第一次运行

只需要修改`input_file`和`output_dir`：

```r
library(asariR)

input_file <- "/absolute/path/to/sample.mzML"
output_dir <- "/absolute/path/to/asariR-results"

stopifnot(file.exists(input_file))

result <- asari_process(
  input = input_file,
  output = output_dir,
  project_name = "my_first_project",
  mode = "pos",
  ppm = 5,
  multicores = 4,
  rt_align = FALSE,
  database_mode = "memory"
)

tables <- asari_read_results(result, table = "both")
View(tables$preferred)
```

单文件使用`rt_align = FALSE`。多样本时，把`input`改成包含mzML的文件夹，
并设置`rt_align = TRUE`。

## 结果在哪里

运行结束后，R会打印项目目录和两张特征表的完整路径。

| 文件 | 含义 |
|---|---|
| `preferred_Feature_table.tsv` | 经过质量筛选的推荐峰，通常优先使用 |
| `export/full_Feature_table.tsv` | 完整峰列表，包含更多低质量候选峰 |
| `project.json` | 输入文件、参数和项目信息 |

也可以读取以前生成的项目：

```r
tables <- asari_read_results(
  "/absolute/path/to/generated-project-folder",
  table = "both"
)

View(tables$preferred)
View(tables$full)
```

`F1`、`F2`等编号只是本次运行的行号。比较Python和R时，
应根据`mz`和`rtime`判断是否为同一个峰，不要只比较编号。

## Commands

Command是在Terminal中告诉软件要执行哪类工作的词。例如：

```text
process --input sample.mzML --ppm 5
```

这里`process`是command，`--input`和`--ppm`是参数。

| Command | 用途 | 对应R函数 |
|---|---|---|
| `process` | 正式处理mzML | `asari_process()` / `asari_process_dims()` |
| `analyze` | 检查一个mzML | `asari_analyze()` |
| `annotate` | 对LC或GC特征表做注释 | `asari_annotate_lc()` / `asari_annotate_gc()` |
| `viz` | 生成离线结果页面 | `asari_dashboard()` |
| `join` | 合并两张或更多特征表 | `asari_merge_feature_tables()` |
| `list_workflows` | 查看支持的工作流 | `asari_available_workflows()` |

第一次使用主要需要`analyze`和`process`。

### Terminal格式

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- COMMAND [OPTIONS]
```

查看帮助：

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- --help
```

检查一个mzML：

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- \
  analyze --input "/absolute/path/to/sample.mzML" --mode pos --ppm 5
```

正式处理：

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- \
  process \
  --input "/absolute/path/to/sample.mzML" \
  --output "/absolute/path/to/results" \
  --project "my_project" \
  --mode pos --ppm 5 --multicores 4
```

查看全部工作流：

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- list_workflows
```

其他commands的参数可以通过`--help`和下方对应R函数的帮助页面查看。

## 常用参数

| 参数 | 含义 | 常见值 |
|---|---|---|
| `input` | 一个mzML或包含mzML的文件夹 | 必须明确填写 |
| `output` | 结果保存位置 | 建议使用绝对路径 |
| `project_name` | 项目名称 | 任意简短名称 |
| `mode` | 离子模式 | `"pos"`或`"neg"` |
| `ppm` | m/z容差 | 常用`5` |
| `multicores` | CPU核心数 | 可先使用`4` |
| `rt_align` | 是否进行保留时间对齐 | 单文件`FALSE`，多文件`TRUE` |
| `database_mode` | 中间数据保存方式 | 普通用户建议`"memory"` |

如果不确定参数，应先使用示例设置处理一个文件。

## 其他功能

| 目的 | R函数 |
|---|---|
| DIMS处理 | `asari_process_dims()` |
| 导出质量轨迹 | `asari_extract_mass_tracks()` |
| 比较特征表 | `asari_compare_features()` |
| 提取目标m/z | `asari_extract_targets()` |
| LC/GC注释 | `asari_annotate_lc()` / `asari_annotate_gc()` |
| 提取、匹配和检索MS2 | `asari_extract_ms2()` / `asari_match_ms2()` / `asari_search_ms2()` |
| QC和统计图 | `asari_feature_qc()` / `asari_qc_report()` / `asari_pca()` |
| 离线结果页面 | `asari_dashboard()` |
| RStudio操作向导 | `asari_gui()` |
| Thermo RAW转换 | `asari_convert_raw()` |

查看某个函数的参数和示例：

```r
?asari_process
?asari_analyze
?asari_annotate_lc
?asari_extract_ms2
?asari_dashboard
```

## Python与网络要求

| 操作 | 要求 |
|---|---|
| 安装R包 | 需要联网 |
| 处理普通mzML | 不需要Python，运行时不需要联网 |
| `database_mode = "ondisk"` | 需要Python写出兼容pickle |
| 读取Python pickle | 需要Python解释器 |
| 首次安装Thermo RAW转换器 | 需要联网；macOS/Linux通常还需要Mono |

已经拥有中心化mzML且使用`database_mode = "memory"`时，
核心处理流程是纯R的。

## 已知限制

1. R版与Python版非常接近，但尚未在大量不同数据上证明完全一致。
2. 峰编号可能不同，版本比较应使用m/z和保留时间。
3. `asari_annotate_lc()`提供本地质量/RT候选匹配，不完全等同于外部JMS/Khipu分组。
4. GC和MS2结果依赖用户提供的校准表、谱库版本和容差。
5. 大型项目应先用少量文件确认模式、ppm和峰检测结果。

## 开发测试

`testthat`只用于检查代码，不参与mzML计算。

```bash
cd "/absolute/path/to/asari-r"
Rscript scripts/run_tests.R
```

测试通过只表示已覆盖的代码没有发现错误，不单独证明所有真实数据与Python完全一致。

## License

MIT License。
