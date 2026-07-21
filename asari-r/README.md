# asariR：安装、运行、参数与结果说明

`asariR`是Python [asari](https://github.com/shuzhao-li-lab/asari) 1.17.0的R实现，
用于高分辨率质谱代谢组学数据预处理。当前公开R包已经接入LC-MS核心流程：

1. 读取中心化（centroided）mzML文件；
2. 构建质量轨道（mass tracks / EIC）；
3. 建立单样本或多样本质量网格；
4. 进行保留时间对齐；
5. 构建复合质量轨道；
6. 检测色谱峰并计算峰质量指标；
7. 导出完整特征表和推荐特征表。

本文档面向第一次使用本软件的人，所有示例都要求用户明确提供输入路径，
程序不会自动寻找隐藏或默认的mzML测试文件。

## 目录

- [1. 当前状态和适用范围](#1-当前状态和适用范围)
- [2. 使用前必须准备什么](#2-使用前必须准备什么)
- [3. 安装方法](#3-安装方法)
- [4. 检查是否安装成功](#4-检查是否安装成功)
- [5. 一分钟快速开始](#5-一分钟快速开始)
- [6. 当前公开指令总览](#6-当前公开指令总览)
- [7. asari_process完整说明](#7-asari_process完整说明)
- [8. 常见运行示例](#8-常见运行示例)
- [9. asari_default_parameters完整说明](#9-asari_default_parameters完整说明)
- [10. 全部高级参数说明](#10-全部高级参数说明)
- [11. 输出目录和文件说明](#11-输出目录和文件说明)
- [12. 特征表每一列的含义](#12-特征表每一列的含义)
- [13. 如何读取和检查结果](#13-如何读取和检查结果)
- [14. 常见错误与解决方法](#14-常见错误与解决方法)
- [15. 更新、重新安装和卸载](#15-更新重新安装和卸载)
- [16. 开发检查和testthat](#16-开发检查和testthat)
- [17. 当前限制](#17-当前限制)

## 1. 当前状态和适用范围

### 1.1 当前最适合处理的数据

- 高分辨率LC-MS数据；
- 已经转换为mzML格式的数据；
- 已经中心化（centroided）的质谱数据；
- 正离子或负离子模式；
- 单个mzML文件；
- 同一实验中的多个mzML文件。

### 1.2 当前已经验证的内容

- R包可以通过`R CMD INSTALL`安装；
- `library(asariR)`可以正常加载；
- `asari_process()`可以处理明确指定的真实mzML文件；
- 已验证样本提取出2992条质量轨道和356个保留时间参考峰；
- 真实单样本运行生成3038个完整特征和2201个推荐特征；
- 与Python asari逐峰比较时，完整表和推荐表的对应率约为99.7%至99.8%；
- 当前自动化测试结果为724项通过、0项失败、0项警告。

### 1.3 当前不应该当作完全交付功能的部分

以下模块已有R对应代码，但当前公开安装版本还没有完成全部真实数据、外部数据库和后端集成验证：

- GC-MS完整工作流；
- MS2谱图处理和谱库检索；
- LC/GC自动注释；
- Dashboard；
- 图形化桌面GUI；
- Python风格的完整命令行子命令；
- 大型真实多样本项目的全面对照验证。

当前最稳妥的公开用法是：通过R中的`asari_process()`运行LC-MS核心流程。

## 2. 使用前必须准备什么

### 2.1 R

本项目已经在R 4.5、macOS ARM64环境中完成安装和真实数据验证。
其他常见R平台理论上也可以使用，但仍应在目标电脑上执行安装检查和小样本测试。

在R中查看当前版本：

```r
R.version.string
```

### 2.2 mzML文件

输入文件必须满足：

- 文件扩展名是严格的`.mzML`；
- 文件真实存在；
- 数据已经centroided；
- 同一批输入中不能有两个同名mzML文件，即使它们位于不同目录。

例如下面的文件名可以被识别：

```text
sample_01.mzML
sample_02.mzML
QC_reference.mzML
```

下面的扩展名不会被当前入口自动识别：

```text
sample_01.mzml
sample_01.MZML
sample_01.raw
```

### 2.3 网络什么时候需要

| 操作 | 是否需要网络 |
|---|---:|
| 第一次安装R依赖 | 需要 |
| 从GitHub安装或更新asariR | 需要 |
| 处理电脑上的mzML文件 | 不需要 |
| 读取已经生成的TSV和JSON结果 | 不需要 |
| 从本地源码重新安装 | 不需要，前提是依赖已经安装 |

## 3. 安装方法

### 3.1 安装mzR和jsonlite

`mzR`负责读取mzML，`jsonlite`负责写出项目参数和记录。这两个包是核心LC流程必需依赖。

在R控制台运行：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("mzR")
install.packages("jsonlite")
```

检查依赖是否安装成功：

```r
requireNamespace("mzR", quietly = TRUE)
requireNamespace("jsonlite", quietly = TRUE)
```

两条命令都应该返回：

```text
TRUE
```

### 3.2 从GitHub安装asariR

第一次使用`remotes`时先安装：

```r
install.packages("remotes")
```

然后安装本项目中的`asari-r`子目录：

```r
remotes::install_github(
  "zidongxu01/asari-python-R",
  subdir = "asari-r",
  upgrade = "never"
)
```

`upgrade = "never"`表示安装asariR时不自动升级电脑上其他无关R包，能够减少环境变化。

### 3.3 从已经下载的本地源码安装

假设项目目录是：

```text
/Users/yourname/Documents/asari-python-R
```

在Terminal中进入项目根目录：

```bash
cd "/Users/yourname/Documents/asari-python-R"
R CMD INSTALL asari-r
```

路径包含空格时必须使用引号。

Windows PowerShell示例：

```powershell
cd "C:\Users\yourname\Documents\asari-python-R"
R CMD INSTALL asari-r
```

### 3.4 安装开发测试工具testthat（普通用户可跳过）

`testthat`只用于开发阶段自动检查代码，不参与mzML计算：

```r
install.packages("testthat")
```

普通用户只处理数据时，不需要为了运行`asari_process()`而安装`testthat`。

## 4. 检查是否安装成功

### 4.1 加载软件包

```r
library(asariR)
```

没有报错就说明R已经找到这个软件包。

### 4.2 查看版本

```r
packageVersion("asariR")
```

当前开发版本应显示类似：

```text
[1] '0.0.0.9000'
```

### 4.3 查看安装位置

```r
find.package("asariR")
```

### 4.4 查看公开函数

```r
getNamespaceExports("asariR")
```

当前主要公开函数是：

```text
asari_process
asari_default_parameters
```

### 4.5 查看帮助

```r
?asari_process
?asari_default_parameters
help(package = "asariR")
```

## 5. 一分钟快速开始

### 5.1 处理一个mzML文件

```r
library(asariR)

result <- asari_process(
  input = "/Users/yourname/data/sample_01.mzML",
  output = "/Users/yourname/results",
  project_name = "sample_01",
  multicores = 1,
  rt_align = FALSE
)
```

单文件不需要跨样本保留时间对齐，因此建议设置：

```r
rt_align = FALSE
```

### 5.2 处理一个目录中的全部mzML文件

```r
result <- asari_process(
  input = "/Users/yourname/data/my_batch",
  output = "/Users/yourname/results",
  project_name = "my_batch",
  multicores = 4,
  rt_align = TRUE
)
```

### 5.3 运行成功后会看到什么

程序会打印类似：

```text
asariR processing completed.
Project directory: /Users/yourname/results/run_my_batch_721175930
Preferred feature table: /Users/yourname/results/run_my_batch_721175930/preferred_Feature_table.tsv
Full feature table: /Users/yourname/results/run_my_batch_721175930/export/full_Feature_table.tsv
```

同时，`result`对象会保存这些路径：

```r
result$project_dir
result$preferred_feature_table
result$full_feature_table
result$input_files
result$parameters
```

## 6. 当前公开指令总览

### 6.1 Terminal指令

| 指令 | 用途 | 普通用户是否常用 |
|---|---|---:|
| `R CMD INSTALL asari-r` | 从本地源码安装 | 是 |
| `R CMD build asari-r` | 构建R源码包压缩文件 | 否，开发者使用 |
| `R CMD check ...` | 检查包结构和文档 | 否，开发者使用 |
| `Rscript scripts/smoke_check.R` | 运行最小源码检查 | 否，开发者使用 |

### 6.2 R指令

| 指令 | 用途 |
|---|---|
| `library(asariR)` | 加载软件包 |
| `packageVersion("asariR")` | 查看版本 |
| `find.package("asariR")` | 查看安装位置 |
| `?asari_process` | 查看处理函数帮助 |
| `asari_process(...)` | 处理mzML并导出特征表 |
| `asari_default_parameters()` | 查看一份完整默认参数副本 |
| `print(result)`或直接输入`result` | 显示结果路径 |
| `read.delim(...)` | 读取TSV特征表 |
| `remove.packages("asariR")` | 卸载软件包 |

### 6.3 为什么没有`asari process ...`命令

Python asari使用命令行子命令。当前公开R包首先交付的是已经验证的R函数入口，
还没有安装一个面向系统Terminal的`asariR process`可执行程序。

因此当前应使用：

```r
asari_process(...)
```

不要把尚未公开的内部`main.R`函数当作稳定用户接口。

## 7. asari_process完整说明

完整函数形式：

```r
asari_process(
  input,
  output = "asari-output",
  project_name = "asari_project",
  mode = NULL,
  ppm = NULL,
  multicores = NULL,
  rt_align = NULL,
  database_mode = NULL,
  recursive = FALSE,
  parameters = list()
)
```

### 7.1 每个参数的含义

| 参数 | 类型 | 默认行为 | 说明 |
|---|---|---|---|
| `input` | 字符串或字符向量 | 必须提供 | 一个mzML文件、多个mzML文件，或一个/多个目录 |
| `output` | 单个字符串 | `"asari-output"` | 存放时间戳项目目录的父目录 |
| `project_name` | 单个字符串 | `"asari_project"` | 项目名称，不能包含`/`或反斜线 |
| `mode` | `"pos"`、`"neg"`或`NULL` | `NULL`时使用默认`"pos"` | 离子模式 |
| `ppm` | 正数或`NULL` | `NULL`时使用5 | m/z匹配容差，映射到`mz_tolerance_ppm` |
| `multicores` | 正整数或`NULL` | `NULL`时使用默认参数4 | 并行工作进程数 |
| `rt_align` | `TRUE`、`FALSE`或`NULL` | `NULL`时使用`TRUE` | 是否进行保留时间对齐 |
| `database_mode` | `"memory"`、`"ondisk"`、`"auto"`或`NULL` | 便捷入口默认`"memory"` | 中间数据存储模式 |
| `recursive` | 逻辑值 | `FALSE` | 是否递归查找输入目录下的子目录 |
| `parameters` | 具名list | `list()` | 高级参数覆盖值 |

### 7.2 input

输入单文件：

```r
input = "/absolute/path/sample.mzML"
```

输入目录：

```r
input = "/absolute/path/batch"
```

输入多个明确文件：

```r
input = c(
  "/absolute/path/sample_01.mzML",
  "/absolute/path/sample_02.mzML",
  "/absolute/path/QC.mzML"
)
```

输入多个目录：

```r
input = c(
  "/absolute/path/batch_1",
  "/absolute/path/batch_2"
)
```

规则：

- 所有路径都会被转换为绝对规范路径；
- 文件会排序，保证输入顺序稳定；
- 重复的完全相同路径会自动去重；
- 不同目录中存在同名mzML时会报错；
- 默认不进入子目录；
- 只识别严格的`.mzML`扩展名。

### 7.3 output

`output`是父目录，不是最终项目目录。每次运行都会在下面创建一个新目录：

```r
output = "/absolute/path/results"
```

可能生成：

```text
/absolute/path/results/run_my_project_721175930
```

如果`output`不存在，程序会递归创建；如果这个路径已经是一个普通文件，程序会报错。

### 7.4 project_name

推荐只使用：

- 字母；
- 数字；
- 下划线；
- 短横线。

推荐示例：

```r
project_name = "MT_batch14_pos"
```

不要使用路径分隔符：

```r
project_name = "batch/14"   # 错误
```

### 7.5 mode

正离子：

```r
mode = "pos"
```

负离子：

```r
mode = "neg"
```

必须根据数据真实采集模式设置，不能为了获得更多峰随意切换。

### 7.6 ppm

```r
ppm = 5
```

`ppm`决定m/z匹配和质量轨道合并容差。默认5适合许多高分辨率质谱数据。
更小的值要求更高质量精度；更大的值会增加不同离子被匹配到一起的风险。

修改前应依据仪器质量精度和QC结果，不建议只为了增加特征数盲目增大。

### 7.7 multicores

```r
multicores = 1
```

含义是允许使用的工作进程数：

- 初次测试建议使用1；
- 多文件正式运行可以使用2、4或其他合适正整数；
- 数值不能是0、负数或小数；
- 增加进程数也会增加内存使用。

### 7.8 rt_align

单文件推荐：

```r
rt_align = FALSE
```

同一实验的多个文件通常使用：

```r
rt_align = TRUE
```

当设置为`TRUE`时，软件会选择参考样本、寻找高质量landmark峰，并建立样本到参考样本的扫描映射。

### 7.9 database_mode

#### memory

```r
database_mode = "memory"
```

- 中间质量轨道保存在R内存中；
- 当前便捷入口默认使用该模式；
- 普通LC处理不需要Python；
- 最适合单样本和较小项目；
- 大项目可能占用较多内存。

#### ondisk

```r
database_mode = "ondisk"
```

- 中间数据写到磁盘；
- 当前实现为了兼容Python asari pickle，需要本机存在Python解释器；
- 可以通过R选项明确指定Python：

```r
options(asariR.python = "/absolute/path/to/python")
```

也可以在启动R前设置环境变量：

```bash
export ASARIR_PYTHON="/absolute/path/to/python"
```

#### auto

```r
database_mode = "auto"
```

根据样本数量和`project_sample_number_small`选择存储模式。
由于它可能选择`ondisk`，当前公开使用中应确保Python可用，或者直接明确使用`memory`。

### 7.10 recursive

目录结构：

```text
data/
├── batch_1/
│   └── a.mzML
└── batch_2/
    └── b.mzML
```

需要递归查找时：

```r
asari_process(
  input = "/absolute/path/data",
  output = "/absolute/path/results",
  recursive = TRUE
)
```

### 7.11 parameters

只覆盖需要修改的项目：

```r
parameters = list(
  min_peak_height = 50000,
  min_intensity_threshold = 500,
  peak_area = "auc"
)
```

推荐传入“短覆盖列表”，不要在没有检查的情况下修改所有参数。

参数优先级是：

```text
asari_process显式参数 > parameters具名列表 > 默认参数
```

例如：

```r
asari_process(
  input = "/path/sample.mzML",
  output = "/path/results",
  ppm = 5,
  parameters = list(mz_tolerance_ppm = 10)
)
```

最终使用`ppm = 5`，因为显式参数优先。

## 8. 常见运行示例

### 8.1 正离子单文件

```r
library(asariR)

result <- asari_process(
  input = "/data/sample_pos.mzML",
  output = "/results",
  project_name = "sample_pos",
  mode = "pos",
  ppm = 5,
  multicores = 1,
  rt_align = FALSE
)
```

### 8.2 负离子单文件

```r
result <- asari_process(
  input = "/data/sample_neg.mzML",
  output = "/results",
  project_name = "sample_neg",
  mode = "neg",
  ppm = 5,
  multicores = 1,
  rt_align = FALSE
)
```

### 8.3 多样本目录并进行RT对齐

```r
result <- asari_process(
  input = "/data/experiment_01",
  output = "/results",
  project_name = "experiment_01",
  mode = "pos",
  ppm = 5,
  multicores = 4,
  rt_align = TRUE,
  database_mode = "memory"
)
```

### 8.4 明确指定参考样本

参考样本可以写完整路径或文件名；内部会按basename匹配：

```r
result <- asari_process(
  input = "/data/experiment_01",
  output = "/results",
  project_name = "experiment_01_refQC",
  rt_align = TRUE,
  parameters = list(
    reference = "QC_reference.mzML"
  )
)
```

如果不指定，软件通常选择anchor m/z配对数量较多的样本作为参考。

### 8.5 降低最低峰高

```r
result <- asari_process(
  input = "/data/sample.mzML",
  output = "/results",
  parameters = list(
    min_peak_height = 50000
  )
)
```

没有同时明确提供派生阈值时，便捷入口会计算：

```text
min_prominence_threshold = as.integer(0.33 × min_peak_height)
cal_min_peak_height = 10 × min_peak_height
```

降低阈值可能增加弱峰，也可能增加噪声峰，必须结合QC和空白样本判断。

### 8.6 修改峰面积算法

简单求和：

```r
parameters = list(peak_area = "sum")
```

局部最大过滤后的近似AUC：

```r
parameters = list(peak_area = "auc")
```

Gaussian拟合面积：

```r
parameters = list(peak_area = "gauss")
```

不要在同一研究的不同批次中随意混用不同峰面积算法。

### 8.7 使用Windows路径

R字符串中推荐使用正斜线：

```r
result <- asari_process(
  input = "C:/Users/yourname/data/sample.mzML",
  output = "C:/Users/yourname/results",
  multicores = 1,
  rt_align = FALSE
)
```

### 8.8 路径包含中文或空格

只要路径写在R字符串引号中即可：

```r
result <- asari_process(
  input = "/Users/yourname/Documents/实验 数据/sample.mzML",
  output = "/Users/yourname/Documents/实验 结果"
)
```

## 9. asari_default_parameters完整说明

查看默认参数：

```r
defaults <- asari_default_parameters()
defaults
```

查看参数名：

```r
names(defaults)
```

查看一个参数：

```r
defaults$min_peak_height
defaults$mz_tolerance_ppm
defaults$rt_align_on
```

每次调用都会返回独立列表：

```r
p1 <- asari_default_parameters()
p2 <- asari_default_parameters()

p1$mode <- "neg"
p2$mode
```

`p2$mode`仍然是`"pos"`，修改`p1`不会污染全局默认值。

如果只改变少数参数，推荐：

```r
result <- asari_process(
  input = "/data/sample.mzML",
  output = "/results",
  parameters = list(
    min_peak_height = 50000,
    peak_area = "sum"
  )
)
```

只有在确实需要保存一份完整配置时，才传入完整默认参数对象：

```r
parameters <- asari_default_parameters()
parameters$mode <- "neg"
parameters$mz_tolerance_ppm <- 5

result <- asari_process(
  input = "/data/sample.mzML",
  output = "/results",
  parameters = parameters
)
```

注意：完整列表中的每个字段都会被视为用户明确提供的值，包括派生阈值。

## 10. 全部高级参数说明

### 10.1 项目和存储参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `project_name` | `"asari_project"` | 项目名；公开入口的同名显式参数会覆盖它 |
| `outdir` | `"output"` | 内部输出前缀；公开入口根据`output`重新设置 |
| `keep_intermediates` | `FALSE` | 保留自Python的配置字段；当前R公开流程尚未完整接线 |
| `reuse_intermediates` | `NULL` | 预留中间文件复用；当前公开流程不建议启用 |
| `database_mode` | `"ondisk"` | 默认参数值；公开便捷入口在未指定时改用`"memory"` |
| `multicores` | `4L` | 并行工作进程数 |

### 10.2 质量和质量轨道参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `mode` | `"pos"` | 离子模式：`"pos"`或`"neg"` |
| `mass_range` | `c(50, 2000)` | 保留的asari配置范围；当前mzML核心提取没有直接使用该字段过滤 |
| `mz_tolerance_ppm` | `5` | m/z聚类、轨道合并和跨样本匹配容差 |
| `correction_tolerance_ppm` | `1` | 是否执行质量校正的ppm阈值 |
| `min_timepoints` | `6L` | 一条有效色谱轨道/峰至少需要的数据点数量 |
| `min_intensity_threshold` | `1000` | 构建质量轨道前过滤低强度点 |
| `min_intensity` | `1000` | 早期R原型兼容别名；真实提取使用`min_intensity_threshold` |

### 10.3 峰检测参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `signal_noise_ratio` | `2` | 峰相对局部噪声的最低倍数 |
| `min_peak_height` | `100000` | 最低峰高 |
| `min_peak_ratio` | `0.001` | 小峰相对同ROI最高峰的最低比例 |
| `wlen` | `25L` | 计算峰prominence使用的保留时间点窗口 |
| `autoheight` | `FALSE` | 自动估计最低峰高；当前公开核心包不建议启用 |
| `gaussian_shape` | `0.5` | 峰检测阶段Gaussian拟合质量下限 |
| `peak_area` | `"sum"` | 峰面积方法：`"sum"`、`"auc"`或`"gauss"` |
| `min_prominence_threshold` | 动态生成 | 默认是`as.integer(0.33 × min_peak_height)` |

### 10.4 autoheight相关参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `min_min_peak_height` | `10` | 自动估计允许返回的最低下限 |
| `dynamic_range` | `1000` | 自动估计动态范围控制参数 |
| `num_files_to_check` | `20L` | 自动估计时抽查的文件数 |

当前GitHub公开核心安装版本可能不包含可选`analyze`模块，因此保持：

```r
autoheight = FALSE
```

需要自动估计时，应先确认完整源码树和真实多样本验证环境都已准备好。

### 10.5 实验和工作流参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `sample_metadata` | `NULL` | 为未来实验设计指导预留，当前核心流程未直接使用 |
| `workflow` | `"LC"` | 工作流类型；当前公开稳定入口应保持`"LC"` |
| `project_sample_number_small` | `10L` | `auto`存储模式和MassGrid策略的小项目分界 |

### 10.6 保留时间对齐参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `reference` | `NULL` | 指定参考样本；可使用完整路径或basename |
| `rt_align_method` | `"lowess"` | 保留的对齐方法名；当前主流程使用LOWESS式校准 |
| `rt_align_on` | `TRUE` | 是否启用RT对齐 |
| `debug_rtime_align` | `FALSE` | 是否输出RT对齐调试信息/图；公开流程尚未全面验证 |
| `drop_unaligned_samples` | `FALSE` | 是否从最终表中删除RT对齐失败样本 |
| `rtime_tolerance` | `50` | 保留的特征RT容差配置；当前核心路径并非所有阶段都直接读取 |
| `cal_min_peak_height` | `100000` | RT校准landmark的最低峰高；便捷入口通常重算为`10 × min_peak_height` |
| `peak_number_rt_calibration` | `20L` | RT校准需要的高质量峰数目标 |
| `max_retention_shift` | `NULL` | RT校准允许的最大峰对位移；`NULL`时内部按数据决定 |
| `num_lowess_iterations` | `3L` | LOWESS迭代次数 |

### 10.7 GC-MS和MS2参数

以下参数已保留并有相应模块代码，但不属于当前公开LC核心入口的完整交付功能：

| 参数 | 默认值 | 预期作用 |
|---|---:|---|
| `ms2_tolerance_in_ppm` | `10` | MS2匹配ppm容差 |
| `ms2_tolerance_in_da` | `0.005` | MS2匹配Da容差 |
| `ri_tolerance` | `30` | 保留指数匹配容差 |
| `score_cutoff_cosine` | `0.4` | 余弦相似度阈值 |
| `score_cutoff_entropy` | `0.4` | 谱图熵相似度阈值 |
| `corr_cutoff` | `0.6` | 特征与base peak跨样本相关性阈值 |
| `max_ri_delta` | `100` | 候选注释最大RI差 |
| `do_mirror_plot` | `TRUE` | 是否生成镜像谱图 |
| `max_core_features` | `20000` | de novo分析特征数上限 |
| `denovo` | `FALSE` | 是否运行经验化合物de novo构建 |

### 10.8 LC注释参数

| 参数 | 默认值 | 当前状态 |
|---|---:|---|
| `anno` | `FALSE` | 当前公开LC核心运行应保持`FALSE` |
| `check_isotope_ratio` | `FALSE` | 注释阶段同位素比例检查开关 |
| `khipu_rtime_tolerance` | `2` | khipu式经验化合物分组RT容差，单位秒 |

完整注释需要额外数据库和适配器。只设置`anno = TRUE`并不能自动安装这些外部后端。

### 10.9 输出和界面参数

| 参数 | 默认值 | 当前作用 |
|---|---:|---|
| `compress` | `FALSE` | 保留的中间文件压缩配置；当前不建议修改 |
| `output_feature_table` | `"Feature_table.tsv"` | 生成`full_`和`preferred_`文件名的基础名称 |
| `mass_grid_mapping` | `"_mass_grid_mapping.csv"` | 质量网格文件名 |
| `annotation_filename` | `"Annotation_table.tsv"` | 预留注释表文件名 |
| `json_empricalCompounds` | `"_empCpd_json.json"` | 保留Python原字段拼写的经验化合物JSON名 |
| `table_for_viz` | `"preferred"` | Dashboard选择推荐表或完整表 |
| `vizualization_max_samples` | `20` | Dashboard显示样本数上限 |

## 11. 输出目录和文件说明

典型LC项目：

```text
results/
└── run_my_experiment_721175930/
    ├── preferred_Feature_table.tsv
    ├── project.json
    ├── README.txt
    └── export/
        ├── full_Feature_table.tsv
        └── _mass_grid_mapping.csv
```

### 11.1 preferred_Feature_table.tsv

推荐优先用于下游统计分析的特征表。

当前默认筛选条件是严格满足：

```text
detection_counts > 0
snr > 2
goodness_fitting > 0.7
cSelectivity > 0.7
```

推荐表不是“绝对真实化合物列表”，仍需要QC、空白过滤、批次检查和后续注释。

### 11.2 export/full_Feature_table.tsv

包含峰检测阶段保留下来的完整特征集合。适合：

- 自己设计不同过滤规则；
- 检查推荐表排除了哪些峰；
- 进行方法学比较；
- 研究弱峰或低选择性区域。

### 11.3 export/_mass_grid_mapping.csv

记录实验级质量网格：

- `mz`列是质量网格代表m/z；
- 每个样本列保存该网格行对应的样本质量轨道ID；
- 缺失对应时可能为空或为缺失值；
- CSV第一列是R写出的行号。

该文件主要用于调试跨样本m/z对齐，不是通常的最终统计输入。

### 11.4 project.json

记录本次运行的有效参数和项目摘要，例如：

- 输入模式；
- ppm；
- 峰高阈值；
- RT对齐设置；
- 参考样本；
- 样本数和扫描数；
- 输出目录；
- 被删除的样本。

为了可重复分析，应把`project.json`与特征表一起保存。

### 11.5 README.txt

每次运行生成的简短结果说明。它来自asari项目级默认说明，不等同于当前这份完整软件手册。

### 11.6 pickle目录

`database_mode = "memory"`时不会依赖磁盘pickle中间数据。

`database_mode = "ondisk"`时，运行过程可能创建`pickle/`目录，并在处理结束后清理中间文件。

## 12. 特征表每一列的含义

| 列名 | 含义 |
|---|---|
| `id_number` | 实验级特征ID，例如`F0`、`F1` |
| `mz` | 特征代表m/z，导出时保留4位小数 |
| `rtime` | 峰顶保留时间，单位秒，保留2位小数 |
| `rtime_left_base` | 峰左边界保留时间，单位秒 |
| `rtime_right_base` | 峰右边界保留时间，单位秒 |
| `parent_masstrack_id` | 产生该峰的实验级复合质量轨道ID，保留Python兼容的0基编号语义 |
| `peak_area` | 复合轨道上的峰面积摘要，算法由`peak_area`参数决定 |
| `cSelectivity` | 色谱选择性，通常越接近1表示峰与邻近峰越容易区分 |
| `goodness_fitting` | Gaussian峰形拟合质量，越高通常表示峰形越接近模型 |
| `snr` | 峰高相对局部噪声的信噪比指标 |
| `detection_counts` | 在多少个有效样本中提取到大于1的面积值 |
| 后续样本列 | 每个样本在该峰边界内计算的峰面积；列名来自mzML文件名 |

注意：

- `peak_area`列是实验级复合轨道的峰面积；
- 每个样本列是该样本自己的对应质量轨道面积；
- 样本列不是简单的峰顶高度；
- 保留时间以秒为单位，不是分钟。

## 13. 如何读取和检查结果

### 13.1 读取推荐表

```r
preferred <- read.delim(
  result$preferred_feature_table,
  check.names = FALSE
)
```

### 13.2 读取完整表

```r
full <- read.delim(
  result$full_feature_table,
  check.names = FALSE
)
```

### 13.3 查看行列数

```r
dim(preferred)
dim(full)
```

### 13.4 查看列名和前几行

```r
names(preferred)
head(preferred)
```

### 13.5 检查m/z和RT范围

```r
range(preferred$mz, na.rm = TRUE)
range(preferred$rtime, na.rm = TRUE)
```

### 13.6 检查质量指标

```r
summary(preferred$snr)
summary(preferred$cSelectivity)
summary(preferred$goodness_fitting)
summary(preferred$detection_counts)
```

### 13.7 识别样本列

固定元数据列共有11个：

```r
metadata_columns <- c(
  "id_number", "mz", "rtime", "rtime_left_base", "rtime_right_base",
  "parent_masstrack_id", "peak_area", "cSelectivity",
  "goodness_fitting", "snr", "detection_counts"
)

sample_columns <- setdiff(names(preferred), metadata_columns)
sample_columns
```

### 13.8 查看本次有效参数

内存中的结果对象：

```r
result$parameters
```

从项目JSON读取：

```r
project <- jsonlite::fromJSON(
  file.path(result$project_dir, "project.json")
)

project$mz_tolerance_ppm
project$min_peak_height
project$reference
```

### 13.9 保存本次R会话环境

```r
sessionInfo()
```

建议把`sessionInfo()`输出和项目结果一起保存，便于记录R版本和依赖版本。

## 14. 常见错误与解决方法

### 14.1 Input path does not exist

原因：输入路径不存在或拼写错误。

检查：

```r
file.exists("/absolute/path/sample.mzML")
dir.exists("/absolute/path/data")
```

### 14.2 No .mzML files were found

可能原因：

- 目录为空；
- 文件位于子目录但`recursive = FALSE`；
- 扩展名不是严格的`.mzML`；
- 输入仍然是RAW等原始厂商格式。

检查：

```r
list.files("/absolute/path/data", recursive = TRUE)
```

### 14.3 Duplicate mzML file name(s) are not supported

原因：不同目录中存在相同basename，例如：

```text
batch_1/sample.mzML
batch_2/sample.mzML
```

解决方法：运行前给文件使用唯一名称，例如：

```text
batch_1_sample.mzML
batch_2_sample.mzML
```

### 14.4 Missing required R package(s): mzR或jsonlite

重新安装依赖：

```r
BiocManager::install("mzR")
install.packages("jsonlite")
```

### 14.5 mode must be either 'pos' or 'neg'

只能使用：

```r
mode = "pos"
```

或者：

```r
mode = "neg"
```

### 14.6 ppm must be one positive number

错误示例：

```r
ppm = 0
ppm = -5
ppm = "five"
```

正确示例：

```r
ppm = 5
```

### 14.7 multicores must be one positive integer

错误示例：

```r
multicores = 0
multicores = 2.5
```

正确示例：

```r
multicores = 1
multicores = 4
```

### 14.8 output exists but is not a directory

原因：`output`指向了一个普通文件，而不是目录。

重新指定父目录：

```r
output = "/absolute/path/results"
```

### 14.9 ondisk模式找不到Python

最简单的解决方法是对普通项目使用：

```r
database_mode = "memory"
```

确实需要磁盘模式时明确设置：

```r
options(asariR.python = "/absolute/path/to/python")
```

### 14.10 autoheight requires the optional analyze module

当前公开核心版本不建议开启`autoheight`。使用明确峰高：

```r
parameters = list(
  autoheight = FALSE,
  min_peak_height = 100000
)
```

### 14.11 No data was processed

常见原因：

- mzML不是centroided；
- mzML损坏或不兼容；
- 强度阈值过高；
- 峰高阈值过高；
- `mzR`无法读取文件；
- 输入实际上没有可用MS1扫描。

建议先用一个明确文件、一个核心运行：

```r
result <- asari_process(
  input = "/absolute/path/one_sample.mzML",
  output = "/absolute/path/test_results",
  multicores = 1,
  rt_align = FALSE,
  database_mode = "memory"
)
```

### 14.12 内存不足或运行很慢

可以依次尝试：

1. 先用一个mzML确认流程；
2. 把`multicores`设为1，减少并发内存；
3. 分批处理数据；
4. 在配置好Python后评估`ondisk`模式；
5. 检查是否把无关目录递归加入输入。

### 14.13 为什么每次输出目录名字不一样

软件会加入时间戳，防止覆盖旧项目。这是预期行为。

应该从返回对象取得准确目录：

```r
result$project_dir
```

不要在脚本中猜测时间戳。

## 15. 更新、重新安装和卸载

### 15.1 从GitHub更新

```r
remotes::install_github(
  "zidongxu01/asari-python-R",
  subdir = "asari-r",
  upgrade = "never",
  force = TRUE
)
```

重新启动R，然后检查：

```r
packageVersion("asariR")
```

### 15.2 从本地源码重新安装

在项目根目录运行：

```bash
R CMD INSTALL asari-r
```

如果R正在加载旧版本，先重新启动R会话。

### 15.3 卸载

```r
remove.packages("asariR")
```

卸载软件包不会自动删除已经生成的结果目录。

### 15.4 查看安装目录

```r
find.package("asariR")
```

## 16. 开发检查和testthat

### 16.1 testthat是什么

`testthat`是R开发测试框架，用来自动检查：

- 函数输入输出；
- 参数验证；
- Python/R行为对应；
- 峰检测；
- 质量轨道和对齐；
- Experiment和Workflow；
- 文件写出和读取；
- 新的公开处理入口。

它不会参与真实mzML计算。

### 16.2 当前测试结果

当前完整测试审计结果：

```text
PASS 724
FAIL 0
WARN 0
SKIP 0
```

### 16.3 最小smoke check

从`asari-r`目录运行：

```bash
Rscript scripts/smoke_check.R
```

期望输出：

```text
Smoke check passed.
```

### 16.4 构建源码包

从项目根目录运行：

```bash
R CMD build asari-r
```

### 16.5 检查包结构

构建后可运行：

```bash
R CMD check --no-manual --no-tests asariR_0.0.0.9000.tar.gz
```

这里的`--no-tests`只检查包结构、依赖、命名空间、文档和示例；
它不替代`testthat`功能测试。

旧测试文件仍在从“逐文件source模式”整理为标准安装包测试布局，
因此普通用户目前不需要直接执行内部完整测试目录。

## 17. 当前限制

1. 当前真实数据逐峰对照主要完成于单样本LC-MS；
2. 多样本m/z和RT逻辑已有自动测试，但仍需要更多真实项目验证；
3. GC-MS、MS2、注释和Dashboard需要额外依赖、数据库或后端；
4. 当前没有安装到Terminal PATH中的完整CLI程序；
5. `autoheight`依赖可选分析模块，公开核心入口默认关闭；
6. `ondisk`兼容pickle模式需要Python；
7. R和Python最终少数峰仍可能因数值库与边界语义产生微小差异；
8. 自动测试通过不代表所有仪器、采集方法和异常文件都已经覆盖。

如果是第一次使用，推荐从以下最小命令开始：

```r
library(asariR)

result <- asari_process(
  input = "/absolute/path/to/one_sample.mzML",
  output = "/absolute/path/to/results",
  project_name = "first_test",
  multicores = 1,
  rt_align = FALSE,
  database_mode = "memory"
)
```

确认单文件结果正常后，再逐步增加样本数量和高级参数。
