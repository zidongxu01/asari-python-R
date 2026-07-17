# asari Python-to-R Mainline / asari Python 转 R 主线

## 中文说明

这个文档的作用是防止我们乱翻译文件。我们的目标不是逐行复制
Python 文件结构，而是用 R 复制 asari 的核心处理流程。

### Python asari 的主线

```text
command_line.py
  -> workflow.py
  -> samples.py / experiment.py
  -> constructors.py
  -> chromatograms.py / mass_functions.py
  -> peaks.py
  -> preferred_Feature_table.tsv
```

这条线的意思是：

```text
接收用户输入和参数
  -> 找到 mzML 文件并创建项目输出目录
  -> 登记样本和实验对象
  -> 构建多样本质量网格和 composite map
  -> 从 mzML 提取 mass tracks / EIC
  -> 在 mass tracks 上找峰
  -> 输出 feature table
```

### 当前 R 原型主线

```text
r-prototype/R/asari_process.R
  -> r-prototype/R/mzml_reader.R
  -> r-prototype/R/mass_tracks.R
  -> r-prototype/R/peak_detection.R
  -> preferred_Feature_table.tsv-style columns
```

### 当前已经完成

| R 文件 | 作用 | 状态 |
| --- | --- | --- |
| `asari_process.R` | 项目外壳：找 mzML、登记样本、创建输出目录、写 `project.json`；显式指定单个 mzML 时可写非空 feature table。 | 单样本原型完成 |
| `mzml_reader.R` | 用 `mzR` 读取 mzML，得到 `scan_id`, `rtime`, `mz`, `intensity`。 | 初版完成 |
| `mass_tracks.R` | 把相近 m/z 的点按 retention time 串成 mass tracks / EIC。 | 粗原型完成 |
| `peak_detection.R` | 在 mass tracks 上找峰，生成 feature 表前 11 列。 | 粗原型完成 |
| `preview_utils.R` | 预览脚本辅助函数，要求用户明确指定 mzML 文件。 | 完成 |

### 当前还没有完成

| Python 模块 | 还缺什么 |
| --- | --- |
| `samples.py` | R 版还没有正式的 sample 对象，目前只是简单 data frame。 |
| `experiment.py` | R 版还没有完整实验对象。 |
| `constructors.py` | 还没有 `MassGrid` / `CompositeMap` / 多样本 feature alignment。 |
| `mass_functions.py` | `asari-r/R/mass_functions.R` 已完成全部顶层函数的独立移植，并完成 Python/SciPy 行为对照；尚未接入 `constructors.py` / MassGrid 完整流程。 |
| `peaks.py` | 只做了粗 peak detection，还没有完整 prominence、Gaussian fitting、cSelectivity 逻辑。 |
| output export | 单样本非空 `preferred_Feature_table.tsv` 已接入；多样本输出还没完成。 |

### 下一步顺序

1. 对照 Python asari 的 `preferred_Feature_table.tsv` 做单样本结果比较。
2. 把已完成的 `mass_functions.py` R 函数接入 `constructors.py` / MassGrid。
3. 验证多样本 m/z 对齐和 reference 更新逻辑。

### 规则

在主线能稳定输出可比较的 feature table 之前，不做 dashboard、
annotation、notebooks、GC annotation 等旁支功能。

---

## English Notes

This document keeps the Python-to-R conversion focused. The goal is not to copy
the Python file structure line by line. The goal is to reproduce asari's core
processing workflow in R.

### Python asari Mainline

```text
command_line.py
  -> workflow.py
  -> samples.py / experiment.py
  -> constructors.py
  -> chromatograms.py / mass_functions.py
  -> peaks.py
  -> preferred_Feature_table.tsv
```

This means:

```text
parse user input and parameters
  -> find mzML files and create project output folders
  -> register samples and experiment state
  -> build mass grids and composite maps across samples
  -> extract mass tracks / EICs from mzML files
  -> detect elution peaks on mass tracks
  -> export feature tables
```

### Current R Prototype Mainline

```text
r-prototype/R/asari_process.R
  -> r-prototype/R/mzml_reader.R
  -> r-prototype/R/mass_tracks.R
  -> r-prototype/R/peak_detection.R
  -> preferred_Feature_table.tsv-style columns
```

### Implemented So Far

| R file | Purpose | Status |
| --- | --- | --- |
| `asari_process.R` | Project shell: find mzML files, register samples, create output folders, write `project.json`; can write a non-empty feature table when one mzML file is explicitly selected. | Single-sample prototype done |
| `mzml_reader.R` | Read mzML with `mzR` and return `scan_id`, `rtime`, `mz`, `intensity`. | Initial version done |
| `mass_tracks.R` | Group nearby m/z points over retention time into mass tracks / EICs. | Rough prototype done |
| `peak_detection.R` | Detect peaks on mass tracks and emit the first 11 feature-table columns. | Rough prototype done |
| `preview_utils.R` | Helper for preview scripts; requires users to explicitly provide the mzML file path. | Done |

### Not Implemented Yet

| Python module | Missing R work |
| --- | --- |
| `samples.py` | No formal R sample object yet; currently represented by a simple data frame. |
| `experiment.py` | No full R experiment object yet. |
| `constructors.py` | No `MassGrid`, `CompositeMap`, or multi-sample feature alignment yet. |
| `mass_functions.py` | All top-level functions are ported in `asari-r/R/mass_functions.R` and checked against Python/SciPy behavior; integration with `constructors.py` / MassGrid remains. |
| `peaks.py` | Only rough peak detection exists; full prominence, Gaussian fitting, and cSelectivity logic are not ported. |
| output export | Single-sample non-empty `preferred_Feature_table.tsv` is wired; multi-sample export is not complete. |

### Next Steps

1. Compare single-sample R output against Python asari's `preferred_Feature_table.tsv`.
2. Connect the completed `mass_functions.py` R functions to `constructors.py` / MassGrid.
3. Validate multi-sample m/z alignment and reference updates.

### Rule

Do not work on dashboard, annotation, notebooks, or GC annotation until the
mainline can produce a comparable feature table.
