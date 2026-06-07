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
| `asari_process.R` | 项目外壳：找 mzML、登记样本、创建输出目录、写 `project.json` 和空 feature table。 | 初版完成 |
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
| `mass_functions.py` | 只做了少量 ppm helper，还没有完整 mass matching 和 clustering。 |
| `peaks.py` | 只做了粗 peak detection，还没有完整 prominence、Gaussian fitting、cSelectivity 逻辑。 |
| output export | 还没有把真实检测结果完整接入 `asari_process()` 写成非空 `preferred_Feature_table.tsv`。 |

### 下一步顺序

1. 把 `mass_tracks.R` 和 `peak_detection.R` 接入 `asari_process()`。
2. 先实现单样本版本，输出真正非空的 feature table。
3. 对照 Python asari 的 `preferred_Feature_table.tsv` 做结果比较。
4. 开始 port `mass_functions.py` 的核心匹配函数。
5. 再做 `constructors.py` 里的多样本对齐逻辑。

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
| `asari_process.R` | Project shell: find mzML files, register samples, create output folders, write `project.json` and an empty feature table. | Initial version done |
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
| `mass_functions.py` | Only minimal ppm helpers exist; full mass matching and clustering are not ported. |
| `peaks.py` | Only rough peak detection exists; full prominence, Gaussian fitting, and cSelectivity logic are not ported. |
| output export | Detected features are not fully wired into `asari_process()` as a non-empty `preferred_Feature_table.tsv`. |

### Next Steps

1. Wire `mass_tracks.R` and `peak_detection.R` into `asari_process()`.
2. Implement a single-sample path that writes a real non-empty feature table.
3. Compare R output against Python asari's `preferred_Feature_table.tsv`.
4. Port the core matching helpers from `mass_functions.py`.
5. Implement the multi-sample alignment logic from `constructors.py`.

### Rule

Do not work on dashboard, annotation, notebooks, or GC annotation until the
mainline can produce a comparable feature table.
