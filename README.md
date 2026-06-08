# asari Python-to-R

中文 | [English](#english)

## 中文说明

这个仓库是
[asari](https://github.com/shuzhao-li/asari) 的 R 语言移植项目。asari
原本是一个用于高分辨率 LC-MS 代谢组数据处理的 Python 程序。

本项目的目标是用 R 复现 asari 的核心科学处理主线，而不是逐行复制
Python 源码目录。当前工作集中在 [`r-prototype/`](r-prototype/) 原型中；
科学主线和移植检查点记录在 [`docs/mainline.md`](docs/mainline.md)。

### 范围

R 版实现围绕 asari 的核心处理流程展开：

1. 读取用户明确指定的 mzML 输入。
2. 登记样本和运行参数。
3. 从 MS1 扫描中提取 mass tracks / extracted ion chromatograms。
4. 在 mass tracks 上检测色谱峰。
5. 写出 asari 风格的 feature table。
6. 从单样本输出扩展到多样本 alignment。

这个仓库不追求逐文件复制 Python 包结构。如果 R 原生组织方式更清晰，
就优先使用 R 原生结构，同时保持科学行为与 Python asari 对齐。

### 当前状态

当前活跃代码位于 [`r-prototype/`](r-prototype/)。

| 模块 | 状态 |
| --- | --- |
| mzML 读取 | 已用 `mzR` 实现初版读取器。 |
| 默认参数 | 已从 Python asari 翻译原型参数。 |
| Mass tracks / EICs | 已实现粗略单样本原型。 |
| Peak detection | 粗略原型可输出 feature-table-style 列。 |
| 单样本流程 | 可从显式 mzML 路径运行，并写出非空 feature table。 |
| 多样本 alignment | 尚未实现。 |
| 完整 Python parity | 尚未声明；下一步是和 Python 输出做对照。 |

详细科学路线和缺失模块见 [`docs/mainline.md`](docs/mainline.md)。

### 快速开始

从仓库根目录运行原型测试：

```sh
Rscript r-prototype/tests/smoke_test.R
Rscript r-prototype/tests/mzml_reader_test.R
Rscript r-prototype/tests/mass_tracks_test.R
Rscript r-prototype/tests/peak_detection_test.R
Rscript r-prototype/tests/single_sample_process_test.R
```

用明确指定的 mzML 文件运行单样本原型：

```sh
Rscript r-prototype/scripts/process_single_sample.R "path/to/file.mzML" "r-output" 20
```

第三个参数是可选的 scan 数量限制，适合快速测试。只有在准备处理完整文件时
才省略它。

预览各个处理阶段：

```sh
Rscript r-prototype/scripts/preview_mzml.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_mass_tracks.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_peaks.R "path/to/file.mzML"
```

预览和处理脚本不应依赖隐藏的默认 mzML 文件；输入文件路径应由用户明确提供。

### 仓库结构

```text
docs/mainline.md      科学主线和移植检查点
r-prototype/R/        当前 R 原型实现
r-prototype/scripts/  需要显式输入的预览和处理脚本
r-prototype/tests/    原型测试
```

### 开发优先级

1. 将当前单样本 R 输出与 Python asari 的 `preferred_Feature_table.tsv`
   做对照。
2. 移植 alignment 所需的核心 mass matching 和 clustering helper。
3. 以 R 原生方式实现多样本 `MassGrid` / `CompositeMap` 行为。

在主线能稳定产出可比较的 feature table 之前，dashboard、annotation、
notebooks、GC annotation 等旁支功能暂不进入范围。

## English

This repository is a working R port of
[asari](https://github.com/shuzhao-li/asari), a Python program for
high-resolution LC-MS metabolomics data processing.

The goal is to reproduce asari's scientific processing mainline in R without
mechanically copying the Python source tree. The current work is organized
around a focused prototype in [`r-prototype/`](r-prototype/) and a mainline map
in [`docs/mainline.md`](docs/mainline.md).

### Scope

The R implementation is being built around the core asari processing path:

1. Read explicit mzML inputs.
2. Register samples and run parameters.
3. Extract mass tracks / extracted ion chromatograms from MS1 scans.
4. Detect chromatographic peaks on mass tracks.
5. Write asari-style feature table output.
6. Extend from single-sample output to multi-sample alignment.

This repository is not trying to clone the Python package layout file by file.
When an R-native structure is clearer, the R port should use it while keeping
the scientific behavior aligned with Python asari.

### Current Status

The active code lives in [`r-prototype/`](r-prototype/).

| Area | Status |
| --- | --- |
| mzML reading | Initial reader implemented with `mzR`. |
| Parameter defaults | Prototype defaults translated from Python asari. |
| Mass tracks / EICs | Rough single-sample prototype implemented. |
| Peak detection | Rough prototype emits feature-table-style columns. |
| Single-sample process | Can run from an explicit mzML path and write a non-empty feature table. |
| Multi-sample alignment | Not implemented yet. |
| Full Python parity | Not claimed yet; comparison against Python output is the next checkpoint. |

For the detailed scientific route and missing pieces, see
[`docs/mainline.md`](docs/mainline.md).

### Quick Start

Run the prototype tests from the repository root:

```sh
Rscript r-prototype/tests/smoke_test.R
Rscript r-prototype/tests/mzml_reader_test.R
Rscript r-prototype/tests/mass_tracks_test.R
Rscript r-prototype/tests/peak_detection_test.R
Rscript r-prototype/tests/single_sample_process_test.R
```

Run a single-sample prototype with an explicit mzML file:

```sh
Rscript r-prototype/scripts/process_single_sample.R "path/to/file.mzML" "r-output" 20
```

The optional third argument limits the number of scans for quick testing. Omit
it only when you are ready to process the full file.

Preview individual stages with explicit mzML input:

```sh
Rscript r-prototype/scripts/preview_mzml.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_mass_tracks.R "path/to/file.mzML"
Rscript r-prototype/scripts/preview_peaks.R "path/to/file.mzML"
```

No preview or processing script should rely on a hidden default mzML file.

### Repository Layout

```text
docs/mainline.md      Scientific mainline and porting checkpoints
r-prototype/R/        Current R prototype implementation
r-prototype/scripts/  Explicit-input preview and processing scripts
r-prototype/tests/    Prototype tests
```

### Development Priorities

1. Compare the current single-sample R output against Python asari's
   `preferred_Feature_table.tsv`.
2. Port the core mass matching and clustering helpers needed for alignment.
3. Implement multi-sample `MassGrid` / `CompositeMap` behavior in an R-native
   shape.

Until the mainline can produce comparable feature tables, dashboard,
annotation, notebooks, and GC annotation work should stay out of scope.
