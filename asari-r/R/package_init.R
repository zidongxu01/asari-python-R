# 对应 Python asari/__init__.py：声明包版本并由R包加载机制统一加载各模块。

# Python 对外提供 __version__ = "1.17.0"；反引号使R保留同名符号。
`__version__` <- "1.17.0"

# 提供符合R命名习惯的同值别名。
asari_version <- `__version__`

# Python __init__ 的 import 会加载主模块；R包在载入时会自动读取R/下所有文件。
ASARI_CORE_MODULES <- c(
  "analyze", "annotate", "chromatograms", "command_line", "constructors",
  "default_parameters", "experiment", "main", "mass_functions", "peaks",
  "samples", "workflow"
)
