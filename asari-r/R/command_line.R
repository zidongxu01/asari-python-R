# 对应 Python asari/command_line.py：该文件只从main模块引入同一入口。

# 使用运行时查找，避免R包按文件名加载时main.R尚未读取。
command_line_main <- function(...) {
  entry <- get("main", mode = "function", inherits = TRUE)
  entry(...)
}
