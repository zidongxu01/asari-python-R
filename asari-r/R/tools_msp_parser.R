# 对应 Python asari/tools/msp_parser.py：解析MSP和MGF质谱库文本。

# Python原版的受控字段同义词表。
MSP_dict <- list(MW = "ExactMass")

# 对应 msp_standarize：只在目标字段不存在时复制旧字段。
msp_standarize <- function(LL, MSP_dict = get("MSP_dict", inherits = TRUE)) {
  lapply(LL, function(entry) {
    for (key in names(MSP_dict)) {
      target <- MSP_dict[[key]]
      if (!is.null(entry[[key]]) && is.null(entry[[target]])) entry[[target]] <- entry[[key]]
    }
    entry
  })
}

# 对应 parse_msp_to_listdict：按空行切分小型MSP文件。
parse_msp_to_listdict <- function(file, field_separator = ": ", return_peaks = TRUE) {
  text <- paste(readLines(file, warn = FALSE), collapse = "\n")
  text <- sub("[\r\n]+$", "", text)
  blocks <- strsplit(text, "\n[[:space:]]*\n", perl = TRUE)[[1L]]
  lapply(blocks, function(block) {
    result <- list()
    peak_lines <- character()
    for (line in strsplit(block, "\n", fixed = TRUE)[[1L]]) {
      position <- regexpr(field_separator, line, fixed = TRUE)[[1L]]
      if (position > 0L) {
        pieces <- strsplit(line, field_separator, fixed = TRUE)[[1L]]
        # Python line.split后只取前两个元素。
        result[[pieces[[1L]]]] <- pieces[[2L]]
      } else if (nzchar(trimws(line))) {
        peak_lines <- c(peak_lines, line)
      }
    }
    if (isTRUE(return_peaks)) {
      peaks <- list()
      for (line in peak_lines) {
        for (token in strsplit(line, ";", fixed = TRUE)[[1L]]) {
          fields <- strsplit(trimws(token), "[[:space:]]+")[[1L]]
          if (length(fields) == 2L) {
            values <- suppressWarnings(as.numeric(fields))
            if (!anyNA(values)) peaks[[length(peaks) + 1L]] <- values
          }
        }
      }
      result$peaks <- peaks
    }
    result
  })
}

# 对应 parse_mgf 内部 def parse_pepmass。
.parse_mgf_pepmass <- function(value) {
  parts <- strsplit(trimws(value), "[[:space:]]+")[[1L]]
  if (length(parts) == 1L) return(list(as.numeric(parts[[1L]]), NULL))
  if (length(parts) >= 2L) return(list(as.numeric(parts[[1L]]), as.numeric(parts[[2L]])))
  list(NULL, NULL)
}

# 对应 parse_mgf 内部 def parse_charge。
.parse_mgf_charge <- function(value) {
  matched <- regmatches(value, regexpr("[0-9]+", value, perl = TRUE))
  if (length(matched) == 0L || !nzchar(matched)) NULL else as.integer(matched)
}

# 对应 parse_mgf：读取BEGIN IONS/END IONS块及参数和峰数组。
parse_mgf <- function(file_path) {
  spectra <- list()
  current <- NULL
  for (raw_line in readLines(file_path, warn = FALSE)) {
    line <- trimws(raw_line)
    if (!nzchar(line) || startsWith(line, "#")) next
    if (toupper(line) == "BEGIN IONS") {
      current <- list(params = list(), mzs = numeric(), intensities = numeric())
    } else if (toupper(line) == "END IONS") {
      if (!is.null(current)) spectra[[length(spectra) + 1L]] <- current
      current <- NULL
    } else if (!is.null(current)) {
      if (grepl("=", line, fixed = TRUE)) {
        position <- regexpr("=", line, fixed = TRUE)[[1L]]
        key <- toupper(trimws(substr(line, 1L, position - 1L)))
        value <- trimws(substr(line, position + 1L, nchar(line)))
        current$params[[key]] <- if (key == "PEPMASS") {
          .parse_mgf_pepmass(value)
        } else if (key == "CHARGE") {
          .parse_mgf_charge(value)
        } else value
      } else {
        parts <- strsplit(line, "[[:space:]]+")[[1L]]
        if (length(parts) >= 2L) {
          values <- suppressWarnings(as.numeric(parts[1:2]))
          if (!anyNA(values)) {
            current$mzs <- c(current$mzs, values[[1L]])
            current$intensities <- c(current$intensities, values[[2L]])
          }
        }
      }
    }
  }
  spectra
}

# 对应 parse_peak_token：解析冒号分隔峰，保留剩余注释文本。
parse_peak_token <- function(token) {
  parts <- strsplit(token, ":", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) return(NULL)
  values <- suppressWarnings(as.numeric(parts[1:2]))
  if (anyNA(values)) return(NULL)
  list(
    mz = values[[1L]],
    intensity = values[[2L]],
    annotation = if (length(parts) > 2L) paste(parts[-(1:2)], collapse = ":") else NULL
  )
}

# 判断Python正则所接受的普通十进制数字，不包含科学记数法。
.msp_plain_number <- function(value) {
  grepl("^-?[0-9]+(\\.[0-9]+)?$", value, perl = TRUE)
}

# 对应 parse_peak_line：兼容空格、逗号和冒号三种峰格式。
parse_peak_line <- function(line) {
  tokens <- strsplit(gsub(",", " ", line, fixed = TRUE), "[[:space:]]+")[[1L]]
  tokens <- tokens[nzchar(tokens)]
  peaks <- list()
  ii <- 1L
  while (ii <= length(tokens)) {
    parsed_pair <- FALSE
    if (ii + 1L <= length(tokens)) {
      values <- suppressWarnings(as.numeric(tokens[ii:(ii + 1L)]))
      if (!anyNA(values)) {
        annotation <- NULL
        if (ii + 2L <= length(tokens) && !.msp_plain_number(tokens[[ii + 2L]])) {
          annotation <- tokens[[ii + 2L]]
          ii <- ii + 3L
        } else {
          ii <- ii + 2L
        }
        peaks[[length(peaks) + 1L]] <- list(
          mz = values[[1L]], intensity = values[[2L]], annotation = annotation
        )
        parsed_pair <- TRUE
      }
    }
    if (!parsed_pair) {
      peak <- parse_peak_token(tokens[[ii]])
      if (!is.null(peak)) peaks[[length(peaks) + 1L]] <- peak
      ii <- ii + 1L
    }
  }
  peaks
}

# 对应 parse_msp：按Name和Num Peaks组织更宽松的MSP记录。
parse_msp <- function(filepath) {
  spectra <- list()
  current <- NULL
  expected_peaks <- NULL
  peak_count <- 0L
  for (raw_line in readLines(filepath, warn = FALSE, encoding = "UTF-8")) {
    line <- trimws(raw_line)
    if (!nzchar(line)) next
    if (grepl(":", line, fixed = TRUE) && !grepl("^[0-9]", line)) {
      position <- regexpr(":", line, fixed = TRUE)[[1L]]
      key <- trimws(substr(line, 1L, position - 1L))
      value <- trimws(substr(line, position + 1L, nchar(line)))
      if (tolower(key) == "name") {
        if (!is.null(current)) spectra[[length(spectra) + 1L]] <- current
        current <- list(Name = value, peaks = list())
        expected_peaks <- NULL
        peak_count <- 0L
        next
      }
      if (is.null(current)) current <- list(peaks = list())
      if (tolower(key) == "num peaks") {
        expected_peaks <- suppressWarnings(as.integer(value))
        if (is.na(expected_peaks)) expected_peaks <- NULL
      } else {
        current[[key]] <- value
      }
      next
    }
    if (is.null(current)) next
    peaks <- parse_peak_line(line)
    current$peaks <- c(current$peaks, peaks)
    peak_count <- peak_count + length(peaks)
    if (!is.null(expected_peaks) && expected_peaks != 0L && peak_count >= expected_peaks) {
      spectra[[length(spectra) + 1L]] <- current
      current <- NULL
      expected_peaks <- NULL
      peak_count <- 0L
    }
  }
  if (!is.null(current)) spectra[[length(spectra) + 1L]] <- current
  spectra
}
