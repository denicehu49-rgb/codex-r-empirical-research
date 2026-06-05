# ------------------------------------------------------------------------------
# File:     R/01_clean/01_clean_copd_charls.R
# Purpose:  清洗 Harmonized CHARLS Version D，生成肺病患者纵向分析样本。
# Inputs:   data/raw/H_CHARLS_D_Data.dta
# Outputs:  data/derived/01_copd_clean.rds
#           output/tables/strobe_flow_counts.csv
# Log:      logs/01_clean_copd_charls.log
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)

source("R/_utils/paths.R")
source("R/_utils/logging.R")

start_log("01_clean_copd_charls")
on.exit(stop_log(), add = TRUE)

required_pkgs <- c("dplyr", "haven", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run scripts/setup_r.R first.")
}

raw_path <- proj_path("data", "raw", "H_CHARLS_D_Data.dta")
derived_dir <- proj_path("data", "derived")
table_dir <- proj_path("output", "tables")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

vars_needed <- c(
  "r3lunge", "r4cesd10", "r3adla_c", "r3shlt", "r3rxlung",
  "r3agey", "ragender", "raeduc_c", "h3rural", "r3smokev",
  "r3drinkev", "r3hibpe", "r3diabe", "r3hearte", "r3cesd10",
  "r3psyche", "r3memrye"
)

cat("读取原始 Stata 数据：", raw_path, "\n", sep = "")
if (!file.exists(raw_path)) {
  stop("Raw data file not found: ", raw_path)
}

raw_data <- haven::read_dta(raw_path)
cat("原始数据行数：", nrow(raw_data), "\n", sep = "")
cat("原始数据列数：", ncol(raw_data), "\n", sep = "")

missing_vars <- setdiff(vars_needed, names(raw_data))
if (length(missing_vars) > 0) {
  stop("Missing required variables in raw data: ", paste(missing_vars, collapse = ", "))
}

# Harmonized CHARLS 常见负值/标签缺失码；本研究变量均不应有负的有效值。
clean_special_missing <- function(x) {
  x <- haven::zap_missing(x)
  x <- haven::zap_labels(x)
  if (is.numeric(x)) x[x < 0] <- NA_real_
  x
}

analysis_vars <- raw_data |>
  dplyr::select(dplyr::all_of(vars_needed)) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), clean_special_missing))

flow <- tibble::tibble(
  step = character(),
  n_remaining = integer(),
  n_excluded = integer(),
  reason = character()
)

add_flow <- function(data, reason, previous_n = NA_integer_) {
  current_n <- nrow(data)
  excluded_n <- if (is.na(previous_n)) 0L else previous_n - current_n
  tibble::add_row(
    flow,
    step = sprintf("%02d", nrow(flow) + 1L),
    n_remaining = current_n,
    n_excluded = excluded_n,
    reason = reason
  )
}

n_previous <- nrow(analysis_vars)
flow <- add_flow(analysis_vars, "原始样本；保留研究所需变量")

analysis_vars <- analysis_vars |>
  dplyr::filter(!is.na(.data$r3lunge), .data$r3lunge == 1)
flow <- add_flow(analysis_vars, "纳入 Wave 3 肺部疾病诊断者（r3lunge = 1）", n_previous)
n_previous <- nrow(analysis_vars)

analysis_vars <- analysis_vars |>
  dplyr::filter(!is.na(.data$r3agey), .data$r3agey >= 45)
flow <- add_flow(analysis_vars, "纳入基线年龄 >= 45 岁", n_previous)
n_previous <- nrow(analysis_vars)

analysis_vars <- analysis_vars |>
  dplyr::filter(!is.na(.data$r4cesd10))
flow <- add_flow(analysis_vars, "纳入 Wave 4 CES-D 总分非缺失者", n_previous)
n_previous <- nrow(analysis_vars)

analysis_vars <- analysis_vars |>
  dplyr::filter(is.na(.data$r3psyche) | .data$r3psyche != 1) |>
  dplyr::filter(is.na(.data$r3memrye) | .data$r3memrye != 1)
flow <- add_flow(analysis_vars, "排除精神疾病或记忆相关疾病报告者（r3psyche = 1 或 r3memrye = 1）", n_previous)

attr(analysis_vars, "strobe_flow") <- flow

readr::write_csv(flow, proj_path("output", "tables", "strobe_flow_counts.csv"))
saveRDS(analysis_vars, proj_path("data", "derived", "01_copd_clean.rds"))

cat("\nSTROBE 纳排人数记录：\n")
print(flow)
cat("\n清洗样本已保存：data/derived/01_copd_clean.rds\n")
cat("STROBE 流程表已保存：output/tables/strobe_flow_counts.csv\n")
