# ------------------------------------------------------------------------------
# File:     R/02_construct/01_construct_copd_variables.R
# Purpose:  构造 CHARLS 肺病患者队列分析变量。
# Inputs:   data/derived/01_copd_clean.rds
# Outputs:  data/derived/02_copd_analysis.rds
#           output/tables/analysis_sample_missingness.csv
# Log:      logs/02_construct_copd_variables.log
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)

source("R/_utils/paths.R")
source("R/_utils/logging.R")

start_log("02_construct_copd_variables")
on.exit(stop_log(), add = TRUE)

required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run scripts/setup_r.R first.")
}

input_path <- proj_path("data", "derived", "01_copd_clean.rds")
output_path <- proj_path("data", "derived", "02_copd_analysis.rds")
table_dir <- proj_path("output", "tables")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Clean data not found: ", input_path, ". Run R/01_clean/01_clean_copd_charls.R first.")
}

cat("读取清洗样本：", input_path, "\n", sep = "")
clean_data <- readRDS(input_path)
cat("清洗样本行数：", nrow(clean_data), "\n", sep = "")

as_binary01 <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_integer_,
    x == 1 ~ 1L,
    x == 0 ~ 0L,
    TRUE ~ NA_integer_
  )
}

analysis_data <- clean_data |>
  dplyr::mutate(
    # 结局变量：CES-D 常用切点 10 分；敏感性分析使用 12 分。
    cesd10_ge10 = dplyr::case_when(
      is.na(.data$r4cesd10) ~ NA_integer_,
      .data$r4cesd10 >= 10 ~ 1L,
      .data$r4cesd10 < 10 ~ 0L
    ),
    cesd10_ge12 = dplyr::case_when(
      is.na(.data$r4cesd10) ~ NA_integer_,
      .data$r4cesd10 >= 12 ~ 1L,
      .data$r4cesd10 < 12 ~ 0L
    ),
    baseline_depressed = dplyr::case_when(
      is.na(.data$r3cesd10) ~ NA_integer_,
      .data$r3cesd10 >= 10 ~ 1L,
      .data$r3cesd10 < 10 ~ 0L
    ),
    incident_depressed = dplyr::case_when(
      .data$baseline_depressed == 0L & .data$cesd10_ge10 == 1L ~ 1L,
      .data$baseline_depressed == 0L & .data$cesd10_ge10 == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # 人口学与行为变量重编码。
    age_group = dplyr::case_when(
      is.na(.data$r3agey) ~ NA_character_,
      .data$r3agey < 60 ~ "45-59",
      .data$r3agey < 75 ~ "60-74",
      TRUE ~ ">=75"
    ),
    female = dplyr::case_when(
      is.na(.data$ragender) ~ NA_integer_,
      .data$ragender == 2 ~ 1L,
      .data$ragender == 1 ~ 0L,
      TRUE ~ NA_integer_
    ),
    gender_group = dplyr::case_when(
      .data$female == 1L ~ "女性",
      .data$female == 0L ~ "男性",
      TRUE ~ NA_character_
    ),
    education_group = dplyr::case_when(
      is.na(.data$raeduc_c) ~ NA_character_,
      .data$raeduc_c <= 1 ~ "小学以下",
      .data$raeduc_c == 2 ~ "小学",
      .data$raeduc_c == 3 ~ "初中",
      .data$raeduc_c >= 4 ~ "高中及以上",
      TRUE ~ NA_character_
    ),
    rural = as_binary01(.data$h3rural),
    rural_group = dplyr::case_when(
      .data$rural == 1L ~ "农村",
      .data$rural == 0L ~ "城镇",
      TRUE ~ NA_character_
    ),
    ever_smoked = as_binary01(.data$r3smokev),
    ever_drank = as_binary01(.data$r3drinkev),
    rxlung = as_binary01(.data$r3rxlung),
    rxlung_group = dplyr::case_when(
      .data$rxlung == 1L ~ "用药",
      .data$rxlung == 0L ~ "未用药",
      TRUE ~ NA_character_
    ),

    # 共病计数：高血压、糖尿病、心脏病。
    hibp = as_binary01(.data$r3hibpe),
    diabetes = as_binary01(.data$r3diabe),
    heart_disease = as_binary01(.data$r3hearte),
    comorbidity_count = rowSums(dplyr::pick(.data$hibp, .data$diabetes, .data$heart_disease), na.rm = TRUE),
    any_comorbidity = dplyr::if_else(.data$comorbidity_count > 0, 1L, 0L),

    # 因子顺序用于表格与作图。
    age_group = factor(.data$age_group, levels = c("45-59", "60-74", ">=75")),
    gender_group = factor(.data$gender_group, levels = c("男性", "女性")),
    education_group = factor(.data$education_group,
                             levels = c("小学以下", "小学", "初中", "高中及以上")),
    rural_group = factor(.data$rural_group, levels = c("城镇", "农村")),
    rxlung_group = factor(.data$rxlung_group, levels = c("未用药", "用药"))
  )

analysis_vars <- c(
  "r3lunge", "r4cesd10", "cesd10_ge10", "cesd10_ge12", "incident_depressed",
  "r3adla_c", "r3shlt", "r3rxlung", "rxlung", "rxlung_group",
  "r3agey", "age_group", "ragender", "female", "gender_group",
  "raeduc_c", "education_group", "h3rural", "rural", "rural_group",
  "r3smokev", "ever_smoked", "r3drinkev", "ever_drank",
  "r3hibpe", "hibp", "r3diabe", "diabetes", "r3hearte", "heart_disease",
  "comorbidity_count", "any_comorbidity", "r3cesd10", "baseline_depressed"
)

missingness <- tibble::tibble(
  variable = analysis_vars,
  n_missing = vapply(analysis_vars, function(v) sum(is.na(analysis_data[[v]])), integer(1)),
  pct_missing = round(100 * .data$n_missing / nrow(analysis_data), 2)
)

readr::write_csv(missingness, proj_path("output", "tables", "analysis_sample_missingness.csv"))
saveRDS(analysis_data, output_path)

cat("\n分析变量缺失情况：\n")
print(missingness)
cat("\n变量构造数据已保存：data/derived/02_copd_analysis.rds\n")
cat("缺失情况表已保存：output/tables/analysis_sample_missingness.csv\n")
