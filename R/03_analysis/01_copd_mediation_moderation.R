# ------------------------------------------------------------------------------
# File:     R/03_analysis/01_copd_mediation_moderation.R
# Purpose:  分析肺病患者 ADL、自评健康、用药行为与随访抑郁症状的关联。
# Inputs:   data/derived/02_copd_analysis.rds
# Outputs:  output/tables/*.csv, output/tables/*.docx
#           output/figures/subgroup_forest_copd.png, .pdf
# Log:      logs/03_analysis_copd_mediation_moderation.log
# ------------------------------------------------------------------------------

if (getRversion() < "4.3.0") stop("Requires R >= 4.3.0; you have ", R.version.string)

source("R/_utils/paths.R")
source("R/_utils/logging.R")
source("R/_utils/theme_journal.R")

start_log("03_analysis_copd_mediation_moderation")
on.exit(stop_log(), add = TRUE)
set.seed(20260605)

required_pkgs <- c(
  "dplyr", "readr", "tibble", "ggplot2", "broom", "gtsummary",
  "flextable", "officer", "lavaan", "EValue"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run scripts/setup_r.R first.")
}

input_path <- proj_path("data", "derived", "02_copd_analysis.rds")
table_dir <- proj_path("output", "tables")
figure_dir <- proj_path("output", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Analysis data not found: ", input_path,
       ". Run R/02_construct/01_construct_copd_variables.R first.")
}

dat <- readRDS(input_path)
cat("分析样本行数：", nrow(dat), "\n", sep = "")

notes <- tibble::tibble(item = character(), note = character())
add_note <- function(item, note) {
  notes <<- tibble::add_row(notes, item = item, note = note)
  cat("[说明] ", item, ": ", note, "\n", sep = "")
}

has_variation <- function(data, var) {
  vals <- stats::na.omit(data[[var]])
  length(unique(vals)) >= 2
}

complete_model_data <- function(data, vars) {
  data |>
    dplyr::select(dplyr::all_of(vars)) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))
}

write_note_table <- function() {
  readr::write_csv(notes, file.path(table_dir, "analysis_notes.csv"))
}

# ------------------------------------------------------------------------------
# 1. Table 1：按肺病用药分组
# ------------------------------------------------------------------------------
cat("\n生成 Table 1：按肺病用药分组\n")
table1_vars <- c(
  "rxlung_group", "r3agey", "age_group", "gender_group", "education_group",
  "rural_group", "ever_smoked", "ever_drank", "hibp", "diabetes",
  "heart_disease", "comorbidity_count", "r3cesd10", "r3adla_c", "r3shlt",
  "r4cesd10", "cesd10_ge10"
)

table1_data <- dat |>
  dplyr::select(dplyr::all_of(table1_vars)) |>
  dplyr::filter(!is.na(.data$rxlung_group))

if (nrow(table1_data) == 0 || !has_variation(table1_data, "rxlung_group")) {
  add_note("Table 1", "r3rxlung/rxlung_group 缺少可分组变异，未生成分组描述表。")
} else {
  tbl1 <- table1_data |>
    gtsummary::tbl_summary(
      by = "rxlung_group",
      statistic = list(
        gtsummary::all_continuous() ~ "{mean} ({sd})",
        gtsummary::all_categorical() ~ "{n} ({p}%)"
      ),
      digits = gtsummary::all_continuous() ~ 2,
      label = list(
        r3agey ~ "年龄（岁）",
        age_group ~ "年龄组",
        gender_group ~ "性别",
        education_group ~ "教育程度",
        rural_group ~ "居住地",
        ever_smoked ~ "曾吸烟",
        ever_drank ~ "曾饮酒",
        hibp ~ "高血压",
        diabetes ~ "糖尿病",
        heart_disease ~ "心脏病",
        comorbidity_count ~ "共病计数",
        r3cesd10 ~ "基线 CES-D 总分",
        r3adla_c ~ "ADL 受限数",
        r3shlt ~ "自评健康",
        r4cesd10 ~ "随访 CES-D 总分",
        cesd10_ge10 ~ "随访抑郁症状（CES-D >= 10）"
      ),
      missing = "ifany"
    ) |>
    gtsummary::add_overall() |>
    gtsummary::add_p() |>
    gtsummary::bold_labels()

  flextable::save_as_docx(
    gtsummary::as_flex_table(tbl1),
    path = file.path(table_dir, "table1_by_rxlung.docx")
  )
  cat("Table 1 已保存：output/tables/table1_by_rxlung.docx\n")
}

# ------------------------------------------------------------------------------
# 2. 逐步调整 Logistic 回归
# ------------------------------------------------------------------------------
cat("\n估计逐步调整 Logistic 回归\n")
base_covars <- c("r3agey", "female", "education_group")
health_covars <- c("rural", "ever_smoked", "ever_drank", "comorbidity_count")
baseline_covars <- c("r3cesd10")

# 严格肺病队列中 r3lunge 通常恒为 1，主效应不可估计；此处自动改用 rxlung 做队列内关联模型。
logit_exposure <- "r3lunge"
if (!has_variation(dat, logit_exposure)) {
  add_note("Logistic 主暴露", "样本已限制为 r3lunge = 1，肺病诊断无变异；逐步 Logistic 改为估计 rxlung 与随访抑郁的队列内关联。")
  logit_exposure <- "rxlung"
}

fit_logit <- function(model_name, covars, outcome = "cesd10_ge10", exposure = logit_exposure) {
  vars <- unique(c(outcome, exposure, covars))
  model_data <- complete_model_data(dat, vars)
  if (nrow(model_data) == 0 || !has_variation(model_data, outcome) || !has_variation(model_data, exposure)) {
    add_note(model_name, "结局或暴露缺少变异，模型未估计。")
    return(NULL)
  }
  fml <- stats::as.formula(paste(outcome, "~", paste(c(exposure, covars), collapse = " + ")))
  fit <- stats::glm(fml, data = model_data, family = stats::binomial())
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(model = model_name, n = nrow(model_data), .before = 1)
}

logit_results <- dplyr::bind_rows(
  fit_logit("Model 1: 未调整", character()),
  fit_logit("Model 2: 人口学调整", base_covars),
  fit_logit("Model 3: 行为与共病调整", c(base_covars, health_covars)),
  fit_logit("Model 4: 加基线抑郁调整", c(base_covars, health_covars, baseline_covars))
)

if (nrow(logit_results) > 0) {
  readr::write_csv(logit_results, file.path(table_dir, "logistic_stepwise_models.csv"))
  cat("逐步 Logistic 结果已保存：output/tables/logistic_stepwise_models.csv\n")
}

# ------------------------------------------------------------------------------
# 3. lavaan 并行中介分析：仅当肺病暴露有变异时估计
# ------------------------------------------------------------------------------
cat("\n估计 lavaan 并行中介模型\n")
mediation_vars <- c("r3lunge", "r3adla_c", "r3shlt", "r4cesd10", base_covars, health_covars, baseline_covars)
mediation_data <- complete_model_data(dat, mediation_vars)

if (!has_variation(mediation_data, "r3lunge")) {
  add_note("lavaan 并行中介", "样本中 r3lunge 无变异，无法估计肺病诊断经 ADL/自评健康到抑郁的中介效应。")
  readr::write_csv(tibble::tibble(result = "not_estimated", reason = "r3lunge has no variation"),
                   file.path(table_dir, "lavaan_parallel_mediation.csv"))
} else {
  med_model <- "
    r3adla_c ~ a1*r3lunge + r3agey + female + education_group + rural + ever_smoked + ever_drank + comorbidity_count + r3cesd10
    r3shlt   ~ a2*r3lunge + r3agey + female + education_group + rural + ever_smoked + ever_drank + comorbidity_count + r3cesd10
    r4cesd10 ~ cprime*r3lunge + b1*r3adla_c + b2*r3shlt + r3agey + female + education_group + rural + ever_smoked + ever_drank + comorbidity_count + r3cesd10
    ind_adl := a1*b1
    ind_shlt := a2*b2
    ind_total := ind_adl + ind_shlt
    total := cprime + ind_total
  "
  med_fit <- lavaan::sem(
    med_model,
    data = mediation_data,
    se = "bootstrap",
    bootstrap = 5000,
    fixed.x = FALSE
  )
  med_results <- lavaan::parameterEstimates(med_fit, boot.ci.type = "perc", standardized = TRUE) |>
    dplyr::as_tibble()
  readr::write_csv(med_results, file.path(table_dir, "lavaan_parallel_mediation.csv"))
  cat("lavaan 中介结果已保存：output/tables/lavaan_parallel_mediation.csv\n")
}

# ------------------------------------------------------------------------------
# 4. 用药调节：肺病主效应交互若不可估计，则估计中介-用药交互
# ------------------------------------------------------------------------------
cat("\n估计用药调节模型\n")
interaction_vars <- c("cesd10_ge10", "rxlung", "r3adla_c", "r3shlt", base_covars, health_covars, baseline_covars)
interaction_data <- complete_model_data(dat, interaction_vars)

if (nrow(interaction_data) == 0 || !has_variation(interaction_data, "rxlung")) {
  add_note("用药调节", "rxlung 缺少变异或完整样本为空，调节模型未估计。")
} else {
  interaction_formula <- stats::as.formula(
    paste(
      "cesd10_ge10 ~ r3adla_c * rxlung + r3shlt * rxlung +",
      paste(c(base_covars, health_covars, baseline_covars), collapse = " + ")
    )
  )
  interaction_fit <- stats::glm(interaction_formula, data = interaction_data, family = stats::binomial())
  interaction_results <- broom::tidy(interaction_fit, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(n = nrow(interaction_data), .before = 1)
  readr::write_csv(interaction_results, file.path(table_dir, "medicator_by_rxlung_interactions.csv"))
  cat("用药调节结果已保存：output/tables/medicator_by_rxlung_interactions.csv\n")
}

# ------------------------------------------------------------------------------
# 5. 亚组森林图
# ------------------------------------------------------------------------------
cat("\n生成亚组森林图\n")
subgroup_specs <- list(
  性别 = "gender_group",
  年龄 = "age_group",
  城乡 = "rural_group",
  用药 = "rxlung_group"
)

estimate_subgroup <- function(group_name, group_var, exposure = logit_exposure) {
  levels_now <- stats::na.omit(unique(dat[[group_var]]))
  dplyr::bind_rows(lapply(levels_now, function(level_value) {
    sub_data <- dat |>
      dplyr::filter(.data[[group_var]] == level_value)
    vars <- unique(c("cesd10_ge10", exposure, base_covars, health_covars, baseline_covars))
    model_data <- complete_model_data(sub_data, vars)
    if (nrow(model_data) == 0 || !has_variation(model_data, "cesd10_ge10") || !has_variation(model_data, exposure)) {
      return(NULL)
    }
    fml <- stats::as.formula(paste("cesd10_ge10 ~", paste(c(exposure, base_covars, health_covars, baseline_covars), collapse = " + ")))
    fit <- stats::glm(fml, data = model_data, family = stats::binomial())
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) |>
      dplyr::filter(.data$term == exposure) |>
      dplyr::transmute(
        subgroup = group_name,
        level = as.character(level_value),
        exposure = exposure,
        n = nrow(model_data),
        estimate = .data$estimate,
        conf.low = .data$conf.low,
        conf.high = .data$conf.high,
        p.value = .data$p.value
      )
  }))
}

forest_results <- dplyr::bind_rows(lapply(names(subgroup_specs), function(nm) {
  estimate_subgroup(nm, subgroup_specs[[nm]])
}))

if (nrow(forest_results) == 0) {
  add_note("亚组森林图", "所有亚组中结局或暴露缺少变异，未能估计森林图点估计。")
  readr::write_csv(tibble::tibble(result = "not_estimated"), file.path(table_dir, "subgroup_forest_results.csv"))
} else {
  readr::write_csv(forest_results, file.path(table_dir, "subgroup_forest_results.csv"))
  forest_plot <- forest_results |>
    dplyr::mutate(label = paste0(.data$subgroup, ": ", .data$level)) |>
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data$label, .data$estimate), y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
                             colour = pal_journal[["navy"]]) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "OR（对数刻度）", title = "亚组分析森林图") +
    theme_journal()
  ggplot2::ggsave(file.path(figure_dir, "subgroup_forest_copd.png"), forest_plot, width = 7, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(figure_dir, "subgroup_forest_copd.pdf"), forest_plot, width = 7, height = 5)
  cat("亚组森林图已保存：output/figures/subgroup_forest_copd.png 和 .pdf\n")
}

# ------------------------------------------------------------------------------
# 6. 敏感性分析：CES-D >= 12、新发抑郁、E-value
# ------------------------------------------------------------------------------
cat("\n估计敏感性分析\n")
sens_ge12 <- fit_logit("Sensitivity: CES-D >= 12", c(base_covars, health_covars, baseline_covars), outcome = "cesd10_ge12")
sens_incident <- fit_logit("Sensitivity: 新发抑郁", c(base_covars, health_covars), outcome = "incident_depressed")
sensitivity_results <- dplyr::bind_rows(sens_ge12, sens_incident)
if (nrow(sensitivity_results) > 0) {
  readr::write_csv(sensitivity_results, file.path(table_dir, "sensitivity_models.csv"))
}

main_effect <- logit_results |>
  dplyr::filter(.data$model == "Model 4: 加基线抑郁调整", .data$term == logit_exposure) |>
  dplyr::slice_head(n = 1)

if (nrow(main_effect) == 1 && all(c("estimate", "conf.low", "conf.high") %in% names(main_effect))) {
  evalue <- EValue::evalues.OR(
    est = main_effect$estimate,
    lo = main_effect$conf.low,
    hi = main_effect$conf.high,
    true = 1
  )
  evalue_tbl <- tibble::as_tibble(as.data.frame(evalue), rownames = "parameter")
  readr::write_csv(evalue_tbl, file.path(table_dir, "evalue_model4.csv"))
  cat("E-value 已保存：output/tables/evalue_model4.csv\n")
} else {
  add_note("E-value", "Model 4 主效应不可用，未计算 E-value。")
}

write_note_table()
cat("\n分析脚本完成。表格目录：output/tables/；图形目录：output/figures/\n")
