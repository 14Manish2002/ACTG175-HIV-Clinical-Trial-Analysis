# =============================================================================
# ACTG 175 Survival Analysis — Kaplan-Meier Estimation
# Script:  03_kaplan_meier.R
# Purpose: KM survival curves — overall, by treatment, and key subgroups
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 03: Kaplan-Meier Analysis ===\n")

# Helper: save survminer plot
save_km_plot <- function(plot_obj, filename, width = 10, height = 8) {
  ggsave(
    filename = file.path(dir_figures, filename),
    plot     = plot_obj,
    width    = width,
    height   = height,
    dpi      = 300,
    bg       = "white"
  )
  cat("Saved:", filename, "\n")
}

# --------------------------------------------------------------------------- #
#  1. Overall survival curve
# --------------------------------------------------------------------------- #
km_overall <- survfit(Surv(time, label) ~ 1, data = actg)

cat("\n--- Overall Survival Summary ---\n")
print(summary(km_overall, times = c(180, 365, 540, 730, 900))$table)

p_km_overall <- ggsurvplot(
  km_overall,
  data           = actg,
  conf.int       = TRUE,
  risk.table     = TRUE,
  risk.table.height = 0.28,
  surv.median.line  = "hv",
  palette        = "#2166AC",
  title          = "Overall Kaplan-Meier Survival Curve",
  subtitle       = "ACTG 175 Trial — Time to AIDS Diagnosis or Death",
  xlab           = "Time (Days)",
  ylab           = "Survival Probability",
  legend.title   = "",
  legend.labs    = "All Patients",
  ggtheme        = theme_actg(),
  risk.table.y.text = FALSE,
  cumevents      = FALSE,
  break.time.by  = 200,
  xlim           = c(0, 1250)
)

save_km_plot(p_km_overall$plot + p_km_overall$table,
             "07_km_overall.png", width = 10, height = 9)

# --------------------------------------------------------------------------- #
#  2. KM by Treatment Arm
# --------------------------------------------------------------------------- #
km_trt <- survfit(Surv(time, label) ~ trt_label, data = actg)

cat("\n--- KM by Treatment Arm ---\n")
print(summary(km_trt, times = c(365, 730, 1000))$table)

# Median survival per group
cat("\nMedian Survival by Treatment:\n")
print(surv_median(km_trt))

p_km_trt <- ggsurvplot(
  km_trt,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.30,
  pval              = TRUE,
  pval.method       = TRUE,
  pval.coord        = c(50, 0.15),
  palette           = trt_colors,
  surv.median.line  = "hv",
  legend.title      = "Treatment",
  legend.labs       = c("ZDV mono", "ZDV + ddI", "ZDV + ddC", "ddI mono"),
  title             = "Kaplan-Meier Survival by Treatment Arm",
  subtitle          = "ACTG 175 Trial — Time to AIDS Diagnosis or Death",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250),
  tables.theme      = theme_cleantable()
)

save_km_plot(p_km_trt$plot + p_km_trt$table,
             "08_km_by_treatment.png", width = 11, height = 10)

# --------------------------------------------------------------------------- #
#  3. KM by Gender
# --------------------------------------------------------------------------- #
km_gender <- survfit(Surv(time, label) ~ gender_label, data = actg)

p_km_gender <- ggsurvplot(
  km_gender,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,
  pval.method       = TRUE,
  palette           = c("#E41A1C", "#377EB8"),
  legend.title      = "Gender",
  legend.labs       = c("Female", "Male"),
  title             = "Kaplan-Meier Survival by Gender",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250)
)

save_km_plot(p_km_gender$plot + p_km_gender$table,
             "09_km_by_gender.png", width = 10, height = 8)

# --------------------------------------------------------------------------- #
#  4. KM by Race
# --------------------------------------------------------------------------- #
km_race <- survfit(Surv(time, label) ~ race_label, data = actg)

p_km_race <- ggsurvplot(
  km_race,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,
  pval.method       = TRUE,
  palette           = c("#FF7F00", "#6A3D9A"),
  legend.title      = "Race",
  legend.labs       = c("White", "Non-White"),
  title             = "Kaplan-Meier Survival by Race",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250)
)

save_km_plot(p_km_race$plot + p_km_race$table,
             "10_km_by_race.png", width = 10, height = 8)

# --------------------------------------------------------------------------- #
#  5. KM by Symptom Status
# --------------------------------------------------------------------------- #
km_symptom <- survfit(Surv(time, label) ~ symptom_label, data = actg)

p_km_symptom <- ggsurvplot(
  km_symptom,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,
  pval.method       = TRUE,
  palette           = c("#1A9850", "#D73027"),
  legend.title      = "Symptom Status",
  legend.labs       = c("Asymptomatic", "Symptomatic"),
  title             = "Kaplan-Meier Survival by Symptom Status at Baseline",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250)
)

save_km_plot(p_km_symptom$plot + p_km_symptom$table,
             "11_km_by_symptom.png", width = 10, height = 8)

# --------------------------------------------------------------------------- #
#  6. KM by Antiretroviral History
# --------------------------------------------------------------------------- #
km_str2 <- survfit(Surv(time, label) ~ str2_label, data = actg)

p_km_str2 <- ggsurvplot(
  km_str2,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,
  pval.method       = TRUE,
  palette           = c("#2166AC", "#B2182B"),
  legend.title      = "ART History",
  legend.labs       = c("Naive", "Experienced"),
  title             = "Kaplan-Meier Survival by Antiretroviral History",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250)
)

save_km_plot(p_km_str2$plot + p_km_str2$table,
             "12_km_by_art_history.png", width = 10, height = 8)

# --------------------------------------------------------------------------- #
#  7. KM by CD4 Response (Improved vs Declined)
# --------------------------------------------------------------------------- #
km_cd4resp <- survfit(Surv(time, label) ~ cd4_response, data = actg)

p_km_cd4resp <- ggsurvplot(
  km_cd4resp,
  data              = actg,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  pval              = TRUE,
  pval.method       = TRUE,
  palette           = c("#1A9850", "#D73027"),
  legend.title      = "CD4 Response",
  legend.labs       = c("Improved/Stable", "Declined"),
  title             = "Kaplan-Meier Survival by CD4 Response at Week 20",
  xlab              = "Time (Days)",
  ylab              = "Survival Probability",
  ggtheme           = theme_actg(),
  risk.table.y.text = FALSE,
  break.time.by     = 200,
  xlim              = c(0, 1250)
)

save_km_plot(p_km_cd4resp$plot + p_km_cd4resp$table,
             "13_km_by_cd4_response.png", width = 10, height = 8)

# --------------------------------------------------------------------------- #
#  8. Save KM summary table
# --------------------------------------------------------------------------- #
km_summary_trt <- data.frame(
  Treatment       = c("ZDV mono", "ZDV + ddI", "ZDV + ddC", "ddI mono"),
  N               = summary(km_trt)$table[, "records"],
  Events          = summary(km_trt)$table[, "events"],
  Event_Rate_Pct  = round(100 * summary(km_trt)$table[, "events"] /
                            summary(km_trt)$table[, "records"], 1),
  Median_Survival = summary(km_trt)$table[, "median"],
  CI_Lower        = summary(km_trt)$table[, "0.95LCL"],
  CI_Upper        = summary(km_trt)$table[, "0.95UCL"]
)

write.csv(km_summary_trt,
          file.path(dir_tables, "km_summary_by_treatment.csv"),
          row.names = FALSE)
cat("\nKM summary table saved.\n")
cat("=== 03: Kaplan-Meier Analysis Complete ===\n\n")
