# =============================================================================
# ACTG 175 Survival Analysis — Publication Tables
# Script:  09_report_tables.R
# Purpose: Generate all publication-quality tables (HTML + CSV)
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 09: Report Tables ===\n")

# --------------------------------------------------------------------------- #
#  1. Table 1 — Baseline Characteristics (gtsummary)
# --------------------------------------------------------------------------- #
cat("\n--- Table 1: Baseline Characteristics ---\n")

tbl_baseline <- actg %>%
  select(
    trt_label,
    age, gender_label, race_label, wtkg,
    karnof, karnof_cat,
    hemo_label, homo_label, drugs_label, symptom_label,
    str2_label, strat_label, preanti,
    cd40, cd420, cd80, cd820
  ) %>%
  tbl_summary(
    by     = trt_label,
    missing = "no",
    label  = list(
      age           ~ "Age (years)",
      gender_label  ~ "Gender",
      race_label    ~ "Race",
      wtkg          ~ "Weight (kg)",
      karnof        ~ "Karnofsky Score",
      karnof_cat    ~ "Karnofsky Category",
      hemo_label    ~ "Hemophilia",
      homo_label    ~ "Homosexual Activity",
      drugs_label   ~ "IV Drug Use",
      symptom_label ~ "Symptomatic at Baseline",
      str2_label    ~ "ART History",
      strat_label   ~ "ART Stratum",
      preanti       ~ "Prior ART (months)",
      cd40          ~ "CD4 Count (baseline, cells/mm³)",
      cd420         ~ "CD4 Count (week 20, cells/mm³)",
      cd80          ~ "CD8 Count (baseline, cells/mm³)",
      cd820         ~ "CD8 Count (week 20, cells/mm³)"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1)
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous()  ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_caption("**Table 1. Baseline Patient Characteristics by Treatment Arm**") %>%
  modify_header(label = "**Variable**")

# Save as HTML
tbl_baseline_html <- as_gt(tbl_baseline)
gt::gtsave(tbl_baseline_html,
           file.path(dir_tables, "table1_baseline_gtsummary.html"))
cat("Table 1 saved as HTML.\n")

# --------------------------------------------------------------------------- #
#  2. Table 2 — Kaplan-Meier Summary
# --------------------------------------------------------------------------- #
cat("\n--- Table 2: KM Survival Summary ---\n")

km_trt <- survfit(Surv(time, label) ~ trt_label, data = actg)

km_tab <- surv_summary(km_trt, data = actg) %>%
  filter(time %in% c(365, 730, 1000)) %>%
  select(trt_label, time, n.risk, n.event, surv, lower, upper) %>%
  mutate(
    Survival_Pct = paste0(round(100 * surv, 1), "% (",
                          round(100 * lower, 1), "–",
                          round(100 * upper, 1), ")"),
    .keep = "unused"
  )

colnames(km_tab) <- c("Treatment", "Time (days)", "At Risk",
                       "Events", "Survival (95% CI)")
print(km_tab)
write.csv(km_tab,
          file.path(dir_tables, "table2_km_survival.csv"),
          row.names = FALSE)

# --------------------------------------------------------------------------- #
#  3. Table 3 — Cox Multivariate Model Results
# --------------------------------------------------------------------------- #
cat("\n--- Table 3: Cox Multivariate Model ---\n")

cox_mod <- readRDS(file.path(dir_data, "cox_model.rds"))

tbl_cox <- cox_mod %>%
  tbl_regression(
    exp        = TRUE,
    label      = list(
      trt_label     ~ "Treatment Arm",
      age           ~ "Age (per year)",
      karnof        ~ "Karnofsky Score (per unit)",
      cd40          ~ "CD4 Baseline (per cell)",
      cd420         ~ "CD4 Week-20 (per cell)",
      cd4_change    ~ "CD4 Change",
      gender_label  ~ "Gender",
      race_label    ~ "Race",
      str2_label    ~ "ART History",
      symptom_label ~ "Symptom Status",
      drugs         ~ "IV Drug Use",
      homo          ~ "Homosexual Activity",
      hemo          ~ "Hemophilia",
      offtrt        ~ "Off Treatment"
    )
  ) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  add_global_p() %>%
  modify_caption("**Table 3. Multivariable Cox Proportional Hazards Model**") %>%
  modify_header(estimate = "**HR**", ci = "**95% CI**", p.value = "**p-value**")

tbl_cox_html <- as_gt(tbl_cox)
gt::gtsave(tbl_cox_html,
           file.path(dir_tables, "table3_cox_model.html"))
cat("Table 3 saved as HTML.\n")

# --------------------------------------------------------------------------- #
#  4. Table 4 — Log-rank Test Results Summary
# --------------------------------------------------------------------------- #
cat("\n--- Table 4: Log-Rank Test Results ---\n")

lr_tab <- read.csv(file.path(dir_tables, "logrank_covariate_tests.csv"))
print(lr_tab)

# --------------------------------------------------------------------------- #
#  5. Table 5 — PH Assumption Test
# --------------------------------------------------------------------------- #
cat("\n--- Table 5: PH Assumption (cox.zph) ---\n")

ph_tab <- read.csv(file.path(dir_tables, "ph_assumption_test.csv"))
print(ph_tab)

# --------------------------------------------------------------------------- #
#  6. Master HTML report (all tables)
# --------------------------------------------------------------------------- #
cat("\n--- Generating combined HTML table report ---\n")

html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>ACTG 175 — Analysis Tables</title>
  <style>
    body  { font-family: Arial, sans-serif; max-width: 1100px; margin: auto;
            padding: 20px; background: #f9f9f9; }
    h1    { color: #1A237E; border-bottom: 3px solid #1A237E; padding-bottom: 8px; }
    h2    { color: #283593; margin-top: 40px; }
    table { border-collapse: collapse; width: 100%%; margin-bottom: 20px;
            background: white; box-shadow: 0 1px 3px rgba(0,0,0,.1); }
    th    { background: #1A237E; color: white; padding: 10px 12px;
            text-align: left; font-size: 13px; }
    td    { padding: 8px 12px; border-bottom: 1px solid #E3E3E3;
            font-size: 12px; }
    tr:nth-child(even) { background: #F5F5F5; }
    tr:hover { background: #EEF2FF; }
    .sig  { color: #B71C1C; font-weight: bold; }
    .caption { font-style: italic; color: #555; margin-bottom: 6px;
               font-size: 13px; }
    footer { text-align: center; color: #aaa; margin-top: 40px;
             font-size: 11px; }
  </style>
</head>
<body>
<h1>ACTG 175 HIV Clinical Trial — Analysis Tables</h1>
<p style="color:#555">Generated: %s | R version: %s</p>

<h2>Table 1. Baseline Characteristics</h2>
<p class="caption">See: output/tables/table1_baseline_gtsummary.html</p>

<h2>Table 2. Kaplan-Meier Survival Summary</h2>
%s

<h2>Table 3. Multivariate Cox PH Model</h2>
<p class="caption">See: output/tables/table3_cox_model.html</p>

<h2>Table 4. Log-Rank Test Results</h2>
%s

<h2>Table 5. PH Assumption Test (Schoenfeld Residuals)</h2>
%s

<footer>ACTG 175 Survival Analysis | HIV Clinical Trials | R Statistical Analysis</footer>
</body></html>',
  Sys.time(),
  R.version.string,
  knitr::kable(km_tab, format = "html", table.attr = 'class="table"'),
  knitr::kable(lr_tab[, c("Covariate", "Chi_Square", "df",
                           "P_value", "Significant")],
               format = "html", table.attr = 'class="table"'),
  knitr::kable(ph_tab[, c("Variable", "Chi_Square", "df",
                           "P_value", "PH_Holds")],
               format = "html", table.attr = 'class="table"')
)

writeLines(html_content,
           file.path(dir_tables, "all_tables_report.html"))
cat("Combined HTML table report saved.\n")

cat("=== 09: Report Tables Complete ===\n\n")
cat("=== ALL ANALYSIS SCRIPTS COMPLETE ===\n")
cat("Output files are in: output/figures/ and output/tables/\n")
