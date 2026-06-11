# =============================================================================
# ACTG 175 Survival Analysis — Exploratory Data Analysis
# Script:  02_eda.R
# Purpose: Descriptive statistics, Table 1, and exploratory visualizations
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 02: Exploratory Data Analysis ===\n")

# --------------------------------------------------------------------------- #
#  1. Table 1 — Baseline Characteristics
# --------------------------------------------------------------------------- #
cat("\n--- Creating Table 1 ---\n")

table1_vars <- c(
  "age", "gender_label", "race_label", "wtkg", "karnof",
  "hemo_label", "homo_label", "drugs_label", "symptom_label",
  "str2_label", "strat_label", "preanti",
  "cd40", "cd420", "cd80", "cd820",
  "offtrt_label", "label"
)

table1_cat_vars <- c(
  "gender_label", "race_label", "hemo_label", "homo_label",
  "drugs_label", "symptom_label", "str2_label", "strat_label",
  "offtrt_label", "label"
)

tbl1 <- CreateTableOne(
  vars     = table1_vars,
  strata   = "trt_label",
  data     = actg,
  factorVars = table1_cat_vars,
  addOverall = TRUE
)

tbl1_print <- print(tbl1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tbl1_print, file = file.path(dir_tables, "table1_baseline.csv"))
cat("Table 1 saved to output/tables/table1_baseline.csv\n")

# --------------------------------------------------------------------------- #
#  2. Distribution Plots
# --------------------------------------------------------------------------- #

# 2a. Time-to-event distribution by treatment
p_time_dist <- ggplot(actg, aes(x = time, fill = trt_label)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~trt_label, nrow = 2) +
  scale_fill_manual(values = trt_colors) +
  scale_x_continuous(breaks = seq(0, 1250, 250)) +
  labs(
    title    = "Distribution of Time-to-Event by Treatment Arm",
    subtitle = "ACTG 175 Trial — 2,139 HIV-Infected Patients",
    x        = "Time to Event or Censoring (Days)",
    y        = "Count",
    fill     = "Treatment",
    caption  = "Shaded = event occurred; white bars = censored patients"
  ) +
  theme(legend.position = "none")

ggsave(file.path(dir_figures, "01_time_distribution.png"),
       p_time_dist, width = 10, height = 7, dpi = 300, bg = "white")

# 2b. CD4 count at baseline and follow-up
cd4_long <- actg %>%
  select(trt_label, cd40, cd420) %>%
  pivot_longer(c(cd40, cd420),
               names_to  = "timepoint",
               values_to = "cd4") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("cd40", "cd420"),
                            labels = c("Baseline", "Week 20")))

p_cd4 <- ggplot(cd4_long, aes(x = timepoint, y = cd4, fill = trt_label)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, alpha = 0.8) +
  facet_wrap(~trt_label) +
  scale_fill_manual(values = trt_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "CD4 Cell Count: Baseline vs. Week 20",
    subtitle = "By Treatment Arm",
    x        = "Time Point",
    y        = "CD4 Cell Count (cells/mm³)",
    fill     = "Treatment"
  ) +
  theme(legend.position = "none")

ggsave(file.path(dir_figures, "02_cd4_boxplot.png"),
       p_cd4, width = 10, height = 6, dpi = 300, bg = "white")

# 2c. CD4 change from baseline
p_cd4_change <- ggplot(actg, aes(x = trt_label, y = cd4_change, fill = trt_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = trt_colors) +
  labs(
    title    = "Change in CD4 Count from Baseline to Week 20",
    subtitle = "Positive values indicate immunological improvement",
    x        = "Treatment Arm",
    y        = "CD4 Change (cells/mm³)",
    fill     = "Treatment"
  ) +
  theme(legend.position = "none")

ggsave(file.path(dir_figures, "03_cd4_change_violin.png"),
       p_cd4_change, width = 9, height = 6, dpi = 300, bg = "white")

# 2d. Age distribution by event
p_age <- ggplot(actg, aes(x = age, fill = factor(label))) +
  geom_density(alpha = 0.6, adjust = 1.2) +
  scale_fill_manual(values  = c("0" = "#2166AC", "1" = "#D73027"),
                    labels  = c("Censored", "Event (AIDS/Death)")) +
  labs(
    title = "Age Distribution by Event Status",
    x     = "Age (years)",
    y     = "Density",
    fill  = "Outcome"
  )

ggsave(file.path(dir_figures, "04_age_density.png"),
       p_age, width = 8, height = 5, dpi = 300, bg = "white")

# 2e. Event rate by treatment
event_rate <- actg %>%
  group_by(trt_label) %>%
  summarise(
    n        = n(),
    events   = sum(label),
    rate_pct = 100 * mean(label),
    .groups  = "drop"
  )

p_event_rate <- ggplot(event_rate, aes(x = trt_label, y = rate_pct, fill = trt_label)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", rate_pct, events)),
            vjust = -0.3, fontface = "bold", size = 4) +
  scale_fill_manual(values = trt_colors) +
  scale_y_continuous(limits = c(0, 35), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Event Rate (AIDS/Death) by Treatment Arm",
    subtitle = "ACTG 175 — Primary Endpoint",
    x        = "Treatment Arm",
    y        = "Event Rate (%)",
    fill     = "Treatment"
  ) +
  theme(legend.position = "none")

ggsave(file.path(dir_figures, "05_event_rate_bar.png"),
       p_event_rate, width = 8, height = 5, dpi = 300, bg = "white")

# 2f. Karnofsky score distribution
p_karnof <- ggplot(actg, aes(x = karnof, fill = trt_label)) +
  geom_bar(position = "dodge", color = "white") +
  scale_fill_manual(values = trt_colors) +
  labs(
    title = "Karnofsky Performance Score Distribution by Treatment",
    x     = "Karnofsky Score",
    y     = "Count",
    fill  = "Treatment"
  )

ggsave(file.path(dir_figures, "06_karnofsky_bar.png"),
       p_karnof, width = 9, height = 5, dpi = 300, bg = "white")

# 2g. Correlation matrix of continuous variables
cat("\n--- Correlation matrix of key continuous variables ---\n")
cont_vars <- actg %>%
  select(time, age, wtkg, karnof, preanti, cd40, cd420, cd80, cd820,
         cd4_change, cd8_change)

cor_mat <- round(cor(cont_vars, use = "complete.obs"), 2)
write.csv(cor_mat, file.path(dir_tables, "correlation_matrix.csv"))

cat("EDA plots saved to output/figures/\n")
cat("=== 02: EDA Complete ===\n\n")
