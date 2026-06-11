# =============================================================================
# ACTG 175 Survival Analysis — Cox Proportional Hazards Model
# Script:  05_cox_model.R
# Purpose: Univariate screening + multivariate Cox PH + forest plot
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 05: Cox Proportional Hazards Model ===\n")

# --------------------------------------------------------------------------- #
#  1. Univariate Cox Models
# --------------------------------------------------------------------------- #
cat("\n--- Univariate Cox Models ---\n")

uni_covariates <- c(
  "trt",        "age",     "wtkg",    "karnof",   "preanti",
  "cd40",       "cd420",   "cd80",    "cd820",    "cd4_change",
  "hemo",       "homo",    "drugs",   "gender",   "race",
  "str2",       "symptom", "offtrt",  "oprior",   "z30"
)

uni_results <- lapply(uni_covariates, function(var) {
  fmla <- as.formula(paste("Surv(time, label) ~", var))
  mod  <- coxph(fmla, data = actg)
  s    <- summary(mod)
  ci   <- exp(confint(mod))

  data.frame(
    Variable  = var,
    HR        = round(exp(coef(mod)), 3),
    CI_Lower  = round(ci[1], 3),
    CI_Upper  = round(ci[2], 3),
    P_value   = round(s$coefficients[, 5], 4),
    Concordance = round(s$concordance[1], 3),
    stringsAsFactors = FALSE
  )
})

uni_df <- do.call(rbind, uni_results)
uni_df$Significant <- ifelse(uni_df$P_value < 0.05, "Yes *", "No")
uni_df <- uni_df[order(uni_df$P_value), ]

cat("\nTop univariate predictors (p < 0.05):\n")
print(uni_df[uni_df$P_value < 0.05, ])

write.csv(uni_df,
          file.path(dir_tables, "cox_univariate_results.csv"),
          row.names = FALSE)

# --------------------------------------------------------------------------- #
#  2. Multivariate Cox Model — Full Model
# --------------------------------------------------------------------------- #
cat("\n--- Multivariate Cox Model (Full) ---\n")

cox_full <- coxph(
  Surv(time, label) ~
    trt_label + age + wtkg + karnof + preanti +
    cd40 + cd420 + cd4_change +
    hemo + homo + drugs + gender_label + race_label +
    str2_label + symptom_label + offtrt,
  data  = actg,
  ties  = "efron",
  x     = TRUE
)

cat("\nFull Cox Model Summary:\n")
print(summary(cox_full))

# --------------------------------------------------------------------------- #
#  3. Multivariate Cox Model — Reduced Model (AIC-based stepwise)
# --------------------------------------------------------------------------- #
cat("\n--- Backward Stepwise Selection (AIC) ---\n")

cox_step <- step(cox_full, direction = "backward", trace = 0)
cat("\nReduced Cox Model Summary:\n")
print(summary(cox_step))

cox_step_summary <- summary(cox_step)

# --------------------------------------------------------------------------- #
#  4. Extract and save model coefficients
# --------------------------------------------------------------------------- #
cox_coef_df <- data.frame(
  Variable   = rownames(cox_step_summary$coefficients),
  HR         = round(exp(cox_step_summary$coefficients[, "coef"]), 3),
  HR_CI_Low  = round(cox_step_summary$conf.int[, "lower .95"], 3),
  HR_CI_High = round(cox_step_summary$conf.int[, "upper .95"], 3),
  Z_score    = round(cox_step_summary$coefficients[, "z"], 3),
  P_value    = round(cox_step_summary$coefficients[, "Pr(>|z|)"], 4),
  stringsAsFactors = FALSE
)
cox_coef_df$Significant <- ifelse(cox_coef_df$P_value < 0.05, "Yes", "No")

write.csv(cox_coef_df,
          file.path(dir_tables, "cox_multivariate_results.csv"),
          row.names = FALSE)
cat("\nCox model results saved.\n")

# --------------------------------------------------------------------------- #
#  5. Forest Plot — Hazard Ratios
# --------------------------------------------------------------------------- #
cat("\n--- Creating Forest Plot ---\n")

# Select meaningful variables for the forest plot
forest_df <- cox_coef_df %>%
  filter(!grepl("Intercept", Variable)) %>%
  mutate(
    Variable = case_when(
      Variable == "trt_labelZDV + ddI"   ~ "ZDV + ddI vs ZDV mono",
      Variable == "trt_labelZDV + ddC"   ~ "ZDV + ddC vs ZDV mono",
      Variable == "trt_labelddI mono"    ~ "ddI mono vs ZDV mono",
      Variable == "age"                  ~ "Age (per year)",
      Variable == "karnof"               ~ "Karnofsky Score (per unit)",
      Variable == "cd40"                 ~ "CD4 Baseline (per 50 cells)",
      Variable == "cd420"                ~ "CD4 Week-20 (per 50 cells)",
      Variable == "cd4_change"           ~ "CD4 Change Baseline→Wk20",
      Variable == "wtkg"                 ~ "Weight (per kg)",
      Variable == "preanti"              ~ "Prior ART (months)",
      Variable == "gender_labelMale"     ~ "Male (vs Female)",
      Variable == "race_labelNon-White"  ~ "Non-White (vs White)",
      Variable == "str2_labelExperienced" ~ "ART Experienced (vs Naive)",
      Variable == "symptom_labelSymptomatic" ~ "Symptomatic at Baseline",
      Variable == "homo"                 ~ "Homosexual Activity",
      Variable == "drugs"                ~ "IV Drug Use",
      Variable == "hemo"                 ~ "Hemophilia",
      Variable == "offtrt"               ~ "Taken Off Treatment",
      TRUE                               ~ Variable
    )
  ) %>%
  arrange(HR)

forest_df$Variable <- factor(forest_df$Variable, levels = forest_df$Variable)

p_forest <- ggplot(forest_df, aes(x = HR, y = Variable)) +
  geom_point(aes(color = Significant), size = 3.5) +
  geom_errorbarh(aes(xmin = HR_CI_Low, xmax = HR_CI_High,
                     color = Significant),
                 height = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("Yes" = "#D73027", "No" = "#92C5DE")) +
  scale_x_log10(breaks = c(0.5, 0.7, 1.0, 1.5, 2.0, 3.0),
                labels = scales::label_number()) +
  labs(
    title    = "Multivariable Cox Proportional Hazards Model",
    subtitle = "Hazard Ratios (95% CI) — Outcome: AIDS/Death",
    x        = "Hazard Ratio (log scale)",
    y        = "",
    color    = "Significant (p<0.05)",
    caption  = "Reference: ZDV monotherapy, Female, White, ART Naive"
  ) +
  theme(
    axis.text.y    = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(file.path(dir_figures, "15_cox_forest_plot.png"),
       p_forest, width = 10, height = 8, dpi = 300, bg = "white")

# --------------------------------------------------------------------------- #
#  6. Predicted survival curves from Cox model
# --------------------------------------------------------------------------- #
cat("\n--- Adjusted Survival Curves from Cox Model ---\n")

# Marginal / adjusted survival by treatment (other covariates at mean/reference)
new_data <- expand.grid(
  trt_label     = levels(actg$trt_label),
  age           = median(actg$age),
  karnof        = median(actg$karnof),
  cd40          = median(actg$cd40),
  cd420         = median(actg$cd420),
  cd4_change    = median(actg$cd4_change),
  gender_label  = "Male",
  race_label    = "White",
  str2_label    = "Naive",
  symptom_label = "Asymptomatic"
)

# Only keep variables that are in the stepwise model
model_vars   <- all.vars(formula(cox_step))[-c(1, 2)]  # remove time, event
new_data_sub <- new_data[, intersect(names(new_data), model_vars), drop = FALSE]

cat("\nCox Model — AIC:", round(AIC(cox_step), 1),
    "| Concordance:", round(summary(cox_step)$concordance[1], 3), "\n")

# Save model object for use in diagnostics
saveRDS(cox_step, file.path(dir_data, "cox_model.rds"))
saveRDS(cox_full, file.path(dir_data, "cox_full_model.rds"))

cat("Cox models saved to data/\n")
cat("=== 05: Cox PH Model Complete ===\n\n")
