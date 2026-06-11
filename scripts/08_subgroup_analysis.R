# =============================================================================
# ACTG 175 Survival Analysis — Subgroup Analysis
# Script:  08_subgroup_analysis.R
# Purpose: Forest plot of HRs across pre-specified patient subgroups
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 08: Subgroup Analysis ===\n")

# --------------------------------------------------------------------------- #
#  1. Define subgroups and run Cox models within each
# --------------------------------------------------------------------------- #
cat("\n--- Running subgroup Cox models ---\n")

# Core Cox formula (treatment effect, adjusted for key confounders)
base_formula <- Surv(time, label) ~ treat + age + karnof + cd40

run_subgroup_cox <- function(data, subgroup_name, subgroup_val) {
  tryCatch({
    mod <- coxph(base_formula, data = data, ties = "efron")
    s   <- summary(mod)
    ci  <- exp(confint(mod))

    trt_row <- grep("treat", rownames(s$coefficients), value = TRUE)[1]
    if (is.na(trt_row)) return(NULL)

    data.frame(
      Subgroup    = subgroup_name,
      Category    = as.character(subgroup_val),
      N           = nrow(data),
      Events      = sum(data$label),
      HR          = round(exp(s$coefficients[trt_row, "coef"]), 3),
      CI_Lower    = round(ci[trt_row, 1], 3),
      CI_Upper    = round(ci[trt_row, 2], 3),
      P_value     = round(s$coefficients[trt_row, "Pr(>|z|)"], 4),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Skipped %s = %s: %s\n", subgroup_name, subgroup_val, e$message))
    NULL
  })
}

# Define subgroup splits
subgroup_list <- list(
  list(var = "gender_label",  name = "Gender"),
  list(var = "race_label",    name = "Race"),
  list(var = "str2_label",    name = "ART History"),
  list(var = "symptom_label", name = "Symptom Status"),
  list(var = "drugs_label",   name = "IV Drug Use"),
  list(var = "hemo_label",    name = "Hemophilia"),
  list(var = "homo_label",    name = "Homosexual Activity"),
  list(var = "age_group",     name = "Age Group"),
  list(var = "karnof_cat",    name = "Karnofsky Category"),
  list(var = "cd4_response",  name = "CD4 Response")
)

subgroup_results <- list()

for (sg in subgroup_list) {
  var   <- sg$var
  sname <- sg$name
  lvls  <- levels(actg[[var]])

  cat(sprintf("Processing subgroup: %s\n", sname))

  for (lv in lvls) {
    sub_data <- actg[actg[[var]] == lv, ]
    if (nrow(sub_data) < 30 || sum(sub_data$label) < 5) {
      cat(sprintf("  Skipped %s = %s (too few observations)\n", sname, lv))
      next
    }
    result <- run_subgroup_cox(sub_data, sname, lv)
    if (!is.null(result)) subgroup_results[[length(subgroup_results) + 1]] <- result
  }
}

subgroup_df <- do.call(rbind, subgroup_results)
write.csv(subgroup_df,
          file.path(dir_tables, "subgroup_analysis.csv"),
          row.names = FALSE)
cat("\nSubgroup results saved.\n")
print(subgroup_df)

# --------------------------------------------------------------------------- #
#  2. Overall treatment effect (all patients)
# --------------------------------------------------------------------------- #
cox_overall <- coxph(base_formula, data = actg, ties = "efron")
cox_sum_ov  <- summary(cox_overall)
trt_row     <- grep("treat", rownames(cox_sum_ov$coefficients), value = TRUE)[1]

overall_row <- data.frame(
  Subgroup = "Overall",
  Category = "All Patients",
  N        = nrow(actg),
  Events   = sum(actg$label),
  HR       = round(exp(cox_sum_ov$coefficients[trt_row, "coef"]), 3),
  CI_Lower = round(exp(confint(cox_overall)[trt_row, 1]), 3),
  CI_Upper = round(exp(confint(cox_overall)[trt_row, 2]), 3),
  P_value  = round(cox_sum_ov$coefficients[trt_row, "Pr(>|z|)"], 4)
)

# Add overall at the top
subgroup_df_all <- rbind(overall_row, subgroup_df)

# --------------------------------------------------------------------------- #
#  3. Forest Plot — Subgroup Analysis
# --------------------------------------------------------------------------- #
cat("\n--- Creating Subgroup Forest Plot ---\n")

subgroup_df_all <- subgroup_df_all %>%
  mutate(
    label_text  = paste0(Subgroup, ": ", Category),
    label_text  = factor(label_text, levels = rev(unique(label_text))),
    is_overall  = Subgroup == "Overall",
    event_label = paste0(Events, "/", N)
  )

p_forest_sg <- ggplot(subgroup_df_all,
                      aes(x = HR, y = label_text)) +
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.8) +

  # Subgroup markers
  geom_point(aes(size = ifelse(is_overall, 5, 3),
                 color = ifelse(is_overall, "Overall", "Subgroup"),
                 shape = ifelse(is_overall, 18, 16)),
             show.legend = TRUE) +

  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper,
                     color = ifelse(is_overall, "Overall", "Subgroup")),
                 height = 0.3, linewidth = 0.7) +

  # Text: sample size
  geom_text(aes(x = max(CI_Upper, na.rm = TRUE) * 1.05,
                label = event_label),
            hjust = 0, size = 3, color = "grey30") +

  scale_color_manual(values = c("Overall" = "#D73027", "Subgroup" = "#2166AC"),
                     guide = "none") +
  scale_size_identity() +
  scale_shape_identity() +
  scale_x_log10(breaks = c(0.4, 0.6, 0.8, 1.0, 1.5, 2.0),
                labels = scales::label_number()) +
  coord_cartesian(xlim = c(0.3, max(subgroup_df_all$CI_Upper, na.rm = TRUE) * 1.5)) +
  labs(
    title    = "Subgroup Analysis — Treatment Effect on AIDS/Death",
    subtitle = "HR < 1 favors combination/alternative therapy vs ZDV monotherapy",
    x        = "Hazard Ratio (log scale) | Events/N",
    y        = "",
    caption  = "Adjusted for age, Karnofsky score, and baseline CD4 count"
  ) +
  theme(
    axis.text.y     = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

ggsave(file.path(dir_figures, "22_subgroup_forest_plot.png"),
       p_forest_sg, width = 12, height = 10, dpi = 300, bg = "white")
cat("Subgroup forest plot saved.\n")

# --------------------------------------------------------------------------- #
#  4. Interaction tests (treatment × subgroup)
# --------------------------------------------------------------------------- #
cat("\n--- Testing Interaction: Treatment × Key Subgroups ---\n")

interaction_vars <- c("gender", "race", "str2", "symptom", "drugs")

interaction_results <- lapply(interaction_vars, function(var) {
  fmla_int  <- as.formula(paste("Surv(time, label) ~ treat *", var,
                                "+ age + karnof + cd40"))
  fmla_main <- as.formula(paste("Surv(time, label) ~ treat +", var,
                                "+ age + karnof + cd40"))
  mod_int   <- coxph(fmla_int,  data = actg, ties = "efron")
  mod_main  <- coxph(fmla_main, data = actg, ties = "efron")
  lrt       <- anova(mod_main, mod_int)

  data.frame(
    Subgroup        = var,
    LRT_Chi_Square  = round(lrt[2, "Chisq"], 3),
    df              = lrt[2, "Df"],
    P_interaction   = round(lrt[2, "P(>|Chi|)"], 4),
    Interaction_Sig = ifelse(lrt[2, "P(>|Chi|)"] < 0.05, "Yes *", "No")
  )
})

interaction_df <- do.call(rbind, interaction_results)
print(interaction_df)

write.csv(interaction_df,
          file.path(dir_tables, "interaction_tests.csv"),
          row.names = FALSE)

cat("=== 08: Subgroup Analysis Complete ===\n\n")
