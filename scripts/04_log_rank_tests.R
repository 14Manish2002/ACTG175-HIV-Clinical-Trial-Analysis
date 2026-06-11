# =============================================================================
# ACTG 175 Survival Analysis — Log-Rank & Statistical Tests
# Script:  04_log_rank_tests.R
# Purpose: Formal group comparisons via log-rank, Wilcoxon, and pairwise tests
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 04: Log-Rank Tests ===\n")

# --------------------------------------------------------------------------- #
#  Helper function: format test results
# --------------------------------------------------------------------------- #
format_test <- function(test_obj, group_name) {
  data.frame(
    Comparison   = group_name,
    Chi_Square   = round(test_obj$chisq, 3),
    df           = length(test_obj$n) - 1,
    P_value      = format.pval(test_obj$pvalue, digits = 4, eps = 0.0001),
    Significant  = ifelse(test_obj$pvalue < 0.05, "Yes ***", "No")
  )
}

# --------------------------------------------------------------------------- #
#  1. Log-rank tests by treatment arm
# --------------------------------------------------------------------------- #
cat("\n--- Log-Rank Test: Treatment Arm ---\n")
lr_trt <- survdiff(Surv(time, label) ~ trt_label, data = actg)
print(lr_trt)

cat(sprintf(
  "\nChi-Square = %.3f, df = %d, p-value = %s\n",
  lr_trt$chisq,
  length(lr_trt$n) - 1,
  format.pval(pchisq(lr_trt$chisq, df = length(lr_trt$n) - 1,
                     lower.tail = FALSE), digits = 4)
))

# --------------------------------------------------------------------------- #
#  2. Pairwise log-rank tests (all treatment pairs)
# --------------------------------------------------------------------------- #
cat("\n--- Pairwise Log-Rank Tests (Treatment Arms) ---\n")
trt_levels <- levels(actg$trt_label)
pairs      <- combn(trt_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(pair) {
  sub_data <- actg %>% filter(trt_label %in% pair) %>%
    mutate(trt_label = droplevels(trt_label))
  test <- survdiff(Surv(time, label) ~ trt_label, data = sub_data)
  pval <- pchisq(test$chisq, df = 1, lower.tail = FALSE)
  data.frame(
    Group1      = pair[1],
    Group2      = pair[2],
    Chi_Square  = round(test$chisq, 3),
    P_value_raw = round(pval, 4),
    P_adj_BH    = NA_real_   # filled after BH correction below
  )
})

pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$P_adj_BH <- round(p.adjust(pairwise_df$P_value_raw, method = "BH"), 4)
pairwise_df$Significant_adj <- ifelse(pairwise_df$P_adj_BH < 0.05, "Yes *", "No")

print(pairwise_df)
write.csv(pairwise_df,
          file.path(dir_tables, "pairwise_logrank_tests.csv"),
          row.names = FALSE)

# --------------------------------------------------------------------------- #
#  3. Log-rank by other key covariates
# --------------------------------------------------------------------------- #
covariates <- list(
  "Gender"        = "gender_label",
  "Race"          = "race_label",
  "Symptom Status" = "symptom_label",
  "ART History"   = "str2_label",
  "IV Drug Use"   = "drugs_label",
  "Hemophilia"    = "hemo_label",
  "Off Treatment" = "offtrt_label",
  "CD4 Response"  = "cd4_response",
  "Age Group"     = "age_group"
)

cat("\n--- Log-Rank Tests: Other Covariates ---\n")

covariate_results <- lapply(names(covariates), function(cov_name) {
  var   <- covariates[[cov_name]]
  fmla  <- as.formula(paste("Surv(time, label) ~", var))
  test  <- survdiff(fmla, data = actg)
  df_val <- length(test$n) - 1
  pval   <- pchisq(test$chisq, df = df_val, lower.tail = FALSE)
  data.frame(
    Covariate   = cov_name,
    Variable    = var,
    Chi_Square  = round(test$chisq, 3),
    df          = df_val,
    P_value     = round(pval, 4),
    Significant = ifelse(pval < 0.05, "Yes", "No")
  )
})

covariate_test_df <- do.call(rbind, covariate_results)
print(covariate_test_df)
write.csv(covariate_test_df,
          file.path(dir_tables, "logrank_covariate_tests.csv"),
          row.names = FALSE)

# --------------------------------------------------------------------------- #
#  4. Wilcoxon (Peto-Peto) test for early differences
# --------------------------------------------------------------------------- #
cat("\n--- Peto-Peto (Wilcoxon) Test: Treatment Arm ---\n")
wilcox_trt <- survdiff(Surv(time, label) ~ trt_label, data = actg, rho = 1)
cat("Peto-Peto chi-square:", round(wilcox_trt$chisq, 3), "\n")
cat("p-value:", format.pval(
  pchisq(wilcox_trt$chisq, df = length(wilcox_trt$n) - 1, lower.tail = FALSE),
  digits = 4
), "\n")

# --------------------------------------------------------------------------- #
#  5. Visualize p-values — covariate significance chart
# --------------------------------------------------------------------------- #
covariate_test_df$neg_log_p <- -log10(covariate_test_df$P_value)
covariate_test_df$Covariate <- factor(
  covariate_test_df$Covariate,
  levels = covariate_test_df$Covariate[order(covariate_test_df$neg_log_p)]
)

p_logrank_bar <- ggplot(covariate_test_df,
                        aes(x = Covariate, y = neg_log_p,
                            fill = Significant)) +
  geom_col(width = 0.7, color = "white") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 0.7, y = -log10(0.05) + 0.1,
           label = "p = 0.05", color = "red", size = 3.5, hjust = 0) +
  scale_fill_manual(values = c("Yes" = "#1B7837", "No" = "#BDBDBD")) +
  coord_flip() +
  labs(
    title    = "Log-Rank Test: -log₁₀(p) for Key Covariates",
    subtitle = "Bars crossing the dashed line indicate statistical significance (p < 0.05)",
    x        = "Covariate",
    y        = "-log₁₀(p-value)",
    fill     = "Significant?"
  )

ggsave(file.path(dir_figures, "14_logrank_significance.png"),
       p_logrank_bar, width = 9, height = 6, dpi = 300, bg = "white")

cat("\nAll log-rank test results saved.\n")
cat("=== 04: Log-Rank Tests Complete ===\n\n")
