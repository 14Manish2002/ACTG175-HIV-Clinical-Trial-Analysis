# =============================================================================
# ACTG 175 Survival Analysis — Cox Model Diagnostics
# Script:  06_model_diagnostics.R
# Purpose: Test PH assumption, residual analysis, influential observations
# =============================================================================

source("R/00_setup.R")
actg    <- readRDS(file.path(dir_data, "actg_processed.rds"))
cox_mod <- readRDS(file.path(dir_data, "cox_model.rds"))
cox_full <- readRDS(file.path(dir_data, "cox_full_model.rds"))

cat("=== 06: Cox Model Diagnostics ===\n")

# --------------------------------------------------------------------------- #
#  1. Schoenfeld Residuals — Test of PH Assumption
# --------------------------------------------------------------------------- #
cat("\n--- Schoenfeld Residual Test (cox.zph) ---\n")

ph_test <- cox.zph(cox_mod)
print(ph_test)

# Save PH test results
ph_df <- data.frame(
  Variable  = rownames(ph_test$table),
  Chi_Square = round(ph_test$table[, "chisq"], 3),
  df         = ph_test$table[, "df"],
  P_value    = round(ph_test$table[, "p"], 4),
  PH_Holds   = ifelse(ph_test$table[, "p"] > 0.05, "Yes", "No (Violation)")
)
print(ph_df)

write.csv(ph_df,
          file.path(dir_tables, "ph_assumption_test.csv"),
          row.names = FALSE)

cat("\nVariables with potential PH violation (p < 0.05):\n")
violations <- ph_df[ph_df$P_value < 0.05 & ph_df$Variable != "GLOBAL", ]
if (nrow(violations) == 0) {
  cat("None — PH assumption appears satisfied for all variables.\n")
} else {
  print(violations)
}

# --------------------------------------------------------------------------- #
#  2. Schoenfeld Residual Plots
# --------------------------------------------------------------------------- #
cat("\n--- Plotting Schoenfeld Residuals ---\n")

png(file.path(dir_figures, "16_schoenfeld_residuals.png"),
    width = 1400, height = 1000, res = 150)
par(mfrow = c(ceiling(length(ph_test$var) / 2), 2),
    mar   = c(4, 4, 3, 1))
plot(ph_test,
     hr   = TRUE,
     lwd  = 2,
     col  = "#2166AC",
     ylab = "Beta(t) Estimate")
dev.off()
cat("Schoenfeld residual plots saved.\n")

# --------------------------------------------------------------------------- #
#  3. Martingale Residuals — Functional form check
# --------------------------------------------------------------------------- #
cat("\n--- Martingale Residuals ---\n")

martingale_res <- residuals(cox_mod, type = "martingale")

# Check functional form for continuous predictors
continuous_vars <- c("age", "karnof", "cd40", "cd4_change")
existing_vars   <- intersect(continuous_vars, names(actg))

if (length(existing_vars) > 0) {
  mart_plots <- lapply(existing_vars, function(var) {
    df_tmp <- data.frame(
      x    = actg[[var]],
      resid = martingale_res
    )
    ggplot(df_tmp, aes(x = x, y = resid)) +
      geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
      geom_smooth(method = "loess", se = TRUE,
                  color = "#D73027", linewidth = 1.2, fill = "#FDAE61") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        title = paste("Martingale Residuals vs.", var),
        x     = var,
        y     = "Martingale Residual"
      )
  })

  mart_combined <- patchwork::wrap_plots(mart_plots, ncol = 2)
  ggsave(file.path(dir_figures, "17_martingale_residuals.png"),
         mart_combined, width = 12, height = 8, dpi = 300, bg = "white")
  cat("Martingale residual plots saved.\n")
}

# --------------------------------------------------------------------------- #
#  4. Deviance Residuals — Outlier detection
# --------------------------------------------------------------------------- #
cat("\n--- Deviance Residuals ---\n")

deviance_res <- residuals(cox_mod, type = "deviance")
actg$deviance_res <- deviance_res

# Flag large residuals (|dev| > 2)
outliers <- actg[abs(deviance_res) > 2, ]
cat(sprintf("Observations with |deviance residual| > 2: %d (%.1f%%)\n",
            nrow(outliers), 100 * nrow(outliers) / nrow(actg)))

p_deviance <- ggplot(data.frame(index = seq_along(deviance_res),
                                 dev   = deviance_res),
                     aes(x = index, y = dev)) +
  geom_point(aes(color = abs(dev) > 2), alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D73027")) +
  labs(
    title   = "Deviance Residuals — Outlier Detection",
    subtitle = "Points beyond ±2 may be influential observations",
    x       = "Observation Index",
    y       = "Deviance Residual",
    color   = "|Deviance| > 2"
  )

ggsave(file.path(dir_figures, "18_deviance_residuals.png"),
       p_deviance, width = 10, height = 5, dpi = 300, bg = "white")

# --------------------------------------------------------------------------- #
#  5. dfbeta Residuals — Influential Observations
# --------------------------------------------------------------------------- #
cat("\n--- dfbeta Residuals (Influential Observations) ---\n")

dfbeta_res <- residuals(cox_mod, type = "dfbeta")

# Plot for treatment coefficients (first few)
n_coef <- min(ncol(dfbeta_res), 4)
dfbeta_df <- as.data.frame(dfbeta_res[, 1:n_coef])
colnames(dfbeta_df) <- paste0("Beta_", 1:n_coef)
dfbeta_df$index <- 1:nrow(dfbeta_df)

dfbeta_long <- dfbeta_df %>%
  pivot_longer(-index, names_to = "Coefficient", values_to = "dfbeta")

p_dfbeta <- ggplot(dfbeta_long, aes(x = index, y = dfbeta)) +
  geom_point(alpha = 0.4, size = 0.6, color = "#4393C3") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Coefficient, scales = "free_y") +
  labs(
    title = "dfbeta Residuals by Coefficient",
    subtitle = "Large values indicate influential observations",
    x = "Observation Index",
    y = "dfbeta"
  )

ggsave(file.path(dir_figures, "19_dfbeta_residuals.png"),
       p_dfbeta, width = 12, height = 6, dpi = 300, bg = "white")

# --------------------------------------------------------------------------- #
#  6. Model fit summary
# --------------------------------------------------------------------------- #
cat("\n--- Model Fit Statistics ---\n")
cox_sum <- summary(cox_mod)
cat(sprintf("Concordance (C-index):  %.4f (SE = %.4f)\n",
            cox_sum$concordance[1], cox_sum$concordance[2]))
cat(sprintf("Log-likelihood:          %.3f\n", cox_sum$logtest["test"]))
cat(sprintf("Likelihood ratio test:   chi² = %.3f, p = %s\n",
            cox_sum$logtest["test"],
            format.pval(cox_sum$logtest["pvalue"], digits = 4)))
cat(sprintf("Wald test:               chi² = %.3f, p = %s\n",
            cox_sum$waldtest["test"],
            format.pval(cox_sum$waldtest["pvalue"], digits = 4)))
cat(sprintf("Score (log-rank) test:   chi² = %.3f, p = %s\n",
            cox_sum$sctest["test"],
            format.pval(cox_sum$sctest["pvalue"], digits = 4)))

model_fit_df <- data.frame(
  Statistic = c("Concordance (C-index)", "Log-rank p-value",
                "Wald p-value", "LR p-value"),
  Value     = c(
    round(cox_sum$concordance[1], 4),
    format.pval(cox_sum$sctest["pvalue"],  digits = 4),
    format.pval(cox_sum$waldtest["pvalue"], digits = 4),
    format.pval(cox_sum$logtest["pvalue"], digits = 4)
  )
)
write.csv(model_fit_df,
          file.path(dir_tables, "cox_model_fit.csv"),
          row.names = FALSE)

cat("=== 06: Model Diagnostics Complete ===\n\n")
