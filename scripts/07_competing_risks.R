# =============================================================================
# ACTG 175 Survival Analysis — Competing Risks Analysis
# Script:  07_competing_risks.R
# Purpose: Cumulative incidence functions & Fine-Gray sub-distribution model
# =============================================================================

source("R/00_setup.R")
actg <- readRDS(file.path(dir_data, "actg_processed.rds"))

cat("=== 07: Competing Risks Analysis ===\n")

# --------------------------------------------------------------------------- #
# Note on ACTG 175 competing risks structure:
#   Event 1 = AIDS diagnosis
#   Event 2 = Death (competing with AIDS)
#   Event 0 = Censored (off-treatment, loss-to-follow-up, end-of-study)
#
# The 'label' variable codes the composite endpoint (AIDS OR death).
# We create synthetic competing event codes for illustration:
#   - label=1 + offtrt=0 → Event type 1 (on-treatment event)
#   - label=1 + offtrt=1 → Event type 2 (off-treatment event, proxy competing)
#   - label=0            → Censored
# --------------------------------------------------------------------------- #

actg <- actg %>%
  mutate(
    event_type = case_when(
      label == 1 & offtrt == 0 ~ 1L,  # Primary event (on-treatment)
      label == 1 & offtrt == 1 ~ 2L,  # Secondary event (competing)
      TRUE                     ~ 0L   # Censored
    ),
    event_type = factor(event_type, levels = 0:2,
                        labels = c("Censored", "On-Tx Event", "Off-Tx Event"))
  )

cat("\nEvent type distribution:\n")
print(table(actg$event_type, actg$trt_label))

# --------------------------------------------------------------------------- #
#  1. Cumulative Incidence Functions — Overall
# --------------------------------------------------------------------------- #
cat("\n--- Cumulative Incidence Functions (Overall) ---\n")

cif_overall <- cuminc(
  ftime   = actg$time,
  fstatus = as.integer(actg$event_type) - 1L,  # 0=cens, 1=primary, 2=competing
  cencode = 0
)

cat("Overall CIF summary:\n")
# Print at key timepoints
print(timepoints(cif_overall, times = c(365, 730, 1000)))

# --------------------------------------------------------------------------- #
#  2. Cumulative Incidence by Treatment Arm
# --------------------------------------------------------------------------- #
cat("\n--- CIF by Treatment Arm ---\n")

cif_trt <- cuminc(
  ftime   = actg$time,
  fstatus = as.integer(actg$event_type) - 1L,
  group   = actg$trt_label,
  cencode = 0
)

cat("\nGray's test for equality of CIFs:\n")
print(cif_trt$Tests)

write.csv(as.data.frame(cif_trt$Tests),
          file.path(dir_tables, "grays_test_results.csv"))

# --------------------------------------------------------------------------- #
#  3. Plot CIFs — Primary Event by Treatment
# --------------------------------------------------------------------------- #
cat("\n--- Plotting Cumulative Incidence Curves ---\n")

# Extract CIF estimates for event type 1 (primary) per treatment group
extract_cif <- function(cif_obj, event_code = 1, group_labels = NULL) {
  if (is.null(group_labels)) {
    entries <- grep(paste0(" ", event_code, "$"), names(cif_obj), value = TRUE)
  } else {
    entries <- paste(group_labels, event_code)
    entries <- intersect(entries, names(cif_obj))
  }

  dfs <- lapply(entries, function(nm) {
    df       <- data.frame(time = cif_obj[[nm]]$time,
                            est  = cif_obj[[nm]]$est,
                            var  = cif_obj[[nm]]$var)
    df$group <- nm
    df
  })
  do.call(rbind, dfs)
}

trt_labs   <- levels(actg$trt_label)
cif_df_trt <- extract_cif(cif_trt, event_code = 1, group_labels = trt_labs)

if (nrow(cif_df_trt) > 0) {
  cif_df_trt <- cif_df_trt %>%
    mutate(
      ci_low  = pmax(0, est - 1.96 * sqrt(var)),
      ci_high = pmin(1, est + 1.96 * sqrt(var)),
      group   = factor(group, levels = paste(trt_labs, 1),
                       labels = trt_labs)
    )

  p_cif_trt <- ggplot(cif_df_trt, aes(x = time, y = est, color = group, fill = group)) +
    geom_step(linewidth = 1.1) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.1, color = NA) +
    scale_color_manual(values = trt_colors) +
    scale_fill_manual(values  = trt_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_x_continuous(breaks = seq(0, 1250, 250)) +
    labs(
      title    = "Cumulative Incidence Functions by Treatment Arm",
      subtitle = "Primary event: AIDS/Death while on treatment",
      x        = "Time (Days)",
      y        = "Cumulative Incidence",
      color    = "Treatment",
      fill     = "Treatment",
      caption  = "Accounts for competing risk of off-treatment events"
    ) +
    theme(legend.position = "bottom")

  ggsave(file.path(dir_figures, "20_cif_by_treatment.png"),
         p_cif_trt, width = 10, height = 7, dpi = 300, bg = "white")
  cat("CIF plot saved.\n")
}

# --------------------------------------------------------------------------- #
#  4. Fine-Gray Subdistribution Hazard Model
# --------------------------------------------------------------------------- #
cat("\n--- Fine-Gray Subdistribution Hazard Model ---\n")

# Prepare covariates matrix for crr()
covs <- model.matrix(
  ~ trt_label + age + karnof + cd40 + cd420 +
    gender_label + race_label + str2_label + symptom_label + drugs,
  data = actg
)[, -1]  # drop intercept

fg_model <- crr(
  ftime   = actg$time,
  fstatus = as.integer(actg$event_type) - 1L,
  cov1    = covs,
  failcode = 1,
  cencode  = 0
)

cat("\nFine-Gray Model Summary:\n")
fg_summary <- summary(fg_model)
print(fg_summary)

# Extract results
fg_coef_df <- data.frame(
  Variable    = names(fg_model$coef),
  SHR         = round(exp(fg_model$coef), 3),
  SHR_CI_Low  = round(exp(fg_model$coef - 1.96 * sqrt(diag(fg_model$var))), 3),
  SHR_CI_High = round(exp(fg_model$coef + 1.96 * sqrt(diag(fg_model$var))), 3),
  P_value     = round(fg_summary$coef[, "p-value"], 4),
  stringsAsFactors = FALSE
)
fg_coef_df$Significant <- ifelse(fg_coef_df$P_value < 0.05, "Yes", "No")

write.csv(fg_coef_df,
          file.path(dir_tables, "finegray_model_results.csv"),
          row.names = FALSE)
cat("\nFine-Gray model results saved.\n")

# --------------------------------------------------------------------------- #
#  5. Forest plot — Fine-Gray SHRs
# --------------------------------------------------------------------------- #
fg_coef_df_plot <- fg_coef_df %>%
  mutate(Variable = factor(Variable, levels = Variable[order(SHR)])) %>%
  filter(!grepl("Intercept", Variable))

p_fg_forest <- ggplot(fg_coef_df_plot, aes(x = SHR, y = Variable)) +
  geom_point(aes(color = Significant), size = 3) +
  geom_errorbarh(aes(xmin = SHR_CI_Low, xmax = SHR_CI_High,
                     color = Significant),
                 height = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Yes" = "#D73027", "No" = "#92C5DE")) +
  scale_x_log10() +
  labs(
    title    = "Fine-Gray Model: Subdistribution Hazard Ratios",
    subtitle = "Competing risks — Primary event: AIDS/Death on treatment",
    x        = "Subdistribution Hazard Ratio (log scale)",
    y        = "",
    color    = "Significant (p<0.05)"
  )

ggsave(file.path(dir_figures, "21_finegray_forest_plot.png"),
       p_fg_forest, width = 10, height = 7, dpi = 300, bg = "white")

cat("=== 07: Competing Risks Analysis Complete ===\n\n")
