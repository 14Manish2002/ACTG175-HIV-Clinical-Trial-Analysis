# =============================================================================
# ACTG 175 Survival Analysis — Data Preprocessing
# Script:  01_data_preprocessing.R
# Purpose: Load raw data, clean, label, and engineer features
# =============================================================================

source("R/00_setup.R")

cat("=== 01: Data Preprocessing ===\n")

# --------------------------------------------------------------------------- #
#  1. Load raw data
# --------------------------------------------------------------------------- #
raw_data <- read.csv(
  file.path(dir_data, "AIDS_ClinicalTrial_GroupStudy175.csv"),
  stringsAsFactors = FALSE
)

cat(sprintf("Raw data loaded: %d rows × %d columns\n", nrow(raw_data), ncol(raw_data)))

# --------------------------------------------------------------------------- #
#  2. Variable labeling & factor creation
# --------------------------------------------------------------------------- #
actg <- raw_data %>%
  mutate(
    # ---- Treatment ----
    trt_label = factor(trt,
      levels = 0:3,
      labels = c("ZDV mono", "ZDV + ddI", "ZDV + ddC", "ddI mono")
    ),
    trt_group = factor(
      ifelse(trt == 0, "Monotherapy (ZDV)", "Combination / Alternative"),
      levels = c("Monotherapy (ZDV)", "Combination / Alternative")
    ),

    # ---- Demographics ----
    gender_label = factor(gender, levels = 0:1, labels = c("Female", "Male")),
    race_label   = factor(race,   levels = 0:1, labels = c("White", "Non-White")),

    # ---- Clinical ----
    hemo_label    = factor(hemo,    levels = 0:1, labels = c("No", "Yes")),
    homo_label    = factor(homo,    levels = 0:1, labels = c("No", "Yes")),
    drugs_label   = factor(drugs,   levels = 0:1, labels = c("No", "Yes")),
    symptom_label = factor(symptom, levels = 0:1, labels = c("Asymptomatic", "Symptomatic")),
    offtrt_label  = factor(offtrt,  levels = 0:1, labels = c("On Treatment", "Off Treatment")),
    str2_label    = factor(str2,    levels = 0:1, labels = c("Naive", "Experienced")),
    strat_label   = factor(strat,
      levels  = 1:3,
      labels  = c("Anti-HIV naive", "≤52 wks ART", ">52 wks ART")
    ),

    # ---- Age groups ----
    age_group = cut(age,
      breaks = c(0, 30, 40, 50, 100),
      labels = c("<30", "30-39", "40-49", "50+"),
      right  = FALSE
    ),

    # ---- CD4 change ----
    cd4_change    = cd420 - cd40,
    cd4_pct_change = ifelse(cd40 > 0, (cd420 - cd40) / cd40 * 100, NA),
    cd4_response   = factor(
      ifelse(cd420 >= cd40, "Improved/Stable", "Declined"),
      levels = c("Improved/Stable", "Declined")
    ),

    # ---- CD8 change ----
    cd8_change = cd820 - cd80,

    # ---- Log transformations for skewed variables ----
    log_cd40  = log1p(cd40),
    log_cd420 = log1p(cd420),
    log_preanti = log1p(preanti),

    # ---- Karnofsky categories ----
    karnof_cat = cut(karnof,
      breaks = c(0, 70, 80, 90, 100),
      labels = c("≤70 (Poor)", "80 (Moderate)", "90 (Good)", "100 (Excellent)"),
      right  = TRUE, include.lowest = TRUE
    )
  )

cat(sprintf("Processed dataset: %d rows × %d columns\n", nrow(actg), ncol(actg)))

# --------------------------------------------------------------------------- #
#  3. Data quality checks
# --------------------------------------------------------------------------- #
cat("\n--- Missing Value Summary ---\n")
missing_summary <- sapply(actg, function(x) sum(is.na(x)))
missing_summary <- missing_summary[missing_summary > 0]
if (length(missing_summary) == 0) {
  cat("No missing values found in any variable.\n")
} else {
  print(missing_summary)
}

cat("\n--- Event Distribution ---\n")
event_tab <- table(Event = actg$label)
cat("Censored (0):", event_tab["0"],
    sprintf("(%.1f%%)\n", 100 * event_tab["0"] / nrow(actg)))
cat("Event    (1):", event_tab["1"],
    sprintf("(%.1f%%)\n", 100 * event_tab["1"] / nrow(actg)))

cat("\n--- Treatment Group Counts ---\n")
print(table(actg$trt_label))

cat("\n--- Time range (days) ---\n")
cat(sprintf("Min: %d  |  Median: %.0f  |  Max: %d\n",
            min(actg$time), median(actg$time), max(actg$time)))

# --------------------------------------------------------------------------- #
#  4. Save processed data
# --------------------------------------------------------------------------- #
saveRDS(actg, file = file.path(dir_data, "actg_processed.rds"))
cat("\nProcessed data saved to data/actg_processed.rds\n")
cat("=== 01: Preprocessing Complete ===\n\n")
