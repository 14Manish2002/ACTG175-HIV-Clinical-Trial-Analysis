# =============================================================================
# ACTG 175 Survival Analysis — Setup & Package Installation
# Script:  00_setup.R
# Author:  Your Name
# Date:    2024
# Purpose: Install and load all required packages; set global options
# =============================================================================

# --------------------------------------------------------------------------- #
#  1. Required packages
# --------------------------------------------------------------------------- #
required_packages <- c(
  # Core survival analysis
  "survival",       # Kaplan-Meier, Cox PH, log-rank
  "survminer",      # Publication-ready survival plots
  "cmprsk",         # Competing risks (Fine-Gray)

  # Data manipulation & display
  "tidyverse",      # dplyr, ggplot2, tidyr, readr, purrr
  "data.table",     # Fast data I/O

  # Tables
  "tableone",       # Baseline characteristics (Table 1)
  "knitr",          # kable() for markdown tables
  "kableExtra",     # Enhanced kable styling
  "gtsummary",      # Clinical summary tables

  # Visualization
  "ggplot2",        # Core plotting (loaded via tidyverse)
  "patchwork",      # Combine ggplot panels
  "RColorBrewer",   # Color palettes
  "scales",         # Axis formatting

  # Forest plots
  "forestplot",     # Forest / meta-analysis style plots

  # Reporting
  "rmarkdown",      # R Markdown report generation
  "here"            # Relative file paths
)

# --------------------------------------------------------------------------- #
#  2. Install missing packages
# --------------------------------------------------------------------------- #
cat("=== Checking and installing required packages ===\n")

missing_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE, repos = "https://cloud.r-project.org")
} else {
  cat("All packages already installed.\n")
}

# --------------------------------------------------------------------------- #
#  3. Load packages
# --------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})

cat("All packages loaded successfully.\n\n")

# --------------------------------------------------------------------------- #
#  4. Global options & theme
# --------------------------------------------------------------------------- #
options(
  scipen       = 999,     # Suppress scientific notation
  digits       = 4,       # Default decimal places
  warn         = 1,       # Warnings as they occur
  stringsAsFactors = FALSE
)

# ggplot2 global theme
theme_actg <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_text(size = base_size, color = "grey40", hjust = 0),
      plot.caption     = element_text(size = base_size - 2, color = "grey50", hjust = 1),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold", size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#f0f4f8"),
      strip.text       = element_text(face = "bold", size = base_size)
    )
}
theme_set(theme_actg())

# Color palette for 4 treatment arms
trt_colors <- c(
  "ZDV mono"   = "#E41A1C",
  "ZDV + ddI"  = "#377EB8",
  "ZDV + ddC"  = "#4DAF4A",
  "ddI mono"   = "#984EA3"
)

# --------------------------------------------------------------------------- #
#  5. Directory paths (using here package)
# --------------------------------------------------------------------------- #
dir_data    <- here::here("data")
dir_output  <- here::here("output")
dir_figures <- here::here("output", "figures")
dir_tables  <- here::here("output", "tables")

# Create output directories if they don't exist
invisible(lapply(
  c(dir_data, dir_output, dir_figures, dir_tables),
  function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
))

cat("=== Setup Complete ===\n")
cat("Project root:     ", here::here(), "\n")
cat("Data directory:   ", dir_data, "\n")
cat("Figures directory:", dir_figures, "\n")
cat("Tables directory: ", dir_tables, "\n\n")
