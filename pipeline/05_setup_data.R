# ==============================================================================
# 05_setup_data.R — Load data only (no model fitting)
#
# Use this instead of source("pipeline/run_pipe.R") when you only need the
# data objects in memory — e.g. before running 05_stratified_analysis.R.
# Skips the slow model fitting (02_model.R) entirely.
#
# Usage:
#   source("pipeline/05_setup_data.R")          # all metros (default)
#   METRO_CODE <- 35620; source("pipeline/05_setup_data.R")  # NY only
# ==============================================================================

# Run from the repository root

# ── Settings (override before sourcing if needed) ──────────────────────────────
if (!exists("METRO_CODE"))       METRO_CODE        <- NULL
if (!exists("METRO_CODES"))      METRO_CODES       <- NULL
if (!exists("YEAR_FILTER"))      YEAR_FILTER       <- NULL
if (!exists("YEAR_FILTERS"))     YEAR_FILTERS      <- NULL
if (!exists("N_WOMEN"))          N_WOMEN           <- 2000L
if (!exists("N_ALTS"))           N_ALTS            <- 10L
if (!exists("SEED"))             SEED              <- 2025L
if (!exists("DDI_FILE"))         DDI_FILE          <- "data/usa_00008.xml"
if (!exists("DO_SPLIT"))         DO_SPLIT          <- TRUE
if (!exists("TRAIN_SHARE"))      TRAIN_SHARE       <- 0.7
if (!exists("MIN_MEN_PER_METRO")) MIN_MEN_PER_METRO <- 100L
if (!exists("NTREE"))            NTREE             <- 50L
if (!exists("MINNODESIZE"))      MINNODESIZE       <- 550L
if (!exists("MINBUCKET"))        MINBUCKET         <- 184L
if (!exists("DEPTH_MAX"))        DEPTH_MAX         <- 10L
if (!exists("NCORES"))           NCORES            <- 2L
if (!exists("TREE_TIMEOUT"))     TREE_TIMEOUT      <- 600L
if (!exists("EPSILON"))          EPSILON           <- 1e-3
if (!exists("ITER_MAX"))         ITER_MAX          <- 5L

# ── Load helpers ───────────────────────────────────────────────────────────────
source("pipeline/00_helpers.R")

# add_spline_cols() is normally defined in 02_model.R — define it here so
# 05_stratified_analysis.R can use it without sourcing 02_model.R
add_spline_cols <- function(df, basis_age, basis_educ) {
  a <- predict(basis_age,  df$age_diff)
  e <- predict(basis_educ, df$educ_diff)
  for (k in seq_len(ncol(a))) df[[paste0("age_ns5_",  k)]] <- a[, k]
  for (k in seq_len(ncol(e))) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

# ── Load data only (01_data.R) ─────────────────────────────────────────────────
cat("==============================================================================\n")
cat(sprintf("DATA SETUP — Metro: %s | N_WOMEN: %d | N_ALTS: %d | SEED: %d\n",
            if (!is.null(METRO_CODE)) as.character(METRO_CODE) else "ALL",
            N_WOMEN, N_ALTS, SEED))
cat("==============================================================================\n")

source("pipeline/01_data.R")

cat("\n==============================================================================\n")
cat("Data loaded. Objects ready: train_data_clogit, test_data_clogit,\n")
cat("                            train_data_v2, test_data_v2\n")
cat("Now run: source('discovery/05_stratified_analysis.R')\n")
cat("==============================================================================\n")
