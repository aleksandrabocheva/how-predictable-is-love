#!/usr/bin/env Rscript
# ==============================================================================
# metro_year.R — Main prediction results: 4 metros × 2 periods
#
# THIS IS THE MAIN DELIVERABLE OF THE PROJECT.
#
# For each of 4 US metros (LA, Atlanta, San Francisco, Boston) in two time
# periods (2016–17 "pre-COVID" and 2022–23 "post-COVID"), fit and evaluate:
#   - M1_v1       : linear conditional logit (main parametric benchmark)
#   - Forest_core : conditional logit forest (main nonparametric comparison)
#   - M1_spline   : spline conditional logit (sensitivity check — shows
#                   a richer parametric form captures what the forest does)
#
# Metric reported: metro_pctile, the mean percentile rank of each woman's
# actual partner among ALL men in her metro-year. Lower = better; 50 = random.
#
# Per-run outputs (written by pipeline/run_pipe.R):
#   results/<config_id>/metrics.csv      — full per-model metrics
#   results/<config_id>/plot_10alt.png
#   results/<config_id>/plot_fullmetro.png
#   results/<config_id>/log.txt
#
# Aggregate output:
#   results/metro_year_summary.csv       — one row per (metro, period, model)
#
# USAGE:
#   cd partner-choice
#   Rscript experiments/metro_year.R
#
# Per-run caching: if a metrics.csv already exists for a config, that run is
# SKIPPED. Delete the file to force a re-run.
#
# EXPECTED RUNTIME: ~2-3 hours (eight pipeline runs, some already cached)
# ==============================================================================

library(tidyverse)

cat("=== metro_year.R — Main prediction results (4 metros × 2 periods) ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

# ------------------------------------------------------------------------------
# Configuration grid
# ------------------------------------------------------------------------------
#
# LA (met 31080) has ~4,400 couples per period; it needs MINNODESIZE=1100
# (≈110 strata/node) to keep trees within the 600s timeout.
#
# Atlanta, SF, Boston are smaller (~1,000–2,100 couples/period) — MINNODESIZE=550
# (≈55 strata/node) is fine and produces ~15–25 splits per tree.

metros <- tribble(
  ~metro_code, ~metro_name,      ~mns,  ~n_couples_approx,
  31080,       "Los Angeles",    1100L, 4400,
  12060,       "Atlanta",         550L, 2100,
  41860,       "San Francisco",   550L, 1500,
  14460,       "Boston",          550L, 1500
)

periods <- tribble(
  ~period_name,   ~year_filters,
  "pre_covid",    c(2016, 2017),
  "post_covid",   c(2022, 2023)
)

cat("Metros:\n")
print(metros)
cat("\nPeriods:\n")
print(periods)
cat("\n")

# ------------------------------------------------------------------------------
# Shared pipeline parameters
# ------------------------------------------------------------------------------
N_WOMEN_SHARED      <- 99999     # use all available couples
N_ALTS_SHARED       <- 10
NTREE_SHARED        <- 30
TREE_TIMEOUT_SHARED <- 600
NCORES_SHARED       <- 4L        # 4 PSOCK workers — safe on 8-core laptop
SEED_SHARED         <- 2025

# ------------------------------------------------------------------------------
# Run pipeline for each (metro, period)
# ------------------------------------------------------------------------------
#
# pipeline/run_pipe.R reads its parameters from the global environment, skips
# the run if metrics.csv already exists, and writes results to a subdirectory
# named by config_id.

all_configs <- expand_grid(metros, periods)

run_one <- function(metro_code, metro_name, mns, n_couples_approx,
                     period_name, year_filters) {
  cat(sprintf("\n========== %s %s (MET=%d, mns=%d, years=%s) ==========\n",
              metro_name, period_name, metro_code, mns,
              paste(year_filters, collapse = "-")))

  # Set pipeline parameters in the global env before sourcing run_pipe.R
  assign("N_WOMEN",      N_WOMEN_SHARED,      envir = .GlobalEnv)
  assign("N_ALTS",       N_ALTS_SHARED,       envir = .GlobalEnv)
  assign("NTREE",        NTREE_SHARED,        envir = .GlobalEnv)
  assign("MINNODESIZE",  mns,                 envir = .GlobalEnv)
  assign("TREE_TIMEOUT", TREE_TIMEOUT_SHARED, envir = .GlobalEnv)
  assign("NCORES",       NCORES_SHARED,       envir = .GlobalEnv)
  assign("SEED",         SEED_SHARED,         envir = .GlobalEnv)
  assign("METRO_CODE",   metro_code,          envir = .GlobalEnv)
  assign("YEAR_FILTERS", year_filters,        envir = .GlobalEnv)
  assign("YEAR_FILTER",  NULL,                envir = .GlobalEnv)
  assign("METRO_CODES",  NULL,                envir = .GlobalEnv)

  # Source the pipeline. It will skip if results already exist for this config.
  source("pipeline/run_pipe.R", local = FALSE)

  # Return the config_id so we can locate results afterwards
  sprintf("w%d_a%d_t%d_mns%d_s%d_met%d_y%s",
          N_WOMEN_SHARED, N_ALTS_SHARED, NTREE_SHARED, mns, SEED_SHARED,
          metro_code, paste(year_filters, collapse = "-"))
}

config_ids <- character(nrow(all_configs))
for (i in seq_len(nrow(all_configs))) {
  row <- all_configs[i, ]
  config_ids[i] <- run_one(
    metro_code       = row$metro_code,
    metro_name       = row$metro_name,
    mns              = row$mns,
    n_couples_approx = row$n_couples_approx,
    period_name      = row$period_name,
    year_filters     = row$year_filters[[1]]
  )
}

# ------------------------------------------------------------------------------
# Aggregate results from all 8 metrics.csv files
# ------------------------------------------------------------------------------
cat("\n========== Aggregating results ==========\n")

summary_rows <- list()
for (i in seq_len(nrow(all_configs))) {
  cfg     <- all_configs[i, ]
  metrics_path <- file.path("results", config_ids[i], "metrics.csv")

  if (!file.exists(metrics_path)) {
    cat(sprintf("  WARN: missing %s\n", metrics_path))
    next
  }

  m <- read.csv(metrics_path, stringsAsFactors = FALSE)
  m$metro_code  <- cfg$metro_code
  m$metro_name  <- cfg$metro_name
  m$period      <- cfg$period_name
  m$years       <- paste(cfg$year_filters[[1]], collapse = "-")
  m$config_id   <- config_ids[i]
  summary_rows[[i]] <- m
}

summary_df <- bind_rows(summary_rows) %>%
  select(metro_code, metro_name, period, years, model,
         metro_pctile, metro_top10, metro_top50, metro_mean_rank, metro_n_men,
         top1, top3, top5, mean_rank, config_id)

out_path <- "results/metro_year_summary.csv"
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("Saved: %s (%d rows)\n", out_path, nrow(summary_df)))

# ------------------------------------------------------------------------------
# Print wide summary table
# ------------------------------------------------------------------------------
cat("\n=== metro_pctile by metro × period × model ===\n")
wide <- summary_df %>%
  select(metro_name, period, model, metro_pctile) %>%
  pivot_wider(names_from = model, values_from = metro_pctile) %>%
  arrange(metro_name, period)

print(wide, n = Inf)

cat(sprintf("\nFinished: %s\n", format(Sys.time())))
