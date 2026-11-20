# ==============================================================================
# run_pipe.R — Master pipeline, single-run entry point
#
# USAGE:
#   Rscript pipeline/run_pipe.R
#   Rscript pipeline/run_pipe.R N_WOMEN=5000 NTREE=50 MINNODESIZE=1100
#   Rscript pipeline/run_pipe.R METRO_CODE=31080 YEAR_FILTER=2019
#
# OUTPUT (all under results/):
#   results/run_<config_id>_metrics.csv  — accuracy metrics
#   results/run_<config_id>_plot.png     — cumulative accuracy plot
#   results/log_<config_id>.txt          — console log
#
# If a metrics file already exists, the run is skipped. Delete it to re-run.
# ==============================================================================

# --- Ensure personal R library is on the path (needed on Tetralith cluster) ---
# On the cluster, R packages are installed to ~/R/packages (not the system library).
# This line tells R to also look there when loading packages with library().
user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

# ==============================================================================
# DEFAULTS
# ==============================================================================
# These are only set if not already defined in the R environment. This allows
# other scripts (like run_scaling.sh) to pre-set values before source()'ing
# this file, while still providing sensible defaults for standalone use.

# Helper: assign a default value to a variable only if it doesn't already exist
d <- function(name, value) { if (!exists(name, envir = .GlobalEnv)) assign(name, value, envir = .GlobalEnv) }

d("N_WOMEN",          2000)
d("N_ALTS",           10)
d("NTREE",            50)
d("MINNODESIZE",      550)       # min obs in node (ROWS, not strata; 550 = 55 strata at N_ALTS=10)
d("MINBUCKET",        184)       # min obs in any leaf
d("DEPTH_MAX",        10)        # max tree depth (0 = unlimited)
d("NCORES",           NULL)      # NULL = auto-detect
d("TREE_TIMEOUT",     600)       # per-tree timeout in seconds
d("MIN_MEN_PER_METRO", 100)      # exclude thin metros
d("METRO_CODE",       NULL)      # single metro filter (e.g. 31080 for LA) 35620 for NY
d("METRO_CODES",      NULL)      # vector of metro codes
d("YEAR_FILTER",      NULL)      # single year filter (e.g. 2019)
d("YEAR_FILTERS",     NULL)      # vector of years (e.g. c(2016, 2017))

d("SEED",             2025)
d("DDI_FILE",         "data/usa_00008.xml")
d("DO_SPLIT",         TRUE)
d("TRAIN_SHARE",      0.7)

d("EPSILON",          1e-3)
d("ITER_MAX",         5)

# Which model families to fit: "clogit" (M1_v1 + M1_spline), "forest"
# (Forest_core + Forest_full), or "both" (the default). Useful when you
# only need the fast clogits (skip the ~minutes-to-hours forest run) or
# only the forests (skip the clogit coefficients).
d("MODELS",           "both")

rm(d)  # clean up helper

# ==============================================================================
# CLI OVERRIDES
# ==============================================================================
# Any KEY=VALUE pairs passed on the command line override the defaults above.
# Example: Rscript pipeline/run_pipe.R N_WOMEN=5000 NTREE=50 METRO_CODE=31080
# This parses each argument, tries to convert to numeric, then boolean/NULL,
# and falls back to string. Only existing variables can be overridden.

cli_args <- commandArgs(trailingOnly = TRUE)
for (arg in cli_args) {
  if (grepl("=", arg)) {
    key <- sub("=.*", "", arg)
    val <- sub(".*=", "", arg)
    parsed <- suppressWarnings(as.numeric(val))
    if (is.na(parsed)) {
      parsed <- switch(toupper(val),
        "TRUE"  = TRUE,
        "FALSE" = FALSE,
        "NULL"  = NULL,
        val
      )
    }
    if (exists(key)) {
      assign(key, parsed)
      cat(sprintf("  CLI override: %s = %s\n", key, val))
    } else {
      warning(sprintf("Unknown parameter: %s (ignored)", key))
    }
  }
}

# ==============================================================================
# AUTO-DETECT CORES (if not set explicitly)
# ==============================================================================
# Mac: 2 cores (Apple Silicon has issues with more PSOCK workers)
# Cluster: reads Slurm's allocation, but NOTE this reports the full node (32),
#          not your --cpus-per-task. For large N, override NCORES explicitly.
# Other: uses all available cores

if (is.null(NCORES)) {
  NCORES <- if (Sys.info()[["sysname"]] == "Darwin") {
    2L  # safe on Apple Silicon
  } else if (nzchar(Sys.getenv("SLURM_CPUS_ON_NODE"))) {
    as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
  } else {
    parallel::detectCores()
  }
}
NCORES <- as.integer(NCORES)

# ==============================================================================
# MODEL SELECTION
# ==============================================================================
# Validate MODELS and derive boolean flags used by 02_model.R / 03_evaluation.R.

if (!MODELS %in% c("clogit", "forest", "both")) {
  stop(sprintf("MODELS must be 'clogit', 'forest', or 'both' (got '%s')", MODELS))
}
FIT_CLOGIT <- MODELS %in% c("clogit", "both")
FIT_FOREST <- MODELS %in% c("forest", "both")

# ==============================================================================
# CONFIG ID
# ==============================================================================

config_id <- sprintf("w%d_a%d_t%d_mns%d_s%d",
                     N_WOMEN, N_ALTS, NTREE, MINNODESIZE, SEED)
if (!is.null(METRO_CODES)) {
  config_id <- paste0(config_id, "_div", length(METRO_CODES))  # e.g. _div6
} else if (!is.null(METRO_CODE)) {
  config_id <- paste0(config_id, "_met", METRO_CODE)
}
if (!is.null(YEAR_FILTERS)) {
  config_id <- paste0(config_id, "_y", paste(YEAR_FILTERS, collapse = "-"))
} else if (!is.null(YEAR_FILTER)) {
  config_id <- paste0(config_id, "_y", YEAR_FILTER)
}
# If only one model family is fitted, mark the output folder so it doesn't
# collide with (or get cache-skipped by) a previous "both" run.
if (MODELS != "both") {
  config_id <- paste0(config_id, "_", MODELS)
}

# Each run gets its own subfolder: results/<config_id>/
results_dir <- file.path("results", config_id)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
results_path <- file.path(results_dir, "metrics.csv")
log_path     <- file.path(results_dir, "log.txt")

# --- Skip if already computed (check BEFORE opening log to avoid overwriting) ---
if (file.exists(results_path)) {
  cat("Result already exists:", results_path, "\n")
  cat("Delete the file to force a re-run.\n")
} else {

  # ==============================================================================
  # CONSOLE LOG — capture all printed output to a file
  # ==============================================================================
  # When running non-interactively (e.g., via Rscript or Slurm), all cat() and
  # message() output is duplicated to a log file alongside the results.
  # sink(split=TRUE) means output goes to BOTH the log file AND the console.
  if (!interactive()) {
    log_con <- file(log_path, open = "wt")
    sink(log_con, type = "output", split = TRUE)
    sink(log_con, type = "message", append = TRUE)
    on.exit({
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"),  silent = TRUE)
      try(close(log_con),         silent = TRUE)
    }, add = TRUE)
  }

  cat("==============================================================================\n")
  cat(paste0("RUN — ", config_id, "\n"))
  cat(paste0("Time: ", Sys.time(), "\n"))
  if (nzchar(Sys.getenv("SLURM_JOB_ID"))) {
    cat(paste0("Slurm job: ", Sys.getenv("SLURM_JOB_ID"),
               "  array task: ", Sys.getenv("SLURM_ARRAY_TASK_ID", "N/A"), "\n"))
  }
  cat("==============================================================================\n")
  cat(sprintf(
    "  N_WOMEN=%d  N_ALTS=%d  NTREE=%d  MINNODESIZE=%d  NCORES=%d  SEED=%d  MODELS=%s\n",
    N_WOMEN, N_ALTS, NTREE, MINNODESIZE, NCORES, SEED, MODELS
  ))

  t_start <- proc.time()

  cat("\n=== 01 DATA ===\n")
  source("pipeline/01_data.R")

  cat("\n=== 02 MODEL ===\n")
  source("pipeline/02_model.R")

  cat("\n=== 03 EVALUATION ===\n")
  source("pipeline/03_evaluation.R")

  t_total <- (proc.time() - t_start)[["elapsed"]]
  cat(sprintf("\n=== DONE === Total time: %.1f min\n", t_total / 60))
  cat("Results saved to:", results_path, "\n")
}
