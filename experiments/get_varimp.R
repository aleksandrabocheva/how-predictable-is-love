# get_varimp.R
# Fits small forests (20 trees, lean=FALSE) and extracts varimp() scores.
# lean=FALSE is required because varimp() needs the stored training data.
# Run from project root: Rscript experiments/get_varimp.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(ipumsr)
  library(CLogitTree)
})

DDI_FILE    <- "data/usa_00008.xml"
SEED        <- 2025L
N_WOMEN     <- 2000L
N_ALTS      <- 10L
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.70
FIT_CLOGIT  <- FALSE
FIT_FOREST  <- FALSE   # we fit manually below
MODELS      <- "none"
NTREE <- MINNODESIZE <- MINBUCKET <- DEPTH_MAX <- NCORES <- TREE_TIMEOUT <- EPSILON <- ITER_MAX <- NULL

set.seed(SEED)
source("pipeline/01_data.R")

# Build training frames (same logic as 02_model.R)
train_core <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>%
  as.data.frame()

train_full <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>%
  as.data.frame()

fit_forest <- function(data, mtry, label) {
  cat(sprintf("\n--- %s (ntree=20, lean=FALSE) ---\n", label))
  t0 <- proc.time()
  rf <- CLogitForest(
    data          = data,
    response      = "choice",
    s             = "strata",
    ntree         = 20L,
    minnodesize   = 550L,
    minbucket     = 183L,
    depth_max     = 10L,
    mtry          = mtry,
    tune.mtry     = FALSE,
    subsample.632 = TRUE,
    ncores        = 4L,
    epsilon       = 1e-3,
    iter.max      = 5L,
    lean          = FALSE,
    tree.timeout  = 120L
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  n_ok <- length(Filter(function(x) !is.null(x) && !is.null(x$tree.b), rf$tree_list))
  cat(sprintf("  Done: %.0fs, %d/20 valid trees\n", elapsed, n_ok))
  rf
}

rf_core <- fit_forest(train_core, mtry = 2L, label = "Forest_core")
rf_full <- fit_forest(train_full, mtry = 5L, label = "Forest_full")

cat("\n=== Forest_core variable importance ===\n")
vi_core <- varimp(rf_core)
print(vi_core)

cat("\n=== Forest_full variable importance ===\n")
vi_full <- varimp(rf_full)
print(vi_full)

# ── save results ───────────────────────────────────────────────────────────────
out_core <- tibble::tibble(
  Variable   = names(vi_core),
  Importance = as.numeric(vi_core),
  Model      = "Forest (core)"
)
out_full <- tibble::tibble(
  Variable   = names(vi_full),
  Importance = as.numeric(vi_full),
  Model      = "Forest (full)"
)
out <- dplyr::bind_rows(out_core, out_full) |>
  dplyr::arrange(Model, dplyr::desc(Importance))

readr::write_csv(out, "results/varimp_forest.csv")
cat("\nSaved: results/varimp_forest.csv\n")
