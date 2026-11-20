#!/usr/bin/env Rscript
# ==============================================================================
# 01_forest_tuning.R — Establish canonical forest configuration
#
# WHAT IT TESTS:
#   Consolidates exp_09/10/11 into one decision script. Before using the
#   conditional random forest as a specification test in later experiments,
#   we establish its optimal configuration on LA 2019 (a single diverse
#   metro in a single year — eliminates temporal trends and metro mixing).
#
# STORY:
#   The forest is our nonparametric reference. If the spline clogit matches
#   it, the spline captures the relevant structure. But "the forest" needs
#   a well-defined configuration. Here we answer three questions:
#
#   Part A: Which feature encoding is better — v1 (abs_diff + direction) or
#           v2 (signed_diff + alter_level)?  → confirms v2.
#
#   Part B: What is optimal mtry for Forest_core (5 v2 predictors)?
#           Scans mtry ∈ {1, 2, 3, 4, 5}.  → expected winner: mtry=2.
#
#   Part C: Does Forest_full (9 v2 predictors, adding race dummies) improve
#           over Forest_core at the mtry recommended by Part B?
#           → expected: no meaningful improvement. Race is captured by same_race.
#
#   Primary metric: full-metro percentile (lower = better, 50 = random).
#   Secondary: top-1 accuracy in 10-alt choice sets (quick sanity check).
#
# HOW TO RUN:
#   cd partner-choice
#   Rscript experiments/01_forest_tuning.R
#
# EXPECTED RUNTIME: ~25-40 min total
#   - Data load: ~2 min
#   - Part A (2 forests × ntree=30): ~8 min
#   - Part B (5 mtry × ntree=30): ~20 min
#   - Part C (1 forest): ~4 min
#   - Full-metro eval: ~2 min
#
# OUTPUT:
#   results/01_forest_tuning.csv — one row per model with:
#     experiment, model, mtry, metro_pctile, top1_10alt, forest_sec
# ==============================================================================

library(tidyverse)
library(survival)
library(CLogitTree)

cat("=== 01_forest_tuning.R — Canonical forest configuration ===\n\n")

# ==============================================================================
# CONFIG
# ==============================================================================

DDI_FILE    <- "data/usa_00008.xml"
METRO_CODE  <- 31080      # Los Angeles
YEAR_FILTER <- 2019
N_WOMEN     <- 2000
N_ALTS      <- 10
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.7
SEED        <- 2025
MIN_MEN_PER_METRO <- 100

NTREE        <- 30
MINNODESIZE  <- 550       # ~1200 couples at LA 2019, ~50 strata/node
TREE_TIMEOUT <- 300       # 5 min/tree — generous at this N
NCORES       <- 2L

mtry_scan_values <- c(1, 2, 3, 4, 5)  # Part B scan

cat(sprintf("Metro: %d (LA)  Year: %d  N_WOMEN: %d  N_ALTS: %d\n",
            METRO_CODE, YEAR_FILTER, N_WOMEN, N_ALTS))
cat(sprintf("ntree: %d  minnodesize: %d  tree_timeout: %d sec  ncores: %d\n\n",
            NTREE, MINNODESIZE, TREE_TIMEOUT, NCORES))

# ==============================================================================
# DATA LOAD (once — shared across all parts)
# ==============================================================================

cat("--- Loading data ---\n")
source("pipeline/00_helpers.R")
t0 <- proc.time()
source("pipeline/01_data.R")
data_sec <- (proc.time() - t0)[["elapsed"]]

n_train <- n_distinct(train_data_v2$strata)
n_test  <- n_distinct(test_data_v2$strata)
cat(sprintf("Data loaded: %.0f sec — train %d strata, test %d strata\n\n",
            data_sec, n_train, n_test))

# ==============================================================================
# FIT M1_v1 — clogit benchmark (fixed reference for all parts)
# ==============================================================================

cat("--- Fitting M1_v1 (clogit benchmark) ---\n")
train_data_clogit$woman_race       <- factor(train_data_clogit$woman_race)
train_data_clogit$race_partner     <- factor(train_data_clogit$race_partner)
train_data_clogit$woman_educ_cat   <- factor(train_data_clogit$woman_educ_cat)
train_data_clogit$educ_partner_cat <- factor(train_data_clogit$educ_partner_cat)

f_v1 <- choice ~ abs_age_diff + older_husband + abs_educ_diff +
  hypergamy_ed + hypogamy_ed + same_race + strata(strata)
m1_v1 <- clogit(f_v1, data = train_data_clogit, method = "efron")

# 10-alt eval for M1_v1
m1_preds <- score_clogit(m1_v1, test_data_clogit)
m1_ranks_10 <- test_data_clogit %>%
  mutate(pred = m1_preds) %>%
  group_by(strata) %>%
  mutate(rank = rank(-pred, ties.method = "average")) %>%
  filter(choice == 1) %>% ungroup()
cat(sprintf("M1_v1: mean_rank=%.2f  top1=%.1f%%\n\n",
            mean(m1_ranks_10$rank), 100 * mean(m1_ranks_10$rank == 1)))

# ==============================================================================
# HELPER: fit one forest
# ==============================================================================

fit_forest <- function(train_df, label, mtry_val,
                       ntree = NTREE, mns = MINNODESIZE) {
  cat(sprintf("  Fitting %s (mtry=%d, ntree=%d, mns=%d)...\n",
              label, mtry_val, ntree, mns))
  set.seed(SEED)
  t0 <- proc.time()
  suppressWarnings({
    rf <- CLogitForest(
      data          = as.data.frame(train_df),
      response      = "choice",
      s             = "strata",
      ntree         = ntree,
      minnodesize   = mns,
      minbucket     = ceiling(mns / 3),
      mtry          = mtry_val,
      tune.mtry     = FALSE,
      subsample.632 = TRUE,
      ncores        = NCORES,
      epsilon       = 1e-3,
      iter.max      = 5,
      lean          = TRUE,
      tree.timeout  = TREE_TIMEOUT
    )
  })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  n_valid <- sum(sapply(rf$tree_list, function(tr) !is.null(tr) && !is.null(tr$tree.b)))
  splits  <- sapply(rf$tree_list, function(tr) {
    if (!is.null(tr) && !is.null(tr$tree.b) && !is.null(tr$tree.b$splits))
      nrow(tr$tree.b$splits) else 0L
  })
  cat(sprintf("    Done: %.0f sec, %d/%d valid trees, mean_splits=%.1f\n",
              elapsed, n_valid, ntree, mean(splits)))
  list(model = rf, elapsed = elapsed, n_valid = n_valid, mean_splits = mean(splits))
}

# ==============================================================================
# HELPER: 10-alt evaluation (quick sanity check)
# ==============================================================================

eval_10alt_forest <- function(rf_obj, test_df) {
  preds <- as.numeric(predict(rf_obj, type = "response", response = "choice",
                               s = "strata", newdata = as.data.frame(test_df)))
  ranks <- test_df %>%
    mutate(pred = preds) %>%
    group_by(strata) %>%
    mutate(rank = rank(-pred, ties.method = "average")) %>%
    filter(choice == 1) %>% ungroup()
  top1 <- mean(ranks$rank == 1)
  top3 <- mean(ranks$rank <= 3)
  cat(sprintf("    10-alt: top1=%.1f%%  top3=%.1f%%  mean_rank=%.2f\n",
              100 * top1, 100 * top3, mean(ranks$rank)))
  top1
}

# ==============================================================================
# FULL-METRO EVALUATION SETUP
# (build cross-join once, reuse across all scoring calls)
# ==============================================================================

build_full_metro <- function() {
  cat("\n--- Building full-metro cross-join ---\n")

  test_women <- sample_women_base %>% filter(strata_id %in% test_strata)

  men_pool_metro <- men_pool %>%
    select(YEAR, MET2013,
           SERIAL_man       = SERIAL,
           PERNUM_man       = PERNUM,
           age_partner      = AGE,
           educ_partner_cat = educ_cat,
           educ_partner_tier = educ_tier,
           educ_yrs_partner = educ_yrs,
           race_partner     = race_refined)

  husbands <- test_women %>%
    select(YEAR, MET2013, SERIAL_man = SERIAL, PERNUM_man = SPLOC) %>%
    mutate(.is_husband = TRUE)

  men_pool_tagged <- men_pool_metro %>%
    left_join(husbands, by = c("YEAR", "MET2013", "SERIAL_man", "PERNUM_man")) %>%
    mutate(.is_husband = coalesce(.is_husband, FALSE)) %>%
    select(-.is_husband)

  full_metro_raw <- test_women %>%
    select(strata_id, YEAR, MET2013, SERIAL, SPLOC,
           woman_age       = AGE,
           woman_educ_cat  = educ_cat,
           woman_educ_tier = educ_tier,
           woman_educ_yrs  = educ_yrs,
           woman_race      = race_refined) %>%
    left_join(men_pool_tagged, by = c("YEAR", "MET2013"),
              relationship = "many-to-many") %>%
    mutate(choice = ifelse(SERIAL == SERIAL_man & SPLOC == PERNUM_man, 1L, 0L))

  cat(sprintf("Full-metro raw: %s rows (%d women x avg %.0f men)\n",
              format(nrow(full_metro_raw), big.mark = ","),
              n_distinct(full_metro_raw$strata_id),
              nrow(full_metro_raw) / n_distinct(full_metro_raw$strata_id)))

  # Engineer v1 features
  fm_v1_eng <- full_metro_raw %>%
    filter(!is.na(woman_educ_tier), !is.na(educ_partner_tier)) %>%
    engineer_features() %>%
    rename(strata = strata_id)
  fm_v1_numeric <- fm_v1_eng %>% select_model_cols()

  # Align v1 columns with training set (add missing one-hot cols as 0)
  train_v1_pred <- setdiff(names(train_data), c("choice", "strata"))
  for (col in train_v1_pred) {
    if (!col %in% names(fm_v1_numeric)) fm_v1_numeric[[col]] <- 0L
  }
  fm_v1_numeric <- fm_v1_numeric[, c("choice", "strata", train_v1_pred)]

  fm_v1_clogit <- bind_cols(
    fm_v1_numeric,
    fm_v1_eng[, c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat"),
              drop = FALSE]
  )
  fm_v1_clogit$woman_race       <- factor(fm_v1_clogit$woman_race)
  fm_v1_clogit$race_partner     <- factor(fm_v1_clogit$race_partner)
  fm_v1_clogit$woman_educ_cat   <- factor(fm_v1_clogit$woman_educ_cat)
  fm_v1_clogit$educ_partner_cat <- factor(fm_v1_clogit$educ_partner_cat)

  # Engineer v2 features
  fm_v2_eng  <- full_metro_raw %>%
    engineer_features_v2() %>%
    rename(strata = strata_id)
  fm_v2_core <- fm_v2_eng %>% select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
  fm_v2_full <- fm_v2_eng %>% select_model_cols_v2(include_race_dummies = TRUE)  %>% as.data.frame()

  # Restrict to strata where husband appears in both feature sets
  keep_v1 <- fm_v1_numeric %>% group_by(strata) %>%
    filter(sum(choice) == 1) %>% pull(strata) %>% unique()
  keep_v2 <- fm_v2_core %>% group_by(strata) %>%
    filter(sum(choice) == 1) %>% pull(strata) %>% unique()
  keep <- intersect(keep_v1, keep_v2)

  fm_v1_numeric <- fm_v1_numeric[fm_v1_numeric$strata %in% keep, ] %>% as.data.frame()
  fm_v1_clogit  <- fm_v1_clogit[fm_v1_clogit$strata   %in% keep, ] %>% as.data.frame()
  fm_v2_core    <- fm_v2_core[fm_v2_core$strata        %in% keep, ] %>% as.data.frame()
  fm_v2_full    <- fm_v2_full[fm_v2_full$strata        %in% keep, ] %>% as.data.frame()

  cat(sprintf("Full-metro: %d strata retained, avg %.0f men/strata\n",
              length(keep),
              nrow(fm_v2_core) / length(keep)))

  list(fm_v1_numeric = fm_v1_numeric,
       fm_v1_clogit  = fm_v1_clogit,
       fm_v2_core    = fm_v2_core,
       fm_v2_full    = fm_v2_full,
       keep_strata   = keep)
}

score_full_metro <- function(rf_obj, fm_df, label) {
  # fm_df must contain 'choice' and 'strata' in addition to predictor columns.
  # Predictions are aligned row-wise with fm_df (same order, same filter).
  preds <- as.numeric(predict(rf_obj, type = "response", response = "choice",
                               s = "strata", newdata = fm_df))
  fm_df %>%
    select(strata, choice) %>%
    mutate(score = preds) %>%
    group_by(strata) %>%
    mutate(rank = rank(-score, ties.method = "random"),
           n_men = n()) %>%
    filter(choice == 1) %>%
    ungroup() %>%
    summarise(
      model        = label,
      metro_pctile = round(100 * mean(rank / n_men), 2),
      metro_top10  = round(mean(rank <= 10), 4),
      metro_n_men  = round(mean(n_men))
    )
}

fm <- build_full_metro()

# M1_v1 full-metro score (fixed reference)
# fm_v1_clogit rows are aligned with fm_v2_core via intersect(keep_strata)
m1_v1_fm_preds <- score_clogit(m1_v1, fm$fm_v1_clogit)
m1_v1_fm <- fm$fm_v1_clogit %>%
  select(strata, choice) %>%
  mutate(score = m1_v1_fm_preds) %>%
  group_by(strata) %>%
  mutate(rank = rank(-score, ties.method = "random"), n_men = n()) %>%
  filter(choice == 1) %>%
  ungroup() %>%
  summarise(
    model        = "M1_v1",
    metro_pctile = round(100 * mean(rank / n_men), 2),
    metro_top10  = round(mean(rank <= 10), 4),
    metro_n_men  = round(mean(n_men))
  )

cat(sprintf("M1_v1 full-metro: pctile=%.1f%%  top10=%.1f%%\n",
            m1_v1_fm$metro_pctile, 100 * m1_v1_fm$metro_top10))

# Prepare train/test column subsets
main6_cols   <- c("abs_age_diff", "older_husband", "abs_educ_diff",
                  "hypergamy_ed", "hypogamy_ed", "same_race")
v2_core_cols <- c("age_diff", "alter_age", "educ_diff", "alter_educ_rank", "same_race")
v2_full_cols <- c(v2_core_cols, "alter_race_white", "alter_race_black",
                  "alter_race_hispanic", "alter_race_asian")

train_v1_sub <- train_data[, c("choice", "strata", main6_cols)]
test_v1_sub  <- test_data[, c("choice", "strata", main6_cols)]

train_core <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
test_core  <- test_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()

train_full <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>% as.data.frame()
test_full  <- test_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>% as.data.frame()

all_results <- list()

# ==============================================================================
# PART A: v1 encoding (main6) vs v2 encoding (core 5-pred), mtry=2
# ==============================================================================

cat("\n======================================\n")
cat("PART A: v1 vs v2 encoding (mtry=2)\n")
cat("======================================\n")

cat("Fitting Forest_v1 (main6, 6 pred, mtry=2)...\n")
rf_v1 <- fit_forest(train_v1_sub, "Forest_v1", mtry_val = 2)
top1_v1 <- eval_10alt_forest(rf_v1$model, test_v1_sub)
fm_v1_res <- score_full_metro(rf_v1$model,
                               fm$fm_v1_numeric[, c("choice", "strata", main6_cols)],
                               "Forest_v1")
cat(sprintf("    Full-metro: pctile=%.1f%%  top10=%.1f%%\n",
            fm_v1_res$metro_pctile, 100 * fm_v1_res$metro_top10))

all_results[["A_v1"]] <- data.frame(
  experiment   = "A",
  model        = "Forest_v1",
  mtry         = 2L,
  metro_pctile = fm_v1_res$metro_pctile,
  top1_10alt   = round(top1_v1, 4),
  forest_sec   = round(rf_v1$elapsed),
  n_valid      = rf_v1$n_valid,
  mean_splits  = round(rf_v1$mean_splits, 1),
  stringsAsFactors = FALSE
)

cat("\nFitting Forest_core (v2, 5 pred, mtry=2)...\n")
rf_core_A <- fit_forest(train_core, "Forest_core_A", mtry_val = 2)
top1_core_A <- eval_10alt_forest(rf_core_A$model, test_core)
fm_core_A_res <- score_full_metro(rf_core_A$model, fm$fm_v2_core, "Forest_core_A")
cat(sprintf("    Full-metro: pctile=%.1f%%  top10=%.1f%%\n",
            fm_core_A_res$metro_pctile, 100 * fm_core_A_res$metro_top10))

all_results[["A_core"]] <- data.frame(
  experiment   = "A",
  model        = "Forest_core",
  mtry         = 2L,
  metro_pctile = fm_core_A_res$metro_pctile,
  top1_10alt   = round(top1_core_A, 4),
  forest_sec   = round(rf_core_A$elapsed),
  n_valid      = rf_core_A$n_valid,
  mean_splits  = round(rf_core_A$mean_splits, 1),
  stringsAsFactors = FALSE
)

cat(sprintf("\nPart A summary: Forest_v1 pctile=%.1f%%  Forest_core pctile=%.1f%%\n",
            fm_v1_res$metro_pctile, fm_core_A_res$metro_pctile))
cat(sprintf("Encoding advantage (v2 - v1): %+.1f pp\n",
            fm_v1_res$metro_pctile - fm_core_A_res$metro_pctile))

# ==============================================================================
# PART B: mtry scan for Forest_core (v2, 5 predictors)
# ==============================================================================

cat("\n======================================\n")
cat("PART B: mtry scan for Forest_core (v2, 5 pred)\n")
cat(sprintf("mtry values: %s\n", paste(mtry_scan_values, collapse = ", ")))
cat("======================================\n")

for (mtry_i in seq_along(mtry_scan_values)) {
  mt <- mtry_scan_values[mtry_i]
  cat(sprintf("\nmtry=%d (%d/%d)\n", mt, mtry_i, length(mtry_scan_values)))

  rf_b <- fit_forest(train_core, paste0("Forest_core_mtry", mt), mtry_val = mt)
  top1_b <- eval_10alt_forest(rf_b$model, test_core)
  fm_b_res <- score_full_metro(rf_b$model, fm$fm_v2_core,
                                paste0("Forest_core_mtry", mt))
  cat(sprintf("    Full-metro: pctile=%.1f%%  top10=%.1f%%\n",
              fm_b_res$metro_pctile, 100 * fm_b_res$metro_top10))

  all_results[[paste0("B_mtry", mt)]] <- data.frame(
    experiment   = "B",
    model        = paste0("Forest_core_mtry", mt),
    mtry         = mt,
    metro_pctile = fm_b_res$metro_pctile,
    top1_10alt   = round(top1_b, 4),
    forest_sec   = round(rf_b$elapsed),
    n_valid      = rf_b$n_valid,
    mean_splits  = round(rf_b$mean_splits, 1),
    stringsAsFactors = FALSE
  )
}

# Identify best mtry from Part B
part_B <- bind_rows(all_results[grepl("^B_", names(all_results))])
best_mtry_idx <- which.min(part_B$metro_pctile)
best_mtry     <- part_B$mtry[best_mtry_idx]
cat(sprintf("\nPart B best mtry: %d (pctile=%.1f%%)\n",
            best_mtry, part_B$metro_pctile[best_mtry_idx]))

# ==============================================================================
# PART C: Forest_full (9 pred, race dummies) vs Forest_core at best mtry
# ==============================================================================

cat("\n======================================\n")
cat(sprintf("PART C: Forest_full (9 pred, mtry=5) vs Forest_core (mtry=%d)\n", best_mtry))
cat("======================================\n")

# Re-run Forest_core with best mtry from Part B for a clean comparison
cat(sprintf("Fitting Forest_core (best mtry=%d)...\n", best_mtry))
rf_core_C <- fit_forest(train_core, "Forest_core_best", mtry_val = best_mtry)
top1_core_C <- eval_10alt_forest(rf_core_C$model, test_core)
fm_core_C_res <- score_full_metro(rf_core_C$model, fm$fm_v2_core, "Forest_core_best")
cat(sprintf("    Full-metro: pctile=%.1f%%  top10=%.1f%%\n",
            fm_core_C_res$metro_pctile, 100 * fm_core_C_res$metro_top10))

all_results[["C_core"]] <- data.frame(
  experiment   = "C",
  model        = "Forest_core",
  mtry         = best_mtry,
  metro_pctile = fm_core_C_res$metro_pctile,
  top1_10alt   = round(top1_core_C, 4),
  forest_sec   = round(rf_core_C$elapsed),
  n_valid      = rf_core_C$n_valid,
  mean_splits  = round(rf_core_C$mean_splits, 1),
  stringsAsFactors = FALSE
)

cat("Fitting Forest_full (9 pred, mtry=5)...\n")
rf_full_C <- fit_forest(train_full, "Forest_full", mtry_val = 5)
top1_full_C <- eval_10alt_forest(rf_full_C$model, test_full)
fm_full_C_res <- score_full_metro(rf_full_C$model, fm$fm_v2_full, "Forest_full")
cat(sprintf("    Full-metro: pctile=%.1f%%  top10=%.1f%%\n",
            fm_full_C_res$metro_pctile, 100 * fm_full_C_res$metro_top10))

all_results[["C_full"]] <- data.frame(
  experiment   = "C",
  model        = "Forest_full",
  mtry         = 5L,
  metro_pctile = fm_full_C_res$metro_pctile,
  top1_10alt   = round(top1_full_C, 4),
  forest_sec   = round(rf_full_C$elapsed),
  n_valid      = rf_full_C$n_valid,
  mean_splits  = round(rf_full_C$mean_splits, 1),
  stringsAsFactors = FALSE
)

cat(sprintf("\nPart C: Forest_core pctile=%.1f%%  Forest_full pctile=%.1f%%\n",
            fm_core_C_res$metro_pctile, fm_full_C_res$metro_pctile))
cat(sprintf("Race dummies advantage: %+.1f pp\n",
            fm_core_C_res$metro_pctile - fm_full_C_res$metro_pctile))

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

res <- bind_rows(all_results)

# Add M1_v1 reference row
m1_row <- data.frame(
  experiment   = "ref",
  model        = "M1_v1",
  mtry         = NA_integer_,
  metro_pctile = m1_v1_fm$metro_pctile,
  top1_10alt   = round(mean(m1_ranks_10$rank == 1), 4),
  forest_sec   = NA_integer_,
  n_valid      = NA_integer_,
  mean_splits  = NA_real_,
  stringsAsFactors = FALSE
)
res <- bind_rows(res, m1_row)

cat("\n\n=== FULL RESULTS TABLE ===\n\n")
print(res[, c("experiment", "model", "mtry", "metro_pctile", "top1_10alt",
              "forest_sec", "n_valid", "mean_splits")],
      row.names = FALSE)

cat("\n=== CONCLUSIONS ===\n")
# Part A
v1_pctile   <- all_results[["A_v1"]]$metro_pctile
core_pctile <- all_results[["A_core"]]$metro_pctile
cat(sprintf("Part A — Encoding: v2 %s v1 by %.1f pp (pctile: %.1f%% vs %.1f%%)\n",
            if (core_pctile < v1_pctile) "beats" else "does NOT beat",
            abs(v1_pctile - core_pctile), core_pctile, v1_pctile))
cat(sprintf("  => Canonical encoding: %s\n",
            if (core_pctile < v1_pctile) "v2 (confirmed)" else "v1 (unexpected)"))

cat(sprintf("Part B — Best mtry for Forest_core: %d (pctile=%.1f%%)\n",
            best_mtry, part_B$metro_pctile[best_mtry_idx]))
cat(sprintf("  => Canonical mtry: %d\n", best_mtry))

full_pctile <- all_results[["C_full"]]$metro_pctile
core_pctile_C <- all_results[["C_core"]]$metro_pctile
cat(sprintf("Part C — Race dummies: Forest_full %s Forest_core by %.1f pp\n",
            if (full_pctile < core_pctile_C) "beats" else "does NOT improve on",
            abs(core_pctile_C - full_pctile)))
cat(sprintf("  => %s\n",
            if (full_pctile < core_pctile_C - 0.5)
              "Consider Forest_full for race-sensitive experiments"
            else
              "Forest_core (v2 encoding, mtry=2) confirmed as canonical reference"))

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(res, "results/01_forest_tuning.csv", row.names = FALSE)
cat("\nSaved results/01_forest_tuning.csv\n")
cat("Done.\n")
