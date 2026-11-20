#!/usr/bin/env Rscript
# ==============================================================================
# 02_forest_scaling.R — Forest runtime and depth scaling with sample size
#
# WHAT IT TESTS:
#   How does CLogitForest runtime and tree depth scale with sample size N,
#   and how does MINNODESIZE control the depth/speed tradeoff?
#
# STORY:
#   Before running forests across the full experiment suite (varying N from
#   small metro-year cells to large pooled samples), we characterise the
#   computational landscape. The key question is: at each N, what MINNODESIZE
#   produces trees that (a) complete within a reasonable timeout and (b) give
#   adequate accuracy? Too small a minnodesize = deep trees that time out.
#   Too large = shallow stumps that underfit.
#
#   We test three regimes:
#     "conservative" (mns ~ 5% of N):  very shallow trees, always fast
#     "standard"     (mns ~ 11% of N): recommended depth (~10-15 strata/node)
#     "permissive"   (mns ~ 22% of N): deeper, slower, risk of timeouts at large N
#
#   This produces the calibration table: given N_WOMEN, pick minnodesize from
#   the "standard" column unless runtime is the binding constraint.
#
# HOW TO RUN:
#   cd partner-choice
#   Rscript experiments/02_forest_scaling.R
#
# EXPECTED RUNTIME: ~30-60 min total
#   3 N_WOMEN values × 3 mns values × ntree=10 ≈ 9 forest fits
#   At N=500: ~1-2 min/fit; at N=2000: ~5-10 min/fit
#
# OUTPUT:
#   results/02_forest_scaling.csv — one row per (N, mns_label, mns) with:
#     n_women, mns_label, mns, elapsed_sec, n_valid_trees, mean_splits_per_tree,
#     metro_pctile, full_metro_n_men
#
#   Prints calibration table at the end (recommended mns per N_WOMEN).
# ==============================================================================

library(tidyverse)
library(CLogitTree)

cat("=== 02_forest_scaling.R — Forest scaling calibration ===\n\n")

# ==============================================================================
# CONFIG
# ==============================================================================

DDI_FILE    <- "data/usa_00008.xml"
METRO_CODE  <- 31080      # Los Angeles (pooled across years)
N_ALTS      <- 10
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.7
SEED        <- 2025
MIN_MEN_PER_METRO <- 100

# No YEAR_FILTER — use all years pooled

N_WOMEN_VALUES <- c(500, 1000, 2000)  # sample sizes to test

# minnodesize as fraction of N_WOMEN — three regimes
# "conservative": mns ~ 5.5% of N  → very shallow trees
# "standard":     mns ~ 11% of N   → recommended 10-15 strata/node
# "permissive":   mns ~ 22% of N   → deeper, slower
MNS_FRACTIONS <- c(
  conservative = 0.055,
  standard     = 0.11,
  permissive   = 0.22
)

NTREE        <- 10        # enough for timing calibration, not for final accuracy
TREE_TIMEOUT <- 600       # 10 min/tree — generous; we want to see what happens
NCORES       <- 2L

v2_core_cols <- c("age_diff", "alter_age", "educ_diff", "alter_educ_rank", "same_race")

cat(sprintf("Metro: %d (LA, all years pooled)\n", METRO_CODE))
cat(sprintf("N_WOMEN values: %s\n", paste(N_WOMEN_VALUES, collapse = ", ")))
cat(sprintf("ntree: %d  ncores: %d  tree_timeout: %d sec\n\n",
            NTREE, NCORES, TREE_TIMEOUT))
cat("mns fractions:\n")
for (nm in names(MNS_FRACTIONS)) {
  cat(sprintf("  %-12s: %.1f%% of N\n", nm, 100 * MNS_FRACTIONS[[nm]]))
}
cat("\n")

# ==============================================================================
# HELPER: build full-metro cross-join for scoring
# ==============================================================================

build_full_metro_v2 <- function() {
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

  fm_v2 <- full_metro_raw %>%
    engineer_features_v2() %>%
    rename(strata = strata_id)

  fm_v2_core <- fm_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()

  # Keep only strata where husband appears
  keep <- fm_v2_core %>%
    group_by(strata) %>%
    filter(sum(choice) == 1) %>%
    pull(strata) %>% unique()

  fm_v2_core[fm_v2_core$strata %in% keep, ]
}

# ==============================================================================
# MAIN LOOP: N_WOMEN × mns_label
# ==============================================================================

all_results <- list()
row_idx     <- 0L

for (N_WOMEN in N_WOMEN_VALUES) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("N_WOMEN = %d\n", N_WOMEN))
  cat(sprintf("============================================================\n"))

  # --- Load data at this N ---
  cat("Loading data...\n")
  source("pipeline/00_helpers.R")
  t0 <- proc.time()
  source("pipeline/01_data.R")
  data_sec <- (proc.time() - t0)[["elapsed"]]

  n_train <- n_distinct(train_data_v2$strata)
  n_test  <- n_distinct(test_data_v2$strata)
  cat(sprintf("Data: %.0f sec — train %d strata, test %d strata\n",
              data_sec, n_train, n_test))

  # Prepare forest training/test data (v2 core, 5 predictors)
  train_core <- train_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()
  test_core <- test_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()

  # Build full-metro cross-join for this N's test set
  cat("Building full-metro cross-join...\n")
  fm_core <- build_full_metro_v2()
  n_fm_strata <- n_distinct(fm_core$strata)
  avg_men <- nrow(fm_core) / n_fm_strata
  cat(sprintf("Full-metro: %d strata, avg %.0f men/strata\n",
              n_fm_strata, avg_men))

  for (mns_label in names(MNS_FRACTIONS)) {
    mns <- round(N_WOMEN * MNS_FRACTIONS[[mns_label]] / 10) * 10L
    mns <- max(mns, 10L)  # floor at 10 to avoid degenerate cases
    mb  <- ceiling(mns / 3L)

    cat(sprintf("\n  mns_label=%-12s  mns=%d  mb=%d\n", mns_label, mns, mb))

    set.seed(SEED)
    t0 <- proc.time()
    suppressWarnings({
      rf <- CLogitForest(
        data          = train_core,
        response      = "choice",
        s             = "strata",
        ntree         = NTREE,
        minnodesize   = mns,
        minbucket     = mb,
        mtry          = 2L,
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

    # Count valid trees and extract splits per tree
    n_valid <- sum(sapply(rf$tree_list, function(tr) !is.null(tr) && !is.null(tr$tree.b)))

    splits_per_tree <- sapply(rf$tree_list, function(tr) {
      if (!is.null(tr) && !is.null(tr$tree.b) && !is.null(tr$tree.b$splits))
        nrow(tr$tree.b$splits)
      else
        0L
    })
    mean_splits <- mean(splits_per_tree)

    cat(sprintf("  Elapsed: %.0f sec, valid trees: %d/%d, mean_splits: %.1f\n",
                elapsed, n_valid, NTREE, mean_splits))

    # Full-metro evaluation
    fm_preds <- as.numeric(
      predict(rf, type = "response", response = "choice",
              s = "strata", newdata = fm_core)
    )

    fm_ranked <- fm_core %>%
      mutate(score = fm_preds) %>%
      group_by(strata) %>%
      mutate(rank = rank(-score, ties.method = "random"),
             n_men = n()) %>%
      filter(choice == 1) %>%
      ungroup()

    metro_pctile <- round(100 * mean(fm_ranked$rank / fm_ranked$n_men), 2)
    metro_top10  <- round(mean(fm_ranked$rank <= 10), 4)

    cat(sprintf("  Full-metro: pctile=%.1f%%  top10=%.1f%%  avg_men=%.0f\n",
                metro_pctile, 100 * metro_top10, mean(fm_ranked$n_men)))

    row_idx <- row_idx + 1L
    all_results[[row_idx]] <- data.frame(
      n_women             = N_WOMEN,
      mns_label           = mns_label,
      mns                 = mns,
      mns_frac_of_n       = round(mns / N_WOMEN, 3),
      elapsed_sec         = round(elapsed),
      n_valid_trees       = n_valid,
      pct_valid_trees     = round(100 * n_valid / NTREE, 1),
      mean_splits_per_tree = round(mean_splits, 1),
      metro_pctile        = metro_pctile,
      metro_top10         = metro_top10,
      full_metro_n_men    = round(mean(fm_ranked$n_men)),
      n_train_strata      = n_train,
      n_test_strata       = n_test,
      stringsAsFactors    = FALSE
    )
  }
}

# ==============================================================================
# ASSEMBLE AND SAVE
# ==============================================================================

res <- bind_rows(all_results)

cat("\n\n=== SCALING CALIBRATION TABLE ===\n\n")
cat(sprintf("%-8s  %-12s  %5s  %4s  %8s  %8s  %8s  %8s\n",
            "N_WOMEN", "mns_label", "mns", "sec", "n_valid", "splits", "pctile", "top10"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  cat(sprintf("%-8d  %-12s  %5d  %4.0f  %4d/%2d  %8.1f  %7.1f%%  %7.1f%%\n",
              r$n_women, r$mns_label, r$mns, r$elapsed_sec,
              r$n_valid_trees, NTREE,
              r$mean_splits_per_tree,
              r$metro_pctile,
              100 * r$metro_top10))
}

# ==============================================================================
# RECOMMENDATION TABLE
# ==============================================================================

cat("\n=== RECOMMENDED MINNODESIZE PER N_WOMEN ===\n\n")
cat("Rule: prefer 'standard' unless >20% of trees timed out at that N;\n")
cat("in that case fall back to 'conservative'.\n\n")

rec <- res %>%
  group_by(n_women) %>%
  group_modify(function(df, grp) {
    # prefer standard; fall back if pct_valid_trees < 80
    std  <- df %>% filter(mns_label == "standard")
    cons <- df %>% filter(mns_label == "conservative")
    perm <- df %>% filter(mns_label == "permissive")

    if (nrow(std) == 0) return(tibble())

    # Recommend standard if >= 80% valid trees; otherwise conservative
    rec_label <- if (std$pct_valid_trees >= 80) "standard" else "conservative"
    rec_row   <- df %>% filter(mns_label == rec_label) %>% head(1)

    tibble(
      recommended   = rec_label,
      rec_mns       = rec_row$mns,
      rec_pctile    = rec_row$metro_pctile,
      std_pctile    = if (nrow(std)  > 0) std$metro_pctile  else NA_real_,
      cons_pctile   = if (nrow(cons) > 0) cons$metro_pctile else NA_real_,
      perm_pctile   = if (nrow(perm) > 0) perm$metro_pctile else NA_real_,
      std_elapsed   = if (nrow(std)  > 0) std$elapsed_sec   else NA_integer_,
      std_pct_valid = if (nrow(std)  > 0) std$pct_valid_trees else NA_real_
    )
  }) %>%
  ungroup()

cat(sprintf("%-8s  %-12s  %6s  %8s  %8s  %8s  %8s\n",
            "N_WOMEN", "recommended", "mns", "pctile", "std_p%", "cons_p%", "perm_p%"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_len(nrow(rec))) {
  r <- rec[i, ]
  cat(sprintf("%-8d  %-12s  %6d  %7.1f%%  %7.1f%%  %7.1f%%  %7.1f%%\n",
              r$n_women, r$recommended, r$rec_mns, r$rec_pctile,
              r$std_pctile, r$cons_pctile, r$perm_pctile))
}

cat("\nNotes:\n")
cat("  pctile: mean percentile rank of true husband among all metro men (lower = better, 50 = random)\n")
cat("  mns: minnodesize used; ~", round(100 * MNS_FRACTIONS["standard"], 0),
    "% of N_WOMEN for 'standard'\n")
cat("  ntree=", NTREE, "(calibration only — use ntree=50+ for final accuracy)\n")
cat("  Forest_core: v2 encoding, 5 predictors, mtry=2\n\n")

cat("Copy the recommended mns column into experiment configs:\n")
for (i in seq_len(nrow(rec))) {
  cat(sprintf("  N_WOMEN=%d -> MINNODESIZE=%d (%s)\n",
              rec$n_women[i], rec$rec_mns[i], rec$recommended[i]))
}

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(res, "results/02_forest_scaling.csv", row.names = FALSE)
cat("\nSaved results/02_forest_scaling.csv\n")
cat("Done.\n")
