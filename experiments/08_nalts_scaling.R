#!/usr/bin/env Rscript
# ==============================================================================
# 08_nalts_scaling.R — Does choice set size affect out-of-sample accuracy?
#
# Story: Trains models on choice sets with N_ALTS = 10, 20, 50 alternatives
# (LA 2019, N_WOMEN=2000). Evaluates both within-the-choice-set (N_ALTS-
# specific top-1 rate) and against the full LA metro (full-metro pctile).
# The full-metro eval is N_ALTS-independent — it always ranks against all
# ~2,600 men in LA 2019 — isolating the effect of training set composition.
#
# Key changes from archive/exp_16:
#   - M1_spline updated to canonical spec: ns(age_diff,5)+ns(educ_diff,3)+same_race
#   - educ encoding uses ordinal rank (v2), not years
#   - Forest_full dropped (race dummies ruled out by 03_race_encoding)
#   - Models: M1_v1, M1_spline, Forest_core only
#
# MINNODESIZE scales with N_ALTS: mns = 50 * (N_ALTS + 1)
#   N_ALTS=10  -> mns=550
#   N_ALTS=20  -> mns=1050
#   N_ALTS=50  -> mns=2550
#
# What it tests:
#   - Does training on larger choice sets improve full-metro generalisation?
#   - Does M1_spline maintain its advantage at all N_ALTS?
#
# USAGE:
#   cd partner-choice
#   Rscript experiments/08_nalts_scaling.R
#
# EXPECTED RUNTIME: ~1-2 hours (N_ALTS=50 forest with mns=2550 is slow)
# ==============================================================================

library(tidyverse)
library(survival)
library(splines)
library(CLogitTree)

cat("=== Experiment 08: N_ALTS Scaling (LA 2019) ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
DDI_FILE    <- "data/usa_00008.xml"
METRO_CODE  <- 31080   # LA
YEAR_FILTER <- 2019
N_WOMEN     <- 2000
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.7
SEED        <- 2025
NCORES      <- 2L
NTREE       <- 30
MIN_MEN_PER_METRO <- 100

# N_ALTS grid: mns = 50 * (N_ALTS + 1)
nalts_grid <- data.frame(
  n_alts       = c(10,   20,   50),
  minnodesize  = c(550,  1050, 2550),
  tree_timeout = c(600,  600,  600)
)

cat("N_ALTS grid:\n")
print(nalts_grid)
cat("\n")

if (!dir.exists("results")) dir.create("results")

source("pipeline/00_helpers.R")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
augment_v1_cols <- function(df) {
  df %>% mutate(
    abs_age_diff  = abs(age_diff),
    older_husband = as.integer(age_diff > 0),
    abs_educ_diff = abs(educ_diff),
    hypergamy_ed  = as.integer(educ_diff > 0),
    hypogamy_ed   = as.integer(educ_diff < 0)
  )
}

add_spline_cols <- function(df, basis_age, basis_educ) {
  a <- predict(basis_age,  df$age_diff)
  e <- predict(basis_educ, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns5_",  k)]] <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

close_psock_connections <- function() {
  open_cons <- showConnections(all = FALSE)
  if (nrow(open_cons) > 0) {
    sock_idx <- which(open_cons[, "class"] == "sockconn")
    for (i in sock_idx) {
      try(close(getConnection(as.integer(rownames(open_cons)[i]))), silent = TRUE)
    }
  }
}

# ------------------------------------------------------------------------------
# Results collector
# ------------------------------------------------------------------------------
results_all <- data.frame()

# ==============================================================================
# MAIN LOOP: one iteration per N_ALTS
# ==============================================================================

for (g in seq_len(nrow(nalts_grid))) {
  N_ALTS      <- nalts_grid$n_alts[g]
  MINNODESIZE <- nalts_grid$minnodesize[g]
  TREE_TIMEOUT <- nalts_grid$tree_timeout[g]

  cat(sprintf("\n========== N_ALTS=%d  mns=%d  timeout=%ds ==========\n",
              N_ALTS, MINNODESIZE, TREE_TIMEOUT))

  # --------------------------------------------------------------------------
  # Load data for this N_ALTS setting
  # --------------------------------------------------------------------------
  t0_data <- proc.time()
  source("pipeline/01_data.R")
  data_sec <- (proc.time() - t0_data)[["elapsed"]]

  n_train <- length(train_strata)
  n_test  <- length(test_strata)
  cat(sprintf("Data: %.0f sec — %d total strata (train %d, test %d)\n",
              data_sec, n_train + n_test, n_train, n_test))

  # --------------------------------------------------------------------------
  # Prepare augmented v2 training data + spline bases
  # --------------------------------------------------------------------------
  train_v2_aug <- augment_v1_cols(train_data_v2)
  basis_age5   <- ns(train_v2_aug$age_diff,  df = 5)
  basis_educ3  <- ns(train_v2_aug$educ_diff, df = 3)
  train_v2_aug <- add_spline_cols(train_v2_aug, basis_age5, basis_educ3)

  # Also augment test data for N_ALTS-specific eval
  test_v2_aug <- augment_v1_cols(test_data_v2)
  test_v2_aug <- add_spline_cols(test_v2_aug, basis_age5, basis_educ3)

  # --------------------------------------------------------------------------
  # Fit M1_v1
  # --------------------------------------------------------------------------
  train_data_clogit$woman_race       <- factor(train_data_clogit$woman_race)
  train_data_clogit$race_partner     <- factor(train_data_clogit$race_partner)
  train_data_clogit$woman_educ_cat   <- factor(train_data_clogit$woman_educ_cat)
  train_data_clogit$educ_partner_cat <- factor(train_data_clogit$educ_partner_cat)

  test_data_clogit$woman_race       <- factor(test_data_clogit$woman_race)
  test_data_clogit$race_partner     <- factor(test_data_clogit$race_partner)
  test_data_clogit$woman_educ_cat   <- factor(test_data_clogit$woman_educ_cat)
  test_data_clogit$educ_partner_cat <- factor(test_data_clogit$educ_partner_cat)

  m1_v1 <- clogit(choice ~ abs_age_diff + older_husband + abs_educ_diff +
                    hypergamy_ed + hypogamy_ed + same_race + strata(strata),
                  data = train_data_clogit, method = "efron")
  cat(sprintf("  M1_v1 same_race coef: %.3f\n", coef(m1_v1)["same_race"]))

  # --------------------------------------------------------------------------
  # Fit M1_spline (canonical: ns(5)+ns(3)+same_race)
  # --------------------------------------------------------------------------
  m1_spline <- clogit(
    choice ~ age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
      educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
      same_race + strata(strata),
    data = train_v2_aug, method = "efron"
  )
  cat(sprintf("  M1_spline same_race coef: %.3f\n", coef(m1_spline)["same_race"]))

  # --------------------------------------------------------------------------
  # Fit Forest_core (v2_core, 5 pred, mtry=2)
  # --------------------------------------------------------------------------
  train_core <- train_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()
  test_core  <- test_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()

  cat(sprintf("  Fitting Forest_core (mns=%d)...\n", MINNODESIZE))
  set.seed(SEED)
  t0 <- proc.time()
  rf_core <- suppressWarnings(
    CLogitForest(
      data = train_core, response = "choice", s = "strata",
      ntree = NTREE, minnodesize = MINNODESIZE,
      minbucket = ceiling(MINNODESIZE / 3),
      mtry = 2, tune.mtry = FALSE,
      subsample.632 = TRUE, ncores = NCORES,
      epsilon = 1e-3, iter.max = 5, lean = TRUE,
      tree.timeout = TREE_TIMEOUT
    )
  )
  forest_sec <- (proc.time() - t0)[["elapsed"]]
  close_psock_connections()

  n_valid <- sum(sapply(rf_core$tree_list,
                        function(x) !is.null(x) && !is.null(x$tree.b)))
  cat(sprintf("  Forest: %.0f sec, %d/%d valid trees\n", forest_sec, n_valid, NTREE))

  # --------------------------------------------------------------------------
  # EVAL 1: N_ALTS-specific (rank among the N_ALTS alternatives in test sets)
  # --------------------------------------------------------------------------
  scores_nalt <- data.frame(
    strata = test_v2_aug$strata,
    choice = test_v2_aug$choice,
    m1_v1    = score_clogit(m1_v1,    test_data_clogit),
    m1_spline = score_clogit(m1_spline, test_v2_aug),
    rf_core  = as.numeric(predict(rf_core, type = "response",
                                  response = "choice", s = "strata",
                                  newdata = test_core))
  )

  ranked_nalt <- scores_nalt %>%
    group_by(strata) %>%
    mutate(
      rank_v1     = rank(-m1_v1,     ties.method = "random"),
      rank_spline = rank(-m1_spline, ties.method = "random"),
      rank_forest = rank(-rf_core,   ties.method = "random")
    ) %>%
    ungroup() %>%
    filter(choice == 1)

  top1_v1     <- mean(ranked_nalt$rank_v1     == 1)
  top1_spline <- mean(ranked_nalt$rank_spline == 1)
  top1_forest <- mean(ranked_nalt$rank_forest == 1)

  cat(sprintf("  %d-alt top1: M1_v1=%.1f%%  M1_spline=%.1f%%  Forest=%.1f%%\n",
              N_ALTS, 100*top1_v1, 100*top1_spline, 100*top1_forest))

  # --------------------------------------------------------------------------
  # EVAL 2: Full-metro (all men in LA 2019)
  # --------------------------------------------------------------------------
  test_women <- sample_women_base %>% filter(strata_id %in% test_strata)

  men_pool_metro <- men_pool %>%
    select(YEAR, MET2013,
           SERIAL_man        = SERIAL,
           PERNUM_man        = PERNUM,
           age_partner       = AGE,
           educ_partner_cat  = educ_cat,
           educ_partner_tier = educ_tier,
           educ_yrs_partner  = educ_yrs,
           race_partner      = race_refined)

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

  # v1 features
  fm_v1 <- full_metro_raw %>%
    filter(!is.na(woman_educ_tier), !is.na(educ_partner_tier)) %>%
    engineer_features() %>%
    rename(strata = strata_id)
  fm_v1_num <- fm_v1 %>% select_model_cols()

  train_v1_pred <- setdiff(names(train_data), c("choice", "strata"))
  for (col in train_v1_pred) {
    if (!col %in% names(fm_v1_num)) fm_v1_num[[col]] <- 0L
  }
  fm_v1_num <- fm_v1_num[, c("choice", "strata", train_v1_pred)]
  fm_v1_clogit <- bind_cols(
    fm_v1_num,
    fm_v1[, c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat"),
          drop = FALSE]
  )
  fm_v1_clogit$woman_race       <- factor(fm_v1_clogit$woman_race)
  fm_v1_clogit$race_partner     <- factor(fm_v1_clogit$race_partner)
  fm_v1_clogit$woman_educ_cat   <- factor(fm_v1_clogit$woman_educ_cat)
  fm_v1_clogit$educ_partner_cat <- factor(fm_v1_clogit$educ_partner_cat)

  # v2 features
  fm_v2 <- full_metro_raw %>%
    engineer_features_v2() %>%
    rename(strata = strata_id)
  fm_v2_core <- fm_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()
  fm_v2_aug <- fm_v2_core %>%
    augment_v1_cols() %>%
    add_spline_cols(basis_age5, basis_educ3)

  # Keep strata where husband survived both tracks
  keep_v1 <- fm_v1_num %>% group_by(strata) %>% filter(sum(choice) == 1) %>%
    pull(strata) %>% unique()
  keep_v2 <- fm_v2_core %>% group_by(strata) %>% filter(sum(choice) == 1) %>%
    pull(strata) %>% unique()
  keep_strata <- intersect(keep_v1, keep_v2)

  fm_v1_clogit <- fm_v1_clogit[fm_v1_clogit$strata %in% keep_strata, ] %>% as.data.frame()
  fm_v2_core   <- fm_v2_core[fm_v2_core$strata %in% keep_strata, ]    %>% as.data.frame()
  fm_v2_aug    <- fm_v2_aug[fm_v2_aug$strata %in% keep_strata, ]      %>% as.data.frame()

  n_fm_strata <- length(keep_strata)
  n_men_per   <- ave(rep(1, nrow(fm_v2_core)), fm_v2_core$strata, FUN = length)

  fm_scores <- data.frame(strata = fm_v2_core$strata, choice = fm_v2_core$choice)
  fm_scores$m1_v1     <- score_clogit(m1_v1,    fm_v1_clogit)
  fm_scores$m1_spline <- score_clogit(m1_spline, fm_v2_aug)
  fm_scores$rf_core   <- as.numeric(predict(rf_core, type = "response",
                                            response = "choice", s = "strata",
                                            newdata = fm_v2_core))

  fm_scores$rank_v1     <- ave(fm_scores$m1_v1,     fm_scores$strata,
                               FUN = function(x) rank(-x, ties.method = "random"))
  fm_scores$rank_spline <- ave(fm_scores$m1_spline, fm_scores$strata,
                               FUN = function(x) rank(-x, ties.method = "random"))
  fm_scores$rank_forest <- ave(fm_scores$rf_core,   fm_scores$strata,
                               FUN = function(x) rank(-x, ties.method = "random"))

  chosen_fm    <- fm_scores[fm_scores$choice == 1, ]
  n_men_chosen <- n_men_per[fm_scores$choice == 1]
  n_men_avg    <- round(mean(n_men_chosen))

  pctile_v1_fm     <- mean(100 * chosen_fm$rank_v1     / n_men_chosen)
  pctile_spline_fm <- mean(100 * chosen_fm$rank_spline / n_men_chosen)
  pctile_forest_fm <- mean(100 * chosen_fm$rank_forest / n_men_chosen)

  cat(sprintf("  Full-metro (%d strata, avg %d men): M1_v1=%.2f%%  M1_spline=%.2f%%  Forest=%.2f%%\n",
              n_fm_strata, n_men_avg,
              pctile_v1_fm, pctile_spline_fm, pctile_forest_fm))

  # --------------------------------------------------------------------------
  # Collect results: one row per (N_ALTS, model)
  # --------------------------------------------------------------------------
  model_rows <- list(
    list(model="M1_v1",     metro_pctile=pctile_v1_fm,     top1_nalt=top1_v1),
    list(model="M1_spline", metro_pctile=pctile_spline_fm, top1_nalt=top1_spline),
    list(model="Forest_core", metro_pctile=pctile_forest_fm, top1_nalt=top1_forest)
  )

  for (mr in model_rows) {
    results_all <- bind_rows(results_all, tibble(
      N_alts       = N_ALTS,
      mns          = MINNODESIZE,
      model        = mr$model,
      metro_pctile = round(mr$metro_pctile, 2),
      top1_nalt    = round(mr$top1_nalt, 4),
      n_fm_strata  = n_fm_strata,
      n_men_avg    = n_men_avg,
      n_valid_trees = if (mr$model == "Forest_core") n_valid else NA_integer_,
      forest_sec   = if (mr$model == "Forest_core") round(forest_sec) else NA_integer_
    ))
  }

  rm(full_metro_raw, fm_v1, fm_v1_num, fm_v1_clogit,
     fm_v2, fm_v2_core, fm_v2_aug, fm_scores,
     train_core, test_core, rf_core,
     scores_nalt, ranked_nalt)
  gc(verbose = FALSE)
}

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n\n=== RESULTS SUMMARY ===\n")
cat(sprintf("%-12s  %5s  %5s  %8s  %8s\n",
            "Model", "NALTS", "mns", "pctile", "top1_nalt"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (j in seq_len(nrow(results_all))) {
  r <- results_all[j, ]
  cat(sprintf("%-12s  %5d  %5d  %7.2f%%  %7.1f%%\n",
              r$model, r$N_alts, r$mns,
              r$metro_pctile, 100 * r$top1_nalt))
}

# Pivot to wide for a comparison table
wide <- results_all %>%
  select(N_alts, model, metro_pctile) %>%
  tidyr::pivot_wider(names_from = model, values_from = metro_pctile)

cat("\n=== Full-metro pctile by N_ALTS ===\n")
print(wide, row.names = FALSE)

# Save
write.csv(results_all, "results/08_nalts_scaling.csv", row.names = FALSE)
cat(sprintf("\nSaved: results/08_nalts_scaling.csv\n"))
cat(sprintf("Finished: %s\n", format(Sys.time())))

