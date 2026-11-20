#!/usr/bin/env Rscript
# ==============================================================================
# 09_n_scaling_no2020.R — Same as 09_n_scaling.R but explicitly excludes 2020
#
# The original 09_n_scaling.R sets YEAR_FILTER <- NULL (scalar, ignored by
# 01_data.R when YEAR_FILTERS is also NULL), which means 2020 ACS data is
# included even though 2020 is officially unreliable (pandemic disruptions,
# low response rates). This script fixes that by setting YEAR_FILTERS to
# c(2016,2017,2018,2019,2021,2022,2023), which 01_data.R reads correctly.
#
# Purpose: check whether excluding 2020 meaningfully changes the n-scaling
# curves before thesis submission. Compare results/09_n_scaling_no2020.csv
# against results/09_n_scaling.csv.
#
# USAGE:
#   cd partner-choice  (project root)
#   Rscript experiments/09_n_scaling_no2020.R
#
# EXPECTED RUNTIME: ~1.5-2.5 hours (dominated by N=5000 forest)
# ==============================================================================

library(tidyverse)
library(survival)
library(splines)
library(CLogitTree)

cat("=== Experiment 09 (no2020): N Scaling (LA 2016-2019, 2021-2023) ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

# ------------------------------------------------------------------------------
# Fixed configuration
# KEY CHANGE vs 09_n_scaling.R: use YEAR_FILTERS (plural, vector) so that
# 01_data.R picks it up and excludes 2020.
# ------------------------------------------------------------------------------
DDI_FILE          <- "data/usa_00008.xml"
METRO_CODE        <- 31080   # LA
YEAR_FILTERS      <- c(2016, 2017, 2018, 2019, 2021, 2022, 2023)   # excludes 2020
YEAR_FILTER       <- NULL    # keep NULL so 01_data.R falls through to YEAR_FILTERS
N_ALTS            <- 10
DO_SPLIT          <- TRUE
TRAIN_SHARE       <- 0.7
SEED              <- 2025
NCORES            <- 2L
MIN_MEN_PER_METRO <- 100

# N_WOMEN scaling grid (identical to 09_n_scaling.R)
n_grid <- data.frame(
  n_women      = c(500,  1000, 2000, 5000),
  minnodesize  = c(110,  220,  550,  1100),
  ntree        = c(30,   30,   30,   20),
  tree_timeout = c(300,  300,  600,  1800)
)

cat("N_WOMEN grid:\n")
print(n_grid)
cat(sprintf("Year filter: %s\n\n", paste(YEAR_FILTERS, collapse = ", ")))

if (!dir.exists("results")) dir.create("results")

source("pipeline/00_helpers.R")

# ------------------------------------------------------------------------------
# Helper functions (identical to 09_n_scaling.R)
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
all_results <- vector("list", nrow(n_grid))

# ==============================================================================
# MAIN LOOP: one iteration per N_WOMEN
# ==============================================================================

for (g in seq_len(nrow(n_grid))) {
  N_WOMEN      <- n_grid$n_women[g]
  MINNODESIZE  <- n_grid$minnodesize[g]
  NTREE        <- n_grid$ntree[g]
  TREE_TIMEOUT <- n_grid$tree_timeout[g]

  cat(sprintf("\n========== N_WOMEN=%d  mns=%d  ntree=%d  timeout=%ds ==========\n",
              N_WOMEN, MINNODESIZE, NTREE, TREE_TIMEOUT))

  # --------------------------------------------------------------------------
  # Load data (YEAR_FILTERS is already set in global env; 01_data.R will use it)
  # --------------------------------------------------------------------------
  t0_data <- proc.time()
  source("pipeline/01_data.R")
  data_sec   <- (proc.time() - t0_data)[["elapsed"]]
  n_train    <- length(train_strata)
  n_test     <- length(test_strata)
  n_couples  <- nrow(sample_women_base)

  cat(sprintf("Data: %.0f sec — %d couples sampled, %d train / %d test strata\n",
              data_sec, n_couples, n_train, n_test))

  # --------------------------------------------------------------------------
  # Prepare augmented v2 training data + spline bases
  # --------------------------------------------------------------------------
  train_v2_aug <- augment_v1_cols(train_data_v2)
  basis_age5   <- ns(train_v2_aug$age_diff,  df = 5)
  basis_educ3  <- ns(train_v2_aug$educ_diff, df = 3)
  train_v2_aug <- add_spline_cols(train_v2_aug, basis_age5, basis_educ3)

  # --------------------------------------------------------------------------
  # Fit M1_v1
  # --------------------------------------------------------------------------
  train_data_clogit$woman_race       <- factor(train_data_clogit$woman_race)
  train_data_clogit$race_partner     <- factor(train_data_clogit$race_partner)
  train_data_clogit$woman_educ_cat   <- factor(train_data_clogit$woman_educ_cat)
  train_data_clogit$educ_partner_cat <- factor(train_data_clogit$educ_partner_cat)

  m1_v1 <- tryCatch(
    clogit(choice ~ abs_age_diff + older_husband + abs_educ_diff +
             hypergamy_ed + hypogamy_ed + same_race + strata(strata),
           data = train_data_clogit, method = "efron"),
    error = function(e) { cat(sprintf("  M1_v1 FAILED: %s\n", e$message)); NULL }
  )

  # --------------------------------------------------------------------------
  # Fit M1_spline (canonical: ns(5)+ns(3)+same_race)
  # --------------------------------------------------------------------------
  m1_spline <- tryCatch(
    clogit(choice ~ age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
             educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
             same_race + strata(strata),
           data = train_v2_aug, method = "efron"),
    error = function(e) { cat(sprintf("  M1_spline FAILED: %s\n", e$message)); NULL }
  )

  if (!is.null(m1_spline))
    cat(sprintf("  M1_spline same_race coef: %.3f\n", coef(m1_spline)["same_race"]))

  # --------------------------------------------------------------------------
  # Fit Forest_core (v2_core, 5 pred, mtry=2)
  # --------------------------------------------------------------------------
  train_core <- train_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()

  cat(sprintf("  Fitting Forest_core (mns=%d, ntree=%d)...\n", MINNODESIZE, NTREE))
  set.seed(SEED)
  t0 <- proc.time()
  rf_core <- tryCatch(
    suppressWarnings(
      CLogitForest(
        data = train_core, response = "choice", s = "strata",
        ntree = NTREE, minnodesize = MINNODESIZE,
        minbucket = ceiling(MINNODESIZE / 3),
        mtry = 2, tune.mtry = FALSE,
        subsample.632 = TRUE, ncores = NCORES,
        epsilon = 1e-3, iter.max = 5, lean = TRUE,
        tree.timeout = TREE_TIMEOUT
      )
    ),
    error = function(e) { cat(sprintf("  Forest FAILED: %s\n", e$message)); NULL }
  )
  forest_sec <- (proc.time() - t0)[["elapsed"]]
  close_psock_connections()

  n_valid <- 0L
  if (!is.null(rf_core)) {
    n_valid <- sum(sapply(rf_core$tree_list,
                          function(x) !is.null(x) && !is.null(x$tree.b)))
    cat(sprintf("  Forest: %.0f sec, %d/%d valid trees\n",
                forest_sec, n_valid, NTREE))
  }

  # --------------------------------------------------------------------------
  # Full-metro evaluation
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

  # Full-metro eval done year-by-year to limit memory
  years_avail <- sort(unique(test_women$YEAR))
  all_chosen_rows <- data.frame()

  train_v1_pred <- setdiff(names(train_data), c("choice", "strata"))

  for (yr in years_avail) {
    tw_yr <- test_women %>% filter(YEAR == yr)
    if (nrow(tw_yr) == 0) next

    mp_yr <- men_pool_tagged %>% filter(YEAR == yr)

    fm_raw <- tw_yr %>%
      select(strata_id, YEAR, MET2013, SERIAL, SPLOC,
             woman_age       = AGE,
             woman_educ_cat  = educ_cat,
             woman_educ_tier = educ_tier,
             woman_educ_yrs  = educ_yrs,
             woman_race      = race_refined) %>%
      left_join(mp_yr, by = c("YEAR", "MET2013"),
                relationship = "many-to-many") %>%
      mutate(choice = ifelse(SERIAL == SERIAL_man & SPLOC == PERNUM_man, 1L, 0L))

    # v1 features
    fv1 <- fm_raw %>%
      filter(!is.na(woman_educ_tier), !is.na(educ_partner_tier)) %>%
      engineer_features() %>%
      rename(strata = strata_id)
    fv1_num <- fv1 %>% select_model_cols()
    for (col in train_v1_pred) {
      if (!col %in% names(fv1_num)) fv1_num[[col]] <- 0L
    }
    fv1_num <- fv1_num[, c("choice", "strata", train_v1_pred)]
    fv1_clogit <- bind_cols(
      fv1_num,
      fv1[, c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat"),
          drop = FALSE]
    )
    fv1_clogit$woman_race       <- factor(fv1_clogit$woman_race)
    fv1_clogit$race_partner     <- factor(fv1_clogit$race_partner)
    fv1_clogit$woman_educ_cat   <- factor(fv1_clogit$woman_educ_cat)
    fv1_clogit$educ_partner_cat <- factor(fv1_clogit$educ_partner_cat)

    # v2 features
    fv2 <- fm_raw %>%
      engineer_features_v2() %>%
      rename(strata = strata_id)
    fv2_core <- fv2 %>%
      select_model_cols_v2(include_race_dummies = FALSE) %>%
      as.data.frame()
    fv2_aug <- fv2_core %>%
      augment_v1_cols() %>%
      add_spline_cols(basis_age5, basis_educ3)

    # Keep strata where husband survived both tracks
    kv1 <- fv1_num %>% group_by(strata) %>% filter(sum(choice) == 1) %>%
      pull(strata) %>% unique()
    kv2 <- fv2_core %>% group_by(strata) %>% filter(sum(choice) == 1) %>%
      pull(strata) %>% unique()
    ks  <- intersect(kv1, kv2)

    fv1_clogit <- fv1_clogit[fv1_clogit$strata %in% ks, ] %>% as.data.frame()
    fv2_core   <- fv2_core[fv2_core$strata %in% ks, ]     %>% as.data.frame()
    fv2_aug    <- fv2_aug[fv2_aug$strata %in% ks, ]       %>% as.data.frame()

    if (nrow(fv1_clogit) == 0) next

    # Score
    sc <- data.frame(strata = fv2_core$strata, choice = fv2_core$choice)

    if (!is.null(m1_v1)) {
      sc$score_v1 <- score_clogit(m1_v1, fv1_clogit)
    } else {
      sc$score_v1 <- NA_real_
    }

    if (!is.null(m1_spline)) {
      sc$score_spline <- score_clogit(m1_spline, fv2_aug)
    } else {
      sc$score_spline <- NA_real_
    }

    if (!is.null(rf_core)) {
      sc$score_forest <- as.numeric(predict(rf_core, type = "response",
                                            response = "choice", s = "strata",
                                            newdata = fv2_core))
    } else {
      sc$score_forest <- NA_real_
    }

    n_men_yr <- ave(rep(1, nrow(sc)), sc$strata, FUN = length)
    sc$rank_v1     <- ave(sc$score_v1,     sc$strata,
                          FUN = function(x) rank(-x, ties.method = "random"))
    sc$rank_spline <- ave(sc$score_spline, sc$strata,
                          FUN = function(x) rank(-x, ties.method = "random"))
    sc$rank_forest <- ave(sc$score_forest, sc$strata,
                          FUN = function(x) rank(-x, ties.method = "random"))
    sc$n_men <- n_men_yr

    all_chosen_rows <- bind_rows(all_chosen_rows, sc[sc$choice == 1, ])
    rm(fm_raw, fv1, fv1_num, fv1_clogit, fv2, fv2_core, fv2_aug, sc)
    gc(verbose = FALSE)
  }

  n_fm_strata   <- nrow(all_chosen_rows)
  n_men_avg     <- round(mean(all_chosen_rows$n_men))
  pctile_v1     <- mean(100 * all_chosen_rows$rank_v1     / all_chosen_rows$n_men, na.rm = TRUE)
  pctile_spline <- mean(100 * all_chosen_rows$rank_spline / all_chosen_rows$n_men, na.rm = TRUE)
  pctile_forest <- mean(100 * all_chosen_rows$rank_forest / all_chosen_rows$n_men, na.rm = TRUE)

  cat(sprintf(
    "  Full-metro (%d strata, avg %d men): M1_v1=%.2f%%  M1_spline=%.2f%%  Forest=%.2f%%\n",
    n_fm_strata, n_men_avg, pctile_v1, pctile_spline, pctile_forest
  ))
  cat(sprintf("  Spline advantage over V1:     %+.2f pp\n", pctile_v1 - pctile_spline))
  cat(sprintf("  Forest-spline gap:            %+.2f pp\n", pctile_forest - pctile_spline))

  all_results[[g]] <- data.frame(
    n_women          = N_WOMEN,
    train_strata     = n_train,
    M1_v1_pctile     = round(pctile_v1,     2),
    M1_spline_pctile = round(pctile_spline,  2),
    Forest_pctile    = round(pctile_forest,  2),
    spline_adv_pp    = round(pctile_v1 - pctile_spline,    2),
    forest_spline_pp = round(pctile_forest - pctile_spline, 2),
    n_fm_strata      = n_fm_strata,
    n_men_avg        = n_men_avg,
    n_valid_trees    = n_valid,
    forest_sec       = round(forest_sec),
    stringsAsFactors = FALSE
  )

  # Free memory before the next N
  rm(rf_core, train_core,
     train_v2_aug, m1_v1, m1_spline,
     train_data, test_data, train_data_clogit, test_data_clogit,
     train_data_v2, test_data_v2,
     all_chosen_rows)
  gc(verbose = FALSE)
}

# ==============================================================================
# SUMMARY
# ==============================================================================
res <- bind_rows(all_results)

cat("\n\n=== SUMMARY TABLE (no2020) ===\n")
cat(sprintf("%-7s  %6s  %8s  %8s  %8s  %8s  %8s\n",
            "N", "strata", "M1_v1", "spline", "Forest",
            "spl-v1", "for-spl"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (j in seq_len(nrow(res))) {
  r <- res[j, ]
  cat(sprintf("%-7d  %6d  %7.2f%%  %7.2f%%  %7.2f%%  %+7.2f  %+7.2f\n",
              r$n_women, r$train_strata,
              r$M1_v1_pctile, r$M1_spline_pctile, r$Forest_pctile,
              -r$spline_adv_pp,
              r$forest_spline_pp))
}
cat("\n  (pctile: lower = better; spline_adv_pp = M1_v1 - M1_spline, positive = spline wins)\n")

# Compare against original results (with 2020) if available
orig_path <- "results/09_n_scaling.csv"
if (file.exists(orig_path)) {
  orig <- read.csv(orig_path)
  cat("\n=== COMPARISON: no2020 vs with-2020 (spline pctile) ===\n")
  cat(sprintf("%-7s  %10s  %10s  %8s\n", "N", "with-2020", "no-2020", "diff"))
  cat(paste(rep("-", 45), collapse = ""), "\n")
  for (j in seq_len(nrow(res))) {
    n_w <- res$n_women[j]
    orig_row <- orig[orig$n_women == n_w, ]
    if (nrow(orig_row) == 1) {
      diff_pp <- res$M1_spline_pctile[j] - orig_row$M1_spline_pctile
      cat(sprintf("%-7d  %9.2f%%  %9.2f%%  %+7.2f pp\n",
                  n_w,
                  orig_row$M1_spline_pctile,
                  res$M1_spline_pctile[j],
                  diff_pp))
    }
  }
  cat("\n  If |diff| < 0.5 pp across all N, removing 2020 does not change conclusions.\n")
}

write.csv(res, "results/09_n_scaling_no2020.csv", row.names = FALSE)
cat(sprintf("\nSaved: results/09_n_scaling_no2020.csv\n"))
cat(sprintf("Finished: %s\n", format(Sys.time())))
