#!/usr/bin/env Rscript
# ==============================================================================
# 04_spec_comparison.R — Parametric specifications vs the forest
#
# Story: M1_v1 (linear clogit) vs Forest_core reveals a gap — the forest knows
# something the linear model doesn't. We then try progressively richer
# parametric specs to close that gap. The key finding: M1_spline with df=5 on
# age and df=3 on education matches the forest on LA 2019 full-metro evaluation.
# The forest served as a nonparametric specification test.
#
# Discovery narrative order:
#   1. M1_v1      (6p)  — linear baseline, misses nonlinear structure
#   2. Forest_core       — nonparametric reference, sets the target
#   3. M1_quad    (8p)  — quadratic terms, partial improvement
#   4. M1_spline_df3 (7p) — df=3 splines, closer but not there yet
#   5. M1_spline  (9p)  — canonical df=5/df=3, matches the forest
#
# Full-metro percentile evaluation for all models.
# Output includes gap_to_forest to trace how each spec closes the gap.
#
# USAGE:
#   cd how-predictable-is-love   (repository root)
#   Rscript experiments/04_spec_comparison.R
#
# RUNTIME: ~20 min (1 forest + 4 clogits + full-metro eval)
# ==============================================================================

library(tidyverse)
library(survival)
library(splines)
library(CLogitTree)

cat("=== Experiment 04: Spec Comparison — Closing the Gap to the Forest ===\n\n")

# --- Settings ---
DDI_FILE    <- "data/usa_00008.xml"
METRO_CODE  <- 31080   # Los Angeles
YEAR_FILTER <- 2019
N_WOMEN     <- 2000
N_ALTS      <- 10
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.7
SEED        <- 2025
NTREE       <- 30
MINNODESIZE <- 550
NCORES      <- 2L
TREE_TIMEOUT <- 300
MIN_MEN_PER_METRO <- 100

cat(sprintf("Config: METRO=%d  YEAR=%d  N_WOMEN=%d  N_ALTS=%d\n",
            METRO_CODE, YEAR_FILTER, N_WOMEN, N_ALTS))

source("pipeline/00_helpers.R")

# --- Load data ---
cat("Loading data...\n")
t0 <- proc.time()
source("pipeline/01_data.R")
cat(sprintf("Data loaded: %.0f sec — %d train / %d test strata\n\n",
            (proc.time() - t0)[["elapsed"]],
            length(train_strata), length(test_strata)))

# ==============================================================================
# PREPARE TRAINING DATA
# ==============================================================================

cat("Preparing features and spline bases...\n")

# v1-style augmentation (needed for M1_v1, M1_quad)
augment_v1_cols <- function(df) {
  df %>% mutate(
    abs_age_diff  = abs(age_diff),
    older_husband = as.integer(age_diff > 0),
    abs_educ_diff = abs(educ_diff),
    hypergamy_ed  = as.integer(educ_diff > 0),
    hypogamy_ed   = as.integer(educ_diff < 0)
  )
}

# Spline bases — fit on training data only; apply same knots to test/full-metro
basis_age3  <- ns(train_data_v2$age_diff,  df = 3)
basis_educ3 <- ns(train_data_v2$educ_diff, df = 3)
basis_age5  <- ns(train_data_v2$age_diff,  df = 5)

# Column-adder for df=3/df=3 splines
add_spline_cols_df3 <- function(df, b_age3, b_educ3) {
  a <- predict(b_age3,  df$age_diff)
  e <- predict(b_educ3, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns3_",  k)]] <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

# Column-adder for df=5/df=3 splines (canonical M1_spline)
add_spline_cols_canonical <- function(df, b_age5, b_educ3) {
  a <- predict(b_age5,  df$age_diff)
  e <- predict(b_educ3, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns5_",  k)]] <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

# Build fully augmented training data
train_aug <- train_data_v2 %>%
  augment_v1_cols() %>%
  add_spline_cols_df3(basis_age3, basis_educ3) %>%
  add_spline_cols_canonical(basis_age5, basis_educ3)

test_aug <- test_data_v2 %>%
  augment_v1_cols() %>%
  add_spline_cols_df3(basis_age3, basis_educ3) %>%
  add_spline_cols_canonical(basis_age5, basis_educ3)

cat("Features ready.\n\n")

# ==============================================================================
# FIT CLOGIT MODELS (following discovery order)
# ==============================================================================

# Step 1: M1_v1 — linear baseline
cat("--- Step 1: Fitting M1_v1 (linear baseline, 6p) ---\n")
f_v1 <- choice ~ abs_age_diff + older_husband + abs_educ_diff +
  hypergamy_ed + hypogamy_ed + same_race + strata(strata)
t0 <- Sys.time()
m1_v1 <- clogit(f_v1, data = train_aug, method = "efron")
cat(sprintf("  Done: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Step 2: Forest_core — fit now so we can show the gap immediately after
cat("\n--- Step 2: Fitting Forest_core (nonparametric reference, 5p, mtry=2) ---\n")
train_core <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
test_core  <- test_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()

set.seed(SEED)
t0 <- proc.time()
suppressWarnings({
  rf_core <- CLogitForest(
    data          = train_core,
    response      = "choice",
    s             = "strata",
    ntree         = NTREE,
    minnodesize   = MINNODESIZE,
    minbucket     = ceiling(MINNODESIZE / 3),
    mtry          = 2,
    tune.mtry     = FALSE,
    subsample.632 = TRUE,
    ncores        = NCORES,
    epsilon       = 1e-3,
    iter.max      = 5,
    lean          = TRUE,
    tree.timeout  = TREE_TIMEOUT
  )
})
n_valid <- sum(sapply(rf_core$tree_list, function(tr) !is.null(tr) && !is.null(tr$tree.b)))
cat(sprintf("  Done: %.0f sec, %d/%d valid trees\n",
            (proc.time() - t0)[["elapsed"]], n_valid, NTREE))

# Step 3: M1_quad — quadratic extension of v1
cat("\n--- Step 3: Fitting M1_quad (quadratic, 8p) ---\n")
f_quad <- choice ~ abs_age_diff + I(abs_age_diff^2) + older_husband +
  abs_educ_diff + I(abs_educ_diff^2) +
  hypergamy_ed + hypogamy_ed + same_race + strata(strata)
t0 <- Sys.time()
m1_quad <- clogit(f_quad, data = train_aug, method = "efron")
cat(sprintf("  Done: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Step 4: M1_spline_df3 — splines with df=3 on both age and education
cat("\n--- Step 4: Fitting M1_spline_df3 (df=3 both, 7p) ---\n")
# basis_age3 and basis_educ3 already computed; columns are age_ns3_1..3, educ_ns3_1..3
f_spline_df3 <- choice ~
  age_ns3_1 + age_ns3_2 + age_ns3_3 +
  educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
  same_race + strata(strata)
t0 <- Sys.time()
m1_spline_df3 <- clogit(f_spline_df3, data = train_aug, method = "efron")
cat(sprintf("  Done: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Step 5: M1_spline — canonical spec df=5 age + df=3 educ
cat("\n--- Step 5: Fitting M1_spline (canonical, 9p) ---\n")
f_spline <- choice ~
  age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
  educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
  same_race + strata(strata)
t0 <- Sys.time()
m1_spline <- clogit(f_spline, data = train_aug, method = "efron")
cat(sprintf("  Done: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ==============================================================================
# FULL-METRO EVALUATION
# ==============================================================================

cat("\n--- Building full-metro choice sets ---\n")

test_women <- sample_women_base %>% filter(strata_id %in% test_strata)

men_pool_metro <- men_pool %>%
  select(
    YEAR, MET2013,
    SERIAL_man        = SERIAL,
    PERNUM_man        = PERNUM,
    age_partner       = AGE,
    educ_partner_cat  = educ_cat,
    educ_partner_tier = educ_tier,
    educ_yrs_partner  = educ_yrs,
    race_partner      = race_refined
  )

husbands <- test_women %>%
  select(YEAR, MET2013, SERIAL_man = SERIAL, PERNUM_man = SPLOC) %>%
  mutate(.is_husband = TRUE)

men_pool_tagged <- men_pool_metro %>%
  left_join(husbands, by = c("YEAR", "MET2013", "SERIAL_man", "PERNUM_man")) %>%
  mutate(.is_husband = coalesce(.is_husband, FALSE)) %>%
  select(-.is_husband)

full_metro_raw <- test_women %>%
  select(
    strata_id, YEAR, MET2013, SERIAL, SPLOC,
    woman_age       = AGE,
    woman_educ_cat  = educ_cat,
    woman_educ_tier = educ_tier,
    woman_educ_yrs  = educ_yrs,
    woman_race      = race_refined
  ) %>%
  left_join(men_pool_tagged, by = c("YEAR", "MET2013"),
            relationship = "many-to-many") %>%
  mutate(choice = ifelse(SERIAL == SERIAL_man & SPLOC == PERNUM_man, 1L, 0L))

cat(sprintf("Full-metro: %s rows\n", format(nrow(full_metro_raw), big.mark = ",")))

# Engineer v2 features and apply spline bases (same knots as training)
fm_v2_raw <- full_metro_raw %>%
  engineer_features_v2() %>%
  rename(strata = strata_id)

keep_strata <- fm_v2_raw %>%
  group_by(strata) %>%
  filter(sum(choice) == 1) %>%
  pull(strata) %>%
  unique()

fm_v2_raw <- fm_v2_raw %>% filter(strata %in% keep_strata)
cat(sprintf("Full-metro: %d strata retained (avg %.0f men)\n",
            length(keep_strata), nrow(fm_v2_raw) / length(keep_strata)))

# Augment full-metro data with all needed columns
fm_aug <- fm_v2_raw %>%
  augment_v1_cols() %>%
  add_spline_cols_df3(basis_age3, basis_educ3) %>%
  add_spline_cols_canonical(basis_age5, basis_educ3)

# Forest needs a subset without extra columns
fm_core <- fm_v2_raw %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>%
  as.data.frame()

# Score all models
cat("\n--- Scoring models on full-metro data ---\n")
fm_scores <- data.frame(strata = fm_aug$strata, choice = fm_aug$choice)
fm_scores$score_v1       <- score_clogit(m1_v1,         fm_aug)
fm_scores$score_quad     <- score_clogit(m1_quad,        fm_aug)
fm_scores$score_spl_df3  <- score_clogit(m1_spline_df3,  fm_aug)
fm_scores$score_spline   <- score_clogit(m1_spline,      fm_aug)
fm_scores$score_forest   <- as.numeric(predict(rf_core,
  type = "response", response = "choice", s = "strata", newdata = fm_core))

# Rank within strata and compute percentile
eval_fm <- function(scores_df, score_col, label, n_params) {
  scores_df %>%
    group_by(strata) %>%
    mutate(
      rank_val = rank(-get(score_col), ties.method = "random"),
      n_men    = n()
    ) %>%
    ungroup() %>%
    filter(choice == 1) %>%
    summarise(
      model        = label,
      n_params     = n_params,
      metro_pctile = mean(100 * rank_val / n_men),
      top10        = mean(rank_val <= 10),
      top50        = mean(rank_val <= 50),
      n_strata     = n()
    )
}

results <- bind_rows(
  eval_fm(fm_scores, "score_v1",       "M1_v1",          6),
  eval_fm(fm_scores, "score_forest",   "Forest_core",    NA_integer_),
  eval_fm(fm_scores, "score_quad",     "M1_quad",        8),
  eval_fm(fm_scores, "score_spl_df3",  "M1_spline_df3",  7),
  eval_fm(fm_scores, "score_spline",   "M1_spline",      9)
)

# Compute gap to forest
forest_pctile <- results$metro_pctile[results$model == "Forest_core"]
results$gap_to_forest <- round(results$metro_pctile - forest_pctile, 2)

# ==============================================================================
# PRINT DISCOVERY NARRATIVE
# ==============================================================================

cat("\n=== Full-metro Results (discovery order) ===\n")
cat(sprintf("%-18s  %5s  %8s  %10s\n", "Model", "param", "pctile", "gap_forest"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  pstr <- if (is.na(r$n_params)) " RF" else sprintf("%5d", r$n_params)
  gstr <- if (is.na(r$gap_to_forest)) "  (ref)" else sprintf("%+7.2f pp", r$gap_to_forest)
  cat(sprintf("%-18s  %5s  %6.2f%%  %10s\n",
              r$model, pstr, r$metro_pctile, gstr))
}

# Interpretation
v1_p  <- results$metro_pctile[results$model == "M1_v1"]
spl_p <- results$metro_pctile[results$model == "M1_spline"]
gap_closed <- v1_p - spl_p
gap_total  <- v1_p - forest_pctile

cat(sprintf("\nGap analysis:\n"))
cat(sprintf("  M1_v1 vs Forest:    %+.2f pp  (linear model misses this)\n",
            v1_p - forest_pctile))
cat(sprintf("  M1_spline vs Forest: %+.2f pp  (splines recover %.0f%% of the gap)\n",
            spl_p - forest_pctile,
            100 * gap_closed / gap_total))
if (spl_p < forest_pctile) {
  cat("=> M1_spline SURPASSES the forest — splines fully close and exceed the nonparametric benchmark.\n")
} else if (abs(spl_p - forest_pctile) < 0.3) {
  cat("=> M1_spline MATCHES the forest — the advantage was purely functional form.\n")
} else {
  cat("=> Gap remains; forest retains some nonparametric advantage.\n")
}

# ==============================================================================
# SAVE
# ==============================================================================

if (!dir.exists("results")) dir.create("results")

results_out <- results %>%
  select(model, n_params, metro_pctile, gap_to_forest, top10, top50, n_strata) %>%
  mutate(across(c(metro_pctile, top10, top50), ~ round(., 4)))

write.csv(results_out, "results/04_spec_comparison.csv", row.names = FALSE)
cat(sprintf("\nSaved: results/04_spec_comparison.csv\n"))
print(results_out, row.names = FALSE)
