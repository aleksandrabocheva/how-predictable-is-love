#!/usr/bin/env Rscript
# ==============================================================================
# 03_race_encoding.R — Does encoding race as dummies add signal beyond same_race?
#
# Story: The forest can in principle learn race-specific preferences beyond a
# same_race binary (e.g., "White women avoid Black men more than other non-White
# men"). We test whether explicit dummies for the partner's race
# (White/Black/Hispanic/Asian) add predictive signal beyond same_race, both in
# the forest and in clogit. The answer: they don't. same_race is sufficient.
#
# Models compared:
#   Forest_core    (5 pred, same_race only, mtry=2) — canonical reference
#   Forest_full    (9 pred, + race dummies,  mtry=5)
#   Forest_norace  (8 pred, race dummies but NO same_race, mtry=5)
#   M1_spline      (9 params: ns(age_diff,5)+ns(educ_diff,3)+same_race)
#   M1_spline_race (13 params: M1_spline + 4 race dummies)
#
# Full-metro percentile evaluation for all 5 models.
# Coefficient table for both clogits showing race dummy significance.
#
# USAGE:
#   cd how-predictable-is-love   (repository root)
#   Rscript experiments/03_race_encoding.R
#
# RUNTIME: ~40 min (3 forests × ~10 min each + full-metro eval)
# ==============================================================================

library(tidyverse)
library(survival)
library(splines)
library(CLogitTree)

cat("=== Experiment 03: Race Encoding Ablation (LA 2019) ===\n\n")

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
cat(sprintf("Forest: ntree=%d  minnodesize=%d  ncores=%d\n\n",
            NTREE, MINNODESIZE, NCORES))

source("pipeline/00_helpers.R")

# --- Load data ---
cat("Loading data...\n")
t0 <- proc.time()
source("pipeline/01_data.R")
cat(sprintf("Data loaded: %.0f sec — %d train / %d test strata\n\n",
            (proc.time() - t0)[["elapsed"]],
            length(train_strata), length(test_strata)))

# ==============================================================================
# PREPARE CLOGIT DATA (spline bases + augmented columns)
# ==============================================================================

cat("Preparing spline bases...\n")

# Augment v2 data with v1-style columns (needed for add_spline_cols and M1_spline)
augment_v2 <- function(df) {
  df %>% mutate(
    abs_age_diff  = abs(age_diff),
    older_husband = as.integer(age_diff > 0),
    abs_educ_diff = abs(educ_diff),
    hypergamy_ed  = as.integer(educ_diff > 0),
    hypogamy_ed   = as.integer(educ_diff < 0)
  )
}

# Spline bases from training distribution
basis_age5  <- ns(train_data_v2$age_diff,  df = 5)
basis_educ3 <- ns(train_data_v2$educ_diff, df = 3)

# Helper: attach spline columns using pre-fitted bases
add_spline_cols <- function(df, basis_age, basis_educ) {
  a <- predict(basis_age,  df$age_diff)
  e <- predict(basis_educ, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns5_",  k)]] <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

train_v2_aug <- train_data_v2 %>% augment_v2() %>% add_spline_cols(basis_age5, basis_educ3)
test_v2_aug  <- test_data_v2  %>% augment_v2() %>% add_spline_cols(basis_age5, basis_educ3)

# ==============================================================================
# CLOGIT MODELS
# ==============================================================================

cat("\n--- Fitting clogit models ---\n")

# M1_spline: canonical 9-param spec
f_spline <- choice ~
  age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
  educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
  same_race + strata(strata)

t0 <- Sys.time()
m1_spline <- clogit(f_spline, data = train_v2_aug, method = "efron")
cat(sprintf("  M1_spline      (9 params): %.1f sec\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# M1_spline_race: add 4 alter race dummies (13 params)
f_spline_race <- choice ~
  age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
  educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
  same_race +
  alter_race_white + alter_race_black + alter_race_hispanic + alter_race_asian +
  strata(strata)

t0 <- Sys.time()
m1_spline_race <- clogit(f_spline_race, data = train_v2_aug, method = "efron")
cat(sprintf("  M1_spline_race (13 params): %.1f sec\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Coefficient table for both clogits
cat("\n--- Coefficient table (M1_spline) ---\n")
cf_spline <- summary(m1_spline)$coefficients
print(round(cf_spline[, c("coef", "exp(coef)", "se(coef)", "Pr(>|z|)")], 4))

cat("\n--- Coefficient table (M1_spline_race) ---\n")
cf_race <- summary(m1_spline_race)$coefficients
print(round(cf_race[, c("coef", "exp(coef)", "se(coef)", "Pr(>|z|)")], 4))

# Highlight race dummy significance
cat("\n--- Race dummy coefficients (M1_spline_race) ---\n")
race_rows <- grep("alter_race_|same_race", rownames(cf_race), value = TRUE)
for (rn in race_rows) {
  sig <- ifelse(cf_race[rn, "Pr(>|z|)"] < 0.05, "*", " ")
  cat(sprintf("  %-25s  coef=%+.4f  p=%6.4f %s\n",
              rn, cf_race[rn, "coef"], cf_race[rn, "Pr(>|z|)"], sig))
}

# ==============================================================================
# FOREST MODELS
# ==============================================================================

cat("\n--- Fitting Forest models ---\n")

fit_forest <- function(train_df, label, mtry_val) {
  cat(sprintf("  Fitting %s (mtry=%d, %d predictors)...\n",
              label, mtry_val, ncol(train_df) - 2))
  set.seed(SEED)
  t0 <- proc.time()
  suppressWarnings({
    rf <- CLogitForest(
      data         = train_df,
      response     = "choice",
      s            = "strata",
      ntree        = NTREE,
      minnodesize  = MINNODESIZE,
      minbucket    = ceiling(MINNODESIZE / 3),
      mtry         = mtry_val,
      tune.mtry    = FALSE,
      subsample.632 = TRUE,
      ncores       = NCORES,
      epsilon      = 1e-3,
      iter.max     = 5,
      lean         = TRUE,
      tree.timeout = TREE_TIMEOUT
    )
  })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  n_valid <- sum(sapply(rf$tree_list, function(tr) !is.null(tr) && !is.null(tr$tree.b)))
  cat(sprintf("    %.0f sec, %d/%d valid trees\n", elapsed, n_valid, NTREE))
  rf
}

# Forest_core: 5 predictors, same_race only, mtry=2
train_core <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
test_core  <- test_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()

rf_core <- fit_forest(train_core, "Forest_core", 2)

# Forest_full: 9 predictors, + race dummies, mtry=5
train_full <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>% as.data.frame()
test_full  <- test_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>% as.data.frame()

set.seed(SEED + 1)  # different seed for full forest
rf_full <- fit_forest(train_full, "Forest_full", 5)

# Forest_norace: race dummies but NO same_race, mtry=5
# Start from full columns then drop same_race
train_norace <- train_full[, setdiff(names(train_full), "same_race")]
test_norace  <- test_full[,  setdiff(names(test_full),  "same_race")]

set.seed(SEED + 2)
rf_norace <- fit_forest(train_norace, "Forest_norace", 5)

# ==============================================================================
# FULL-METRO EVALUATION
# ==============================================================================

cat("\n--- Building full-metro choice sets ---\n")

test_women <- sample_women_base %>% filter(strata_id %in% test_strata)

men_pool_metro <- men_pool %>%
  select(
    YEAR, MET2013,
    SERIAL_man    = SERIAL,
    PERNUM_man    = PERNUM,
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

cat(sprintf("Full-metro raw: %s rows (%d women × avg %d men)\n",
            format(nrow(full_metro_raw), big.mark = ","),
            n_distinct(full_metro_raw$strata_id),
            round(nrow(full_metro_raw) / n_distinct(full_metro_raw$strata_id))))

# Engineer v2 features
fm_v2_eng <- full_metro_raw %>%
  engineer_features_v2() %>%
  rename(strata = strata_id)

# Keep only strata where the husband is present (choice == 1 for exactly 1 row)
keep_strata <- fm_v2_eng %>%
  group_by(strata) %>%
  filter(sum(choice) == 1) %>%
  pull(strata) %>%
  unique()

fm_v2_eng <- fm_v2_eng %>% filter(strata %in% keep_strata)
cat(sprintf("Full-metro: %d strata retained, avg %.0f men per strata\n",
            length(keep_strata),
            nrow(fm_v2_eng) / length(keep_strata)))

# Prepare full-metro data subsets for each model type
fm_core   <- fm_v2_eng %>% select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
fm_full   <- fm_v2_eng %>% select_model_cols_v2(include_race_dummies = TRUE)  %>% as.data.frame()
fm_norace <- fm_full[, setdiff(names(fm_full), "same_race")]

# Add spline columns for clogit evaluation
fm_clogit <- fm_v2_eng %>%
  rename(strata = strata) %>%  # already renamed above
  augment_v2() %>%
  add_spline_cols(basis_age5, basis_educ3) %>%
  filter(strata %in% keep_strata)

# Score all models
cat("\n--- Scoring models on full-metro data ---\n")

fm_scores <- data.frame(
  strata = fm_clogit$strata,
  choice = fm_clogit$choice
)

fm_scores$score_spline <- score_clogit(m1_spline,      fm_clogit)
fm_scores$score_spline_race <- score_clogit(m1_spline_race, fm_clogit)
fm_scores$score_forest_core  <- as.numeric(predict(rf_core,  type = "response",
  response = "choice", s = "strata", newdata = fm_core[fm_core$strata %in% keep_strata, ]))
fm_scores$score_forest_full  <- as.numeric(predict(rf_full,  type = "response",
  response = "choice", s = "strata", newdata = fm_full[fm_full$strata %in% keep_strata, ]))
fm_scores$score_forest_norace <- as.numeric(predict(rf_norace, type = "response",
  response = "choice", s = "strata", newdata = fm_norace[fm_norace$strata %in% keep_strata, ]))

# Rank within strata and compute metro percentile
rank_and_pctile <- function(scores_df, score_col, label) {
  scores_df %>%
    group_by(strata) %>%
    mutate(
      rank_val = rank(-get(score_col), ties.method = "random"),
      n_men    = n()
    ) %>%
    ungroup() %>%
    filter(choice == 1) %>%
    mutate(
      metro_pctile = 100 * rank_val / n_men,
      model        = label
    ) %>%
    summarise(
      model        = first(model),
      metro_pctile = mean(metro_pctile),
      top10        = mean(rank_val <= 10),
      top50        = mean(rank_val <= 50),
      n_strata     = n()
    )
}

cat("\n=== Full-metro Results ===\n")

model_scores <- list(
  list(col = "score_forest_core",   label = "Forest_core",    n_pred = 5,  model_type = "forest"),
  list(col = "score_forest_full",   label = "Forest_full",    n_pred = 9,  model_type = "forest"),
  list(col = "score_forest_norace", label = "Forest_norace",  n_pred = 8,  model_type = "forest"),
  list(col = "score_spline",        label = "M1_spline",      n_pred = 9,  model_type = "clogit"),
  list(col = "score_spline_race",   label = "M1_spline_race", n_pred = 13, model_type = "clogit")
)

results <- data.frame()
for (ms in model_scores) {
  r <- rank_and_pctile(fm_scores, ms$col, ms$label)
  r$n_pred    <- ms$n_pred
  r$model_type <- ms$model_type
  results <- bind_rows(results, r)

  cat(sprintf("  %-18s  pctile=%5.2f%%  top10=%5.1f%%  top50=%5.1f%%\n",
              ms$label, r$metro_pctile, 100 * r$top10, 100 * r$top50))
}

# Key deltas
core_p  <- results$metro_pctile[results$model == "Forest_core"]
full_p  <- results$metro_pctile[results$model == "Forest_full"]
nor_p   <- results$metro_pctile[results$model == "Forest_norace"]
spl_p   <- results$metro_pctile[results$model == "M1_spline"]
splr_p  <- results$metro_pctile[results$model == "M1_spline_race"]

cat(sprintf("\n--- Race dummy deltas ---\n"))
cat(sprintf("  Forest: full vs core:        %+.2f pp\n", full_p - core_p))
cat(sprintf("  Forest: norace vs core:      %+.2f pp  (dummies alone, no same_race)\n", nor_p - core_p))
cat(sprintf("  Clogit: spline_race vs spl:  %+.2f pp\n", splr_p - spl_p))
cat(sprintf("\nConclusion: same_race is %s; race dummies add %s signal.\n",
            if (abs(full_p - nor_p) < 0.5) "not redundant vs dummies" else "mostly sufficient",
            if (abs(splr_p - spl_p) < 0.3 && abs(full_p - core_p) < 0.3) "negligible" else "some"))

# ==============================================================================
# SAVE
# ==============================================================================

if (!dir.exists("results")) dir.create("results")

results_out <- results %>%
  select(model, model_type, n_pred, metro_pctile, top10, top50, n_strata) %>%
  mutate(across(c(metro_pctile, top10, top50), ~ round(., 4)))

write.csv(results_out, "results/03_race_encoding.csv", row.names = FALSE)
cat(sprintf("\nSaved: results/03_race_encoding.csv\n"))
print(results_out, row.names = FALSE)
