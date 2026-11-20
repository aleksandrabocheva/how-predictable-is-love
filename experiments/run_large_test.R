# run_large_test.R
# Addresses reviewer concern about small test set.
#
# Design:
#   - Train on 2,000 women (same seed/sample as main analysis)
#   - Test  on 10,000 FRESH women drawn from the remaining eligible pool
#   - All four models (clogit + forest); forest trains on same 2,000 so
#     fitting time is identical to the original run (~12 min)
#
# Run from project root:
#   Rscript experiments/run_large_test.R
#
# Output: results/large_test/metrics_large_test.csv

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(ipumsr)
library(survival)
library(splines)
library(CLogitTree)
source("pipeline/00_helpers.R")

# --- Parameters ---------------------------------------------------------------
SEED            <- 2025
N_TRAIN         <- 2000
N_TEST          <- 10000
N_ALTS          <- 10
MIN_MEN_PER_METRO <- 100
DDI_FILE        <- "data/usa_00008.xml"
OUT_DIR         <- "results/large_test"
NTREE           <- 50
MINNODESIZE     <- 550
MINBUCKET       <- 184
DEPTH_MAX       <- 10
NCORES          <- 2L
TREE_TIMEOUT    <- 600
EPSILON         <- 1e-3
ITER_MAX        <- 5
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Large held-out test ===\n")
cat(sprintf("Train: %d women  |  Test: %d fresh women  |  Alts: %d\n",
            N_TRAIN, N_TEST, N_ALTS))

# --- Load and filter data (same as 01_data.R) ---------------------------------
cat("\nLoading data...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(ddi,
  vars = c("YEAR","SERIAL","PERNUM","SPLOC","SEX","AGE","EDUC",
           "RACE","HISPAN","MET2013","YRMARR","MARST")
) %>% zap_labels()

recode_race <- function(r, h) {
  case_when(h > 0 ~ "Hispanic", r == 1 ~ "White", r == 2 ~ "Black",
            r %in% 4:6 ~ "Asian", TRUE ~ "Other")
}

data_clean <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0) %>%
  mutate(marriage_duration = YEAR - YRMARR) %>%
  filter(marriage_duration <= 3) %>%
  mutate(
    race_refined = recode_race(RACE, HISPAN),
    educ_cat = case_when(
      EDUC < 6      ~ "1_LessHS",
      EDUC == 6     ~ "2_HS_Grad",
      EDUC %in% 7:9 ~ "3_SomeColl",
      EDUC == 10    ~ "4_Bachelor",
      EDUC >= 11    ~ "5_Advanced",
      TRUE          ~ "Unknown"
    ),
    educ_tier = case_when(
      educ_cat %in% c("1_LessHS","2_HS_Grad")      ~ "Low",
      educ_cat == "3_SomeColl"                      ~ "Mid",
      educ_cat %in% c("4_Bachelor","5_Advanced")    ~ "High",
      TRUE                                          ~ NA_character_
    ),
    educ_yrs = case_when(
      EDUC == 0 ~ 0,  EDUC == 1 ~ 4,  EDUC == 2 ~ 8,
      EDUC == 3 ~ 9,  EDUC == 4 ~ 10, EDUC == 5 ~ 11,
      EDUC == 6 ~ 12, EDUC == 7 ~ 13, EDUC == 8 ~ 14,
      EDUC == 9 ~ 15, EDUC == 10 ~ 16, EDUC >= 11 ~ 18,
      TRUE ~ NA_real_
    )
  )

men_pool   <- data_clean %>% filter(SEX == 1)
women_pool <- data_clean %>% filter(SEX == 2)

metro_pool_sizes <- men_pool %>%
  count(MET2013, YEAR, name = "n_men_in_market")

# Join each woman to her husband
couples_women <- women_pool %>%
  inner_join(
    men_pool %>% select(YEAR, SERIAL, PERNUM,
                        age_partner      = AGE,
                        educ_partner_cat = educ_cat,
                        educ_partner_tier = educ_tier,
                        educ_yrs_partner = educ_yrs,
                        race_partner     = race_refined),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  ) %>%
  left_join(metro_pool_sizes, by = c("MET2013", "YEAR")) %>%
  filter(!is.na(n_men_in_market), n_men_in_market >= MIN_MEN_PER_METRO) %>%
  select(-n_men_in_market)

cat(sprintf("Eligible women after all filters: %s\n",
            format(nrow(couples_women), big.mark = ",")))

# --- Sample train and test WITHOUT overlap ------------------------------------
set.seed(SEED)
train_idx <- sample(nrow(couples_women), size = N_TRAIN)
train_women <- couples_women[train_idx, ] %>% mutate(strata_id = row_number())

remaining   <- couples_women[-train_idx, ]
n_test_actual <- min(N_TEST, nrow(remaining))
test_idx    <- sample(nrow(remaining), size = n_test_actual)
test_women  <- remaining[test_idx, ] %>% mutate(strata_id = row_number())

cat(sprintf("Training women: %d  |  Test women: %d  |  Unused: %s\n",
            nrow(train_women), nrow(test_women),
            format(nrow(remaining) - n_test_actual, big.mark = ",")))

# --- Build choice sets --------------------------------------------------------
generate_choice_sets <- function(input_women, input_men, n_alts) {
  n_random <- n_alts - 1
  build_set <- function(w) {
    pool <- input_men %>%
      filter(YEAR == w$YEAR, MET2013 == w$MET2013) %>%
      filter(!(SERIAL == w$SERIAL & PERNUM == w$SPLOC))
    if (nrow(pool) < n_random) return(NULL)
    alts <- pool %>%
      slice_sample(n = n_random) %>%
      mutate(choice = 0) %>%
      rename(age_partner      = AGE,
             educ_partner_cat = educ_cat,
             educ_partner_tier = educ_tier,
             educ_yrs_partner = educ_yrs,
             race_partner     = race_refined)
    chosen <- tibble(
      age_partner      = w$age_partner,
      educ_partner_cat = w$educ_partner_cat,
      educ_partner_tier = w$educ_partner_tier,
      educ_yrs_partner = w$educ_yrs_partner,
      race_partner     = w$race_partner,
      choice = 1
    )
    bind_rows(chosen, alts) %>%
      mutate(strata_id   = w$strata_id,
             woman_age   = w$AGE,
             woman_educ_cat  = w$educ_cat,
             woman_educ_tier = w$educ_tier,
             woman_educ_yrs  = w$educ_yrs,
             woman_race  = w$race_refined)
  }
  purrr::map_dfr(split(input_women, seq_len(nrow(input_women))), build_set)
}

men_pool_small <- men_pool %>%
  semi_join(
    bind_rows(
      train_women %>% distinct(YEAR, MET2013),
      test_women  %>% distinct(YEAR, MET2013)
    ),
    by = c("YEAR", "MET2013")
  )

cat("Building training choice sets...\n")
raw_train <- generate_choice_sets(train_women, men_pool_small, N_ALTS)

cat("Building test choice sets...\n")
raw_test  <- generate_choice_sets(test_women,  men_pool_small, N_ALTS)

# --- Feature engineering ------------------------------------------------------
clean_v1 <- function(raw) {
  raw %>%
    engineer_features() %>%
    rename(strata = strata_id) %>%
    filter(if_all(c("choice","strata","abs_age_diff","older_husband",
                    "abs_educ_diff","hypergamy_ed","hypogamy_ed","same_race"),
                  ~ !is.na(.))) %>%
    group_by(strata) %>%
    filter(sum(choice) == 1, n() >= 2) %>%
    ungroup()
}

clean_v2 <- function(raw) {
  raw %>%
    engineer_features_v2() %>%
    rename(strata = strata_id) %>%
    filter(if_all(c("choice","strata","age_diff","alter_age",
                    "educ_diff","alter_educ_rank","same_race"), ~ !is.na(.))) %>%
    group_by(strata) %>%
    filter(sum(choice) == 1, n() >= 2) %>%
    ungroup()
}

cat("Engineering features...\n")
train_v1 <- clean_v1(raw_train)
test_v1  <- clean_v1(raw_test)
train_v2 <- clean_v2(raw_train)
test_v2  <- clean_v2(raw_test)

# Keep only strata that survived both encodings
keep_train <- intersect(unique(train_v1$strata), unique(train_v2$strata))
keep_test  <- intersect(unique(test_v1$strata),  unique(test_v2$strata))
train_v1 <- train_v1 %>% filter(strata %in% keep_train)
train_v2 <- train_v2 %>% filter(strata %in% keep_train)
test_v1  <- test_v1  %>% filter(strata %in% keep_test)
test_v2  <- test_v2  %>% filter(strata %in% keep_test)

cat(sprintf("Training strata: %d  |  Test strata: %d\n",
            n_distinct(train_v1$strata), n_distinct(test_v1$strata)))

# --- Fit models ---------------------------------------------------------------
cat("\nFitting M1_v1 (linear clogit)...\n")
m1_v1 <- clogit(
  choice ~ abs_age_diff + older_husband + abs_educ_diff +
    hypergamy_ed + hypogamy_ed + same_race + strata(strata),
  data   = train_v1 %>% as.data.frame(),
  method = "efron"
)

cat("Fitting M1_spline (spline clogit)...\n")
train_v2_aug <- train_v2 %>%
  mutate(abs_age_diff  = abs(age_diff),
         older_husband = as.integer(age_diff > 0),
         abs_educ_diff = abs(educ_diff),
         hypergamy_ed  = as.integer(educ_diff > 0),
         hypogamy_ed   = as.integer(educ_diff < 0))

basis_age5  <- ns(train_v2_aug$age_diff,  df = 5)
basis_educ3 <- ns(train_v2_aug$educ_diff, df = 3)

add_splines <- function(df) {
  a <- predict(basis_age5,  df$age_diff)
  e <- predict(basis_educ3, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns5_",  k)]] <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}
train_v2_aug <- add_splines(train_v2_aug)

m1_spline <- clogit(
  choice ~ age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
    educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
    same_race + strata(strata),
  data   = train_v2_aug %>% as.data.frame(),
  method = "efron"
)

# --- Fit forests --------------------------------------------------------------
cat("\nFitting Forest_core (5 predictors, v2 encoding)...\n")
train_core <- train_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>%
  as.data.frame()

set.seed(SEED)
t0 <- proc.time()
rf_core <- CLogitForest(
  data = train_core, response = "choice", s = "strata",
  ntree = NTREE, minnodesize = MINNODESIZE, minbucket = MINBUCKET,
  depth_max = DEPTH_MAX, mtry = 2, tune.mtry = FALSE,
  subsample.632 = TRUE, ncores = NCORES,
  epsilon = EPSILON, iter.max = ITER_MAX,
  lean = TRUE, tree.timeout = TREE_TIMEOUT
)
cat(sprintf("  Done: %.0fs\n", (proc.time() - t0)[["elapsed"]]))

cat("Fitting Forest_full (9 predictors, v2 encoding)...\n")
train_full <- train_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>%
  as.data.frame()

set.seed(SEED + 1)
t0 <- proc.time()
rf_full <- CLogitForest(
  data = train_full, response = "choice", s = "strata",
  ntree = NTREE, minnodesize = MINNODESIZE, minbucket = MINBUCKET,
  depth_max = DEPTH_MAX, mtry = 5, tune.mtry = FALSE,
  subsample.632 = TRUE, ncores = NCORES,
  epsilon = EPSILON, iter.max = ITER_MAX,
  lean = TRUE, tree.timeout = TREE_TIMEOUT
)
cat(sprintf("  Done: %.0fs\n", (proc.time() - t0)[["elapsed"]]))

# --- Evaluate on large test set -----------------------------------------------
top_k_accuracy <- function(scores, choice, strata, k) {
  tibble(scores = scores, choice = choice, strata = strata) %>%
    group_by(strata) %>%
    summarise(hit = rank(-scores, ties.method = "first")[choice == 1] <= k,
              .groups = "drop") %>%
    pull(hit) %>% mean()
}

evaluate <- function(model, test_df, label) {
  scores <- score_clogit(model, test_df %>% as.data.frame())
  tibble(
    model   = label,
    n_strata = n_distinct(test_df$strata),
    top1    = top_k_accuracy(scores, test_df$choice, test_df$strata, 1),
    top3    = top_k_accuracy(scores, test_df$choice, test_df$strata, 3),
    top5    = top_k_accuracy(scores, test_df$choice, test_df$strata, 5),
    random_top1 = 1 / N_ALTS
  )
}

cat("\nEvaluating on large test set...\n")
test_v2_aug <- test_v2 %>%
  mutate(abs_age_diff  = abs(age_diff),
         older_husband = as.integer(age_diff > 0),
         abs_educ_diff = abs(educ_diff),
         hypergamy_ed  = as.integer(educ_diff > 0),
         hypogamy_ed   = as.integer(educ_diff < 0)) %>%
  add_splines()

# Forest scoring uses predict() directly
score_forest <- function(rf, test_df) {
  as.numeric(predict(rf, newdata = test_df %>% as.data.frame()))
}

evaluate_forest <- function(rf, test_df, label) {
  scores <- score_forest(rf, test_df)
  tibble(
    model    = label,
    n_strata = n_distinct(test_df$strata),
    top1     = top_k_accuracy(scores, test_df$choice, test_df$strata, 1),
    top3     = top_k_accuracy(scores, test_df$choice, test_df$strata, 3),
    top5     = top_k_accuracy(scores, test_df$choice, test_df$strata, 5),
    random_top1 = 1 / N_ALTS
  )
}

test_core <- test_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>%
  as.data.frame()
test_full <- test_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>%
  as.data.frame()

results <- bind_rows(
  evaluate(m1_v1,     test_v1,     "M1_v1"),
  evaluate(m1_spline, test_v2_aug, "M1_spline"),
  evaluate_forest(rf_core, test_core, "Forest_core"),
  evaluate_forest(rf_full, test_full, "Forest_full")
)

# --- Report -------------------------------------------------------------------
cat("\n=== RESULTS (large held-out test) ===\n")
cat(sprintf("  Training: %d women  |  Test: %d women\n",
            n_distinct(train_v1$strata), n_distinct(test_v1$strata)))
cat(sprintf("  Random baseline (top-1 among %d): %.1f%%\n\n",
            N_ALTS, 100 / N_ALTS))

results %>%
  mutate(across(c(top1, top3, top5), ~ paste0(round(. * 100, 1), "%"))) %>%
  select(model, n_strata, top1, top3, top5) %>%
  print()

out_path <- file.path(OUT_DIR, "metrics_large_test.csv")
write_csv(results, out_path)
cat(sprintf("\nSaved: %s\n", out_path))
