# gen_ablation_h3.R
# Variable-importance ablation for H3: which dimension contributes most?
#
# Fits 7 linear conditional logit models on the same 2,000-woman training split,
# each using a different subset of predictors (age, education, race).
# Evaluates every model on the full metropolitan market (metro percentile rank)
# and on the 10-alternative test (top-1 accuracy).
#
# Addresses supervisor feedback: not enough to drop only race — must also drop
# age alone and education alone to establish which contributes most.
#
# Run from project root: Rscript experiments/gen_ablation_h3.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(ipumsr)
  library(CLogitTree)
})

# ── pipeline parameters (must match main run) ─────────────────────────────────
DDI_FILE    <- "data/usa_00008.xml"
SEED        <- 2025L
N_WOMEN     <- 2000L
N_ALTS      <- 10L
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.70
FIT_CLOGIT  <- TRUE
FIT_FOREST  <- FALSE
MODELS      <- "clogit"
NTREE <- MINNODESIZE <- MINBUCKET <- DEPTH_MAX <- NCORES <- TREE_TIMEOUT <- EPSILON <- ITER_MAX <- NULL

set.seed(SEED)
source("pipeline/01_data.R")

# ── prepare training / test clogit data frames ────────────────────────────────
fix_factors <- function(df, ref = NULL) {
  for (col in c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat")) {
    lvls <- if (!is.null(ref)) levels(ref[[col]]) else NULL
    df[[col]] <- factor(df[[col]], levels = lvls)
  }
  df
}
train_df <- fix_factors(train_data_clogit)
test_df  <- fix_factors(test_data_clogit, ref = train_df)

# ── build full-metro cross-join (identical to 03_evaluation.R) ────────────────
cat("\nBuilding full-metro cross-join...\n")

test_women <- sample_women_base %>% filter(strata_id %in% test_strata)

men_pool_metro <- men_pool %>%
  select(YEAR, MET2013,
         SERIAL_man       = SERIAL,
         PERNUM_man       = PERNUM,
         age_partner      = AGE,
         educ_partner_cat = educ_cat,
         educ_partner_tier= educ_tier,
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

cat(sprintf("Full-metro raw: %s rows\n", format(nrow(full_metro_raw), big.mark = ",")))

fm_v1 <- full_metro_raw %>%
  filter(!is.na(woman_educ_tier), !is.na(educ_partner_tier)) %>%
  engineer_features() %>%
  rename(strata = strata_id)

fm_v1 <- fix_factors(fm_v1, ref = train_df)

# keep strata with exactly one true husband
keep_strata <- fm_v1 %>%
  group_by(strata) %>%
  summarise(n_chosen = sum(choice), .groups = "drop") %>%
  filter(n_chosen == 1) %>%
  pull(strata)

fm_v1 <- fm_v1 %>% filter(strata %in% keep_strata)

n_men_vec <- ave(rep(1L, nrow(fm_v1)), fm_v1$strata, FUN = length)
fm_choice  <- fm_v1$choice

cat(sprintf("Full-metro: %d strata retained, avg %.0f men/strata\n",
            length(keep_strata), mean(n_men_vec[fm_choice == 1])))

# ── helper: fit one clogit and evaluate ───────────────────────────────────────
eval_ablation <- function(label, formula_rhs) {
  f <- as.formula(paste0("choice ~ ", formula_rhs, " + strata(strata)"))

  m <- tryCatch(
    clogit(f, data = train_df, method = "efron"),
    error = function(e) { cat("  ERROR in", label, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)

  # 10-alt accuracy
  sc10 <- data.frame(
    strata = test_df$strata,
    choice = test_df$choice,
    score  = score_clogit(m, test_df)
  )
  rank10 <- sc10 %>%
    group_by(strata) %>%
    mutate(rank = rank(-score, ties.method = "random")) %>%
    filter(choice == 1)
  top1 <- mean(rank10$rank == 1)

  # full-metro evaluation
  fm_sc <- score_clogit(m, fm_v1)
  fm_df <- data.frame(strata = fm_v1$strata, choice = fm_choice,
                      score = fm_sc, n_men = n_men_vec)
  chosen_fm <- fm_df %>%
    group_by(strata) %>%
    mutate(rank = rank(-score, ties.method = "random")) %>%
    filter(choice == 1) %>%
    mutate(pctile = rank / n_men * 100)
  metro_pctile <- mean(chosen_fm$pctile)

  cat(sprintf("  %-45s  top-1=%4.1f%%  metro=%.2f%%\n",
              label, top1 * 100, metro_pctile))

  tibble(Model = label, `Top-1 (%)` = round(top1 * 100, 1),
         `Metro rank (%)` = round(metro_pctile, 2))
}

# ── run all ablation models ───────────────────────────────────────────────────
cat("\nFitting ablation models...\n")

specs <- list(
  c("All three: age + education + race",
    "abs_age_diff + older_husband + abs_educ_diff + hypergamy_ed + hypogamy_ed + same_race"),
  c("Drop race  (age + education only)",
    "abs_age_diff + older_husband + abs_educ_diff + hypergamy_ed + hypogamy_ed"),
  c("Drop age   (education + race only)",
    "abs_educ_diff + hypergamy_ed + hypogamy_ed + same_race"),
  c("Drop educ  (age + race only)",
    "abs_age_diff + older_husband + same_race"),
  c("Race alone",    "same_race"),
  c("Age alone",     "abs_age_diff + older_husband"),
  c("Education alone","abs_educ_diff + hypergamy_ed + hypogamy_ed")
)

results <- bind_rows(lapply(specs, function(s) eval_ablation(s[1], s[2])))

# add random baseline and compute losses
full_metro_val <- results$`Metro rank (%)`[1]
full_top1_val  <- results$`Top-1 (%)`[1]

results <- results %>%
  mutate(
    `Loss in metro rank (pp)` = round(`Metro rank (%)` - full_metro_val, 2),
    `Loss in top-1 (pp)`      = round(`Top-1 (%)` - full_top1_val, 1)
  )
results$`Loss in metro rank (pp)`[1] <- NA
results$`Loss in top-1 (pp)`[1]      <- NA

random_row <- tibble(
  Model = "Random baseline",
  `Top-1 (%)` = 10.0, `Metro rank (%)` = 50.0,
  `Loss in metro rank (pp)` = NA_real_,
  `Loss in top-1 (pp)` = NA_real_
)
results <- bind_rows(random_row, results)

# ── print and save ────────────────────────────────────────────────────────────
cat("\n=== ABLATION RESULTS (H3: variable importance) ===\n")
print(as.data.frame(results), row.names = FALSE)

out_path <- "results/ablation_h3.csv"
write_csv(results, out_path)
cat(sprintf("\nSaved: %s\n", out_path))
