# ==============================================================================
# 03_evaluation.R — Rank-based prediction accuracy
#
# Evaluates all 4 models on two tasks:
#
#   1. N-ALT EVALUATION (quick sanity check)
#      Within each test choice set (the same N_ALTS alternatives used in training),
#      rank the true partner. Report top-1 accuracy (% of times the model ranks
#      the true partner #1 out of N_ALTS alternatives).
#
#   2. FULL-METRO EVALUATION (the real test)
#      For each test woman, rank her true partner among ALL men in her metro area
#      (typically 2,000-3,000 men). This is much harder and more realistic than
#      the 10-alt task. The key metric is the percentile rank: how far up the
#      ranking is the true partner? Lower = better, 50 = random.
#
# WHY TWO EVALUATIONS?
#   The N-alt evaluation uses the same random alternatives the model trained on.
#   It measures how well the model distinguishes the true partner from a small
#   random sample. But in reality, a woman's "choice" is from ALL available men,
#   not 10 random ones. The full-metro evaluation tests whether the model's
#   learned preferences generalize to ranking among the entire market.
#
# EXPECTS from 01_data.R:
#   test_data, test_data_clogit, test_data_v2 — test set in both encodings
#   sample_women_base — focal women with demographics
#   men_pool — all men in the data (for full-metro cross-join)
#   test_strata, train_strata — which choice sets are train vs test
# EXPECTS from 02_model.R:
#   m1_v1, m1_spline5 — fitted clogit model objects
#   rf_core, rf_full — fitted CLogitForest objects
#   basis_age5, basis_educ3 — spline basis objects (for consistent knot placement)
#   add_spline_cols — helper function to apply spline bases
# ==============================================================================

library(splines)

eval_label <- if (DO_SPLIT) "held-out test" else "in-sample"
cat(sprintf("Evaluating on %s data: %d strata\n", eval_label,
            length(test_strata)))

# Which models are actually fitted? Each of these is NULL if its family
# was skipped in 02_model.R (via MODELS="clogit" or MODELS="forest").
has_m1_v1   <- !is.null(m1_v1)
has_spline  <- !is.null(m1_spline5)
has_rf_core <- !is.null(rf_core)
has_rf_full <- !is.null(rf_full)

if (!any(c(has_m1_v1, has_spline, has_rf_core, has_rf_full))) {
  stop("No fitted models found — nothing to evaluate.")
}


# ==============================================================================
# EVALUATION 1: N-alt (rank within training choice sets)
# ==============================================================================

cat(sprintf("\n--- %d-alt Evaluation ---\n", N_ALTS))

# Prepare test data for scoring. The four models use different predictor encodings
# (see 01_data.R STEP 6 for what "v1" and "v2" mean):
#   - M1_v1 needs v1 columns (abs_age_diff, older_husband, etc.)
#   - M1_spline5 needs v2 columns (age_diff, educ_diff) + spline basis columns
#   - Forest_core needs v2 core columns (5 predictors)
#   - Forest_full needs v2 full columns (9 predictors, including race dummies)
# We build "test_v2_aug" (only needed if any clogit family was fit: it adds both
# v1 columns and spline basis columns).
if (has_spline) {
  test_v2_aug <- test_data_v2 %>%
    mutate(
      abs_age_diff  = abs(age_diff),
      older_husband = as.integer(age_diff > 0),
      abs_educ_diff = abs(educ_diff),
      hypergamy_ed  = as.integer(educ_diff > 0),
      hypogamy_ed   = as.integer(educ_diff < 0)
    ) %>%
    add_spline_cols(basis_age5, basis_educ3)
}

# Forest test data (only if a forest was fit)
if (has_rf_core) {
  test_core <- test_data_v2 %>%
    select_model_cols_v2(include_race_dummies = FALSE) %>%
    as.data.frame()
}
if (has_rf_full) {
  test_full <- test_data_v2 %>%
    select_model_cols_v2(include_race_dummies = TRUE) %>%
    as.data.frame()
}

# SCORE each fitted model. score_clogit() returns exp(X*beta). predict() on
# CLogitForest averages predictions across all trees.
scores_10alt <- data.frame(
  strata = test_data_v2$strata,
  choice = test_data_v2$choice
)
if (has_m1_v1)   scores_10alt$m1_v1   <- score_clogit(m1_v1, test_data_clogit)
if (has_spline)  scores_10alt$m1_spl5 <- score_clogit(m1_spline5, test_v2_aug)
if (has_rf_core) scores_10alt$rf_core <- as.numeric(predict(rf_core,
  type = "response", response = "choice", s = "strata", newdata = test_core))
if (has_rf_full) scores_10alt$rf_full <- as.numeric(predict(rf_full,
  type = "response", response = "choice", s = "strata", newdata = test_full))

# RANK each fitted model within each choice set (highest score = rank 1)
ranked_10alt <- scores_10alt %>% group_by(strata)
if (has_m1_v1)   ranked_10alt <- ranked_10alt %>% mutate(rank_m1_v1   = rank(-m1_v1,   ties.method = "random"))
if (has_spline)  ranked_10alt <- ranked_10alt %>% mutate(rank_m1_spl5 = rank(-m1_spl5, ties.method = "random"))
if (has_rf_core) ranked_10alt <- ranked_10alt %>% mutate(rank_rf_core = rank(-rf_core, ties.method = "random"))
if (has_rf_full) ranked_10alt <- ranked_10alt %>% mutate(rank_rf_full = rank(-rf_full, ties.method = "random"))
ranked_10alt <- ranked_10alt %>% ungroup()

# Keep only the true partner's rank (the row where choice == 1)
chosen_10alt <- ranked_10alt %>% filter(choice == 1)

# Report accuracy: what fraction of true partners were ranked #1, top-3, top-5?
eval_10alt <- function(rank_col, label) {
  ranks <- chosen_10alt[[rank_col]]
  top1 <- mean(ranks == 1)    # % ranked first
  top3 <- mean(ranks <= 3)    # % in top 3
  top5 <- mean(ranks <= 5)    # % in top 5
  mr   <- mean(ranks)         # average rank (1 = perfect, N_ALTS/2 = random)
  cat(sprintf("  %-15s  top1=%5.1f%%  top3=%5.1f%%  top5=%5.1f%%  mean_rank=%.2f\n",
              label, 100*top1, 100*top3, 100*top5, mr))
  data.frame(model = label, top1 = round(top1, 4), top3 = round(top3, 4),
             top5 = round(top5, 4), mean_rank = round(mr, 3),
             stringsAsFactors = FALSE)
}

res_10_rows <- list()
if (has_m1_v1)   res_10_rows[["M1_v1"]]       <- eval_10alt("rank_m1_v1",   "M1_v1")
if (has_spline)  res_10_rows[["M1_spline5"]]  <- eval_10alt("rank_m1_spl5", "M1_spline5")
if (has_rf_core) res_10_rows[["Forest_core"]] <- eval_10alt("rank_rf_core", "Forest_core")
if (has_rf_full) res_10_rows[["Forest_full"]] <- eval_10alt("rank_rf_full", "Forest_full")
res_10 <- do.call(rbind, res_10_rows)


# ==============================================================================
# EVALUATION 2: Full-metro (rank among ALL men in the metro)
# ==============================================================================

cat("\n--- Full-metro Evaluation ---\n")

# Get the test women (we need their demographics + metro/year for the cross-join)
test_women <- sample_women_base %>%
  filter(strata_id %in% test_strata)

# Prepare the men pool with renamed columns matching the choice data format
men_pool_metro <- men_pool %>%
  select(YEAR, MET2013,
         SERIAL_man  = SERIAL,
         PERNUM_man  = PERNUM,
         age_partner       = AGE,
         educ_partner_cat  = educ_cat,
         educ_partner_tier = educ_tier,
         educ_yrs_partner  = educ_yrs,
         race_partner      = race_refined)

# Tag which men are actual husbands of test women (for setting choice = 1)
husbands <- test_women %>%
  select(YEAR, MET2013, SERIAL_man = SERIAL, PERNUM_man = SPLOC) %>%
  mutate(.is_husband = TRUE)

men_pool_tagged <- men_pool_metro %>%
  left_join(husbands, by = c("YEAR", "MET2013", "SERIAL_man", "PERNUM_man")) %>%
  mutate(.is_husband = coalesce(.is_husband, FALSE)) %>%
  select(-.is_husband)

cat(sprintf("Full-metro: %d test women x %d men in pool\n",
            nrow(test_women), nrow(men_pool_tagged)))

# BUILD THE FULL-METRO DATA: cross-join each test woman with ALL men in her
# metro-year. This is a many-to-many join that produces a very large data frame
# (e.g., 735 women x 2,608 men = 1.9 million rows for LA 2019).
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

cat(sprintf("Full-metro raw: %s rows\n",
            format(nrow(full_metro_raw), big.mark = ",")))

# ENGINEER FEATURES on the full-metro data (same functions as training data).
# Only prepare the encodings needed for the fitted models.

# v1 features — needed only for M1_v1 scoring
if (has_m1_v1) {
  fm_v1 <- full_metro_raw %>%
    filter(!is.na(woman_educ_tier), !is.na(educ_partner_tier)) %>%
    engineer_features() %>%
    rename(strata = strata_id)

  fm_v1_numeric <- fm_v1 %>% select_model_cols()
  # Add any one-hot columns that exist in training but not in this metro's data
  train_v1_pred <- setdiff(names(train_data), c("choice", "strata"))
  for (col in train_v1_pred) {
    if (!col %in% names(fm_v1_numeric)) fm_v1_numeric[[col]] <- 0L
  }
  fm_v1_numeric <- fm_v1_numeric[, c("choice", "strata", train_v1_pred)]

  # Add factor columns for clogit scoring
  fm_v1_clogit <- bind_cols(
    fm_v1_numeric,
    fm_v1[, c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat"), drop = FALSE]
  )
  fm_v1_clogit$woman_race       <- factor(fm_v1_clogit$woman_race)
  fm_v1_clogit$race_partner     <- factor(fm_v1_clogit$race_partner)
  fm_v1_clogit$woman_educ_cat   <- factor(fm_v1_clogit$woman_educ_cat)
  fm_v1_clogit$educ_partner_cat <- factor(fm_v1_clogit$educ_partner_cat)
}

# v2 features — needed for M1_spline, Forest_core, and Forest_full
fm_v2 <- full_metro_raw %>%
  engineer_features_v2() %>%
  rename(strata = strata_id)

fm_v2_core <- fm_v2 %>% select_model_cols_v2(include_race_dummies = FALSE) %>% as.data.frame()
if (has_rf_full) {
  fm_v2_full <- fm_v2 %>% select_model_cols_v2(include_race_dummies = TRUE) %>% as.data.frame()
}
if (has_spline) {
  fm_v2_aug <- fm_v2_core %>%
    mutate(
      abs_age_diff  = abs(age_diff),
      older_husband = as.integer(age_diff > 0),
      abs_educ_diff = abs(educ_diff),
      hypergamy_ed  = as.integer(educ_diff > 0),
      hypogamy_ed   = as.integer(educ_diff < 0)
    ) %>%
    add_spline_cols(basis_age5, basis_educ3)
}

# Some strata lose their true partner during feature engineering (e.g., husband's
# education was "Unknown"). Drop these — a choice set without a correct answer
# can't be scored.
keep_v2 <- fm_v2_core %>% group_by(strata) %>%
  filter(sum(choice) == 1) %>% pull(strata) %>% unique()

if (has_m1_v1) {
  keep_v1 <- fm_v1_numeric %>% group_by(strata) %>%
    filter(sum(choice) == 1) %>% pull(strata) %>% unique()
  keep_strata <- intersect(keep_v1, keep_v2)
} else {
  keep_strata <- keep_v2
}

if (has_m1_v1) {
  fm_v1_numeric <- fm_v1_numeric[fm_v1_numeric$strata %in% keep_strata, ] %>% as.data.frame()
  fm_v1_clogit  <- fm_v1_clogit[fm_v1_clogit$strata %in% keep_strata, ] %>% as.data.frame()
}
fm_v2_core    <- fm_v2_core[fm_v2_core$strata %in% keep_strata, ] %>% as.data.frame()
if (has_rf_full) fm_v2_full <- fm_v2_full[fm_v2_full$strata %in% keep_strata, ] %>% as.data.frame()
if (has_spline)  fm_v2_aug  <- fm_v2_aug[fm_v2_aug$strata %in% keep_strata, ] %>% as.data.frame()

n_fm_strata <- length(keep_strata)
cat(sprintf("Full-metro: %d strata retained\n", n_fm_strata))

# SCORE each fitted model on the full-metro data
fm_scores <- data.frame(
  strata = fm_v2_core$strata,
  choice = fm_v2_core$choice
)
if (has_m1_v1)   fm_scores$m1_v1   <- score_clogit(m1_v1, fm_v1_clogit)
if (has_spline)  fm_scores$m1_spl5 <- score_clogit(m1_spline5, fm_v2_aug)
if (has_rf_core) fm_scores$rf_core <- as.numeric(predict(rf_core,
  type = "response", response = "choice", s = "strata", newdata = fm_v2_core))
if (has_rf_full) fm_scores$rf_full <- as.numeric(predict(rf_full,
  type = "response", response = "choice", s = "strata", newdata = fm_v2_full))

# RANK within each strata (each test woman's full metro market)
score_cols_fm <- c()
if (has_m1_v1)   score_cols_fm <- c(score_cols_fm, "m1_v1")
if (has_spline)  score_cols_fm <- c(score_cols_fm, "m1_spl5")
if (has_rf_core) score_cols_fm <- c(score_cols_fm, "rf_core")
if (has_rf_full) score_cols_fm <- c(score_cols_fm, "rf_full")
for (mc in score_cols_fm) {
  fm_scores[[paste0("rank_", mc)]] <- ave(
    fm_scores[[mc]], fm_scores$strata,
    FUN = function(x) rank(-x, ties.method = "random")
  )
}

# How many men in each strata's market? (varies by metro size)
n_men_per_strata <- ave(rep(1, nrow(fm_scores)), fm_scores$strata, FUN = length)
# Keep only the true partner's rank
chosen_fm <- fm_scores[fm_scores$choice == 1, ]
n_men_chosen <- n_men_per_strata[fm_scores$choice == 1]

# Report full-metro accuracy
res_fm <- data.frame()

eval_fm <- function(rank_col, label) {
  ranks <- chosen_fm[[rank_col]]
  # Percentile: rank / total men * 100 (lower = better, 50 = random)
  pctile <- 100 * ranks / n_men_chosen
  mr  <- mean(ranks)          # mean absolute rank
  mp  <- mean(pctile)         # mean percentile (the primary metric)
  t10 <- mean(ranks <= 10)    # % of true partners in top 10 men
  t50 <- mean(ranks <= 50)    # % of true partners in top 50 men
  cat(sprintf("  %-15s  mean_rank=%5.0f  pctile=%5.2f%%  top10=%5.1f%%  top50=%5.1f%%\n",
              label, mr, mp, 100*t10, 100*t50))
  data.frame(model = label, metro_mean_rank = round(mr, 1),
             metro_pctile = round(mp, 2), metro_top10 = round(t10, 4),
             metro_top50 = round(t50, 4), metro_n_men = round(mean(n_men_chosen)),
             stringsAsFactors = FALSE)
}

cat(sprintf("\nFull-metro results (%d strata, avg %d men):\n",
            n_fm_strata, round(mean(n_men_chosen))))

res_fm_rows <- list()
if (has_m1_v1)   res_fm_rows[["M1_v1"]]       <- eval_fm("rank_m1_v1",   "M1_v1")
if (has_spline)  res_fm_rows[["M1_spline5"]]  <- eval_fm("rank_m1_spl5", "M1_spline5")
if (has_rf_core) res_fm_rows[["Forest_core"]] <- eval_fm("rank_rf_core", "Forest_core")
if (has_rf_full) res_fm_rows[["Forest_full"]] <- eval_fm("rank_rf_full", "Forest_full")
res_fm <- do.call(rbind, res_fm_rows)

# Show random baseline for reference
cat(sprintf("  %-15s  mean_rank=%5.0f  pctile=%5.2f%%\n",
            "Random", round(mean(n_men_chosen) / 2), 50.0))


# ==============================================================================
# COMBINED RESULTS TABLE
# ==============================================================================

# Merge 10-alt and full-metro results into one table
results <- merge(res_10, res_fm, by = "model", all = TRUE)

cat("\n=== SUMMARY ===\n")
cat(sprintf("%-15s  %6s  %6s  %8s  %8s\n",
            "Model", "top1", "top3", "pctile", "top10"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  cat(sprintf("%-15s  %5.1f%%  %5.1f%%  %7.2f%%  %7.1f%%\n",
              r$model, 100*r$top1, 100*r$top3,
              r$metro_pctile, 100*r$metro_top10))
}

# Save results as CSV
if (exists("results_path")) {
  write.csv(results, results_path, row.names = FALSE)
  cat(sprintf("\nResults saved to: %s\n", results_path))
}


# ==============================================================================
# PLOTS
# ==============================================================================
#
# Each plot renders only lines for the models that were actually fitted.
# Palette + linetypes are shared so visuals are stable across partial-MODELS runs.

palette <- c(
  "Forest_core" = "#2563eb",
  "Forest_full" = "#9333ea",
  "M1_spline5"  = "#059669",
  "M1_v1"       = "#f97316",
  "Random"      = "grey50"
)
linetypes <- c(
  "Forest_core" = "solid",
  "Forest_full" = "solid",
  "M1_spline5"  = "solid",
  "M1_v1"       = "dotted",
  "Random"      = "dashed"
)

# Mapping from fitted-model slot -> (label, rank column in chosen_* data)
plot_models <- list(
  list(label = "Forest_core", rank = "rank_rf_core", present = has_rf_core),
  list(label = "Forest_full", rank = "rank_rf_full", present = has_rf_full),
  list(label = "M1_spline5",  rank = "rank_m1_spl5", present = has_spline),
  list(label = "M1_v1",       rank = "rank_m1_v1",   present = has_m1_v1)
)
plot_models <- Filter(function(x) x$present, plot_models)
model_labels <- vapply(plot_models, function(x) x$label, character(1))

# Cumulative accuracy at rank cutoff k = 1..K, in long form (one row per
# (model, k)). Then a single geom_line handles all lines.
build_cum_long <- function(chosen_df, K) {
  # Random baseline: k / n_men averaged across women; for 10-alt it's k/K
  if (K == N_ALTS) {
    rnd <- (1:K) / K
  } else {
    rnd <- sapply(1:K, function(k) mean(k / n_men_chosen))
  }
  rows <- list(data.frame(k = 1:K, model = "Random", pct = rnd, stringsAsFactors = FALSE))
  for (m in plot_models) {
    pct <- sapply(1:K, function(k) mean(chosen_df[[m$rank]] <= k))
    rows[[m$label]] <- data.frame(k = 1:K, model = m$label, pct = pct, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$model <- factor(out$model, levels = c(model_labels, "Random"))
  out
}

# --- Plot 1: Cumulative accuracy (10-alt) ---
MAX_K <- N_ALTS
cum_long <- build_cum_long(chosen_10alt, MAX_K)

p <- ggplot(cum_long, aes(x = k, y = pct, colour = model, linetype = model)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 1:MAX_K) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_manual(values = palette) +
  scale_linetype_manual(values = linetypes) +
  labs(
    title = paste0("Cumulative rank accuracy (", eval_label, ")"),
    subtitle = sprintf("Strata: %d | Choice set size: %d", nrow(chosen_10alt), N_ALTS),
    x = "Top k", y = "% of true partners ranked within top k",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p)
if (exists("config_id")) {
  ggsave(file.path(results_dir, "plot_10alt.png"),
         plot = p, width = 10, height = 6, dpi = 150)
}

# --- Plot 2: Cumulative accuracy (full-metro) ---
FM_MAX_K <- 200
fm_long <- build_cum_long(chosen_fm, FM_MAX_K)

# y-axis upper limit = max attained accuracy across fitted models at k=200
ymax <- max(0.1, fm_long$pct[fm_long$k == FM_MAX_K & fm_long$model != "Random"])

p_fm <- ggplot(fm_long, aes(x = k, y = pct, colour = model, linetype = model)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = c(1, 5, 10, 25, 50, 100, 150, 200)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, ymax)) +
  scale_colour_manual(values = palette) +
  scale_linetype_manual(values = linetypes) +
  labs(
    title = "Full-metro cumulative rank accuracy",
    subtitle = sprintf("Strata: %d | Avg men per metro: %d", nrow(chosen_fm), round(mean(n_men_chosen))),
    x = "Top k (among all men in metro)",
    y = "% of true partners ranked within top k",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p_fm)
if (exists("config_id")) {
  ggsave(file.path(results_dir, "plot_fullmetro.png"),
         plot = p_fm, width = 10, height = 6, dpi = 150)
}
