# ==============================================================================
# 05_stratified_analysis.R — Stratified models (race + education tier)
#
# Must be run AFTER source("pipeline/05_setup_data.R") — all data objects and
# helper functions must already be in the environment.
#
# CONTROL WHAT RUNS:
#   RUN_FOREST <- FALSE   # default — run clogit only (A + B), inspect results
#   RUN_FOREST <- TRUE    # then set this and re-source to also run Forest_core (C)
#
# Typical workflow:
#   1. source("pipeline/05_setup_data.R")
#   2. RUN_FOREST <- FALSE; source("discovery/05_stratified_analysis.R")
#      → inspect clogit results
#   3. RUN_FOREST <- TRUE;  source("discovery/05_stratified_analysis.R")
#      → now adds forest results
#
# What this does:
#   For each subgroup (race group / education tier), it:
#     (A) Fits M1_v1        — baseline clogit, v1 features (6 params)
#     (B) Fits M1_spline5   — spline clogit, v2 features (9 params)
#     (C) Fits Forest_core  — CLogitForest, v2 core features (5 predictors)
#                             Only runs if RUN_FOREST = TRUE AND clogit converged.
#
# Research questions:
#   • Does predictability vary across racial / education groups?
#   • Does the forest gain an advantage within a more homogeneous group?
#
# Outputs:
#   results/stratified_results.rds   — tidy results table (one row per group)
#   results/stratified_top1.png      — bar chart: top-1 accuracy by group × model
# ==============================================================================

library(tidyverse)
library(survival)
library(splines)
library(CLogitTree)
library(ggplot2)

# score_clogit() lives in 00_helpers.R (already sourced by 05_setup_data.R)
# add_spline_cols() was defined in 02_model.R / 05_setup_data.R

# ── RUN_FOREST flag ────────────────────────────────────────────────────────────
# Set to FALSE (default) to run clogit only and inspect results first.
# Set to TRUE to also fit Forest_core per group (adds ~5-30 min depending on group size).
if (!exists("RUN_FOREST")) RUN_FOREST <- FALSE
cat(sprintf("\nRUN_FOREST = %s  (%s)\n",
            RUN_FOREST,
            if (RUN_FOREST) "will fit Forest_core after each clogit"
            else            "clogit only — set RUN_FOREST <- TRUE to add forest"))


# ── 1. Verify the environment ──────────────────────────────────────────────────
# These objects must exist from run_pipe.R → 01_data.R / 02_model.R
required_objs <- c(
  "train_data_clogit", "test_data_clogit",  # v1 clogit track (has woman_race)
  "train_data_v2",     "test_data_v2",       # v2 track (for spline + forest)
  "N_ALTS", "SEED", "NTREE", "NCORES"
)
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing_objs) > 0)
  stop("Missing objects — run source('pipeline/run_pipe.R') first:\n  ",
       paste(missing_objs, collapse = ", "))

# Forest hyperparameters — use run_pipe.R values if present, else safe defaults
if (!exists("MINNODESIZE"))  MINNODESIZE  <- 550L
if (!exists("DEPTH_MAX"))    DEPTH_MAX    <- 10L
if (!exists("TREE_TIMEOUT")) TREE_TIMEOUT <- 600L
if (!exists("EPSILON"))      EPSILON      <- 1e-3
if (!exists("ITER_MAX"))     ITER_MAX     <- 5L


# ── 2. Reconstruct full (pre-split) data ──────────────────────────────────────
# run_pipe.R split into train/test. Binding them back gives the full data,
# which we then re-split PER GROUP so each group has its own 70/30 split.

all_clogit <- bind_rows(train_data_clogit, test_data_clogit)  # v1 features + woman_race
all_v2     <- bind_rows(train_data_v2,     test_data_v2)      # v2 features (signed diffs)

# Derive education tier from woman_educ_cat (not stored separately in clogit track)
all_clogit <- all_clogit %>%
  mutate(woman_educ_tier = case_when(
    woman_educ_cat %in% c("1_LessHS", "2_HS_Grad")    ~ "Low",
    woman_educ_cat == "3_SomeColl"                     ~ "Mid",
    woman_educ_cat %in% c("4_Bachelor", "5_Advanced")  ~ "High",
    TRUE ~ NA_character_
  ))

# One row per stratum — use the chosen row (choice==1) to identify each woman
strata_attrs <- all_clogit %>%
  filter(choice == 1) %>%
  select(strata, woman_race, woman_educ_tier) %>%
  mutate(woman_race = as.character(woman_race)) %>%
  distinct()


# ── 3. Sample size check — BEFORE any modelling ───────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("STEP 0: Sample size check\n")
cat("══════════════════════════════════════════════════════\n")

cat(sprintf("\nTotal strata in data: %d  |  N_ALTS: %d  |  Random top-1: %.1f%%\n",
            nrow(strata_attrs), N_ALTS, 100 / N_ALTS))

cat("\nRace group sizes (n strata):\n")
race_counts <- strata_attrs %>%
  count(woman_race, name = "n_strata") %>%
  arrange(desc(n_strata)) %>%
  mutate(
    pct         = round(n_strata / sum(n_strata) * 100, 1),
    train_70pct = floor(n_strata * 0.7),
    test_30pct  = n_strata - floor(n_strata * 0.7),
    feasible    = ifelse(n_strata >= 30, "YES", "SKIP (<30)")
  )
print(race_counts)

cat("\nEducation tier sizes (n strata):\n")
educ_counts <- strata_attrs %>%
  count(woman_educ_tier, name = "n_strata") %>%
  arrange(desc(n_strata)) %>%
  mutate(
    pct         = round(n_strata / sum(n_strata) * 100, 1),
    train_70pct = floor(n_strata * 0.7),
    test_30pct  = n_strata - floor(n_strata * 0.7),
    feasible    = ifelse(n_strata >= 30, "YES", "SKIP (<30)")
  )
print(educ_counts)

cat("\nGroups with < 30 strata will be skipped.\n")
cat("Groups with < 100 strata may produce noisy forest results — interpret carefully.\n")


# ── 4. Core fitting + evaluation function ─────────────────────────────────────
#
# fit_and_eval_group(group_label, strata_ids)
#
# Runs steps A → B → C in sequence:
#   A: M1_v1 (baseline clogit, v1 features)
#   B: M1_spline5 (spline clogit, v2 features — spline basis re-fitted per group)
#   C: Forest_core (CLogitForest, v2 core 5 predictors — only if A succeeded)
#
# Returns a one-row tibble with all metrics and group metadata.
# ------------------------------------------------------------------------------
fit_and_eval_group <- function(group_label, strata_ids) {

  cat(sprintf("\n══════ %s ══════\n", group_label))

  # Filter both tracks to this group's strata
  clogit_g <- all_clogit %>% filter(strata %in% strata_ids) %>% as.data.frame()
  v2_g     <- all_v2     %>% filter(strata %in% strata_ids) %>% as.data.frame()

  n_total <- n_distinct(clogit_g$strata)
  cat(sprintf("  Strata: %d  |  Rows (v1): %d  |  Rows (v2): %d\n",
              n_total, nrow(clogit_g), nrow(v2_g)))

  # Minimum strata for clogit to converge reliably
  MIN_STRATA <- 30L
  if (n_total < MIN_STRATA) {
    cat(sprintf("  [SKIP] Too few strata (need >= %d)\n", MIN_STRATA))
    return(tibble(
      group = group_label, n_strata = n_total, n_test = NA_integer_,
      skipped = TRUE,
      m1_top1 = NA_real_,     m1_top3 = NA_real_,     m1_mean_rank = NA_real_,
      spline_top1 = NA_real_, spline_top3 = NA_real_, spline_mean_rank = NA_real_,
      forest_top1 = NA_real_, forest_top3 = NA_real_, forest_mean_rank = NA_real_
    ))
  }

  # ── Re-split 70/30 within this group ──────────────────────────────────────
  all_ids_g <- unique(clogit_g$strata)
  set.seed(SEED)
  tr_ids <- sample(all_ids_g, floor(0.7 * length(all_ids_g)))
  te_ids <- setdiff(all_ids_g, tr_ids)

  tr_clogit <- clogit_g %>% filter(strata %in% tr_ids) %>% as.data.frame()
  te_clogit <- clogit_g %>% filter(strata %in% te_ids) %>% as.data.frame()
  tr_v2     <- v2_g     %>% filter(strata %in% tr_ids) %>% as.data.frame()
  te_v2     <- v2_g     %>% filter(strata %in% te_ids) %>% as.data.frame()

  cat(sprintf("  Split: %d train | %d test strata\n", length(tr_ids), length(te_ids)))

  # Ensure factor levels for clogit
  for (col in c("woman_race", "race_partner", "woman_educ_cat", "educ_partner_cat")) {
    if (col %in% names(tr_clogit)) {
      tr_clogit[[col]] <- factor(tr_clogit[[col]])
      te_clogit[[col]] <- factor(te_clogit[[col]])
    }
  }

  # Adaptive forest node size: ~4% of training rows (same relative constraint
  # as the main run: MINNODESIZE=550 on 14000 rows ≈ 3.9%)
  n_train_rows    <- nrow(tr_clogit)
  grp_minnodesize <- max(N_ALTS * 5L, round(n_train_rows * 0.04))
  grp_minbucket   <- ceiling(grp_minnodesize / 3L)

  # Results container — pre-filled with NA so skipped steps leave NA not error
  results <- tibble(
    group = group_label, n_strata = n_total, n_test = length(te_ids), skipped = FALSE,
    m1_top1 = NA_real_,     m1_top3 = NA_real_,     m1_mean_rank = NA_real_,
    spline_top1 = NA_real_, spline_top3 = NA_real_, spline_mean_rank = NA_real_,
    forest_top1 = NA_real_, forest_top3 = NA_real_, forest_mean_rank = NA_real_
  )
  clogit_ok <- FALSE  # gate for forest: only run if clogit converged


  # ── STEP A: M1_v1 — Baseline CLogit on v1 features ────────────────────────
  cat("  [A] M1_v1 (baseline clogit, 6 params)...\n")

  m1 <- tryCatch(
    clogit(
      choice ~ abs_age_diff + older_husband + abs_educ_diff +
               hypergamy_ed + hypogamy_ed + same_race + strata(strata),
      data   = tr_clogit,
      method = "efron"
    ),
    error = function(e) { cat("      Error:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(m1)) {
    te_clogit$pred_m1 <- score_clogit(m1, te_clogit)
    ranks_m1 <- te_clogit %>%
      group_by(strata) %>%
      mutate(r = rank(-pred_m1, ties.method = "average")) %>%
      filter(choice == 1) %>%
      pull(r)
    results$m1_top1      <- mean(ranks_m1 == 1)
    results$m1_top3      <- mean(ranks_m1 <= 3)
    results$m1_mean_rank <- mean(ranks_m1)
    cat(sprintf("      top-1=%.1f%%  top-3=%.1f%%  mean_rank=%.2f\n",
                100 * results$m1_top1, 100 * results$m1_top3, results$m1_mean_rank))
    clogit_ok <- TRUE
  }


  # ── STEP B: M1_spline5 — Spline CLogit on v2 features ────────────────────
  # Re-fit spline bases on this group's training data so knots reflect
  # the group's own age/education distribution (not the full-sample knots).
  cat("  [B] M1_spline5 (spline clogit, 9 params)...\n")

  grp_basis_age  <- NULL
  grp_basis_educ <- NULL

  m1_sp <- tryCatch({
    grp_basis_age  <- ns(tr_v2$age_diff,  df = 5)
    grp_basis_educ <- ns(tr_v2$educ_diff, df = 3)

    # add_spline_cols() was defined in 02_model.R and is in the environment
    tr_v2_sp <- add_spline_cols(tr_v2, grp_basis_age, grp_basis_educ)

    clogit(
      choice ~ age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
               educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
               same_race + strata(strata),
      data   = tr_v2_sp,
      method = "efron"
    )
  }, error = function(e) { cat("      Error:", conditionMessage(e), "\n"); NULL })

  if (!is.null(m1_sp) && !is.null(grp_basis_age)) {
    # Apply the SAME knots to the test data (crucial — don't refit on test)
    te_v2_sp <- add_spline_cols(te_v2, grp_basis_age, grp_basis_educ)

    te_v2_sp$pred_spline <- score_clogit(m1_sp, te_v2_sp)
    ranks_sp <- te_v2_sp %>%
      group_by(strata) %>%
      mutate(r = rank(-pred_spline, ties.method = "average")) %>%
      filter(choice == 1) %>%
      pull(r)
    results$spline_top1      <- mean(ranks_sp == 1)
    results$spline_top3      <- mean(ranks_sp <= 3)
    results$spline_mean_rank <- mean(ranks_sp)
    cat(sprintf("      top-1=%.1f%%  top-3=%.1f%%  mean_rank=%.2f\n",
                100 * results$spline_top1, 100 * results$spline_top3,
                results$spline_mean_rank))
  }


  # ── STEP C: Forest_core — only if RUN_FOREST=TRUE AND clogit converged ────
  if (!RUN_FOREST) {
    cat("  [C] Forest_core skipped (RUN_FOREST=FALSE).\n")
    return(results)
  }
  if (!clogit_ok) {
    cat("  [C] Forest_core skipped — clogit did not converge in step A.\n")
    return(results)
  }

  cat(sprintf("  [C] Forest_core (5 v2 predictors, mtry=2)...\n"))
  cat(sprintf("      ntree=%d  minnodesize=%d  minbucket=%d  ncores=%d\n",
              NTREE, grp_minnodesize, grp_minbucket, NCORES))

  # v2 core: 5 predictors only (same as main Forest_core)
  core_cols <- c("choice", "strata",
                 "age_diff", "alter_age", "educ_diff", "alter_educ_rank", "same_race")
  tr_core <- tr_v2[, core_cols, drop = FALSE]
  te_core <- te_v2[, core_cols, drop = FALSE]

  rf <- tryCatch({
    set.seed(SEED)
    t0 <- proc.time()
    f  <- CLogitForest(
      data          = tr_core,
      response      = "choice",
      s             = "strata",
      ntree         = NTREE,
      minnodesize   = grp_minnodesize,
      minbucket     = grp_minbucket,
      depth_max     = if (DEPTH_MAX == 0L) NULL else DEPTH_MAX,
      mtry          = 2L,        # sqrt(5) ≈ 2, same as main Forest_core
      tune.mtry     = FALSE,
      subsample.632 = TRUE,
      ncores        = NCORES,
      epsilon       = EPSILON,
      iter.max      = ITER_MAX,
      lean          = TRUE,
      tree.timeout  = TREE_TIMEOUT
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]
    n_valid <- length(Filter(function(x) !is.null(x) && !is.null(x$tree.b), f$tree_list))
    cat(sprintf("      Done: %.0fs | Valid trees: %d/%d\n", elapsed, n_valid, NTREE))
    if (n_valid == 0L) stop("All trees failed or timed out")
    f
  }, error = function(e) { cat("      Forest error:", conditionMessage(e), "\n"); NULL })

  # Close any stray PSOCK sockets CLogitForest opened
  open_cons <- showConnections(all = FALSE)
  if (nrow(open_cons) > 0) {
    sock_idx <- which(open_cons[, "class"] == "sockconn")
    for (i in sock_idx)
      try(close(getConnection(as.integer(rownames(open_cons)[i]))), silent = TRUE)
  }

  if (!is.null(rf)) {
    te_core$pred_forest <- as.numeric(
      predict(rf, type = "response", response = "choice",
              s = "strata", newdata = te_core)
    )
    ranks_f <- te_core %>%
      group_by(strata) %>%
      mutate(r = rank(-pred_forest, ties.method = "average")) %>%
      filter(choice == 1) %>%
      pull(r)
    results$forest_top1      <- mean(ranks_f == 1)
    results$forest_top3      <- mean(ranks_f <= 3)
    results$forest_mean_rank <- mean(ranks_f)
    cat(sprintf("      top-1=%.1f%%  top-3=%.1f%%  mean_rank=%.2f\n",
                100 * results$forest_top1, 100 * results$forest_top3,
                results$forest_mean_rank))
  }

  results
}


# ── 5. PART 1: Race-stratified models ─────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════════════\n")
cat("PART 1: Race-Stratified Models\n")
cat("══════════════════════════════════════════════════════\n")
cat("Each race group gets its own train/test split and independently fitted models.\n")
cat("Question: Does predictability vary by race? Does the forest gain an edge?\n\n")

race_groups <- c("White", "Black", "Hispanic", "Asian", "Other")

race_results <- purrr::map_dfr(race_groups, function(grp) {
  ids <- strata_attrs %>% filter(woman_race == grp) %>% pull(strata)
  fit_and_eval_group(paste0("Race: ", grp), ids)
}) %>%
  mutate(stratification = "Race")

cat("\n─── Race Results ───────────────────────────────────\n")
race_results %>%
  filter(!skipped) %>%
  transmute(
    group,
    n_strata,
    n_test,
    M1_v1_top1     = paste0(round(m1_top1     * 100, 1), "%"),
    Spline_top1    = paste0(round(spline_top1 * 100, 1), "%"),
    Forest_top1    = paste0(round(forest_top1 * 100, 1), "%"),
    Forest_vs_M1   = paste0(round((forest_top1 - m1_top1) * 100, 1), " pp"),
    M1_mean_rank   = round(m1_mean_rank,     2),
    Forest_rank    = round(forest_mean_rank, 2)
  ) %>%
  print(n = 100)


# ── 6. PART 2: Education-tier stratified models ───────────────────────────────
cat("\n\n══════════════════════════════════════════════════════\n")
cat("PART 2: Education-Tier Stratified Models\n")
cat("Low = <HS + HS grad  |  Mid = some college  |  High = bachelor+\n")
cat("══════════════════════════════════════════════════════\n")
cat("Question: Is partner choice more predictable among the highly educated?\n\n")

educ_groups <- c("Low", "Mid", "High")

educ_results <- purrr::map_dfr(educ_groups, function(grp) {
  ids <- strata_attrs %>% filter(woman_educ_tier == grp) %>% pull(strata)
  fit_and_eval_group(paste0("Educ: ", grp), ids)
}) %>%
  mutate(stratification = "Education")

cat("\n─── Education-Tier Results ─────────────────────────\n")
educ_results %>%
  filter(!skipped) %>%
  transmute(
    group,
    n_strata,
    n_test,
    M1_v1_top1   = paste0(round(m1_top1     * 100, 1), "%"),
    Spline_top1  = paste0(round(spline_top1 * 100, 1), "%"),
    Forest_top1  = paste0(round(forest_top1 * 100, 1), "%"),
    Forest_vs_M1 = paste0(round((forest_top1 - m1_top1) * 100, 1), " pp"),
    M1_mean_rank = round(m1_mean_rank,     2),
    Forest_rank  = round(forest_mean_rank, 2)
  ) %>%
  print(n = 100)


# ── 7. Combined summary + plot ─────────────────────────────────────────────────
all_results <- bind_rows(race_results, educ_results)
random_top1 <- 100 / N_ALTS

cat("\n\n══════════════════════════════════════════════════════\n")
cat("COMBINED SUMMARY TABLE\n")
cat(sprintf("Random baseline: %.1f%%  |  N_ALTS=%d  |  NTREE=%d\n",
            random_top1, N_ALTS, NTREE))
cat("══════════════════════════════════════════════════════\n")

all_results %>%
  filter(!skipped) %>%
  transmute(
    stratification, group, n_strata, n_test,
    M1_pct      = round(m1_top1     * 100, 1),
    Spline_pct  = round(spline_top1 * 100, 1),
    Forest_pct  = round(forest_top1 * 100, 1),
    Forest_vs_M1_pp = round((forest_top1 - m1_top1) * 100, 1)
  ) %>%
  print(n = 100)

# Bar chart: top-1 accuracy by group and model
plot_data <- all_results %>%
  filter(!skipped) %>%
  select(stratification, group, m1_top1, spline_top1, forest_top1) %>%
  pivot_longer(
    cols      = c(m1_top1, spline_top1, forest_top1),
    names_to  = "model",
    values_to = "top1"
  ) %>%
  mutate(
    model = recode(model,
                   "m1_top1"     = "M1_v1 (CLogit)",
                   "spline_top1" = "M1_spline5",
                   "forest_top1" = "Forest_core"),
    model    = factor(model, levels = c("M1_v1 (CLogit)", "M1_spline5", "Forest_core")),
    top1_pct = top1 * 100,
    group    = factor(group, levels = rev(unique(group)))
  )

p <- ggplot(plot_data, aes(x = group, y = top1_pct, fill = model)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = random_top1, linetype = "dashed",
             colour = "grey50", linewidth = 0.8) +
  annotate("text", x = Inf, y = random_top1 + 0.8,
           label = sprintf("Random (%.0f%%)", random_top1),
           hjust = 1.05, size = 3.2, colour = "grey40") +
  facet_wrap(~ stratification, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(
    "M1_v1 (CLogit)" = "#e11d48",
    "M1_spline5"     = "#f97316",
    "Forest_core"    = "#2563eb"
  )) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand  = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Top-1 Accuracy by Subgroup and Model",
    subtitle = sprintf(
      "Each group fitted independently — N_ALTS=%d | NTREE=%d | Random baseline=%.0f%%",
      N_ALTS, NTREE, random_top1
    ),
    x    = NULL,
    y    = "Top-1 accuracy (%)",
    fill = "Model"
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text      = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(p)

# Save results next to main pipeline results
out_dir <- if (exists("results_dir")) results_dir else "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "stratified_top1.png"), p, width = 9, height = 8, dpi = 150)
saveRDS(all_results, file.path(out_dir, "stratified_results.rds"))

cat(sprintf("\nPlot saved : %s\n", file.path(out_dir, "stratified_top1.png")))
cat(sprintf("Results saved: %s\n", file.path(out_dir, "stratified_results.rds")))

cat("\n══════════════════════════════════════════════════════\n")
cat("05_stratified_analysis.R complete.\n")
cat("══════════════════════════════════════════════════════\n")
