# ==============================================================================
# 02_model.R — Fit all models: two clogit benchmarks + two forests
#
# This script fits four models on the training data:
#
#   1. M1_v1 (6 params) — conditional logit using the v1 encoding (absolute
#      differences + direction indicators): abs_age_diff, older_husband,
#      abs_educ_diff, hypergamy_ed, hypogamy_ed, same_race.
#      This is the standard linear benchmark.
#
#   2. M1_spline (9 params) — conditional logit using the v2 encoding (signed
#      differences) with natural splines: df=5 on age_diff (continuous) and
#      df=3 on educ_diff (ordinal 1-5 scale), plus same_race. Natural splines
#      are smooth curves that capture nonlinear preferences (e.g., "small age
#      gaps are fine but large ones are penalized") that the linear M1_v1 misses.
#
#   3. Forest_core (5 predictors, mtry=2) — conditional logit random forest using
#      v2 encoding "core" features: age_diff, alter_age, educ_diff,
#      alter_educ_yrs, same_race. Each tree partitions the data and fits a
#      conditional logit in each leaf. The forest averages across 50 trees.
#
#   4. Forest_full (9 predictors, mtry=5) — same as Forest_core but adds 4
#      alter race dummies from the v2 encoding (alter_race_white/black/hispanic/
#      asian). Tests whether the forest can learn race-specific preferences
#      beyond the same_race binary.
#
# The "v1" and "v2" labels refer to two different predictor encodings created
# in 01_data.R. See the comments in STEP 6 of that file for a full explanation.
#
# WHAT IS A CONDITIONAL LOGIT (clogit)?
#   A regression model for discrete choice data. Within each choice set (strata),
#   it estimates which characteristics make an alternative more likely to be chosen.
#   The key feature: it conditions on the chooser's fixed traits (absorbed by the
#   strata fixed effect), so only *alternative-varying* characteristics matter.
#
# WHAT IS CLogitForest?
#   A random forest where each tree is a conditional logit tree. At each internal
#   node, the tree finds the best binary split (e.g., age_diff > 5) that maximizes
#   the conditional logit likelihood in the child nodes. The forest averages
#   predictions across many such trees, each fitted on a bootstrap subsample.
#
# EXPECTS from 01_data.R:
#   train_data, train_data_clogit, train_data_v2, train_strata
# EXPECTS from run_pipe.R:
#   SEED, NTREE, MINNODESIZE, MINBUCKET, NCORES, DEPTH_MAX,
#   TREE_TIMEOUT, EPSILON, ITER_MAX
# ==============================================================================

library(splines)  # for ns() — natural spline basis functions

cat(sprintf("Training strata: %d (v1: %d rows, v2: %d rows)\n",
            length(train_strata), nrow(train_data), nrow(train_data_v2)))

# FIT_CLOGIT and FIT_FOREST are set by run_pipe.R from the MODELS parameter.
# If 02_model.R is sourced standalone (e.g. during development), default to
# fitting both model families.
if (!exists("FIT_CLOGIT")) FIT_CLOGIT <- TRUE
if (!exists("FIT_FOREST")) FIT_FOREST <- TRUE

# Initialise fitted model objects to NULL so 03_evaluation.R can check
# existence (`!is.null(m1_v1)`) regardless of which families were fit.
m1_v1        <- NULL
m1_spline5   <- NULL
rf_core      <- NULL
rf_full      <- NULL
basis_age5   <- NULL
basis_educ3  <- NULL
add_spline_cols <- function(df, ...) df   # no-op; redefined if clogit fit

cat(sprintf("\nMODELS: %s (clogit=%s, forest=%s)\n",
            if (FIT_CLOGIT && FIT_FOREST) "both" else if (FIT_CLOGIT) "clogit" else "forest",
            FIT_CLOGIT, FIT_FOREST))


# ==============================================================================
# CLOGIT FAMILY: M1_v1 + M1_spline
# ==============================================================================

if (FIT_CLOGIT) {

# ==============================================================================
# MODEL 1: M1_v1 — Linear conditional logit (6 parameters)
# ==============================================================================

cat("\n--- M1_v1: Baseline CLogit (6 params) ---\n")

# The 6 predictor columns used in the baseline model
main6_cols <- c("abs_age_diff", "older_husband", "abs_educ_diff",
                "hypergamy_ed", "hypogamy_ed", "same_race")

# clogit() requires factor columns to have proper levels (not just character strings)
train_data_clogit$woman_race       <- factor(train_data_clogit$woman_race)
train_data_clogit$race_partner     <- factor(train_data_clogit$race_partner)
train_data_clogit$woman_educ_cat   <- factor(train_data_clogit$woman_educ_cat)
train_data_clogit$educ_partner_cat <- factor(train_data_clogit$educ_partner_cat)

# The formula:
#   choice ~ predictors + strata(strata)
# - choice is the outcome (1 = chosen, 0 = not chosen)
# - strata(strata) tells clogit which rows belong to the same choice set
# - method="efron" is a correction for tied event times (standard choice)
f_m1_v1 <- choice ~ abs_age_diff + older_husband + abs_educ_diff +
  hypergamy_ed + hypogamy_ed + same_race + strata(strata)
m1_v1 <- clogit(f_m1_v1, data = train_data_clogit, method = "efron")

cat("Coefficients:\n")
print(round(coef(m1_v1), 4))
# Expected signs: abs_age_diff negative (prefer similar age),
# older_husband positive (slight preference for older men),
# same_race positive (strong preference for same race)


# ==============================================================================
# MODEL 2: M1_spline5 — Spline conditional logit (11 parameters)
# ==============================================================================

cat("\n--- M1_spline: Spline CLogit (9 params) ---\n")

# We need the v2 data (which has signed age_diff and ordinal educ_diff) for
# the splines. Also add the v1 columns (abs_age_diff etc.) so the same data
# frame can be used for scoring M1_v1 during evaluation.
train_v2_aug <- train_data_v2 %>%
  mutate(
    abs_age_diff  = abs(age_diff),
    older_husband = as.integer(age_diff > 0),
    abs_educ_diff = abs(educ_diff),
    hypergamy_ed  = as.integer(educ_diff > 0),
    hypogamy_ed   = as.integer(educ_diff < 0)
  )

# PRE-COMPUTE SPLINE BASES on the training data.
# ns() (natural splines) creates smooth basis functions. The spline is linear
# beyond the boundary knots (no wild extrapolation). These basis objects are
# saved and reused on test/full-metro data via predict() to ensure the SAME
# knot positions are used everywhere.
#
# Age:       df=5 (continuous variable, many unique values, rich shape to learn)
# Education: df=3 (ordinal 1-5, only 9 possible difference values -4 to +4,
#            so df=3 is appropriate — df=5 would nearly saturate the scale)
basis_age5  <- ns(train_v2_aug$age_diff,  df = 5)
basis_educ3 <- ns(train_v2_aug$educ_diff, df = 3)

# Helper: apply saved spline bases to any data frame
# predict(basis_obj, new_values) evaluates the spline at new_values
# using the knots determined from the original training data
add_spline_cols <- function(df, basis_age, basis_educ) {
  a <- predict(basis_age,  df$age_diff)
  e <- predict(basis_educ, df$educ_diff)
  for (k in 1:ncol(a)) df[[paste0("age_ns5_", k)]]  <- a[, k]
  for (k in 1:ncol(e)) df[[paste0("educ_ns3_", k)]] <- e[, k]
  df
}

train_v2_aug <- add_spline_cols(train_v2_aug, basis_age5, basis_educ3)

# Formula: 5 age spline columns + 3 educ spline columns + same_race = 9 params
f_spline <- choice ~ age_ns5_1 + age_ns5_2 + age_ns5_3 + age_ns5_4 + age_ns5_5 +
  educ_ns3_1 + educ_ns3_2 + educ_ns3_3 +
  same_race + strata(strata)

m1_spline5 <- clogit(f_spline, data = train_v2_aug, method = "efron")

cat("Coefficients:\n")
print(round(coef(m1_spline5), 4))
# Spline coefficients don't have intuitive interpretations individually —
# they jointly define a smooth curve. The same_race coefficient should be
# similar to M1_v1's (~1.9).

}  # end FIT_CLOGIT


# ==============================================================================
# FOREST FAMILY: Forest_core + Forest_full
# ==============================================================================

if (FIT_FOREST) {

# ==============================================================================
# MODEL 3: Forest v2_core — 5 predictors, mtry=2
# ==============================================================================

cat("\n--- Forest v2_core (5 pred, mtry=2) ---\n")

MTRY_CORE <- 2  # at each split, randomly select 2 of 5 predictors to consider

# Prepare training data: only the 5 core v2 predictors + choice + strata
train_core <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = FALSE) %>%
  as.data.frame()

# Set defaults for any hyperparameters not set by run_pipe.R
# (this allows 02_model.R to be sourced standalone for testing)
if (!exists("MINNODESIZE"))  MINNODESIZE  <- 550
if (!exists("MINBUCKET"))    MINBUCKET    <- ceiling(MINNODESIZE / 3)
if (!exists("DEPTH_MAX"))    DEPTH_MAX    <- 10
if (!exists("NCORES"))       NCORES       <- 2L
if (!exists("TREE_TIMEOUT")) TREE_TIMEOUT <- 600
if (!exists("EPSILON"))      EPSILON      <- 1e-3
if (!exists("ITER_MAX"))     ITER_MAX     <- 5

cat(sprintf("  ntree=%d  mns=%d  mtry=%d  ncores=%d  timeout=%ds\n",
            NTREE, MINNODESIZE, MTRY_CORE, NCORES, TREE_TIMEOUT))

set.seed(SEED)
t0 <- proc.time()

# CLogitForest parameters explained:
#   data:          training data (choice + strata + predictors)
#   response:      name of the outcome column ("choice")
#   s:             name of the strata column ("strata")
#   ntree:         number of trees in the forest (50)
#   minnodesize:   minimum number of ROWS (not strata) in a node before it
#                  can be split further. Controls tree depth. With N_ALTS=10,
#                  each strata has 11 rows, so mns=550 means ~50 strata/node.
#   minbucket:     minimum ROWS in a leaf node (usually mns/3)
#   depth_max:     maximum tree depth (10 = hard ceiling, rarely binding)
#   mtry:          number of randomly selected predictors to consider at each
#                  split. Lower = more diverse trees = better ensemble.
#   tune.mtry:     if TRUE, automatically select mtry (we set it manually)
#   subsample.632: each tree uses a 63.2% bootstrap subsample (standard RF)
#   ncores:        parallel workers for tree fitting (each gets a full data copy)
#   epsilon:       convergence tolerance for the Newton-Raphson optimizer inside
#                  each split evaluation (looser = faster, fine for split selection)
#   iter.max:      max Newton-Raphson iterations per split evaluation
#   lean:          if TRUE, don't store full data in tree objects (saves memory)
#   tree.timeout:  max seconds per tree (kills stuck trees)
rf_core <- CLogitForest(
  data          = train_core,
  response      = "choice",
  s             = "strata",
  ntree         = NTREE,
  minnodesize   = MINNODESIZE,
  minbucket     = MINBUCKET,
  depth_max     = if (DEPTH_MAX == 0) NULL else DEPTH_MAX,
  mtry          = MTRY_CORE,
  tune.mtry     = FALSE,
  subsample.632 = TRUE,
  ncores        = NCORES,
  epsilon       = EPSILON,
  iter.max      = ITER_MAX,
  lean          = TRUE,
  tree.timeout  = TREE_TIMEOUT
)

elapsed_core <- (proc.time() - t0)[["elapsed"]]
# Count valid trees (some may have timed out or failed, returning tree.b = NULL)
n_valid_core <- length(Filter(function(x) !is.null(x) && !is.null(x$tree.b), rf_core$tree_list))
cat(sprintf("  Done: %.0fs, %d/%d valid trees\n", elapsed_core, n_valid_core, NTREE))

if (n_valid_core == 0) stop("Forest_core: all trees failed or timed out. Check MINNODESIZE and TREE_TIMEOUT.")

# Clean up PSOCK parallel worker connections (prevents zombie processes)
open_cons <- showConnections(all = FALSE)
if (nrow(open_cons) > 0) {
  sock_idx <- which(open_cons[, "class"] == "sockconn")
  for (i in sock_idx) try(close(getConnection(as.integer(rownames(open_cons)[i]))), silent = TRUE)
}


# ==============================================================================
# MODEL 4: Forest v2_full — 9 predictors, mtry=5
# ==============================================================================

cat("\n--- Forest v2_full (9 pred, mtry=5) ---\n")

MTRY_FULL <- 5  # higher mtry because there are more predictors (9 vs 5)

# Same as core but includes alter_race_* dummies
train_full <- train_data_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>%
  as.data.frame()

cat(sprintf("  ntree=%d  mns=%d  mtry=%d  ncores=%d  timeout=%ds\n",
            NTREE, MINNODESIZE, MTRY_FULL, NCORES, TREE_TIMEOUT))

# Use a different seed so the two forests have independent bootstrap samples
set.seed(SEED + 1)
t0 <- proc.time()

rf_full <- CLogitForest(
  data          = train_full,
  response      = "choice",
  s             = "strata",
  ntree         = NTREE,
  minnodesize   = MINNODESIZE,
  minbucket     = MINBUCKET,
  depth_max     = if (DEPTH_MAX == 0) NULL else DEPTH_MAX,
  mtry          = MTRY_FULL,
  tune.mtry     = FALSE,
  subsample.632 = TRUE,
  ncores        = NCORES,
  epsilon       = EPSILON,
  iter.max      = ITER_MAX,
  lean          = TRUE,
  tree.timeout  = TREE_TIMEOUT
)

elapsed_full <- (proc.time() - t0)[["elapsed"]]
n_valid_full <- length(Filter(function(x) !is.null(x) && !is.null(x$tree.b), rf_full$tree_list))
cat(sprintf("  Done: %.0fs, %d/%d valid trees\n", elapsed_full, n_valid_full, NTREE))

if (n_valid_full == 0) stop("Forest_full: all trees failed or timed out. Check MINNODESIZE and TREE_TIMEOUT.")

# Clean up parallel connections
open_cons <- showConnections(all = FALSE)
if (nrow(open_cons) > 0) {
  sock_idx <- which(open_cons[, "class"] == "sockconn")
  for (i in sock_idx) try(close(getConnection(as.integer(rownames(open_cons)[i]))), silent = TRUE)
}

cat(sprintf("\nTotal forest fitting time: %.1f min\n",
            (elapsed_core + elapsed_full) / 60))

}  # end FIT_FOREST
