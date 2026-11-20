# ==============================================================================
# 00_helpers.R — Shared helper functions used across the pipeline
#
# This file defines three things:
#   1. Feature engineering functions that turn raw census data into model-ready
#      predictor variables (two encodings: v1 and v2)
#   2. Column selection functions that keep only the variables models need
#   3. A scoring function that generates predictions from fitted clogit models
#
# Source this file at the top of any script that needs these functions.
# Currently sourced by: 01_data.R (which is sourced by everything else)
# ==============================================================================


# ==============================================================================
# FEATURE ENGINEERING
#
# Two encodings of the same raw data — "v1" and "v2" — used by different models.
# Both describe how similar a (woman, man) pair is, but represent it differently:
#
#   v1: absolute differences + direction indicators
#       e.g., abs_age_diff=5, older_husband=1 → "5 years apart, he's older"
#       Best for linear models (M1_v1 clogit)
#
#   v2: signed differences + partner's absolute values
#       e.g., age_diff=+5, alter_age=35 → "he's 5 years older and is 35"
#       Best for nonlinear models (spline clogit, forest)
#
# The choice between v1 and v2 matters because linear models can only learn
# "bigger gaps are worse" from abs_diff, while splines and trees can learn
# "gaps of 3-5 years are preferred but gaps over 15 are strongly penalized"
# from the signed difference.
# ==============================================================================

# ==============================================================================
# v1 encoding (absolute differences + direction indicators)
# ==============================================================================

#' Build predictor variables from raw choice-set data (v1 encoding)
#'
#' Takes a data frame where each row is a (woman, potential partner) pair
#' and creates variables that describe how similar or different they are.
#'
#' The v1 encoding splits each difference into magnitude + direction:
#'   - abs_age_diff:  how many years apart (unsigned)
#'   - older_husband: 1 if the man is older, 0 if not
#'   - abs_educ_diff: how far apart in education (unsigned, ordinal scale)
#'   - hypergamy_ed:  1 if the man is more educated ("marrying up")
#'   - hypogamy_ed:   1 if the man is less educated ("marrying down")
#'   - same_race:     1 if both are the same race/ethnicity
#'
#' It also creates one-hot interaction dummies (rxr_* and exe_*) that give
#' the forest a separate binary variable for every specific race pairing
#' (e.g., rxr_White_Black = 1 if the woman is White and the man is Black)
#' and every education pairing. These turned out to add no predictive value
#' beyond same_race, but they're kept for completeness.
#'
#' @param df Data frame with columns: woman_age, age_partner, woman_educ_cat,
#'   educ_partner_cat, woman_educ_tier, educ_partner_tier, woman_race, race_partner
#' @return Data frame with all original columns plus engineered features
engineer_features <- function(df) {

  # --- Continuous and binary features ---
  # mutate() adds new columns to the data frame without removing existing ones
  out <- df %>%
    mutate(
      # Convert education categories to numeric for computing differences
      # factor() with ordered levels turns "1_LessHS" -> 1, "2_HS_Grad" -> 2, etc.
      woman_educ_num   = as.integer(factor(
        woman_educ_cat,
        levels = c("1_LessHS", "2_HS_Grad", "3_SomeColl", "4_Bachelor", "5_Advanced")
      )),
      educ_partner_num = as.integer(factor(
        educ_partner_cat,
        levels = c("1_LessHS", "2_HS_Grad", "3_SomeColl", "4_Bachelor", "5_Advanced")
      )),

      # Age difference: absolute value (how far apart, ignoring who is older)
      abs_age_diff = abs(age_partner - woman_age),
      # Direction: 1 if husband is older, 0 otherwise
      older_husband = as.integer(age_partner > woman_age),

      # Education difference: absolute value (how far apart on the 1-5 scale)
      abs_educ_diff = abs(educ_partner_num - woman_educ_num),
      # Direction: "hypergamy" = woman marries up (man more educated)
      hypergamy_ed = as.integer(educ_partner_num > woman_educ_num),
      # "hypogamy" = woman marries down (man less educated)
      hypogamy_ed  = as.integer(educ_partner_num < woman_educ_num),

      # Race match: 1 if both partners are the same race/ethnicity
      same_race    = as.integer(woman_race == race_partner),

      # Interaction terms as factors — these will be expanded into
      # one-hot dummies below (one column per unique pairing)
      race_x_race = interaction(woman_race, race_partner, sep = "_", drop = TRUE),
      educ_x_educ = interaction(woman_educ_cat, educ_partner_cat, sep = "_", drop = TRUE)
    )

  # Drop rows where education couldn't be coded (e.g., "Unknown" category)
  out <- out %>%
    filter(
      !is.na(woman_educ_num),
      !is.na(educ_partner_num),
      woman_educ_cat  != "Unknown",
      educ_partner_cat != "Unknown"
    )

  # --- One-hot encode the interaction terms ---
  # model.matrix() with ~ 0 + ... creates a binary (0/1) column for each
  # unique level of the factor, WITHOUT an intercept column.
  # This gives the forest individual split points for each race pairing
  # (e.g., rxr_White_Black, rxr_Hispanic_Hispanic, etc.)

  X_race <- model.matrix(~ 0 + race_x_race, data = out)
  # Rename columns from "race_x_raceWhite_Black" to "rxr_White_Black"
  colnames(X_race) <- sub("^race_x_race", "rxr_", colnames(X_race))

  X_educ <- model.matrix(~ 0 + educ_x_educ, data = out)
  colnames(X_educ) <- sub("^educ_x_educ", "exe_", colnames(X_educ))

  # Append the dummy columns to the data frame
  # bind_cols() is a tidyverse function that pastes columns side by side
  out <- bind_cols(out, as.data.frame(X_race), as.data.frame(X_educ))

  out
}


#' Keep only the columns that go into models (v1 encoding)
#'
#' Drops raw identifiers (SERIAL, PERNUM, etc.) and intermediate variables,
#' keeping only: choice (outcome), strata (group ID), and predictor columns.
#' This ensures the same columns are used in training and prediction.
#'
#' @param df Data frame output from engineer_features()
#' @param strata_col Name of the strata column (default "strata")
#' @return Data frame with choice, strata, and predictor columns only
select_model_cols <- function(df, strata_col = "strata") {
  # The 6 "main6" predictors used in the baseline clogit
  predictor_cols <- c(
    "abs_age_diff", "older_husband",
    "abs_educ_diff", "hypergamy_ed", "hypogamy_ed",
    "same_race"
  )
  # Also include all one-hot dummy columns (names starting with rxr_ or exe_)
  # grep() with value=TRUE returns the matching column names as strings
  dummy_cols <- grep("^(rxr_|exe_)", names(df), value = TRUE)

  keep <- c("choice", strata_col, predictor_cols, dummy_cols)
  df[, keep, drop = FALSE]
}


# ==============================================================================
# v2 encoding (signed differences + partner's absolute values)
# ==============================================================================

#' Build predictor variables from raw choice-set data (v2 encoding)
#'
#' The v2 encoding uses signed differences and the partner's absolute values:
#'   - age_diff:        husband_age - wife_age (positive = husband older)
#'   - alter_age:       husband's age (absolute, not relative to wife)
#'   - educ_diff:       husband_educ_rank - wife_educ_rank (positive = husband
#'                      more educated). Uses the ordinal rank (1-5) of education
#'                      categories, NOT years of schooling — we don't assume
#'                      equal spacing between education levels.
#'   - alter_educ_rank: husband's education rank (1-5 ordinal)
#'   - same_race:       1 if both same race
#'   + optional alter race dummies (alter_race_white, alter_race_black, etc.)
#'
#' Why two encodings? The v1 encoding (abs_diff + direction) was designed for
#' the linear clogit, where abs_age_diff captures "people prefer similar ages"
#' and older_husband captures the directional preference. But a tree can learn
#' more from signed differences + absolute levels: it can split on "man is
#' over 35" (alter_age > 35) which groups men by absolute age, regardless of
#' the woman's age. The clogit can't use alter_age without overfitting, but
#' the forest can.
#'
#' @param df Data frame with: woman_age, age_partner, woman_educ_cat,
#'   educ_partner_cat, woman_race, race_partner
#' @return Data frame with all original columns plus v2 features
engineer_features_v2 <- function(df) {

  # Education ordinal rank: 1 = less than HS, 2 = HS grad, ..., 5 = advanced degree
  # This is the same 1-5 scale used in v1 (woman_educ_num / educ_partner_num).
  # We use ordinal rank rather than "years of schooling" because education
  # categories are ordinal, not cardinal — the gap between "some college" and
  # "bachelor's" isn't necessarily equal to "HS" and "some college".
  educ_levels <- c("1_LessHS", "2_HS_Grad", "3_SomeColl", "4_Bachelor", "5_Advanced")

  out <- df %>%
    mutate(
      # Signed age difference: positive means husband is older
      # (e.g., +5 = husband is 5 years older, -3 = husband is 3 years younger)
      age_diff  = age_partner - woman_age,
      # The man's actual age (not relative to the woman)
      alter_age = age_partner,

      # Education ordinal ranks (1-5)
      woman_educ_rank   = as.integer(factor(woman_educ_cat, levels = educ_levels)),
      partner_educ_rank = as.integer(factor(educ_partner_cat, levels = educ_levels)),

      # Signed education difference on the ordinal scale
      # (e.g., +1 = husband one category higher, -2 = husband two categories lower)
      educ_diff       = partner_educ_rank - woman_educ_rank,
      # The man's education rank (1-5)
      alter_educ_rank = partner_educ_rank,

      # Race match: 1 if both are the same race/ethnicity
      same_race = as.integer(woman_race == race_partner),

      # Alter (= partner's) race as individual dummies
      # "Other" is the reference category (no dummy for it)
      alter_race_white    = as.integer(race_partner == "White"),
      alter_race_black    = as.integer(race_partner == "Black"),
      alter_race_hispanic = as.integer(race_partner == "Hispanic"),
      alter_race_asian    = as.integer(race_partner == "Asian")
    )

  # Remove rows where education couldn't be ranked (Unknown category → NA)
  out <- out %>%
    filter(!is.na(educ_diff), !is.na(alter_educ_rank)) %>%
    select(-woman_educ_rank, -partner_educ_rank)  # drop intermediate columns

  out
}


#' Keep only the columns that go into models (v2 encoding)
#'
#' @param df Data frame output from engineer_features_v2()
#' @param strata_col Name of the strata column (default "strata")
#' @param include_race_dummies If TRUE (default), include alter_race_* columns.
#'   Set to FALSE for the "v2_core" feature set (5 predictors only).
#' @return Data frame with choice, strata, and v2 predictor columns only
select_model_cols_v2 <- function(df, strata_col = "strata",
                                  include_race_dummies = TRUE) {
  core_cols <- c("age_diff", "alter_age", "educ_diff", "alter_educ_rank",
                 "same_race")
  if (include_race_dummies) {
    core_cols <- c(core_cols, "alter_race_white", "alter_race_black",
                   "alter_race_hispanic", "alter_race_asian")
  }
  keep <- c("choice", strata_col, core_cols)
  df[, keep, drop = FALSE]
}


# ==============================================================================
# SCORING: predict from a fitted conditional logit model
# ==============================================================================

#' Score a conditional logit model on new data
#'
#' Given a fitted clogit model and new data, returns exp(X * beta) — the
#' "relative utility" of each alternative. Higher values = model thinks this
#' alternative is more likely to be chosen. Within each choice set (strata),
#' the alternative with the highest score is the model's top prediction.
#'
#' WHY NOT JUST USE predict()? The standard predict.coxph() function calls
#' model.frame() internally, which checks that all factor levels in the new
#' data match the training data — INCLUDING the strata variable. But test
#' data has different strata IDs than training data, so predict() always
#' rejects it with "new levels" errors. Since the strata term doesn't
#' actually affect the prediction (it only partitions the likelihood during
#' model fitting), we simply strip it out and compute X*beta manually.
#'
#' HOW IT WORKS:
#'   1. Get the model's coefficients (beta vector)
#'   2. Extract the formula, remove the strata() term
#'   3. Build a design matrix X from the new data using model.matrix()
#'      (this handles factor expansion, polynomial terms, spline bases, etc.)
#'   4. Compute exp(X %*% beta) — matrix multiplication of design matrix
#'      by coefficient vector, then exponentiate
#'
#' @param model A fitted clogit (survival::clogit) model object
#' @param newdata Data frame with the same predictor columns as training data
#' @return Numeric vector of predicted relative utilities (one per row)
score_clogit <- function(model, newdata) {

  # Extract the fitted coefficients (one number per predictor)
  coefs <- coef(model)

  # Align factor levels in newdata to match what the model expects.
  # If the training data had race levels c("White", "Black", "Hispanic", "Asian", "Other"),
  # the new data must have the same levels in the same order, even if some
  # levels don't appear in the new data. model$xlevels stores this info.
  if (!is.null(model$xlevels)) {
    for (nm in names(model$xlevels)) {
      # Skip the strata term — we're going to remove it anyway
      if (grepl("^strata\\(", nm)) next
      if (nm %in% names(newdata)) {
        newdata[[nm]] <- factor(newdata[[nm]], levels = model$xlevels[[nm]])
      }
    }
  }

  # Extract the right-hand side of the model formula, dropping strata()
  # terms() parses the formula; attr(, "term.labels") gives the variable names
  tt <- terms(formula(model), specials = "strata")
  rhs_labels <- attr(tt, "term.labels")
  rhs_labels <- rhs_labels[!grepl("^strata\\(", rhs_labels)]
  # reformulate() rebuilds the formula from the remaining terms
  f_clean <- reformulate(rhs_labels)

  # Build the design matrix: this expands factors into dummies,
  # evaluates I() terms (like I(age_diff^2)), and applies spline
  # bases (ns(), bs()). The result is a numeric matrix where each
  # column corresponds to one model coefficient.
  mm <- model.matrix(f_clean, data = newdata)
  # Remove the intercept column if present (clogit has no intercept)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }

  # Matrix multiplication: X %*% beta gives the linear predictor for each row.
  # exp() converts to relative utilities (proportional to choice probabilities).
  # We select only the columns that match coefficient names, in case the design
  # matrix has extra columns the model didn't use.
  as.numeric(exp(mm[, names(coefs), drop = FALSE] %*% coefs))
}

