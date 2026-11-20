# ==============================================================================
# 01_data.R — Data loading, cleaning, and choice set generation
#
# This is the data preparation step of the pipeline. It:
#   1. Loads US census microdata from an IPUMS extract
#   2. Identifies recently married couples in metro areas
#   3. Builds "choice sets" for each woman: her actual partner + random
#      alternative men from the same metro area and year
#   4. Engineers predictor variables (two encodings: v1 and v2)
#   5. Splits into training and test sets
#
# WHAT IS A CHOICE SET?
#   For each focal woman, we create a group of N_ALTS men: her true partner
#   (choice = 1) and N_ALTS - 1 randomly sampled men from the same marriage
#   market (choice = 0). Each group is called a "stratum" or "choice set".
#   The model's job is to identify the true partner from the alternatives.
#
# WHAT IS A MARRIAGE MARKET?
#   We define it as a metro area × year combination. A woman in LA in 2019
#   "chooses" from men in LA in 2019. This is a simplification — real markets
#   overlap and people move — but it's standard in the literature.
#
# EXPECTS these variables to be set before sourcing (set by run_pipe.R):
#   SEED, N_WOMEN, N_ALTS, DDI_FILE, DO_SPLIT, TRAIN_SHARE
#   Optional: METRO_CODE, METRO_CODES, YEAR_FILTER, MIN_MEN_PER_METRO
#
# PRODUCES these objects in the R environment (used by 02_model.R):
#   train_data / test_data                — v1 numeric data (for forest)
#   train_data_clogit / test_data_clogit  — v1 + factor columns (for clogit)
#   train_data_v2 / test_data_v2          — v2 data (for spline clogit + forest)
#   train_strata / test_strata            — which strata IDs are in each set
#   sample_women_base                     — the focal women (needed for full-metro eval)
#   men_pool                              — all men in the data (needed for full-metro eval)
# ==============================================================================

# --- Load required packages ---
library(tidyverse)   # data manipulation (dplyr, tidyr, ggplot2, etc.)
library(ipumsr)      # reading IPUMS microdata extracts
library(survival)    # clogit() for conditional logistic regression
library(CLogitTree)  # conditional logit forests (custom modified package)

# Load shared helper functions (feature engineering, scoring)
source("pipeline/00_helpers.R")


# ==============================================================================
# STEP 1: Load IPUMS census microdata
# ==============================================================================

# read_ipums_ddi() reads the XML codebook that describes the data structure
# read_ipums_micro() reads the actual data file using that codebook
# zap_labels() removes IPUMS label attributes (we just want plain numbers)
ddi <- read_ipums_ddi(DDI_FILE)

data <- read_ipums_micro(
  ddi,
  vars = c("YEAR",    # survey year (2016-2023)
           "SERIAL",  # household ID
           "PERNUM",  # person number within household
           "SPLOC",   # spouse's person number (0 = no spouse present)
           "SEX",     # 1 = male, 2 = female
           "AGE",     # age in years
           "EDUC",    # education level (IPUMS codes 0-11)
           "RACE",    # race (IPUMS codes)
           "HISPAN",  # Hispanic origin (0 = not Hispanic, 1+ = Hispanic subtypes)
           "METRO",   # metropolitan status
           "MET2013", # metro area code (2013 delineation)
           "YRMARR",  # year of marriage
           "MARST")   # marital status (1-2 = married spouse present)
) %>% zap_labels()


# ==============================================================================
# STEP 2: Clean and recode
# ==============================================================================

# Helper: collapse detailed IPUMS race + Hispanic codes into 5 categories
recode_race_pattern <- function(r, h) {
  case_when(
    h > 0       ~ "Hispanic",     # Hispanic of any race
    r == 1      ~ "White",        # non-Hispanic White
    r == 2      ~ "Black",        # non-Hispanic Black
    r %in% 4:6  ~ "Asian",        # non-Hispanic Asian/Pacific Islander
    TRUE        ~ "Other"         # all other combinations
  )
}

# Print sample size diagnostics
cat("Total Rows (Observations):", format(nrow(data), big.mark = ","), "\n")

n_people <- data %>% distinct(SERIAL, PERNUM) %>% nrow()
cat("Total Unique People:", format(n_people, big.mark = ","), "\n")

# Count couples: each couple appears twice (once per spouse), so divide by 2
n_couples_raw <- data %>%
  filter(SPLOC > 0, MARST %in% c(1, 2)) %>%
  nrow() / 2
cat("Total Married Co-residing Couples:", format(n_couples_raw, big.mark = ","), "\n")

# Apply sample restrictions:
# - YRMARR > 0: valid year of marriage recorded
# - SPLOC > 0: spouse is present in the household
# - MET2013 > 0: lives in an identified metropolitan area
# - marriage_duration <= 3: married within the last 3 years (recent matches)
data_clean <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0) %>%
  mutate(marriage_duration = YEAR - YRMARR) %>%
  filter(marriage_duration <= 3) %>%
  mutate(
    # Collapse race into 5 categories
    race_refined = recode_race_pattern(RACE, HISPAN),

    # Education: 5 ordered categories from IPUMS codes
    educ_cat = case_when(
      EDUC < 6      ~ "1_LessHS",      # less than high school
      EDUC == 6     ~ "2_HS_Grad",      # high school graduate
      EDUC %in% 7:9 ~ "3_SomeColl",    # some college / associate's
      EDUC == 10    ~ "4_Bachelor",     # bachelor's degree
      EDUC >= 11    ~ "5_Advanced",     # master's, professional, doctorate
      TRUE          ~ "Unknown"
    ),

    # Broader tiers (used in some interaction dummies)
    educ_tier = case_when(
      educ_cat %in% c("1_LessHS", "2_HS_Grad") ~ "Low",
      educ_cat == "3_SomeColl"                  ~ "Mid",
      educ_cat %in% c("4_Bachelor", "5_Advanced") ~ "High",
      TRUE                                       ~ NA_character_
    ),

    # Approximate years of education (used in descriptive statistics only —
    # the v2 encoding uses ordinal rank 1-5, not years; see 00_helpers.R)
    educ_yrs = case_when(
      EDUC == 0  ~ 0,   EDUC == 1  ~ 4,   EDUC == 2  ~ 8,
      EDUC == 3  ~ 9,   EDUC == 4  ~ 10,  EDUC == 5  ~ 11,
      EDUC == 6  ~ 12,  EDUC == 7  ~ 13,  EDUC == 8  ~ 14,
      EDUC == 9  ~ 15,  EDUC == 10 ~ 16,  EDUC >= 11 ~ 18,
      TRUE       ~ NA_real_
    )
  )


# ==============================================================================
# STEP 3: Optional geographic/temporal filters
# ==============================================================================

# For case studies (e.g., LA only), set METRO_CODE or METRO_CODES before sourcing.
# METRO_CODES (a vector like c(31080, 35620)) takes precedence over METRO_CODE.
if (!exists("METRO_CODE")) METRO_CODE <- NULL
if (!exists("METRO_CODES")) METRO_CODES <- NULL
if (!exists("YEAR_FILTER")) YEAR_FILTER <- NULL
if (!exists("YEAR_FILTERS")) YEAR_FILTERS <- NULL

if (!is.null(METRO_CODES)) {
  data_clean <- data_clean %>% filter(MET2013 %in% METRO_CODES)
  cat(sprintf("Metro filter: %d metros — %s rows remain\n",
              length(METRO_CODES), format(nrow(data_clean), big.mark = ",")))
} else if (!is.null(METRO_CODE)) {
  data_clean <- data_clean %>% filter(MET2013 == METRO_CODE)
  cat(sprintf("Metro filter: MET2013=%d — %s rows remain\n",
              METRO_CODE, format(nrow(data_clean), big.mark = ",")))
}
# YEAR_FILTERS (vector) takes precedence over YEAR_FILTER (scalar)
if (!is.null(YEAR_FILTERS)) {
  data_clean <- data_clean %>% filter(YEAR %in% YEAR_FILTERS)
  cat(sprintf("Year filter: %d years (%s) — %s rows remain\n",
              length(YEAR_FILTERS), paste(YEAR_FILTERS, collapse = ","),
              format(nrow(data_clean), big.mark = ",")))
} else if (!is.null(YEAR_FILTER)) {
  data_clean <- data_clean %>% filter(YEAR == YEAR_FILTER)
  cat(sprintf("Year filter: YEAR=%d — %s rows remain\n",
              YEAR_FILTER, format(nrow(data_clean), big.mark = ",")))
}

# Separate into men and women pools
men_pool   <- data_clean %>% filter(SEX == 1)
women_pool <- data_clean %>% filter(SEX == 2)


# ==============================================================================
# STEP 4: Build the sample of focal women and their true partners
# ==============================================================================

set.seed(SEED)  # for reproducibility of random sampling

# Join each woman to her actual husband's characteristics
# inner_join matches on: same year, same household (SERIAL), and the woman's
# SPLOC (spouse location) equals the man's PERNUM (person number)
couples_women <- women_pool %>%
  inner_join(
    men_pool %>% select(YEAR, SERIAL, PERNUM,
                        age_partner       = AGE,
                        educ_partner_cat  = educ_cat,
                        educ_partner_tier = educ_tier,
                        educ_yrs_partner  = educ_yrs,
                        race_partner      = race_refined),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  )

# Remove women from very small metros where there aren't enough men to form
# meaningful choice sets (default threshold: 100 men)
if (!exists("MIN_MEN_PER_METRO")) MIN_MEN_PER_METRO <- 100

metro_pool_sizes <- men_pool %>%
  filter(MET2013 > 0) %>%
  count(MET2013, YEAR, name = "n_men_in_market")

n_before <- nrow(couples_women)
couples_women <- couples_women %>%
  left_join(metro_pool_sizes, by = c("MET2013", "YEAR")) %>%
  filter(!is.na(n_men_in_market), n_men_in_market >= MIN_MEN_PER_METRO) %>%
  select(-n_men_in_market)
cat(paste0("Thin-market filter (>= ", MIN_MEN_PER_METRO, " men): ",
           n_before - nrow(couples_women), " women removed, ",
           nrow(couples_women), " eligible women remain\n"))

if (nrow(couples_women) == 0) {
  stop("No couples found after filtering. Check METRO_CODE/YEAR_FILTER settings.")
}

# Sample N_WOMEN focal women (or use all if N_WOMEN exceeds available)
# Each woman gets a unique strata_id that will identify her choice set
if (N_WOMEN >= nrow(couples_women)) {
  sample_women_base <- couples_women %>%
    mutate(strata_id = row_number())
  cat(sprintf("Using ALL %d couples (N_WOMEN=%d >= available)\n",
              nrow(sample_women_base), N_WOMEN))
} else {
  sample_women_base <- couples_women %>%
    slice_sample(n = N_WOMEN) %>%    # random sample without replacement
    mutate(strata_id = row_number())
}


# ==============================================================================
# STEP 5: Generate choice sets
# ==============================================================================

#' Build choice sets for all focal women
#'
#' For each woman, creates one row for her true partner (choice = 1) and
#' n_alts - 1 rows for randomly sampled alternative men from her marriage
#' market (choice = 0). All rows in a choice set share the same strata_id.
#'
#' @param input_women Data frame of focal women (with strata_id)
#' @param input_men Data frame of all available men
#' @param n_alts Total alternatives per choice set (including the true partner)
#' @return Data frame with n_alts rows per woman, containing both woman and
#'   partner characteristics
generate_choice_data_women <- function(input_women, input_men, n_alts) {
  n_random <- n_alts - 1  # number of random alternatives to sample

  # Build one choice set for a single woman
  build_set <- function(w) {
    # Find all men in this woman's metro area and year, excluding her husband
    market_pool <- input_men %>%
      filter(YEAR == w$YEAR, MET2013 == w$MET2013) %>%
      filter(!(SERIAL == w$SERIAL & PERNUM == w$SPLOC))

    # If there aren't enough men, skip this woman
    if (nrow(market_pool) < n_random) return(NULL)

    # Sample n_random alternative men (these are the "unchosen" options)
    market <- market_pool %>%
      slice_sample(n = n_random) %>%
      mutate(choice = 0) %>%
      rename(age_partner       = AGE,
             educ_partner_cat  = educ_cat,
             educ_partner_tier = educ_tier,
             educ_yrs_partner  = educ_yrs,
             race_partner      = race_refined)

    # The true partner (the man she actually married)
    chosen <- tibble(
      age_partner       = w$age_partner,
      educ_partner_cat  = w$educ_partner_cat,
      educ_partner_tier = w$educ_partner_tier,
      educ_yrs_partner  = w$educ_yrs_partner,
      race_partner      = w$race_partner,
      choice = 1
    )

    # Combine: true partner + alternatives, with the woman's characteristics
    bind_rows(chosen, market) %>%
      mutate(
        strata_id       = w$strata_id,
        woman_age       = w$AGE,
        woman_educ_cat  = w$educ_cat,
        woman_educ_tier = w$educ_tier,
        woman_educ_yrs  = w$educ_yrs,
        woman_race      = w$race_refined
      )
  }

  # Apply build_set to each woman and combine all results into one data frame
  # map_dfr splits the data frame into rows, applies the function, and row-binds
  purrr::map_dfr(split(input_women, seq_len(nrow(input_women))), build_set)
}

# Restrict men pool to only the metro-years we're using (saves memory)
men_pool_small <- men_pool %>%
  semi_join(
    sample_women_base %>% distinct(YEAR, MET2013),
    by = c("YEAR", "MET2013")
  )

cat("Generating choice data...\n")
raw_choice_data <- generate_choice_data_women(
  sample_women_base,
  men_pool_small,
  n_alts = N_ALTS
)


# ==============================================================================
# STEP 6: Feature engineering — two encodings
# ==============================================================================

# We create two versions ("encodings") of the predictor variables from the
# same raw data. Both describe the same thing — how similar two people are —
# but represent it differently:
#
# v1 encoding: absolute differences + direction indicators
#   Example: a couple where the man is 5 years older and more educated becomes
#     abs_age_diff=5, older_husband=1, abs_educ_diff=2, hypergamy_ed=1, same_race=1
#   Also includes one-hot dummies for every race×race and educ×educ pairing.
#   Used by: M1_v1 (the linear clogit benchmark)
#   Why: splitting into magnitude + direction is natural for a linear model that
#     can only learn "bigger gaps are worse" + "older husbands are slightly preferred"
#
# v2 encoding: signed differences + partner's absolute values
#   Example: the same couple becomes
#     age_diff=+5, alter_age=35, educ_diff=+1, alter_educ_rank=4, same_race=1
#   Education uses ordinal rank (1=LessHS ... 5=Advanced), not years of schooling.
#   Used by: M1_spline5 (spline clogit), Forest_core, Forest_full
#   Why: the signed difference lets splines capture nonlinear preferences (e.g.,
#     small gaps are fine but large gaps are penalized — a U-shape). The partner's
#     absolute level (alter_age, alter_educ_yrs) gives the forest extra information
#     to split on (e.g., "men over 40" regardless of the woman's age).
#
# The variables named with "_v1" or "_v2" throughout the code refer to these
# two encodings. For example, train_data_v2 is the training data with v2 columns.

cat("Processing variables...\n")

# --- v1 track ---
# engineer_features() creates the abs_diff/direction/one-hot encoding
engineered <- raw_choice_data %>%
  engineer_features() %>%
  rename(strata = strata_id)

# Remove rows with NAs, but ONLY in model-relevant columns.
# (The raw data has NAs in identifier columns for the "chosen" row,
# which is fine — we only care about predictor completeness.)
model_cols <- c("choice", "strata",
                "abs_age_diff", "older_husband",
                "abs_educ_diff", "hypergamy_ed", "hypogamy_ed",
                "same_race",
                grep("^(rxr_|exe_)", names(engineered), value = TRUE),
                "woman_race", "race_partner",
                "woman_educ_cat", "educ_partner_cat")

engineered <- engineered %>%
  # if_all() checks that ALL specified columns are non-NA for each row
  filter(if_all(all_of(model_cols), ~ !is.na(.))) %>%
  group_by(strata) %>%
  filter(
    sum(choice) == 1,  # the true partner must still be in the set
    n() >= 2           # need at least 1 alternative
  ) %>%
  ungroup()

# Report how many choice sets survived the filters
n_complete <- sum(table(engineered$strata) == N_ALTS)
n_partial  <- sum(table(engineered$strata) < N_ALTS & table(engineered$strata) >= 2)
cat(sprintf("Strata: %d complete (%d alts), %d partial (lost some alts to Unknown/NA filter)\n",
            n_complete, N_ALTS, n_partial))

# Numeric-only version for the forest (drops factor columns the forest can't use)
choice_data_final <- engineered %>%
  select_model_cols() %>%
  as.data.frame()

# Version with factor columns for the clogit (it needs factors for its formula)
clogit_factor_cols <- c("choice", "strata",
                        "woman_race", "race_partner",
                        "woman_educ_cat", "educ_partner_cat")
clogit_extra <- engineered[, intersect(clogit_factor_cols, names(engineered)), drop = FALSE]
choice_data_clogit <- bind_cols(choice_data_final, clogit_extra[, setdiff(names(clogit_extra),
                                                                          names(choice_data_final))])

n_predictors <- ncol(choice_data_final) - 2  # minus choice and strata
cat(paste0("Data ready (v1). Rows: ", nrow(choice_data_final),
           ", Strata: ", n_distinct(choice_data_final$strata),
           ", Predictors: ", n_predictors, "\n"))

# --- v2 track ---
# engineer_features_v2() creates the signed diff + alter level encoding
engineered_v2 <- raw_choice_data %>%
  engineer_features_v2() %>%
  rename(strata = strata_id)

v2_model_cols <- c("choice", "strata",
                   "age_diff", "alter_age", "educ_diff", "alter_educ_rank",
                   "same_race",
                   "alter_race_white", "alter_race_black",
                   "alter_race_hispanic", "alter_race_asian")
engineered_v2 <- engineered_v2 %>%
  filter(if_all(all_of(v2_model_cols), ~ !is.na(.))) %>%
  group_by(strata) %>%
  filter(sum(choice) == 1, n() >= 2) %>%
  ungroup()

choice_data_final_v2 <- engineered_v2 %>%
  select_model_cols_v2(include_race_dummies = TRUE) %>%
  as.data.frame()

n_pred_v2 <- ncol(choice_data_final_v2) - 2
cat(paste0("Data ready (v2). Rows: ", nrow(choice_data_final_v2),
           ", Strata: ", n_distinct(choice_data_final_v2$strata),
           ", Predictors: ", n_pred_v2, "\n"))


# ==============================================================================
# STEP 7: Train / test split
# ==============================================================================

# Split at the STRATA level (not the row level). All rows in a choice set
# go to either training or test — never split across.
# The split is 70/30 by default (set by TRAIN_SHARE in run_pipe.R).

if (DO_SPLIT) {
  # Use only strata that survived BOTH v1 and v2 feature engineering.
  # The two pipelines can drop different rows (v1 drops Unknown educ_cat,
  # v2 drops NA educ_yrs), so a stratum might exist in one but not the other.
  # Using the intersection ensures v1 and v2 data are perfectly aligned.
  all_strata <- intersect(
    unique(choice_data_final$strata),
    unique(choice_data_final_v2$strata)
  )
  train_strata <- sample(all_strata, size = floor(TRAIN_SHARE * length(all_strata)))
  test_strata  <- setdiff(all_strata, train_strata)

  # v1 data
  train_data        <- choice_data_final  %>% filter(strata %in% train_strata) %>% as.data.frame()
  test_data         <- choice_data_final  %>% filter(strata %in% test_strata)  %>% as.data.frame()
  train_data_clogit <- choice_data_clogit %>% filter(strata %in% train_strata) %>% as.data.frame()
  test_data_clogit  <- choice_data_clogit %>% filter(strata %in% test_strata)  %>% as.data.frame()

  # v2 data (same strata assignment — both encodings see the same women)
  train_data_v2 <- choice_data_final_v2 %>% filter(strata %in% train_strata) %>% as.data.frame()
  test_data_v2  <- choice_data_final_v2 %>% filter(strata %in% test_strata)  %>% as.data.frame()

  cat(paste0("Split: ", length(train_strata), " train / ",
             length(test_strata), " test strata\n"))
} else {
  # No split: use all data for both training and evaluation (for diagnostics only)
  train_strata <- unique(choice_data_final$strata)
  test_strata  <- train_strata
  train_data        <- choice_data_final
  test_data         <- choice_data_final
  train_data_clogit <- choice_data_clogit
  test_data_clogit  <- choice_data_clogit
  train_data_v2     <- choice_data_final_v2
  test_data_v2      <- choice_data_final_v2
  cat("No split: training and evaluating on the same data\n")
}
