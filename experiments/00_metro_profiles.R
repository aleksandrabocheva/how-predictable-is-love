#!/usr/bin/env Rscript
# ==============================================================================
# 00_metro_profiles.R — Metro area composition profiles
#
# WHAT IT DOES:
#   Loads IPUMS ACS data, applies the same filters as pipeline/01_data.R
#   (recently married couples, identified metro areas), and profiles the 20
#   largest US metros by race composition, diversity, couple counts, and
#   education distribution.
#
# STORY:
#   Before running predictive models in specific metros, we need to understand
#   what markets look like: how many couples are available, how racially diverse
#   each metro is, and which metros have enough data to produce stable estimates.
#   This script produces the reference table used to select metros for targeted
#   experiments (e.g., LA = most diverse large metro) and to interpret why
#   same_race has different importance in different markets.
#
# HOW TO RUN:
#   cd partner-choice
#   Rscript experiments/00_metro_profiles.R
#
# EXPECTED RUNTIME: ~3-5 min (data loading dominates)
#
# OUTPUT:
#   results/00_metro_profiles.csv — compact reference table (20 largest metros)
#   (Also prints full distributions to console)
#
# USED BY:
#   - experiments/01_forest_tuning.R  (to confirm LA is the right test metro)
#   - experiments/04_spec_comparison.R (to select multi-metro test set)
#   - paper/tables: Table 1 metro descriptives
# ==============================================================================

library(tidyverse)
library(ipumsr)

cat("=== 00_metro_profiles.R — Metro area composition profiles ===\n\n")

# --- Config ---
DDI_FILE <- "data/usa_00008.xml"
N_TOP    <- 20   # number of metros to include in output CSV

# --- Load raw IPUMS data ---
cat("Loading IPUMS data from", DDI_FILE, "...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(
  ddi,
  vars = c("YEAR", "SERIAL", "PERNUM", "SPLOC", "SEX", "AGE",
           "EDUC", "RACE", "HISPAN", "METRO", "MET2013", "YRMARR", "MARST")
) %>% zap_labels()

cat("Loaded:", format(nrow(data), big.mark = ","), "rows\n")

# ==============================================================================
# RECODE — same logic as pipeline/01_data.R
# ==============================================================================

recode_race <- function(r, h) {
  case_when(
    h > 0      ~ "Hispanic",
    r == 1     ~ "White",
    r == 2     ~ "Black",
    r %in% 4:6 ~ "Asian",
    TRUE       ~ "Other"
  )
}

# Education: ordinal 1-5 rank (used throughout the modelling pipeline)
# Note: educ_yrs is also retained for descriptive context (approximate only —
# the ordinal rank is the canonical representation in models).
educ_levels <- c("1_LessHS", "2_HS_Grad", "3_SomeColl", "4_Bachelor", "5_Advanced")

data_clean <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0) %>%
  mutate(
    marriage_duration = YEAR - YRMARR,
    race = recode_race(RACE, HISPAN),
    # Ordinal education rank (1-5) — canonical representation
    educ_cat = case_when(
      EDUC < 6      ~ "1_LessHS",
      EDUC == 6     ~ "2_HS_Grad",
      EDUC %in% 7:9 ~ "3_SomeColl",
      EDUC == 10    ~ "4_Bachelor",
      EDUC >= 11    ~ "5_Advanced",
      TRUE          ~ "Unknown"
    ),
    educ_rank = as.integer(factor(educ_cat, levels = educ_levels)),
    # Approximate years of schooling — for descriptive stats only, NOT for models
    educ_yrs_approx = case_when(
      EDUC == 0  ~ 0,
      EDUC == 1  ~ 4,
      EDUC == 2  ~ 8,
      EDUC == 3  ~ 9,
      EDUC == 4  ~ 10,
      EDUC == 5  ~ 11,
      EDUC == 6  ~ 12,
      EDUC == 7  ~ 13,
      EDUC == 8  ~ 14,
      EDUC == 9  ~ 15,
      EDUC == 10 ~ 16,
      EDUC == 11 ~ 18,
      TRUE       ~ NA_real_
    )
  ) %>%
  filter(marriage_duration <= 3)

cat("After filter (recently married, metro identified):",
    format(nrow(data_clean), big.mark = ","), "rows\n")

men   <- data_clean %>% filter(SEX == 1)
women <- data_clean %>% filter(SEX == 2)

cat(sprintf("Men: %s   Women: %s\n",
            format(nrow(men), big.mark = ","),
            format(nrow(women), big.mark = ",")))

# ==============================================================================
# BUILD COUPLES
# ==============================================================================

cat("\nBuilding couples (household join)...\n")
couples <- women %>%
  inner_join(
    men %>% select(YEAR, SERIAL, PERNUM,
                   age_m      = AGE,
                   race_m     = race,
                   educ_cat_m = educ_cat,
                   educ_rank_m = educ_rank),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  )
cat("Couples:", format(nrow(couples), big.mark = ","), "\n")

# ==============================================================================
# GET METRO LABELS FROM DDI
# ==============================================================================

met_info <- ddi$var_info %>%
  filter(var_name == "MET2013") %>%
  pull(val_labels) %>%
  .[[1]]

metro_labels <- tibble(
  MET2013    = as.integer(met_info$val),
  metro_name = met_info$lbl
) %>% filter(MET2013 > 0)

# ==============================================================================
# RACE COMPOSITION (men and women separately)
# ==============================================================================

cat("\nComputing race composition...\n")

race_wide <- function(df, sex_label) {
  df %>%
    count(MET2013, race) %>%
    group_by(MET2013) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = race, values_from = pct,
                values_fill = 0, names_prefix = paste0(sex_label, "_pct_"))
}

race_men   <- race_wide(men,   "m")
race_women <- race_wide(women, "w")

# ==============================================================================
# AGE AND EDUCATION STATS (men and women separately)
# ==============================================================================

age_educ_stats <- function(df, sex_label) {
  df %>%
    group_by(MET2013) %>%
    summarise(
      "{sex_label}_n"              := n(),
      "{sex_label}_age_mean"       := round(mean(AGE, na.rm = TRUE), 1),
      "{sex_label}_age_sd"         := round(sd(AGE, na.rm = TRUE), 1),
      # educ_rank_mean: mean ordinal rank (1-5) — primary representation
      "{sex_label}_educ_rank_mean" := round(mean(educ_rank, na.rm = TRUE), 2),
      # educ_yrs_mean: approximate years (descriptive only, NOT used in models)
      "{sex_label}_educ_yrs_approx_mean" := round(mean(educ_yrs_approx, na.rm = TRUE), 1),
      .groups = "drop"
    )
}

stats_men   <- age_educ_stats(men,   "m")
stats_women <- age_educ_stats(women, "w")

# ==============================================================================
# COUPLE-LEVEL STATS
# ==============================================================================

cat("Computing couple-level statistics...\n")

couple_stats <- couples %>%
  group_by(MET2013) %>%
  summarise(
    n_couples      = n(),
    same_race_pct  = round(mean(race == race_m) * 100, 1),
    mean_age_diff  = round(mean(age_m - AGE), 2),        # positive = husband older
    mean_educ_diff = round(mean(educ_rank_m - educ_rank, na.rm = TRUE), 3),
    # Fraction at each education pairing type
    pct_hypergamy  = round(100 * mean(educ_rank_m > educ_rank, na.rm = TRUE), 1),
    pct_homogamy   = round(100 * mean(educ_rank_m == educ_rank, na.rm = TRUE), 1),
    pct_hypogamy   = round(100 * mean(educ_rank_m < educ_rank, na.rm = TRUE), 1),
    .groups = "drop"
  )

# ==============================================================================
# DIVERSITY INDEX
#
# Herfindahl-Hirschman Index (HHI) of race proportions: sum(p^2)
#   HHI = 1  → perfectly homogeneous (one race = 100%)
#   HHI = 0  → perfectly diverse (infinitely many equal groups)
# We report 1 - HHI as the "diversity index" so higher = more diverse.
# Primary use: select diverse metros for race-sensitivity analyses.
# ==============================================================================

cat("Computing diversity indices (HHI)...\n")

diversity <- data_clean %>%
  count(MET2013, race) %>%
  group_by(MET2013) %>%
  mutate(share = n / sum(n)) %>%
  summarise(
    hhi             = round(sum(share^2), 4),
    diversity_index = round(1 - sum(share^2), 4),  # 0 = homogeneous, 1 = fully diverse
    dominant_race   = race[which.max(share)],
    dominant_pct    = round(max(share) * 100, 1),
    n_races_gt5pct  = sum(share > 0.05),
    .groups = "drop"
  )

# ==============================================================================
# MARKET SIZES (metro × year) — used for MIN_MEN_PER_METRO filter
# ==============================================================================

cat("Computing per-year market sizes...\n")

market_sizes <- men %>%
  count(MET2013, YEAR, name = "n_men") %>%
  group_by(MET2013) %>%
  summarise(
    n_years_in_data    = n(),
    min_men_per_year   = min(n_men),
    max_men_per_year   = max(n_men),
    mean_men_per_year  = round(mean(n_men)),
    total_men          = sum(n_men),
    .groups = "drop"
  )

# ==============================================================================
# ASSEMBLE FULL SUMMARY
# ==============================================================================

cat("Assembling metro summary table...\n")

metro_summary <- metro_labels %>%
  inner_join(couple_stats,  by = "MET2013") %>%
  left_join(stats_men,      by = "MET2013") %>%
  left_join(stats_women,    by = "MET2013") %>%
  left_join(race_men,       by = "MET2013") %>%
  left_join(race_women,     by = "MET2013") %>%
  left_join(diversity,      by = "MET2013") %>%
  left_join(market_sizes,   by = "MET2013") %>%
  arrange(desc(n_couples))

# ==============================================================================
# CONSOLE OUTPUT
# ==============================================================================

cat("\n=== TOP", N_TOP, "METROS BY N COUPLES ===\n")
metro_summary %>%
  select(metro_name, n_couples, m_n, w_n,
         diversity_index, dominant_race, dominant_pct,
         same_race_pct, n_races_gt5pct) %>%
  head(N_TOP) %>%
  print(n = N_TOP, width = 200)

cat("\n=== METROS WITH DIVERSITY INDEX > 0.5 (most racially diverse) ===\n")
metro_summary %>%
  filter(diversity_index > 0.5) %>%
  select(metro_name, n_couples, diversity_index, dominant_race,
         dominant_pct, same_race_pct, n_races_gt5pct) %>%
  arrange(desc(diversity_index)) %>%
  print(n = 50, width = 200)

cat("\n=== METROS >90% ONE RACE (most homogeneous — limited race variation) ===\n")
metro_summary %>%
  filter(dominant_pct > 90) %>%
  select(metro_name, n_couples, dominant_race, dominant_pct, diversity_index) %>%
  arrange(desc(n_couples)) %>%
  print(n = 30, width = 200)

cat("\n=== MARKET SIZE FILTER IMPACT ===\n")
cat("A metro is excluded from an experiment if ANY year has < threshold men.\n\n")
for (threshold in c(50, 100, 200, 500)) {
  excluded <- market_sizes %>% filter(min_men_per_year < threshold)
  retained <- market_sizes %>% filter(min_men_per_year >= threshold)
  n_excluded_couples <- couples %>%
    semi_join(excluded, by = "MET2013") %>% nrow()
  cat(sprintf("  Threshold >= %3d men/year: %2d metros excluded (%s couples lost), %2d retained\n",
              threshold,
              nrow(excluded),
              format(n_excluded_couples, big.mark = ","),
              nrow(retained)))
}

cat("\n=== DATA BY YEAR ===\n")
data_clean %>%
  count(YEAR, SEX) %>%
  pivot_wider(names_from = SEX, values_from = n, names_prefix = "sex_") %>%
  rename(men = sex_1, women = sex_2) %>%
  print()

cat("\n=== EDUCATION DISTRIBUTION IN TOP 5 METROS (women) ===\n")
top5 <- metro_summary %>% head(5) %>% pull(MET2013)
women %>%
  filter(MET2013 %in% top5) %>%
  left_join(metro_labels, by = "MET2013") %>%
  filter(educ_cat != "Unknown") %>%
  count(metro_name, educ_cat) %>%
  group_by(metro_name) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = educ_cat, values_from = pct, values_fill = 0) %>%
  print(width = 200)

# ==============================================================================
# SAVE
# ==============================================================================

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Compact output: top N_TOP metros, key columns only
metro_compact <- metro_summary %>%
  head(N_TOP) %>%
  select(
    MET2013, metro_name, n_couples, m_n, w_n,
    hhi, diversity_index, dominant_race, dominant_pct, n_races_gt5pct,
    same_race_pct, pct_hypergamy, pct_homogamy, pct_hypogamy,
    mean_age_diff, mean_educ_diff,
    m_age_mean, w_age_mean,
    m_educ_rank_mean, w_educ_rank_mean,
    m_educ_yrs_approx_mean, w_educ_yrs_approx_mean,
    min_men_per_year, mean_men_per_year, n_years_in_data
  )

write_csv(metro_compact, "results/00_metro_profiles.csv")
cat(sprintf("\nSaved results/00_metro_profiles.csv (%d metros, %d columns)\n",
            nrow(metro_compact), ncol(metro_compact)))

# Also save the full summary as RDS for interactive exploration
saveRDS(metro_summary, "results/00_metro_profiles_full.rds")
cat("Saved results/00_metro_profiles_full.rds\n")

cat("\n=== QUICK REFERENCE: LA vs NYC ===\n")
metro_summary %>%
  filter(MET2013 %in% c(31080, 35620)) %>%  # LA=31080, NYC=35620
  select(metro_name, n_couples, diversity_index, dominant_race,
         dominant_pct, same_race_pct, m_n, min_men_per_year) %>%
  print(width = 200)

cat("\nDone.\n")
