# gen_fig_age_dist.R
# Generates figures/01_age_distribution.png — boxplot of age gap by husband age group.
# Run from the project root:  Rscript experiments/gen_fig_age_dist.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(ipumsr)
library(viridis)

DDI_FILE          <- "data/usa_00008.xml"
MIN_MEN_PER_METRO <- 100
OUT_FILE          <- "text_thesis/figures/01_age_distribution.png"

pal_cb <- viridis::viridis(8, begin = 0.05, end = 0.78)

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.3, colour = "grey88"),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      strip.text        = element_text(face = "bold", size = rel(0.9)),
      axis.text         = element_text(colour = "grey20"),
      legend.key        = element_blank(),
      legend.background = element_rect(fill = "white", colour = "grey80"),
      plot.margin       = margin(8, 8, 8, 8)
    )
}
theme_set(theme_pub())

recode_race <- function(r, h) {
  case_when(h > 0 ~ "Hispanic", r == 1 ~ "White", r == 2 ~ "Black",
            r %in% 4:6 ~ "Asian", TRUE ~ "Other")
}

cat("Loading IPUMS data...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(ddi,
  vars = c("YEAR","SERIAL","PERNUM","SPLOC","SEX","AGE","EDUC",
           "RACE","HISPAN","MET2013","YRMARR","MARST")
) %>% zap_labels()

cat("Filtering and joining couples...\n")
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
    )
  )

men_pool   <- data_clean %>% filter(SEX == 1)
women_pool <- data_clean %>% filter(SEX == 2)

metro_pool_sizes <- men_pool %>%
  count(MET2013, YEAR, name = "n_men_in_market")

couples_clean <- women_pool %>%
  inner_join(
    men_pool %>% select(YEAR, SERIAL, PERNUM,
                        man_age      = AGE,
                        man_educ_cat = educ_cat,
                        man_race     = race_refined),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  ) %>%
  left_join(metro_pool_sizes, by = c("MET2013", "YEAR")) %>%
  filter(!is.na(n_men_in_market), n_men_in_market >= MIN_MEN_PER_METRO) %>%
  select(-n_men_in_market) %>%
  mutate(
    woman_age    = AGE,
    woman_race   = race_refined,
    woman_educ_cat = educ_cat,
    age_gap      = man_age - woman_age
  )

cat(sprintf("Couples available: %s\n", format(nrow(couples_clean), big.mark = ",")))

p <- couples_clean %>%
  mutate(
    age_group = cut(
      man_age,
      breaks = c(18, 23, 28, 33, 38, 43, 48, 53, 58, 65),
      right  = FALSE,
      labels = c("18-22", "23-27", "28-32", "33-37",
                 "38-42", "43-47", "48-52", "53-57", "58-64")
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  ggplot(aes(x = age_group, y = age_gap)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_boxplot(
    fill          = pal_cb[1], alpha        = 0.55,
    colour        = "grey30",  linewidth    = 0.4,
    outlier.size  = 0.5,       outlier.alpha = 0.2
  ) +
  coord_cartesian(ylim = c(-20, 30)) +
  labs(
    x = "Husband's Age Group (years)",
    y = "Age Gap (husband - wife, years)"
  )

ggsave(OUT_FILE, plot = p, width = 7, height = 4, dpi = 300)
cat("Saved:", OUT_FILE, "\n")
