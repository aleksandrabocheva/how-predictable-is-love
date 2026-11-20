# gen_fig_race_homogamy.R
# Generates text_thesis/figures/04_race_homogamy.png
# Racial endogamy matrix: rows = wife's race, columns = husband's race,
# cells = % of women in that row who married a man in the column group.
# Run from the project root:  Rscript experiments/gen_fig_race_homogamy.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(ipumsr)

DDI_FILE          <- "data/usa_00008.xml"
MIN_MEN_PER_METRO <- 100
OUT_FILE          <- "text_thesis/figures/04_race_homogamy.png"

recode_race <- function(r, h) {
  case_when(
    h > 0        ~ "Hispanic",
    r == 1       ~ "White",
    r == 2       ~ "Black",
    r %in% 4:6   ~ "Asian",
    TRUE         ~ "Other"
  )
}

cat("Loading IPUMS data...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(
  ddi,
  vars = c("YEAR", "SERIAL", "PERNUM", "SPLOC",
           "SEX", "AGE", "RACE", "HISPAN", "MET2013", "YRMARR")
) %>% zap_labels()

cat("Filtering and joining couples...\n")
data_clean <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0,
         (YEAR - YRMARR) <= 3) %>%
  mutate(race = recode_race(RACE, HISPAN))

men_pool   <- data_clean %>% filter(SEX == 1)
women_pool <- data_clean %>% filter(SEX == 2)

metro_sizes <- men_pool %>% count(MET2013, YEAR, name = "n_men")

couples <- women_pool %>%
  inner_join(
    men_pool %>% select(YEAR, SERIAL, PERNUM, man_race = race),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  ) %>%
  left_join(metro_sizes, by = c("MET2013", "YEAR")) %>%
  filter(!is.na(n_men), n_men >= MIN_MEN_PER_METRO) %>%
  rename(wife_race = race)

cat(sprintf("Couples for matrix: %s\n", format(nrow(couples), big.mark = ",")))

# Ordered factor levels for consistent row/column ordering
race_levels <- c("White", "Black", "Hispanic", "Asian", "Other")

mat <- couples %>%
  mutate(
    wife_race = factor(wife_race, levels = race_levels),
    man_race  = factor(man_race,  levels = race_levels)
  ) %>%
  count(wife_race, man_race) %>%
  group_by(wife_race) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p <- mat %>%
  ggplot(aes(x = man_race, y = fct_rev(wife_race), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            size = 3.2, colour = "grey15") +
  scale_fill_gradient(
    low  = "#fef0d9",
    high = "#b30000",
    name = "% of row",
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(position = "bottom") +
  labs(
    x = "Husband's Race",
    y = "Wife's Race"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "right"
  )

ggsave(OUT_FILE, plot = p, width = 6.5, height = 4.5, dpi = 300)
cat("Saved:", OUT_FILE, "\n")
