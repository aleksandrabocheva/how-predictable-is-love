# gen_fig_education_homogamy.R
# Produces text_thesis/figures/03_education_homogamy.png
# Educational homogamy matrix: rows = wife's education, columns = husband's education,
# cells = % of women in that row who married a man in the column education level.
# Run from the project root:  Rscript experiments/gen_fig_education_homogamy.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ipumsr)
})

DDI_FILE          <- "data/usa_00008.xml"
MIN_MEN_PER_METRO <- 100
OUT_FILE          <- "text_thesis/figures/03_education_homogamy.png"

pal <- viridis::viridis(8, begin = 0.05, end = 0.78)

recode_educ <- function(educ) {
  case_when(
    educ < 6      ~ "1_LessHS",
    educ == 6     ~ "2_HS_Grad",
    educ %in% 7:9 ~ "3_SomeColl",
    educ == 10    ~ "4_Bachelor",
    educ >= 11    ~ "5_Advanced",
    TRUE          ~ NA_character_
  )
}

educ_levels <- c("1_LessHS", "2_HS_Grad", "3_SomeColl", "4_Bachelor", "5_Advanced")
educ_labels <- c("< HS", "HS Grad", "Some College", "Bachelor's", "Advanced")

cat("Loading IPUMS data...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(
  ddi,
  vars = c("YEAR", "SERIAL", "PERNUM", "SPLOC",
           "SEX", "EDUC", "MET2013", "YRMARR")
) %>% zap_labels()

cat("Filtering and joining couples...\n")
data_clean <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0,
         (YEAR - YRMARR) <= 3) %>%
  mutate(educ_cat = recode_educ(EDUC)) %>%
  filter(!is.na(educ_cat))

men_pool   <- data_clean %>% filter(SEX == 1)
women_pool <- data_clean %>% filter(SEX == 2)

metro_sizes <- men_pool %>% count(MET2013, YEAR, name = "n_men")

couples <- women_pool %>%
  inner_join(
    men_pool %>% select(YEAR, SERIAL, PERNUM, man_educ_cat = educ_cat),
    by = c("YEAR", "SERIAL", "SPLOC" = "PERNUM")
  ) %>%
  left_join(metro_sizes, by = c("MET2013", "YEAR")) %>%
  filter(!is.na(n_men), n_men >= MIN_MEN_PER_METRO) %>%
  rename(woman_educ_cat = educ_cat)

cat(sprintf("Couples for matrix: %s\n", format(nrow(couples), big.mark = ",")))

mat <- couples %>%
  mutate(
    woman_educ_cat = factor(woman_educ_cat, levels = educ_levels, labels = educ_labels),
    man_educ_cat   = factor(man_educ_cat,   levels = educ_levels, labels = educ_labels)
  ) %>%
  count(woman_educ_cat, man_educ_cat) %>%
  group_by(woman_educ_cat) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p <- mat %>%
  ggplot(aes(x = man_educ_cat, y = fct_rev(woman_educ_cat), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            size = 3.2, colour = "grey15") +
  scale_fill_gradient(
    low  = "white",
    high = pal[2],
    name = "% of row",
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(position = "bottom") +
  labs(
    x = "Husband's Education",
    y = "Wife's Education"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave(OUT_FILE, plot = p, width = 6.5, height = 4.5, dpi = 300, bg = "white")
cat("Saved:", OUT_FILE, "\n")
