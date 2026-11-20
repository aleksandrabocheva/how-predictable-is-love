# gen_fig_metro_diversity.R
# Produces text_thesis/figures/R4_metro_diversity.png
#
# Scatter plot: racial diversity index vs same-race marriage rate
# across the 20 largest US metros. Used to explain why predictability
# varies across cities in the multi-metro robustness analysis.
#
# Run from project root: Rscript experiments/gen_fig_metro_diversity.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

pal <- viridis::viridis(8, begin = 0.05, end = 0.78)

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.3, colour = "grey88"),
      axis.ticks        = element_blank(),
      axis.title        = element_text(size = 10),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      legend.position   = "none"
    )
}

# ── load data ─────────────────────────────────────────────────────────────────
metros <- read_csv("results/00_metro_profiles.csv", show_col_types = FALSE) |>
  mutate(
    # short label: first city name only
    city_label = gsub("[-,].*", "", metro_name),
    city_label = trimws(city_label),
    # flag the four focal metros used in the multi-metro analysis
    focal = MET2013 %in% c(31080, 12060, 41860, 14460),
    focal_name = case_when(
      MET2013 == 31080 ~ "Los Angeles",
      MET2013 == 12060 ~ "Atlanta",
      MET2013 == 41860 ~ "San Francisco",
      MET2013 == 14460 ~ "Boston",
      TRUE             ~ NA_character_
    )
  )

col_focal <- pal[6]   # teal/green for the four highlighted metros
col_other <- "grey60"

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(metros, aes(x = diversity_index, y = same_race_pct)) +

  # smooth trend line across all metros
  geom_smooth(method = "lm", se = TRUE, colour = "grey80",
              fill = "grey90", linewidth = 0.7, linetype = "solid") +

  # all background metros (grey)
  geom_point(data = filter(metros, !focal),
             colour = col_other, size = 3, alpha = 0.75) +

  # four focal metros (coloured, larger)
  geom_point(data = filter(metros, focal),
             colour = col_focal, size = 4.5) +

  # labels for ALL cities
  geom_text_repel(
    aes(label = city_label),
    data          = filter(metros, !focal),
    colour        = "grey50",
    size          = 2.8,
    max.overlaps  = 20,
    seed          = 42,
    segment.color = "grey75",
    segment.size  = 0.3,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +

  # labels for focal cities (darker, bold-weight via fontface)
  geom_text_repel(
    aes(label = city_label),
    data          = filter(metros, focal),
    colour        = col_focal,
    size          = 3.3,
    fontface      = "bold",
    max.overlaps  = 20,
    seed          = 42,
    segment.color = col_focal,
    segment.size  = 0.4,
    box.padding   = 0.4,
    point.padding = 0.3,
    nudge_y       = 0.5
  ) +

  scale_x_continuous(
    limits = c(0.3, 0.75),
    labels = function(x) round(x, 2),
    name   = "Racial diversity index (0 = homogeneous, 1 = fully diverse)"
  ) +
  scale_y_continuous(
    limits = c(65, 92),
    labels = function(y) paste0(y, "%"),
    name   = "Same-race marriage rate (%)"
  ) +
  labs() +
  theme_pub()

ggsave("text_thesis/figures/R4_metro_diversity.png",
       plot = p, width = 8, height = 5.5, dpi = 300, bg = "white")
cat("Saved: text_thesis/figures/R4_metro_diversity.png\n")
