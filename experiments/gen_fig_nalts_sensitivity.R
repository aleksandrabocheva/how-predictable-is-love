# gen_fig_nalts_sensitivity.R
# Generates figures/A1_nalts_sensitivity.png — N_ALTS robustness check.
# Run from the project root:  Rscript experiments/gen_fig_nalts_sensitivity.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(viridis)

OUT_FILE <- "text_thesis/figures/A1_nalts_sensitivity.png"

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

model_labels <- c(
  "M1_v1"      = "Linear Clogit",
  "M1_spline"  = "Spline Clogit",
  "Forest_core" = "Forest"
)

d <- read_csv("results/08_nalts_scaling.csv", show_col_types = FALSE) |>
  mutate(
    model       = factor(model, levels = names(model_labels), labels = model_labels),
    top1_pct    = top1_nalt * 100,
    random_top1 = 100 / N_alts
  )

# Panel A: Top-1 accuracy — decreases mechanically as task gets harder
p_top1 <- d |>
  ggplot(aes(x = N_alts, y = top1_pct, colour = model, shape = model)) +
  # random baseline (varies with N_alts)
  geom_line(aes(x = N_alts, y = random_top1),
            data = d |> distinct(N_alts, random_top1),
            colour = "grey60", linetype = "dashed", linewidth = 0.6,
            inherit.aes = FALSE) +
  annotate("text", x = 47, y = 3.5, label = "Random\nbaseline",
           colour = "grey50", size = 3, hjust = 1) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = pal_cb[c(1, 3, 5)], name = NULL) +
  scale_shape_manual(values  = c(16, 17, 15),       name = NULL) +
  scale_x_continuous(breaks = c(10, 20, 50)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 60)) +
  labs(
    x = "Choice set size (N alternatives)",
    y = "Top-1 accuracy (%)"
  ) +
  theme(legend.position = "bottom")

# Panel B: Metro percentile rank — stable across N_ALTS (key robustness finding)
p_metro <- d |>
  ggplot(aes(x = N_alts, y = metro_pctile, colour = model, shape = model)) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey60",
             linewidth = 0.6) +
  annotate("text", x = 47, y = 52, label = "Random\nbaseline (50%)",
           colour = "grey50", size = 3, hjust = 1) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = pal_cb[c(1, 3, 5)], name = NULL) +
  scale_shape_manual(values  = c(16, 17, 15),       name = NULL) +
  scale_x_continuous(breaks = c(10, 20, 50)) +
  scale_y_continuous(limits = c(0, 55),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x = "Choice set size (N alternatives)",
    y = "Percentile rank of true husband (lower = better)"
  ) +
  theme(legend.position = "bottom")

library(patchwork)
p <- (p_top1 | p_metro) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(OUT_FILE, plot = p, width = 9, height = 4, dpi = 300)
cat("Saved:", OUT_FILE, "\n")
