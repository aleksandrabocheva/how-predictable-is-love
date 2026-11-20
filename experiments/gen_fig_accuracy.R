# gen_fig_accuracy.R
# Produces text_thesis/figures/R3_accuracy.png — graphical model comparison
# (replaces the plain accuracy table with a Cleveland dot plot).
# Run from project root: Rscript experiments/gen_fig_accuracy.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(patchwork)

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
      strip.background  = element_rect(fill = "grey96", colour = NA),
      strip.text        = element_text(size = 10),
      legend.position   = "bottom",
      legend.title      = element_blank()
    )
}

# ── load results ──────────────────────────────────────────────────────────────
metrics <- read_csv("results/w2000_a10_t50_mns550_s2025/metrics.csv",
                    show_col_types = FALSE) |>
  mutate(
    Model = case_when(
      model == "M1_v1"       ~ "Linear Clogit",
      model == "M1_spline5"  ~ "Spline Clogit",
      model == "Forest_core" ~ "Forest (core)",
      model == "Forest_full" ~ "Forest (full)"
    ),
    Model = factor(Model, levels = c("Linear Clogit", "Spline Clogit",
                                     "Forest (core)", "Forest (full)"))
  )

model_colours <- setNames(pal[c(2, 4, 6, 7)],
                          c("Linear Clogit", "Spline Clogit",
                            "Forest (core)", "Forest (full)"))

# =============================================================================
# PANEL A — N-alternative accuracy (top-1 / top-3 / top-5)
# =============================================================================
acc_long <- metrics |>
  select(Model, top1, top3, top5) |>
  pivot_longer(cols = c(top1, top3, top5),
               names_to  = "metric",
               values_to = "accuracy") |>
  mutate(
    metric = factor(metric,
                    levels = c("top1", "top3", "top5"),
                    labels = c("Top-1", "Top-3", "Top-5")),
    accuracy_pct = accuracy * 100
  )

# random baselines for 10-alternative choice set
random_baselines <- tibble(
  metric   = factor(c("Top-1", "Top-3", "Top-5"),
                    levels = c("Top-1", "Top-3", "Top-5")),
  baseline = c(10, 30, 50)
)

p_acc <- ggplot(acc_long, aes(x = accuracy_pct, y = Model, colour = Model)) +
  geom_vline(data = random_baselines,
             aes(xintercept = baseline),
             linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = model_colours) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  labs(
    x = "Accuracy (%)",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "none")

# =============================================================================
# PANEL B — Full metropolitan market: percentile rank
# =============================================================================
p_metro <- ggplot(metrics,
                  aes(x = metro_pctile, y = Model, colour = Model)) +
  geom_vline(xintercept = 50, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_point(size = 4) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(limits = c(0, 55),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x = "Mean percentile rank of true husband (lower = better)",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"))

# =============================================================================
# Combine and save
# =============================================================================
p_combined <- p_acc / p_metro +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("text_thesis/figures/R3_accuracy.png",
       plot = p_combined, width = 9, height = 7, dpi = 300, bg = "white")
cat("Saved: text_thesis/figures/R3_accuracy.png\n")
