# gen_fig_coef_and_splines.R
# Produces two figures for the Results section:
#   text_thesis/figures/R1_coef_m1v1.png   — coefficient plot for linear clogit
#   text_thesis/figures/R2_spline_curves.png — marginal-effect curves for spline clogit
# Run from project root: Rscript experiments/gen_fig_coef_and_splines.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(tidyverse)
library(survival)
library(splines)
library(ipumsr)
library(CLogitTree)
library(patchwork)
library(broom)

# ── pipeline parameters (must match main run) ─────────────────────────────────
DDI_FILE    <- "data/usa_00008.xml"
SEED        <- 2025L
N_WOMEN     <- 2000L
N_ALTS      <- 10L
DO_SPLIT    <- TRUE
TRAIN_SHARE <- 0.70
FIT_CLOGIT  <- TRUE
FIT_FOREST  <- FALSE
MODELS      <- "clogit"
NTREE <- MINNODESIZE <- MINBUCKET <- DEPTH_MAX <- NCORES <- TREE_TIMEOUT <- EPSILON <- ITER_MAX <- NULL

set.seed(SEED)

source("pipeline/01_data.R")
source("pipeline/02_model.R")

# ── colour palette (thesis viridis) ──────────────────────────────────────────
pal <- viridis::viridis(8, begin = 0.05, end = 0.78)
col_v1    <- pal[2]
col_spline <- pal[5]

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.3, colour = "grey88"),
      axis.ticks        = element_blank(),
      axis.title        = element_text(size = 10),
      legend.position   = "none",
      plot.title        = element_blank(),
      plot.subtitle     = element_blank()
    )
}

# =============================================================================
# FIGURE 1 — Coefficient plot for M1_v1
# =============================================================================

tidy_m1 <- broom::tidy(m1_v1, conf.int = TRUE, exponentiate = FALSE) |>
  mutate(
    term_label = case_when(
      term == "abs_age_diff"  ~ "Age gap (absolute, years)",
      term == "older_husband" ~ "Man is older",
      term == "abs_educ_diff" ~ "Education gap (absolute, levels)",
      term == "hypergamy_ed"  ~ "Man more educated (hypergamy)",
      term == "hypogamy_ed"   ~ "Man less educated (hypogamy)",
      term == "same_race"     ~ "Same race / ethnicity",
      TRUE ~ term
    ),
    term_label = factor(term_label, levels = rev(c(
      "Same race / ethnicity",
      "Man more educated (hypergamy)",
      "Man less educated (hypogamy)",
      "Education gap (absolute, levels)",
      "Man is older",
      "Age gap (absolute, years)"
    ))),
    significant = p.value < 0.05
  )

p_coef <- ggplot(tidy_m1, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.25, colour = col_v1, linewidth = 0.7,
                orientation = "y") +
  geom_point(aes(colour = significant), size = 3) +
  scale_colour_manual(values = c("TRUE" = col_v1, "FALSE" = "grey60")) +
  labs(
    x = "Coefficient",
    y = NULL
  ) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 10))

ggsave("text_thesis/figures/R1_coef_m1v1.png",
       plot = p_coef, width = 7, height = 4, dpi = 300, bg = "white")
cat("Saved: text_thesis/figures/R1_coef_m1v1.png\n")

# =============================================================================
# FIGURE 2 — Spline marginal-effect curves (age gap + education gap)
# =============================================================================
# Strategy: vary one predictor across its observed range, hold other predictors
# at their median/reference values, predict V_ij (linear predictor), and plot.

# clip age grid to 99th percentile to avoid wild spline extrapolation at extremes
age_p99_lo <- quantile(train_v2_aug$age_diff, 0.01, na.rm = TRUE)  # -39
age_p99_hi <- quantile(train_v2_aug$age_diff, 0.99, na.rm = TRUE)  #  44

age_grid <- seq(age_p99_lo, age_p99_hi, length.out = 200)
educ_med <- median(train_v2_aug$educ_diff, na.rm = TRUE)

newdata_age <- tibble(
  age_diff  = age_grid,
  educ_diff = educ_med,
  same_race = 0L
) |>
  add_spline_cols(basis_age = basis_age5, basis_educ = basis_educ3)

spline_cols <- c(paste0("age_ns5_",  1:5),
                 paste0("educ_ns3_", 1:3),
                 "same_race")
coefs <- coef(m1_spline5)

lp_age <- as.matrix(newdata_age[, spline_cols]) %*% coefs
lp_age <- lp_age - mean(lp_age)

df_age <- tibble(x = age_grid, lp = as.vector(lp_age))

# --- Education gap curve ---
educ_grid <- seq(
  min(train_v2_aug$educ_diff, na.rm = TRUE),
  max(train_v2_aug$educ_diff, na.rm = TRUE),
  length.out = 100
)
age_med <- median(train_v2_aug$age_diff, na.rm = TRUE)

newdata_educ <- tibble(
  age_diff  = age_med,
  educ_diff = educ_grid,
  same_race = 0L
) |>
  add_spline_cols(basis_age = basis_age5, basis_educ = basis_educ3)

lp_educ <- as.matrix(newdata_educ[, spline_cols]) %*% coefs
lp_educ <- lp_educ - mean(lp_educ)

df_educ <- tibble(x = educ_grid, lp = lp_educ)

# --- Plot ---
p_age <- ggplot(df_age, aes(x = x, y = lp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_line(colour = col_spline, linewidth = 1.1) +
  scale_x_continuous(
    breaks = c(-30, -20, -10, 0, 10, 20, 30, 40),
    labels = c("-30", "-20", "-10", "0", "10", "20", "30", "40")
  ) +
  labs(
    x = "Age gap (husband minus wife, years; negative = wife older)",
    y = "Deviation from mean (log-probability scale)"
  ) +
  theme_pub()

p_educ <- ggplot(df_educ, aes(x = x, y = lp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_line(colour = col_spline, linewidth = 1.1) +
  scale_x_continuous(
    breaks = -4:4,
    labels = c("-4", "-3", "-2", "-1", "0\n(same)", "+1", "+2", "+3", "+4")
  ) +
  labs(
    x = "Education gap (husband minus wife, levels; negative = wife more educated)",
    y = "Deviation from mean (log-probability scale)"
  ) +
  theme_pub()

p_splines <- p_age | p_educ

ggsave("text_thesis/figures/R2_spline_curves.png",
       plot = p_splines, width = 10, height = 4.5, dpi = 300, bg = "white")
cat("Saved: text_thesis/figures/R2_spline_curves.png\n")
