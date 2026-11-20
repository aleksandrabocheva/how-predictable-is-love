# gen_pipeline_diagram.R
# Produces figures/pipeline_diagram.png — two-panel flowchart of the full
# data processing and modelling pipeline.
# Run from project root: Rscript experiments/gen_pipeline_diagram.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(ggplot2)
library(patchwork)

OUT_FILE <- "text_thesis/figures/pipeline_diagram.png"

# ── colour palette (matches thesis viridis palette) ─────────────────────────
col_data   <- "#3B0F70"   # dark purple  — raw data / filtering steps
col_feat   <- "#8C2981"   # mid purple   — feature engineering
col_v1     <- "#DE4968"   # pink-red     — v1 encoding branch
col_v2     <- "#FE9F6D"   # orange       — v2 encoding branch
col_model  <- "#FCFDBF"   # pale yellow  — model boxes (light, dark text)
col_eval   <- "#000004"   # near-black   — evaluation
txt_light  <- "white"
txt_dark   <- "grey20"
arrow_col  <- "grey40"

# ── helper: draw a rounded rectangle label ───────────────────────────────────
box <- function(x, y, w, h, label, fill, colour = "white",
                size = 3, fontface = "plain", text_colour = NULL) {
  tc <- if (is.null(text_colour)) {
    if (fill %in% c(col_model)) txt_dark else txt_light
  } else text_colour

  list(
    annotate("rect",
             xmin = x - w/2, xmax = x + w/2,
             ymin = y - h/2, ymax = y + h/2,
             fill = fill, colour = colour, linewidth = 0.4),
    annotate("text",
             x = x, y = y, label = label,
             colour = tc, size = size, fontface = fontface,
             lineheight = 0.95, hjust = 0.5, vjust = 0.5)
  )
}

# ── helper: straight arrow ───────────────────────────────────────────────────
arr <- function(x1, y1, x2, y2) {
  annotate("segment",
           x = x1, xend = x2, y = y1, yend = y2,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"))
}

# ── helper: label on a line ──────────────────────────────────────────────────
lbl <- function(x, y, label, size = 2.5, colour = "grey50") {
  annotate("text", x = x, y = y, label = label,
           colour = colour, size = size, hjust = 0.5)
}

# =============================================================================
# PANEL 1 — Data Pipeline
# =============================================================================
bw <- 3.6   # box width
bh <- 0.65  # box height
cx <- 5     # centre x

p_data <- ggplot() +
  # ── boxes ──────────────────────────────────────────────────────────────────
  box(cx, 13.0, bw, bh,
      "IPUMS ACS 2016–2023\n22.8 M individuals",
      col_data, size = 2.8, fontface = "bold") +

  box(cx, 11.2, bw, bh,
      "Basic filters\n(recently married · metro · co-residing)",
      col_data, size = 2.7) +

  box(cx,  9.4, bw, bh,
      "Match each woman to her husband\n(same household · spouse identifier)",
      col_data, size = 2.7) +

  box(cx,  7.6, bw, bh,
      "Remove thin markets\n(< 100 men per metro–year)",
      col_data, size = 2.7) +

  box(cx,  5.8, bw, bh,
      "Sample 2,000 focal women",
      col_data, size = 2.7) +

  box(cx,  4.0, bw, bh,
      "Build choice sets\ntrue husband  +  9 random men",
      col_data, size = 2.7) +

  box(cx,  2.2, bw, bh,
      "Feature engineering\n(V1 encoding  ·  V2 encoding)",
      col_feat, size = 2.7, fontface = "bold") +

  # ── arrows ─────────────────────────────────────────────────────────────────
  arr(cx, 12.67, cx, 11.55) +
  lbl(cx, 12.10, "325,808 women remain") +

  arr(cx, 10.87, cx, 9.72) +
  lbl(cx, 10.30, "307,511 matched") +

  arr(cx,  9.07, cx, 7.93) +
  lbl(cx,  8.50, "258,861 eligible women") +

  arr(cx,  7.27, cx, 6.13) +

  arr(cx,  5.47, cx, 4.33) +
  lbl(cx,  4.90, "20,000 rows  (2,000 × 10)") +

  arr(cx,  3.67, cx, 2.53) +

  # ── theme ──────────────────────────────────────────────────────────────────
  scale_x_continuous(limits = c(2.5, 7.5)) +
  scale_y_continuous(limits = c(1.5, 13.8)) +
  labs(title = "Data Pipeline") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11,
                                  hjust = 0.5, margin = margin(b = 4)))


# =============================================================================
# PANEL 2 — Modelling Pipeline
# =============================================================================
# Two branches: V1 (left, x=3) and V2 (right, x=7)
# Then all four models → evaluation

bw2 <- 2.8

p_model <- ggplot() +

  # ── Feature engineering fork ───────────────────────────────────────────────
  box(5, 13.0, bw2, bh,
      "Choice set data\n(from data pipeline)",
      col_feat, size = 2.7, fontface = "bold") +

  # V1 branch
  box(3,  11.0, bw2, bh,
      "V1 Encoding\nabs. differences + direction",
      col_v1, size = 2.7) +

  # V2 branch
  box(7,  11.0, bw2, bh,
      "V2 Encoding\nsigned diffs + absolute levels",
      col_v2, size = 2.7, text_colour = txt_dark) +

  # Models
  box(3,   8.8, bw2, bh,
      "M1 — Linear Clogit\n6 parameters",
      col_v1, size = 2.6) +

  box(6,   8.8, bw2, bh,
      "M2 — Spline Clogit\n9 parameters",
      col_v2, size = 2.6, text_colour = txt_dark) +

  box(8.2, 7.5, bw2, bh,
      "Forest (core)\n5 predictors",
      col_v2, size = 2.6, text_colour = txt_dark) +

  box(8.2, 6.2, bw2, bh,
      "Forest (full)\n9 predictors",
      col_v2, size = 2.6, text_colour = txt_dark) +

  # Train / test split
  box(5,   5.0, bw2, bh,
      "Train / Test Split\n70% train  ·  30% test (stratum level)",
      col_data, size = 2.6) +

  # Evaluation boxes
  box(3,   3.0, bw2, bh,
      "N-alternative test\ntop-1 / top-3 / top-5 accuracy",
      col_eval, size = 2.6) +

  box(7,   3.0, bw2, bh,
      "Full metro market\npercentile rank of true husband",
      col_eval, size = 2.6) +

  # ── arrows ─────────────────────────────────────────────────────────────────
  # Fork from choice set data
  annotate("segment", x = 5, xend = 3, y = 12.67, yend = 11.33,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("segment", x = 5, xend = 7, y = 12.67, yend = 11.33,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +

  # V1 → M1_v1
  arr(3, 10.67, 3, 9.13) +

  # V2 → M1_spline
  annotate("segment", x = 7, xend = 6, y = 10.67, yend = 9.13,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +

  # V2 → Forest_core
  annotate("segment", x = 7, xend = 8.2, y = 10.67, yend = 7.83,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +

  # Forest_core → Forest_full
  arr(8.2, 7.17, 8.2, 6.53) +
  lbl(8.8, 6.85, "+ race\ndummies", size = 2.2) +

  # All models → train/test split
  annotate("segment", x = 3,   xend = 4.5, y = 8.47, yend = 5.33,
           colour = arrow_col, linewidth = 0.4, linetype = "dashed") +
  annotate("segment", x = 6,   xend = 5.3, y = 8.47, yend = 5.33,
           colour = arrow_col, linewidth = 0.4, linetype = "dashed") +
  annotate("segment", x = 8.2, xend = 5.5, y = 5.87, yend = 5.33,
           colour = arrow_col, linewidth = 0.4, linetype = "dashed") +

  # Train/test → evaluation
  annotate("segment", x = 5, xend = 3, y = 4.67, yend = 3.33,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("segment", x = 5, xend = 7, y = 4.67, yend = 3.33,
           colour = arrow_col, linewidth = 0.5,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +

  # ── legend ─────────────────────────────────────────────────────────────────
  annotate("rect",  xmin = 1.0, xmax = 2.2, ymin = 1.4, ymax = 1.8,
           fill = col_v1, colour = "white", linewidth = 0.3) +
  annotate("text",  x = 2.4, y = 1.6, label = "V1 branch",
           hjust = 0, size = 2.4, colour = col_v1) +
  annotate("rect",  xmin = 4.2, xmax = 5.4, ymin = 1.4, ymax = 1.8,
           fill = col_v2, colour = "white", linewidth = 0.3) +
  annotate("text",  x = 5.6, y = 1.6, label = "V2 branch",
           hjust = 0, size = 2.4, colour = txt_dark) +

  scale_x_continuous(limits = c(0.8, 10.0)) +
  scale_y_continuous(limits = c(1.2, 13.8)) +
  labs(title = "Modelling Pipeline") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11,
                                  hjust = 0.5, margin = margin(b = 4)))


# =============================================================================
# Combine and save
# =============================================================================
p_combined <- p_data | p_model

ggsave(OUT_FILE, plot = p_combined, width = 11, height = 7, dpi = 300,
       bg = "white")
cat("Saved:", OUT_FILE, "\n")
