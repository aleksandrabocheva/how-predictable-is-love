# gen_fig_market_sizes.R
# Produces text_thesis/figures/07_market_sizes.png
# Histogram of metropolitan marriage market sizes (eligible men per metro x year cell).
# Run from the project root:  Rscript experiments/gen_fig_market_sizes.R

user_lib <- file.path(Sys.getenv("HOME"), "R", "packages")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ipumsr)
  library(scales)
})

DDI_FILE          <- "data/usa_00008.xml"
MIN_MEN_PER_METRO <- 100
OUT_FILE          <- "text_thesis/figures/07_market_sizes.png"

pal <- viridis::viridis(8, begin = 0.05, end = 0.78)

cat("Loading IPUMS data...\n")
ddi  <- read_ipums_ddi(DDI_FILE)
data <- read_ipums_micro(
  ddi,
  vars = c("YEAR", "SERIAL", "PERNUM", "SPLOC",
           "SEX", "MET2013", "YRMARR")
) %>% zap_labels()

cat("Building men pool...\n")
men_pool <- data %>%
  filter(YRMARR > 0, SPLOC > 0, MET2013 > 0,
         (YEAR - YRMARR) <= 3,
         SEX == 1)

market_sizes <- men_pool %>%
  count(MET2013, YEAR, name = "n_men") %>%
  filter(n_men >= MIN_MEN_PER_METRO)

cat(sprintf("Metro-year cells (>= %d men): %s\n",
            MIN_MEN_PER_METRO, format(nrow(market_sizes), big.mark = ",")))

p <- ggplot(market_sizes, aes(x = n_men)) +
  geom_histogram(bins = 30, fill = pal[3], colour = "white",
                 alpha = 0.85, linewidth = 0.3) +
  scale_x_log10(
    labels = comma,
    breaks = c(100, 300, 1000, 3000),
    name   = "Eligible Men per Metro-Year Cell (log scale)"
  ) +
  scale_y_continuous(labels = comma,
                     name   = "Number of Metro-Year Cells") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, colour = "grey88"),
    axis.ticks        = element_blank()
  )

ggsave(OUT_FILE, plot = p, width = 7, height = 4.5, dpi = 300, bg = "white")
cat("Saved:", OUT_FILE, "\n")
