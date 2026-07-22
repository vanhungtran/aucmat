# ===========================================================================
# aucmat hex logo — stylized ROC curve with biomarker ranking
# Design: ROC curve (discrimination) + ranked dots (biomarker screening)
# ===========================================================================
library(hexSticker)
library(ggplot2)
library(showtext)

font_add_google("Righteous", "righteous")
showtext_auto()

set.seed(42)

# ── Generate a clean ROC-style curve + colorful biomarker dots ─────
# Background: smooth ROC curve spanning the hex space
roc_x <- seq(0, 1, length.out = 200)
roc_y <- 1 - exp(-3.5 * roc_x)  # realistic concave ROC shape

roc_df <- data.frame(x = roc_x, y = roc_y)

# Diagonal reference line
diag_df <- data.frame(x = c(0, 1), y = c(0, 1))

# Biomarker dots: scattered along/near the ROC curve
# Use the aucmat accent palette for colorful dots
accent_cols <- c("#2c7da0", "#EE7733", "#009988", "#CC3311",
                 "#33BBEE", "#AA3377", "#228833", "#4477AA",
                 "#DDCC77", "#66CCEE", "#88CCAA", "#CC6677")

n_dots <- 20
dots <- data.frame(
  x = sort(runif(n_dots, 0.06, 0.90)),
  y = 1 - exp(-3.4 * sort(runif(n_dots, 0.06, 0.90)))
)
# Push some dots above the curve (strong biomarkers)
dots$y[1:8]  <- dots$y[1:8] + runif(8, 0, 0.15)
dots$y <- pmin(dots$y, 0.92)
# Assign colors cycling through the palette
dots$colour <- rep(accent_cols, length.out = n_dots)
# Vary size: stronger biomarkers = larger
dots$size <- seq(3.5, 1.8, length.out = n_dots)

p <- ggplot() +
  # Diagonal reference
  geom_line(data = diag_df, aes(x = x, y = y),
            linetype = "dashed", linewidth = 0.25,
            colour = "#8899AA", alpha = 0.7) +
  # ROC curve
  geom_line(data = roc_df, aes(x = x, y = y),
            linewidth = 1.8, colour = "#2c7da0") +
  # Colorful biomarker dots
  geom_point(data = dots, aes(x = x, y = y, size = size, colour = I(colour)),
             alpha = 0.92) +
  scale_size_identity() +
  coord_fixed(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

sticker(p,
  package        = "aucmat",
  p_size         = 22,
  p_color        = "#1a5276",
  p_family       = "righteous",
  p_y            = 1.48,
  s_x            = 1,
  s_y            = 0.73,
  s_width        = 1.3,
  s_height       = 1.05,
  h_fill         = "#F8F9FA",
  h_color        = "#2c7da0",
  h_size         = 1.8,
  filename       = "man/figures/logo.png",
  dpi            = 300,
  white_around_sticker = FALSE
)

cat("Logo saved to man/figures/logo.png\n")
