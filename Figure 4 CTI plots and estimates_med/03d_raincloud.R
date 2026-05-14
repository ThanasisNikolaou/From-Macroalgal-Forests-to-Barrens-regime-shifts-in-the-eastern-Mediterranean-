library(tidyverse)
library(ggdist)
library(cowplot)
library(ggtext)
library(here)

df <- read.csv(
  here("data", "med", "dataset_with_CTI_med.csv"),
  header = TRUE,
  sep = ","
)

# Ensure region is a factor with defined levels
df$region <- factor(df$region, levels = c("NAS", "SAL"))


# Colour Palette & Theme
colors     <- c("NAS" = "#2166ac", "SAL" = "#b2182b")
fill_alpha <- 0.55

# plot
theme_raincloud <- function() {
  theme_classic(base_size = 13) +
    theme(
      # Axes
      axis.line        = element_line(linewidth = 0.5, colour = "grey30"),
      axis.ticks       = element_line(linewidth = 0.4, colour = "grey30"),
      axis.text        = element_text(size = 12, colour = "grey15"),
      axis.title.x     = element_text(size = 13, face = "bold",
                                      colour = "grey10", margin = margin(t = 8)),
      axis.title.y     = element_text(size = 13, face = "bold",
                                      colour = "grey10", margin = margin(r = 8)),
      # Panel
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      # Legend
      legend.position  = "none",
      # Plot margins
      plot.margin      = margin(10, 15, 10, 10)
    )
}

# Statistical Tests
# --- CTI (exclude NAs — sites where all species were absent) ---
df_CTI <- df %>% filter(!is.na(CTI_med))

shapiro_CTI <- by(df_CTI$CTI_med, df_CTI$region, shapiro.test)
wilcox_CTI  <- wilcox.test(CTI_med ~ region, data = df_CTI, exact = FALSE)

# --- CHi ---
shapiro_CHi <- by(df$CHi, df$region, shapiro.test)
wilcox_CHi  <- wilcox.test(CHi ~ region, data = df, exact = FALSE)

# --- Descriptive Statistics ---
desc_stats <- function(data, var, group_var = "region") {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n      = n(),
      mean   = mean(.data[[var]], na.rm = TRUE),
      sd     = sd(.data[[var]], na.rm = TRUE),
      median = median(.data[[var]], na.rm = TRUE),
      IQR    = IQR(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = var)
}

stats_CTI <- desc_stats(df_CTI, "CTI_med")
stats_CHi <- desc_stats(df, "CHi")



# Build Significance Label Helper
sig_label <- function(p) {
  if (p < 0.001) "p < 0.001 ***"
  else if (p < 0.01)  paste0("p = ", format(round(p, 3), nsmall = 3), " **")
  else if (p < 0.05)  paste0("p = ", format(round(p, 3), nsmall = 3), " *")
  else paste0("p = ", format(round(p, 3), nsmall = 3), " ns")
}

label_CTI <- sig_label(wilcox_CTI$p.value)
label_CHi <- sig_label(wilcox_CHi$p.value)


# Raincloud Plot Function
raincloud_horizontal <- function(data,
                                 yvar,
                                 xlab_title,
                                 stat_label,
                                 test_type = "Wilcoxon rank-sum test") {

  # Numeric y positions: NAS = 1, SAL = 0  (gap = 1 unit)
  data <- data %>%
    mutate(
      region_y = ifelse(region == "NAS", 1, 0),
      region_f = factor(region, levels = c("NAS", "SAL"))
    )

  # x-axis range for positioning annotations
  x_range  <- range(data[[yvar]], na.rm = TRUE)
  x_span   <- diff(x_range)
  x_annot  <- x_range[1] + x_span * 0.02   # left-justify annotation

  p <- ggplot(data,
              aes(y     = region_y,
                  x     = .data[[yvar]],
                  fill  = region_f,
                  color = region_f,
                  group = region_f)) +

    # Density (half-eye) 
    stat_halfeye(
      orientation   = "horizontal",
      adjust        = 1.4,        
      width         = 0.55,
      .width        = 0,
      trim          = T,
      expand        = TRUE,     
      justification = -0.22,
      point_colour  = NA,
      alpha         = fill_alpha
    ) +

    # Boxplot 
    geom_boxplot(
      aes(group = region_y),
      width         = 0.12,
      outlier.shape = NA,
      alpha         = 0.70,
      linewidth     = 0.55
    ) +

    # Raw data points
    geom_point(
      position = position_jitter(height = 0.06, seed = 42),
      size     = 1.8,
      alpha    = 0.55,
      shape    = 16
    ) +

    # Scales 
    scale_fill_manual(values  = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("SAL", "NAS"),
      expand = expansion(mult = c(0.05, 0.45))   # room for densities above
    ) +

    #  Significance bracket
    # Horizontal line
    annotate("segment",
             x = x_range[1] + x_span * 0.60,
             xend = x_range[1] + x_span * 0.60,
             y = 0, yend = 1,
             colour = "grey20", linewidth = 0.5) +
    annotate("segment",
             x = x_range[1] + x_span * 0.57,
             xend = x_range[1] + x_span * 0.60,
             y = 0, yend = 0,
             colour = "grey20", linewidth = 0.5) +
    annotate("segment",
             x = x_range[1] + x_span * 0.57,
             xend = x_range[1] + x_span * 0.60,
             y = 1, yend = 1,
             colour = "grey20", linewidth = 0.5) +
    annotate("text",
             x     = x_range[1] + x_span * 0.625,
             y     = 0.5,
             label = stat_label,
             size  = 3.5,
             hjust = 0,
             colour = "grey20",
             fontface = "italic") +

    # Test type label (bottom-left) 
    annotate("text",
             x     = x_annot,
             y     = -0.38,
             label = test_type,
             hjust = 0,
             size  = 4,
             colour = "grey20",
             fontface = "italic") +

    # Axes & labels
    labs(x = xlab_title, y = "Region") +

    theme_raincloud()

  return(p)
}


# Create Individual Plots
plot_CTI <- raincloud_horizontal(
  data       = df_CTI,
  yvar       = "CTI_med",
  xlab_title = expression("CTI (°C)"),
  stat_label = label_CTI,
  test_type  = "Wilcoxon rank-sum test"
)

plot_CHi <- raincloud_horizontal(
  data       = df,
  yvar       = "CHi",
  xlab_title = expression(bold(paste("CHI (g consumed/day/m"^2, ")"))),
  stat_label = label_CHi,
  test_type  = "Wilcoxon rank-sum test"
)


# Combine Plots
final_plot <- cowplot::plot_grid(
  plot_CTI, plot_CHi,
  ncol   = 1,
  nrow   = 2,
  labels = c("a", "b"),
  label_size = 14,
  label_fontface = "bold",
  rel_heights = c(1, 1)
)

# Add a shared title
title_row <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Community Temperature & Heat Indices: NAS vs SAL",
    fontface  = "bold",
    size      = 15,
    colour    = "grey10",
    x         = 0.5,
    hjust     = 0.5
  )

final_plot_titled <- cowplot::plot_grid(
  title_row, final_plot,
  ncol        = 1,
  rel_heights = c(0.06, 1)
)

print(final_plot_titled)


# Export
out_dir <- here("data", "med")

# High-resolution PNG
ggsave(
  filename = file.path(out_dir, "Raincloud_CTI_CHi_NAS_SAL.png"),
  plot     = final_plot_titled,
  width    = 10,
  height   = 9,
  dpi      = 300,
  bg       = "white"
)

