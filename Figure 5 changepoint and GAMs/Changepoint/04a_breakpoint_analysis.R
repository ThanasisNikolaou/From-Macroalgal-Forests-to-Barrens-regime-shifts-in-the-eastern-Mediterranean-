library(strucchange)
library(ggplot2)
library(dplyr)


# function to run the full breakpoint analysis for a given response variable
run_breakpoint_analysis <- function(data, y_col, x_col, y_label, file_prefix) {
  
  y <- data[[y_col]]
  x <- data[[x_col]]
  n <- length(y)
  
  canopy_ts <- ts(y)
  
  # structural change analysis
  fstat <- Fstats(canopy_ts ~ 1)
  sc_result <- sctest(fstat)
  print(sc_result)
  
  set.seed(123)
  bp <- breakpoints(canopy_ts ~ 1, breaks = 1)
  print(summary(bp))
  ci_bp <- confint(bp)
  
  # extract breakpoint location and CI bounds
  bp_index <- bp$breakpoints[1]
  breakpoint_value <- x[bp_index]
  ci_lower <- x[ci_bp$confint[1]]
  ci_upper <- x[ci_bp$confint[3]]
  
  # group means (breakpoint row goes to upper group)
  mean_lower <- mean(y[x < breakpoint_value],  na.rm = TRUE)
  mean_upper <- mean(y[x >= breakpoint_value], na.rm = TRUE)
  
  # align F-stats to the full data length
  plot_fstat <- c(NA, as.numeric(fstat$Fstats), NA, NA)
  fstat_df   <- data.frame(index = x, fstat = plot_fstat)
  
  # --- F-statistic plot ---
  p_fstat <- ggplot(fstat_df, aes(x = index, y = fstat)) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    geom_vline(xintercept = breakpoint_value, linetype = "dotted",
               color = "black", linewidth = 1) +
    labs(x = "Grazing index", y = "F Statistic") +
    theme_classic(base_size = 16) +
    theme(
      axis.text   = element_text(color = "black"),
      axis.line   = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  ggsave(paste0(file_prefix, "_fstat_plot.tiff"),
         p_fstat, width = 6, height = 5, dpi = 300)
  
  # --- Canopy cover plot ---
  plot_df <- data.frame(x = x, y = y)
  
  p_canopy <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(size = 2.5) +
    annotate("rect",
             xmin = ci_lower, xmax = ci_upper,
             ymin = -Inf, ymax = Inf,
             fill = "lightgrey", alpha = 0.4) +
    geom_vline(xintercept = breakpoint_value, linetype = "dotted",
               color = "black", linewidth = 1) +
    annotate("segment",
             x = min(x), xend = breakpoint_value,
             y = mean_lower, yend = mean_lower,
             color = "red", linewidth = 1.2) +
    annotate("segment",
             x = breakpoint_value, xend = max(x),
             y = mean_upper, yend = mean_upper,
             color = "red", linewidth = 1.2) +
    labs(x = "Grazing index", y = y_label) +
    ylim(0, 60) +
    theme_classic(base_size = 16) +
    theme(
      axis.text    = element_text(color = "black"),
      axis.line    = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  ggsave(paste0(file_prefix, "_breakpoint_plot.tiff"),
         p_canopy, width = 6, height = 5, dpi = 300)
  
  # return key values in case they're needed later
  invisible(list(
    breakpoint   = breakpoint_value,
    ci           = c(ci_lower, ci_upper),
    mean_lower   = mean_lower,
    mean_upper   = mean_upper,
    sctest       = sc_result
  ))
}


# --- average canopy cover ---
data <- read.delim("clipboard") # read (copy) the first two columns corresponding to the average coverage per bin
run_breakpoint_analysis(data,   y_col = "canopy",     x_col = "CH",
                        y_label = "Average canopy cover (%)",
                        file_prefix = "average_spp")

# --- maximum canopy cover ---
data_2 <- read.delim("clipboard") # read (copy) the latter two columns corresponding to the maximum coverage per bin
run_breakpoint_analysis(data_2, y_col = "canopy_max", x_col = "CH",
                        y_label = "Maximum canopy cover (%)",
                        file_prefix = "max_spp")
