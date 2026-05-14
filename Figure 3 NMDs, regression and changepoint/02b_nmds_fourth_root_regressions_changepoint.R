# Script 02b — NMDS1 regressions and changepoint analysis
#
# Loads the NMDS object saved by Script 02a -> regresses NMDS1 scores against annual average, minimum, and maximum
# temperature, and runs changepoint analysis (PELT / BIC) on NMDS1 sorted by each temperature variable.
#
# Input:  data/nmds_fourth_root.rds   (NMDS object from Script 02a)
#         data/dataset_eastMed_reefs_abundance_density.csv
# Output: figures/NMDS1_vs_temperature.png
#         figures/NMDS1_vs_temperature_changepoints.png

library(vegan)
library(ggplot2)
library(patchwork)
library(changepoint)
library(here)


# --- Load data -------------------------------------------------------

data <- read.csv(
  here("dataset_eastMed_reefs_abundance_density.csv"),
  header = TRUE
)


# --- Reproduce the same empty-row removal as Script 02a -------------
#

species_cols <- c("pl", "al", "ss", "sl", "sr", "sc")
species_data <- data[, species_cols]

species_transformed <- species_data^0.25

empty_rows <- rowSums(species_transformed) == 0
if (any(empty_rows)) {
  cat("Removing", sum(empty_rows), "empty row(s):", which(empty_rows), "\n")
  data <- data[!empty_rows, ]
}


# --- Load NMDS object saved by Script 02a ---------------------------
#

nmds <- readRDS(here("data", "nmds_fourth_root.rds"))
cat("Stress value:", nmds$stress, "\n")

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
plot_data   <- cbind(nmds_scores, data)

region_colors <- c("NAS" = "#7aa3cd", "SAL" = "#d17480")

fig_dir <- here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# --- Changepoint analysis (mean-only, PELT / BIC) -------------------
#
# NMDS1 scores are sorted by each temperature variable before analysis.
#

temp_vars   <- c("avg_temp", "min_temp", "max_temp")
temp_labels <- c("Average Temperature", "Minimum Temperature", "Maximum Temperature")

# cp_results stores the mean-only changepoint temperature value for each variable (NULL if none detected). 
cp_results <- list()

for (i in seq_along(temp_vars)) {

  var   <- temp_vars[i]
  label <- temp_labels[i]

  sorted_data  <- plot_data[order(plot_data[[var]]), ]
  nmds1_sorted <- sorted_data$NMDS1

  # Run all three variants for diagnostic completeness
  cp_meanvar <- cpt.meanvar(nmds1_sorted, method = "PELT", penalty = "BIC")
  cp_mean    <- cpt.mean(nmds1_sorted,    method = "PELT", penalty = "BIC")
  cp_var     <- cpt.var(nmds1_sorted,     method = "PELT", penalty = "BIC")

  cat("\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n")
  cat("NMDS1 vs", label, "\n")
  cat("\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\u2501\n")

  # Mean & Variance (diagnostic only — not plotted)
  cp_mv_locs <- cpts(cp_meanvar)
  if (length(cp_mv_locs) == 0) {
    cat("Mean & Variance [diagnostic]: No changepoint detected\n")
  } else {
    cat("Mean & Variance [diagnostic]: changepoint(s) at", var, "=",
        round(sorted_data[[var]][cp_mv_locs], 3), "\n")
  }

  # Mean only — used for plotting
  cp_m_locs <- cpts(cp_mean)
  if (length(cp_m_locs) == 0) {
    cat("Mean only [plotted]:          No changepoint detected\n")
    cp_results[[var]] <- NULL
  } else {
    cp_m_temps <- sorted_data[[var]][cp_m_locs]
    cat("Mean only [plotted]:          changepoint(s) at", var, "=",
        round(cp_m_temps, 3), "\n")
    # Store the first changepoint temperature value for this variable
    cp_results[[var]] <- cp_m_temps[1]
  }

  # Variance only (diagnostic only — not plotted)
  cp_v_locs <- cpts(cp_var)
  if (length(cp_v_locs) == 0) {
    cat("Variance only [diagnostic]:   No changepoint detected\n")
  } else {
    cat("Variance only [diagnostic]:   changepoint(s) at", var, "=",
        round(sorted_data[[var]][cp_v_locs], 3), "\n")
  }

  cat("\n")
}


# --- Plot function: NMDS1 vs temperature (with changepoint line) ----
#

make_plot <- function(x_var, x_label) {

  cp_value <- cp_results[[x_var]]

  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = NMDS1, color = region)) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey50", linewidth = 0.7) +
    geom_point(size = 3.5, alpha = 0.85) +
    # Single pooled regression line across both regions combined
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                linetype = "dashed", linewidth = 0.7) +
    scale_color_manual(values = region_colors) +
    labs(x = x_label, y = "NMDS1", color = "Region",
         title = paste("NMDS1 vs", x_label)) +
    theme_classic(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )

  if (!is.null(cp_value)) {
    p <- p +
      geom_vline(xintercept = cp_value, linetype = "dotted",
                 color = "black", linewidth = 0.8) +
      annotate("text", x = cp_value, y = Inf,
               label = paste0("Changepoint at ", round(cp_value, 2)),
               hjust = -0.1, vjust = 1.5, size = 5,
               color = "black", fontface = "italic")
  } else {
    p <- p +
      annotate("text", x = Inf, y = Inf,
               label = "No changepoint detected",
               hjust = 1.1, vjust = 1.5, size = 5,
               color = "black", fontface = "italic")
  }

  return(p)
}


# --- Produce and save figures ----------------------------------------

p1 <- make_plot("avg_temp", "Average Temperature")
p2 <- make_plot("min_temp", "Minimum Temperature")
p3 <- make_plot("max_temp", "Maximum Temperature")

combined_plot <- (p1 / p2 / p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "NMDS1 vs Temperature Variables",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

print(combined_plot)

ggsave(
  file.path(fig_dir, "NMDS1_vs_temperature_changepoints.png"),
  plot = combined_plot,
  dpi = 300, width = 10, height = 12, units = "in"
)

cat("Saved: figures/NMDS1_vs_temperature_changepoints.png\n")
