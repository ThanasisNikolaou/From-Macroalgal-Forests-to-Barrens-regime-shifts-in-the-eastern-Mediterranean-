# Script 02a — NMDS ordination (fourth-root transformation, Bray-Curtis)
#
# Runs NMDS on fourth-root transformed herbivore abundance data,  plots the ordination coloured by region (NAS / SAL), and saves
# the NMDS object for use in Script 02b (regressions & changepoints).
#
#
# Input:  dataset_eastMed_reefs_abundance_density.csv
# Output: nmds_fourth_root.rds
#         NMDS_fourth_root.png

library(vegan)
library(ggplot2)
library(here)


# --- Load data -------------------------------------------------------

data <- read.csv(
  here("dataset_eastMed_reefs_abundance_density.csv"),
  header = TRUE
)


# --- Fourth-root transformation --------------------------------------

species_cols <- c("pl", "al", "ss", "sl", "sr", "sc")
species_data <- data[, species_cols]

# Fourth-root trasnformation
species_transformed <- species_data^0.25

# Drop rows that are all zeros after transformation 
empty_rows <- rowSums(species_transformed) == 0
if (any(empty_rows)) {
  cat("Removing", sum(empty_rows), "empty row(s):", which(empty_rows), "\n")
  species_transformed <- species_transformed[!empty_rows, ]
  data               <- data[!empty_rows, ]
}


# --- Run NMDS --------------------------------------------------------

set.seed(123)
nmds <- metaMDS(species_transformed, distance = "bray", k = 2, trymax = 100)

cat("Stress value:", nmds$stress, "\n")


# --- Save NMDS object for Script 02b ---------------------------------
#

out_dir <- here("data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(nmds, file.path(out_dir, "nmds_fourth_root.rds"))
cat("NMDS object saved to: data/nmds_fourth_root.rds\n")


# --- Plot ordination -------------------------------------------------

nmds_scores         <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$region  <- data$region

region_colors <- c("NAS" = "#7aa3cd", "SAL" = "#d17480")

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = region)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = region_colors) +
  stat_ellipse(aes(group = region), type = "t",
               linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Stress = ", round(nmds$stress, 3)),
           hjust = 1.1, vjust = 1.5, size = 3.8, color = "grey30") +
  labs(title = "NMDS",
       x = "NMDS1", y = "NMDS2", color = "Region") +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

fig_dir <- here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(fig_dir, "NMDS_fourth_root.png"),
  dpi = 300, width = 8, height = 10, units = "in"
)

