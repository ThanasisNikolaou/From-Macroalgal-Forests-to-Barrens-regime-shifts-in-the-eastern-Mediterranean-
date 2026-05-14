# Script 3 of 3 — Mediterranean CTI pipeline
# Compute Species Thermal Index (Ti) and Community Thermal Index (CTI)
#
# This script brings together the SST-matched Mediterranean occurrence records from Script 2 and the survey dataset to 
# produce the final CTI values. Ti is estimated from Mediterranean occurrences only.
# Steps:
#   1. Load SST-matched occurrence records from Script 2
#   2. Remove non-Mediterranean records (Atlantic Portugal + Bay of Biscay)
#   3. Plot the SST distribution for each species (1 °C bins) with Q5, Q95, mean and Ti marked
#   4. Compute Ti = (Q5 + Q95) / 2 for each species
#   5. Load the survey dataset
#   6. Map occurrence records for each species
#   7. Compute CTI for each site
#   8. Save the Ti table and the CTI-enriched survey dataset
#
# Input:  data/med/occurrences_med_sst.rds
#         data/dataset_eastMed_reefs_abundance_density.csv
#           (place the survey CSV inside the project's data/ folder)
# Output: data/med/Ti_med.csv
#         data/med/dataset_with_CTI_med.csv
#         data/med/maps/map_med_<species>.png


# --- Libraries -------------------------------------------------------

library(dplyr)            
library(ggplot2)           
library(here)              
library(sf)   
library(rnaturalearth)     
library(rnaturalearthdata)

# --- Species-to-column mapping ---------------------------------------
#
# Redefined here rather than read from the saved RDS.

species_col_map <- c(
  "Paracentrotus lividus" = "pl",
  "Arbacia lixula"        = "al",
  "Sarpa salpa"           = "ss",
  "Siganus luridus"       = "sl",
  "Siganus rivulatus"     = "sr",
  "Sparisoma cretense"    = "sc"
)


# --- Load data from Script 2 -----------------------------------------

cat("=== Loading SST-matched Mediterranean occurrences ===\n")

saved        <- readRDS(here("data", "med",
                             "occurrences_med_sst.rds"))
all_occs_sst <- saved$occurrences_sst   # named list: one df per species
species_list <- saved$species_list      # character vector of species names

# Recover the bounding box so maps use the same extent
MED_LON_MIN <- saved$MED_LON_MIN
MED_LON_MAX <- saved$MED_LON_MAX
MED_LAT_MIN <- saved$MED_LAT_MIN
MED_LAT_MAX <- saved$MED_LAT_MAX


# --- Remove non-Mediterranean records and known erroneous points -----
#
# The Mediterranean bounding box used in 03a extends west to -18.1°, which inadvertently captures records from coastal 
# Portugal and the Bay of Biscay. These are removed here without re-running the time-consuming fetch and SST-matching steps.
#
# Rules applied:
#   - Exclude records west of the Gibraltar Strait (lon < -6)
#   - Exclude Bay of Biscay records (lon < 0 AND lat > 42)
#   - Exclude a known on-land OBIS record for Siganus luridus in Cyprus
#     (lon 32–35, lat 34–36, source OBIS) identified during visual inspection

cat("\n=== Removing non-Mediterranean and erroneous records ===\n")

all_occs_sst <- lapply(all_occs_sst, function(df) {
  if (is.null(df)) return(NULL)
  n_before <- nrow(df)
  df <- df %>%
    filter(
      decimalLongitude >= -6,                           # exclude Atlantic/Portugal
      !(decimalLongitude < 0 & decimalLatitude > 42),  # exclude Bay of Biscay
      # Exclude known on-land OBIS record for Siganus luridus in Cyprus. 
      # There is only one OBIS record for this species in that window, so the broad bounds are unambiguous.
      !(species == "Siganus luridus" &
        decimalLongitude < 33        &
        decimalLatitude  > 34        &
        decimalLatitude  < 35        &
        source == "OBIS")
    )
  n_removed <- n_before - nrow(df)
  if (n_removed > 0)
    cat(sprintf("  %-25s  removed %d record(s)\n",
                df$species[1], n_removed))
  df
})

#  check —> should print 0
n_leaked <- lapply(all_occs_sst, function(df) {
  if (is.null(df)) return(0L)
  nrow(filter(df, decimalLongitude < -6 |
                  (decimalLongitude < 0 & decimalLatitude > 42)))
}) %>% unlist() %>% sum()
cat(sprintf("  Records still outside Mediterranean area: %d\n", n_leaked))


# Verify species_col_map covers all species
if (!all(species_list %in% names(species_col_map))) {
  missing <- species_list[!species_list %in% names(species_col_map)]
  stop(paste("species_col_map missing entries for:",
             paste(missing, collapse = ", ")))
}

cat("  Records per species (Mediterranean, with valid GODAS SST):\n")
for (sp in species_list) {
  n_ok <- sum(!is.na(all_occs_sst[[sp]]$godas_sst))
  cat(sprintf("  %-25s  %d\n", sp, n_ok))
}


# --- SST distribution plots ------------------------------------------
#
# One histogram per species showing the GODAS SST distribution VS all Mediterranean occurrence records with a valid SST value.

plot_sst_hist <- function(sst_vals, sp_name,
                          q5, q95, sst_mean, Ti,
                          bar_col = "lightgray") {
  
  if (length(sst_vals) < 5) return(invisible(NULL))
  
  hist(
    sst_vals,
    breaks = seq(floor(min(sst_vals,   na.rm = TRUE)),
                 ceiling(max(sst_vals, na.rm = TRUE)),
                 by = 0.5),                     # 0.5°C bins
    main   = paste(sp_name,
                   sprintf("— n = %d  |  Ti = %.2f°C",
                           length(sst_vals), Ti)),
    xlab   = "GODAS SST (°C)",
    ylab   = "Number of records",
    col    = bar_col,
    border = "white"
  )
  
  abline(v = c(q5, q95), col = "red",       lwd = 2, lty = 2)
  abline(v = sst_mean,   col = "blue",      lwd = 2, lty = 1)
  abline(v = Ti,         col = "darkgreen", lwd = 2, lty = 3)
  
  legend("topright",
         legend = c("Q5 / Q95", "Mean SST", "Ti"),
         col    = c("red", "blue", "darkgreen"),
         lty    = c(2, 1, 3), lwd = 2, cex = 0.8)
}

cat("\n=== Plotting and saving SST distributions ===\n")

hist_dir <- here("data", "med", "histograms")
dir.create(hist_dir, showWarnings = FALSE, recursive = TRUE)

for (sp in species_list) {
  
  sst_vals <- all_occs_sst[[sp]]$godas_sst
  sst_vals <- sst_vals[!is.na(sst_vals)]
  
  if (length(sst_vals) > 20) {
    
    sp_clean <- gsub(" ", "_", sp)
    png_path <- file.path(hist_dir, sprintf("hist_med_%s.png", sp_clean))
    
    png(png_path, width = 8, height = 5, units = "in", res = 300)
    
    plot_sst_hist(
      sst_vals = sst_vals,
      sp_name  = sp,
      q5       = as.numeric(quantile(sst_vals, 0.05)),
      q95      = as.numeric(quantile(sst_vals, 0.95)),
      sst_mean = mean(sst_vals),
      Ti       = as.numeric(
        (quantile(sst_vals, 0.05) + quantile(sst_vals, 0.95)) / 2
      ),
      bar_col  = "lightgray"
    )
    
    dev.off()
    
    cat(sprintf("  %-25s  saved: %s\n", sp, basename(png_path)))
  }
}

# --- Compute Species Thermal Index (Ti) ------------------------------
#
# Ti = (Q5 + Q95) / 2

compute_Ti <- function(sp_name, occs_list) {

  occs <- occs_list[[sp_name]]
  if (is.null(occs)) return(NULL)

  sst_vals <- as.numeric(occs$godas_sst[!is.na(occs$godas_sst)])

  cat(sprintf("  %-25s  n SST records: %d\n", sp_name, length(sst_vals)))

  if (length(sst_vals) < 20) {
    warning(paste("Too few SST records for", sp_name, "— Ti set to NA"))
    return(NULL)
  }

  q5  <- as.numeric(unname(quantile(sst_vals, probs = 0.05)))
  q95 <- as.numeric(unname(quantile(sst_vals, probs = 0.95)))
  Ti  <- as.numeric((q5 + q95) / 2)

  # returns 0 if a category is absent from the table, avoiding zero-length NA entries in the result
  src     <- table(occs$source[!is.na(occs$godas_sst)])
  get_src <- function(label) {
    v <- src[label]
    if (length(v) == 0 || is.na(v)) 0L else as.integer(v)
  }

  # Build as a list and verify all elements are length 1 before calling as.data.frame() to prevent row mismatch errors
  result <- list(
    Species    = sp_name,
    AbundCol   = as.character(unname(species_col_map[sp_name])),
    n_total    = as.integer(length(sst_vals)),
    n_OBIS     = get_src("OBIS"),
    n_GBIF     = get_src("GBIF"),
    n_BOTH     = get_src("BOTH"),
    SST_mean   = round(mean(sst_vals),   4),
    SST_median = round(median(sst_vals), 4),
    SST_q5     = round(q5,               4),
    SST_q95    = round(q95,              4),
    Ti         = round(Ti,               4)
  )

  lens <- sapply(result, length)
  if (any(lens != 1)) {
    cat("  DEBUG — element lengths for", sp_name, ":\n")
    print(lens[lens != 1])
    stop(paste("Non-scalar element in compute_Ti for", sp_name))
  }

  as.data.frame(result, stringsAsFactors = FALSE)
}

Ti_results <- bind_rows(
  lapply(species_list, compute_Ti, occs_list = all_occs_sst)
)

cat("\n=== Ti results (Mediterranean) ===\n")
print(Ti_results)

missing_sp <- setdiff(species_list, Ti_results$Species)
if (length(missing_sp) > 0) {
  cat("\n\u26a0 WARNING: No Ti computed for:\n")
  cat(paste(" ", missing_sp, collapse = "\n"), "\n")
  cat("  These species will be excluded from CTI.\n")
}


# --- Load survey dataset ---------------------------------------------

dataset <- read.csv(
  here("data", "dataset_eastMed_reefs_abundance_density.csv"),
  header = TRUE,
  sep    = ","
)

# Verify all required abundance columns are present
required_cols <- unname(species_col_map)   # "pl", "al", "ss", "sl", "sr", "sc"
missing_cols  <- required_cols[!required_cols %in% names(dataset)]
if (length(missing_cols) > 0) {
  stop(paste("Missing abundance columns in dataset:",
             paste(missing_cols, collapse = ", ")))
}

cat(sprintf("\n  Survey dataset loaded: %d sites\n", nrow(dataset)))


# --- Occurrence maps — one map per species, Mediterranean extent -----
#
map_dir <- here("data", "med", "maps")
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

# Load world land polygons
world_land <- rnaturalearth::ne_countries(
  scale       = "medium",   # 1:50m resolution
  returnclass = "sf"
)

for (sp in species_list) {

  occ_map <- all_occs_sst[[sp]] %>%
    filter(
      !is.na(decimalLongitude),
      !is.na(decimalLatitude),
      !is.na(godas_sst)
    )

  if (nrow(occ_map) == 0) {
    cat(sprintf("  %-25s  no valid records to map — skipping\n", sp))
    next
  }

  Ti_row    <- Ti_results %>% filter(Species == sp)
  Ti_val    <- if (nrow(Ti_row) > 0) round(Ti_row$Ti, 2) else NA
  src_tab   <- table(occ_map$source)
  src_str   <- paste(names(src_tab), src_tab, sep = " = ", collapse = "   ")
  map_title <- sprintf("%s   (n = %d  |  Ti = %.2f°C)",
                       sp, nrow(occ_map), Ti_val)
  sst_min   <- floor(min(occ_map$godas_sst,   na.rm = TRUE))
  sst_max   <- ceiling(max(occ_map$godas_sst, na.rm = TRUE))

  p <- ggplot() +

    theme_bw() +
    theme(
      panel.background  = element_rect(fill = "#D6EAF8"),  # light blue ocean
      # Fonts sized for a 6-panel figure (3 col x 2 row): each panel is
      # displayed at ~1/3 width and ~1/2 height, so all text is scaled
      # up to remain readable after assembly in Photoshop.
      plot.title        = element_text(face = "italic", size = 20),
      plot.subtitle     = element_text(size = 15, colour = "gray40"),
      legend.position   = "right",
      legend.title      = element_text(size = 16),
      legend.text       = element_text(size = 14),
      legend.key.height = unit(2, "cm"),
      panel.grid.major  = element_line(colour = "white", linewidth = 0.3),
      axis.text         = element_text(size = 14),
      axis.title        = element_text(size = 16)
    ) +

    # Land layer drawn first so occurrence points appear on top
    geom_sf(data = world_land, fill = "grey70", colour = NA) +

    # Occurrence points coloured by GODAS SST
    # Point size increased so dots remain visible at panel scale
    geom_point(
      data  = occ_map,
      aes(x = decimalLongitude, y = decimalLatitude, colour = godas_sst),
      size  = 2.0,
      alpha = 0.7
    ) +

    # Plasma palette: no white, perceptually uniform, colour-blind safe
    scale_colour_viridis_c(
      option = "plasma",
      direction = 1,
      limits = c(sst_min, sst_max),
      name   = "GODAS SST (°C)",
      guide  = guide_colourbar(
        barheight   = unit(8, "cm"),
        label.theme = element_text(size = 14),
        title.theme = element_text(size = 16)
      )
    ) +

    # Zoom to the Mediterranean bounding box (western edge clipped to -6)
    coord_sf(
      xlim   = c(max(MED_LON_MIN, -6), MED_LON_MAX),
      ylim   = c(MED_LAT_MIN, MED_LAT_MAX),
      expand = FALSE
    ) +

    labs(
      title    = map_title,
      subtitle = paste("Source:", src_str),
      x        = "Longitude",
      y        = "Latitude"
    )

  print(p)

  sp_clean <- gsub(" ", "_", sp)
  png_path <- file.path(map_dir, sprintf("map_med_%s.png", sp_clean))

  ggsave(
    filename = png_path,
    plot     = p,
    width    = 12,
    height   = 7,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )

  cat(sprintf("  %-25s  mapped: %d  |  saved: %s\n",
              sp, nrow(occ_map), basename(png_path)))
}


# --- Build Ti vector aligned to abundance columns -------------------
#
# CTI formula:
#   CTI = sum(log10(N_i + 1) * Ti) / sum(log10(N_i + 1))
#

Ti_clean <- Ti_results %>% filter(!is.na(Ti))
Ti_vec   <- setNames(Ti_clean$Ti, Ti_clean$AbundCol) 

print(Ti_vec)

abund_cols <- names(Ti_vec)
abund_mat  <- dataset[, abund_cols, drop = FALSE]

# Log10-transform matrix
abund_log  <- log10(abund_mat + 1)


# --- Calculate CTI for each survey site ------------------------------
#

dataset$CTI_med <- apply(abund_log, 1, function(x) {
  denom <- sum(x, na.rm = TRUE)
  if (denom == 0) return(NA_real_)
  sum(x * Ti_vec, na.rm = TRUE) / denom
})


# --- Inspect and summarise results -----------------------------------

cat("\n=== CTI summary (Mediterranean GODAS Ti) ===\n")
print(summary(dataset$CTI_med))

# Mean CTI and CHi by region
cat("\n=== Mean CTI and CHi by region ===\n")
region_means <- dataset %>%
  group_by(region) %>%
  summarise(
    n        = n(),
    mean_CTI = round(mean(CTI_med, na.rm = TRUE), 3),
    sd_CTI   = round(sd(CTI_med,   na.rm = TRUE), 3),
    mean_CHi = round(mean(CHi,     na.rm = TRUE), 3),
    .groups  = "drop"
  )
print(region_means)


# --- Save outputs ----------------------------------------------------

out_dir <- here("data", "med")

write.csv(
  Ti_results,
  file.path(out_dir, "Ti_med.csv"),
  row.names = FALSE
)

write.csv(
  dataset,
  file.path(out_dir, "dataset_with_CTI_med.csv"),
  row.names = FALSE
)

cat("\n=== Script 3 (Mediterranean) complete ===\n")
cat("Ti table    : data/med/Ti_med.csv\n")
cat("CTI dataset : data/med/dataset_with_CTI_med.csv\n")
cat("Maps        : data/med/maps/\n")
