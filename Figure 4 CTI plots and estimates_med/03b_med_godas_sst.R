# Script 2 of 3 — Mediterranean CTI pipeline
# Download GODAS SST and match to occurrence records
#
# This script downloads annual GODAS potential temperature files (1980–2021) from NOAA PSL, skipping any files already on disk.
# For each occurrence record from Script 1, it extracts the annual mean SST at the record's coordinates in its observation year.
#
# GODAS longitude is stored as 0–360°, so negative longitudes are
# converted before extraction.
#
# Input:  data/med/occurrences_med_raw.rds
# Output: data/med/occurrences_med_sst.rds
#
# GODAS product details:
#   Variable  : pottmp (potential temperature)
#   Depth     : surface layer (depth <= 5 m)
#   Resolution: 1/3° latitude x 1° longitude
#   Units     : Kelvin -> converted to Celsius here
#   Period    : 1980–2021, one file per year (~120 MB each)
#   Source    : https://psl.noaa.gov/data/gridded/data.godas.html


# --- Libraries -------------------------------------------------------

library(tidync)
library(dplyr)
library(here)


# --- Settings --------------------------------------------------------

# NOAA servers can be slow and files are ~120 MB each; 6000 seconds 
options(timeout = 6000)

YEAR_MIN <- 1980
YEAR_MAX <- 2021
years    <- YEAR_MIN:YEAR_MAX

# Rather than extracting a single grid cell, we average all GODAS cells within a small window around each occurrence point.
EXTRACT_RANGE_LON <- 0.5   # +/- 0.5° = one full 1° cell in longitude
EXTRACT_RANGE_LAT <- 1/6   # +/- 1/6° = one full 1/3° cell in latitude # matches one GODAS latitude cell

# Directory to store downloaded GODAS NetCDF files
GODAS_DIR <- here("data", "godas_nc")
dir.create(GODAS_DIR, showWarnings = FALSE, recursive = TRUE)

# Output directory
out_dir <- here("data", "med")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- Load occurrences from Script 1 ----------------------------------

cat("=== Loading Mediterranean occurrences from Script 1 ===\n")

saved        <- readRDS(file.path(out_dir, "occurrences_med_raw.rds"))
all_occs     <- saved$occurrences
species_list <- saved$species_list

cat("  Records per species (Mediterranean only):\n")
print(sapply(all_occs, nrow))


# --- Download GODAS NetCDF files -------------------------------------
#
# One file per year from NOAA PSL. Files are saved locally and reused on subsequent runs, so this step only needs to run once in full.
# Total disk space required: approximately 42 files × 120 MB ≈ 5 GB.
# If you already ran the global pipeline, these files are already present and will be skipped.

cat("\n=== Downloading GODAS NetCDF files (1980–2021) ===\n")
cat("  Saving to:", GODAS_DIR, "\n\n")

for (yr in years) {

  dest <- file.path(GODAS_DIR, sprintf("pottmp.%d.nc", yr))

  # Skip if already downloaded and file size is plausible (> 50 MB)
  if (file.exists(dest) && file.size(dest) > 50e6) {
    cat(sprintf("  %d  already on disk (%.0f MB) — skipping\n",
                yr, file.size(dest) / 1e6))
    next
  }

  url <- sprintf(
    "https://downloads.psl.noaa.gov/Datasets/godas/pottmp.%d.nc", yr
  )

  cat(sprintf("  %d  downloading ... ", yr))

  tryCatch({
    # mode = "wb" is required for binary files on Windows
    download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    cat(sprintf("done (%.0f MB)\n", file.size(dest) / 1e6))
  }, error = function(e) {
    cat("FAILED:", e$message, "\n")
    if (file.exists(dest)) file.remove(dest)   # remove partial download
  })
}

# Record which years are available for SST extraction
downloaded_years <- years[sapply(years, function(yr) {
  dest <- file.path(GODAS_DIR, sprintf("pottmp.%d.nc", yr))
  file.exists(dest) && file.size(dest) > 50e6
})]

cat(sprintf("\n  Available GODAS files: %d / %d years\n",
            length(downloaded_years), length(years)))

failed_years <- setdiff(years, downloaded_years)
if (length(failed_years) > 0) {
  cat("  Missing years:", paste(failed_years, collapse = ", "), "\n")
  cat("  Records from these years will receive NA for SST.\n")
}


# --- Function: extract GODAS SST at one occurrence point ------------
# Defines a small spatial box around the occurrence point using the EXTRACT_RANGE constants
#
# GODAS longitude convention — GODAS stores longitude as 0–360°. Negative longitudes must be converted before the 
# spatial filter is applied (e.g. -5° W becomes 355° in GODAS coordinates). 
#
# Only the surface layer (depth <= 5 m)
#
# Annual mean 

extract_godas_sst <- function(lon_180, lat, yr,
                               range_lon = EXTRACT_RANGE_LON,
                               range_lat = EXTRACT_RANGE_LAT) {

  nc_path <- file.path(GODAS_DIR, sprintf("pottmp.%d.nc", yr))
  if (!file.exists(nc_path)) return(NA_real_)

  # Convert -180:180 longitude to 0:360 as used by GODAS
  lon_360 <- ifelse(lon_180 < 0, 360 + lon_180, lon_180)

  # Define extraction window in GODAS coordinates
  lon_min <- lon_360 - range_lon
  lon_max <- lon_360 + range_lon
  lat_min <- lat    - range_lat
  lat_max <- lat    + range_lat

  tryCatch({

    # Open NetCDF file
    nc <- tidync::tidync(nc_path)

    # hyper_filter() applies spatial and depth filters before reading any data into memory
    sub <- nc %>%
      tidync::hyper_filter(
        lon   = lon > lon_min & lon <= lon_max,
        lat   = lat > lat_min & lat <= lat_max,
        level = level <= 5     # surface layer only
      ) %>%
      tidync::hyper_tibble()   # read the small filtered subset into R

    if (nrow(sub) == 0) return(NA_real_)

    # Convert Kelvin to Celsius and average all months and cells
    mean(sub$pottmp - 273.15, na.rm = TRUE)

  }, error = function(e) NA_real_)
}

# --- Match GODAS SST to all occurrence records ----------------------
#
# Each occurrence record is matched to the GODAS annual mean SST 
# Records for which no valid ocean cell exists within the extraction  window receive NA and will be excluded from Ti calculations in Script 3.

match_sst_for_species <- function(sp_name, occs_list) {

  occs <- occs_list[[sp_name]]
  cat(sprintf("\n  Processing: %-25s  (%d records)\n", sp_name, nrow(occs)))

  # Separate records within and outside GODAS temporal coverage
  occs_valid   <- occs %>% filter(year %in% downloaded_years)
  occs_missing <- occs %>%
    filter(!year %in% downloaded_years) %>%
    mutate(godas_sst = NA_real_)

  n_skip <- nrow(occs) - nrow(occs_valid)
  if (n_skip > 0)
    cat(sprintf("    \u26a0 %d records outside GODAS year range \u2192 NA\n", n_skip))

  # Group records by year —> each annual file is opened once
  occs_by_year <- split(occs_valid, occs_valid$year)

  matched_by_year <- lapply(names(occs_by_year), function(yr_chr) {

    chunk <- occs_by_year[[yr_chr]]
    yr    <- as.integer(yr_chr)
    cat(sprintf("    Year %d: %d records\n", yr, nrow(chunk)))

    # Extract SST for each record in this year
    chunk$godas_sst <- sapply(seq_len(nrow(chunk)), function(i) {
      extract_godas_sst(
        lon_180 = chunk$decimalLongitude[i],
        lat     = chunk$decimalLatitude[i],
        yr      = yr
      )
    })

    return(chunk)
  })

  bind_rows(c(matched_by_year, list(occs_missing)))
}

all_occs_sst <- lapply(species_list, match_sst_for_species,
                        occs_list = all_occs)
names(all_occs_sst) <- species_list


# --- Summary diagnostics ---------------------------------------------

cat("\n=== SST matching summary ===\n")

for (sp in species_list) {
  df   <- all_occs_sst[[sp]]
  n_ok <- sum(!is.na(df$godas_sst))
  cat(sprintf(
    "  %-25s  total: %5d  |  with SST: %5d  |  missing: %5d\n",
    sp, nrow(df), n_ok, nrow(df) - n_ok
  ))
}


# --- Save output -----------------------------------------------------

saveRDS(
  list(
    occurrences_sst = all_occs_sst,
    species_list    = species_list,
    species_col_map = saved$species_col_map,
    MED_LON_MIN     = saved$MED_LON_MIN,
    MED_LON_MAX     = saved$MED_LON_MAX,
    MED_LAT_MIN     = saved$MED_LAT_MIN,
    MED_LAT_MAX     = saved$MED_LAT_MAX
  ),
  file.path(out_dir, "occurrences_med_sst.rds")
)

cat("\n=== Script 2 (Mediterranean) complete ===\n")
cat("Saved to: data/med/occurrences_med_sst.rds\n")
cat("Next: run 03_med_CTI.R\n")

