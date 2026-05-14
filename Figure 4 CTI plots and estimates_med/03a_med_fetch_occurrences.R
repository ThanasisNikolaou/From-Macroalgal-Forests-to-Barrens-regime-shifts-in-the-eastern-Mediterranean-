# Script 1 of 3 — Mediterranean CTI pipeline
# Fetch species occurrence records from OBIS and GBIF
#
# This script downloads occurrence data for all herbivore species, restricts records to the Mediterranean bounding box, 
# filters to valid coordinates and years within the GODAS SST period (1980–2021), merges the two databases, removes duplicate
# observations, and saves the cleaned dataset for use in Script 2.
# Records without a year are dropped because they cannot be matched to a year-specific SST value.
#
# Mediterranean bounding box:
#   Longitude: -18.1° to 36.3° E
#   Latitude :  31.0° to 46.0° N
#
# Output: data/med/occurrences_med_raw.rds

# --- Libraries -------------------------------------------------------

library(robis)   
library(rgbif)   
library(dplyr)   
library(here)

# --- Settings --------------------------------------------------------

# Species names must match the dataset
species_list <- c(
  "Paracentrotus lividus",
  "Arbacia lixula",
  "Sarpa salpa",
  "Siganus luridus",
  "Siganus rivulatus",
  "Sparisoma cretense"
)

# Map each species name to its abundance column in the CSV
species_col_map <- c(
  "Paracentrotus lividus" = "pl",
  "Arbacia lixula"        = "al",
  "Sarpa salpa"           = "ss",
  "Siganus luridus"       = "sl",
  "Siganus rivulatus"     = "sr",
  "Sparisoma cretense"    = "sc"
)

# Temporal bounds of the GODAS SST product 
YEAR_MIN <- 1980
YEAR_MAX <- 2021

# Mediterranean bounding box (all occurrence records outside this area are excluded before any further processing)
MED_LON_MIN <- -18.1
MED_LON_MAX <-  36.3
MED_LAT_MIN <-  31.0
MED_LAT_MAX <-  46.0

# Output directory
out_dir <- here("data", "med")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- GBIF credentials ------------------------------------------------
#
# GBIF requires authentication for bulk downloads. Store your
# credentials as environment variables so they are never hard-coded:
#   usethis::edit_r_environ()
# Add the following three lines, then restart R:
#   GBIF_USER=your_username
#   GBIF_PWD=your_password
#   GBIF_EMAIL=your_email

if (Sys.getenv("GBIF_USER") == "") {
  stop(paste(
    "GBIF credentials not found.\n",
    "Run usethis::edit_r_environ(), add your credentials, restart R."
  ))
}


# --- standardise columns from OBIS or GBIF ------------------
#
# OBIS and GBIF use different column names and data types. 
# This function harmonises bot into a single consisntent four-column format
# species / long / lat / year / source

standardise_occs <- function(df, source_label, sp_name) {

  # Lowercase all column names to handle any OBIS/GBIF inconsistencies
  names(df) <- tolower(names(df))

  df %>%
    transmute(
      species          = sp_name,                        # canonical species name
      decimalLongitude = as.numeric(decimallongitude),   # longitude
      decimalLatitude  = as.numeric(decimallatitude),    # latitude
      year             = as.integer(year),               # observation year
      source           = source_label                    # "OBIS" or "GBIF"
    ) %>%
    # Remove records with missing coordinates or year,
    # records outside the GODAS temporal coverage,
    # and records outside the Mediterranean bounding box
    filter(
      !is.na(decimalLongitude),
      !is.na(decimalLatitude),
      !is.na(year),
      year >= YEAR_MIN,
      year <= YEAR_MAX,
      decimalLongitude >= MED_LON_MIN,
      decimalLongitude <= MED_LON_MAX,
      decimalLatitude  >= MED_LAT_MIN,
      decimalLatitude  <= MED_LAT_MAX,
      # Exclude Atlantic records west of Gibraltar Strait
      !(decimalLongitude < -6),
      # Exclude Bay of Biscay records
      !(decimalLongitude < 0 & decimalLatitude > 42)
    )
}


# --- Fetch OBIS occurrences (Mediterranean bounding box) -------------
#
# the bounding box is passed directly to the OBIS API via the geometry argument, which reduces the volume of data transferred. 

cat("=== Fetching OBIS occurrences (Mediterranean) ===\n")

obis_raw <- lapply(species_list, function(sp) {
  cat("  Querying OBIS:", sp, "\n")
  tryCatch(
    robis::occurrence(
      scientificname = sp,
      geometry = sprintf(
        "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
        MED_LON_MIN, MED_LAT_MIN,
        MED_LON_MAX, MED_LAT_MIN,
        MED_LON_MAX, MED_LAT_MAX,
        MED_LON_MIN, MED_LAT_MAX,
        MED_LON_MIN, MED_LAT_MIN # close WKT polygon
      )
    ),
    error = function(e) {
      warning(paste("OBIS failed for", sp, ":", e$message))
      return(NULL)
    }
  )
})
names(obis_raw) <- species_list

# Standardise columns and apply year + bounding box filter
obis_std <- lapply(species_list, function(sp) {
  if (is.null(obis_raw[[sp]])) return(NULL)
  standardise_occs(obis_raw[[sp]], "OBIS", sp)
})
names(obis_std) <- species_list

cat("\n  OBIS record counts (after year and bounding box filter):\n")
print(sapply(obis_std, function(x) if (is.null(x)) 0L else nrow(x)))

# how many records are lost to the year filte?
cat("\n  Records lost to year filter (NA year or outside 1980-2021):\n")
for (sp in species_list) {
  raw  <- if (is.null(obis_raw[[sp]])) 0L else nrow(obis_raw[[sp]])
  kept <- if (is.null(obis_std[[sp]])) 0L else nrow(obis_std[[sp]])
  lost <- raw - kept
  pct  <- if (raw > 0) round(100 * lost / raw, 1) else 0
  cat(sprintf("  %-25s  raw: %6d  |  kept: %5d  |  lost: %5d (%.1f%%)\n",
              sp, raw, kept, lost, pct))
}


# --- Fetch GBIF occurrences via the download API ---------------------
#
# The Mediterranean bounding box is passed as a server-side predicate alongside the year filter,
# so only relevant records are downloaded. 

cat("\n=== Fetching GBIF taxon keys ===\n")

# Each species needs a GBIF numeric taxon key for filtering
gbif_keys <- sapply(species_list, function(sp) {
  res <- rgbif::name_backbone(name = sp, rank = "species")
  if (!is.null(res$usageKey)) {
    cat(sprintf("  %-25s  key: %s\n", sp, res$usageKey))
    return(res$usageKey)
  } else {
    warning(paste("No GBIF key found for", sp))
    return(NA_integer_)
  }
})

cat("\n=== Submitting GBIF downloads (one at a time) ===\n")

gbif_data <- lapply(species_list, function(sp) {

  key <- gbif_keys[sp]
  if (is.na(key)) {
    warning(paste("Skipping", sp, "— no GBIF key"))
    return(NULL)
  }

  cat("\n  Submitting download for:", sp, "\n")

  dl <- rgbif::occ_download(
    rgbif::pred("taxonKey",           key),
    rgbif::pred("hasCoordinate",      TRUE),       # must have coordinates
    rgbif::pred("hasGeospatialIssue", FALSE),      # exclude flagged records
    rgbif::pred_gte("year",           YEAR_MIN),   # within GODAS period
    rgbif::pred_lte("year",           YEAR_MAX),
    rgbif::pred_within(sprintf(       # Mediterranean bounding box
      "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
      MED_LON_MIN, MED_LAT_MIN,
      MED_LON_MAX, MED_LAT_MIN,
      MED_LON_MAX, MED_LAT_MAX,
      MED_LON_MIN, MED_LAT_MAX,
      MED_LON_MIN, MED_LAT_MIN
    )),
    format = "SIMPLE_CSV"
  )

  # Poll GBIF every 30 seconds until the download is ready
  cat("  Waiting for GBIF to prepare download:", sp, "\n")
  repeat {
    status <- rgbif::occ_download_meta(dl)$status
    cat(sprintf("    [%s] %s\n", format(Sys.time(), "%H:%M:%S"), status))
    if (status == "SUCCEEDED") break
    if (status == "KILLED")
      stop(paste("GBIF download killed for", sp))
    Sys.sleep(30)
  }

  dat <- rgbif::occ_download_get(dl) %>%
    rgbif::occ_download_import()

  cat(sprintf("  %-25s  records: %d\n", sp, nrow(dat)))
  return(dat)
})
names(gbif_data) <- species_list

# Standardise GBIF columns and apply year + bounding box filter
gbif_std <- lapply(species_list, function(sp) {
  if (is.null(gbif_data[[sp]])) return(NULL)
  standardise_occs(gbif_data[[sp]], "GBIF", sp)
})
names(gbif_std) <- species_list

cat("\n  GBIF record counts (after year and bounding box filter):\n")
print(sapply(gbif_std, function(x) if (is.null(x)) 0L else nrow(x)))


# --- Pool OBIS + GBIF, label source, and deduplicate ----------------
#
# Records sharing the same decimal coordinates and observation year are treated as duplicates and reduced to a single entry
#   "OBIS" — found only in OBIS
#   "GBIF" — found only in GBIF
#   "BOTH" — the same observation was found to both databases

cat("\n=== Pooling OBIS + GBIF, flagging sources, deduplicating ===\n")

deduplicate_pooled <- function(sp) {

  obis_df <- obis_std[[sp]]
  gbif_df <- gbif_std[[sp]]

  # Combine all records from both databases into one data frame
  combined <- bind_rows(obis_df, gbif_df)

  if (is.null(combined) || nrow(combined) == 0) {
    warning(paste("No records for", sp))
    return(NULL)
  }

  n_before <- nrow(combined)

  # For each unique (lon, lat, year) group, determine which databases contributed and assign the appropriate source label
  coord_sources <- combined %>%
    group_by(decimalLongitude, decimalLatitude, year) %>%
    summarise(
      sources = paste(sort(unique(source)), collapse = "+"),
      .groups = "drop"
    ) %>%
    mutate(source_final = case_when(
      sources == "GBIF" ~ "GBIF",   # only in GBIF
      sources == "OBIS" ~ "OBIS",   # only in OBIS
      TRUE              ~ "BOTH"    # same observation in both databases
    ))

  # Keep one random record per (lon, lat, year) combination
  deduped <- combined %>%
    group_by(decimalLongitude, decimalLatitude, year) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-source) %>%
    left_join(
      coord_sources %>%
        select(decimalLongitude, decimalLatitude, year, source_final),
      by = c("decimalLongitude", "decimalLatitude", "year")
    ) %>%
    rename(source = source_final)

  n_after <- nrow(deduped)

  cat(sprintf(
    "  %-25s  pooled: %5d  |  deduped: %5d  |  removed: %5d\n",
    sp, n_before, n_after, n_before - n_after
  ))

  src_tab <- table(deduped$source)
  cat("    Sources: ")
  cat(paste(names(src_tab), src_tab, sep = "=", collapse = "  "), "\n")

  return(deduped)
}

all_occs <- lapply(species_list, deduplicate_pooled)
names(all_occs) <- species_list


# --- Save output -----------------------------------------------------

saveRDS(
  list(
    occurrences     = all_occs,
    species_list    = species_list,
    species_col_map = species_col_map,
    YEAR_MIN        = YEAR_MIN,
    YEAR_MAX        = YEAR_MAX,
    MED_LON_MIN     = MED_LON_MIN,
    MED_LON_MAX     = MED_LON_MAX,
    MED_LAT_MIN     = MED_LAT_MIN,
    MED_LAT_MAX     = MED_LAT_MAX
  ),
  file.path(out_dir, "occurrences_med_raw.rds")
)

cat("\n=== Script 1 (Mediterranean) complete ===\n")
cat("Saved to: data/med/occurrences_med_raw.rds\n")
cat("Next: run 02_med_godas_sst.R\n")
