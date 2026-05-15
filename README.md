# From Macroalgal Forests to Barrens: High Herbivory Pressure Drives Regime Shifts in eastern Mediterranean Rocky Reef Ecosystems

**Athanasios Nikolaou¹, Deevesh Ashley Hemraj², Jacob Carstensen², George Tsirtsis¹, Stelios Katsanevakis¹**

¹ Department of Marine Sciences, University of the Aegean, University Hill 81100, Mytilene, Lesvos, Greece  
² Department of Ecoscience, Aarhus University, Frederiksborgvej 399, Roskilde, DK-4000 Denmark  

📧 [thanasisnik072@gmail.com](mailto:thanasisnik072@gmail.com)

---

## Overview

This repository contains all data and R scripts necessary to reproduce the analyses and figures reported in the associated manuscript. The study quantifies macroalgal forest state and herbivore community structure along a tropicalisation gradient spanning the North Aegean Sea to the Levantine Sea, develops a Cumulative Herbivory Index (CHI) integrating all herbivorous species, and examines regime shifts in macroalgal forests driven by cumulative herbivory pressure.

## Abstract

Macroalgal forests are key structural components of Mediterranean rocky reefs, forming three-dimensional habitats that support biodiversity, ecosystem functioning, and essential ecosystem services. However, these forests are in decline, often due to overgrazing, which can trigger regime shifts from productive forests to impoverished barren or turf-dominated states. Here, we quantify the state of macroalgal forests and herbivore communities along a tropicalisation gradient, spanning the colder North Aegean to the warmer, highly bio-invaded Levantine Sea. We developed a cumulative herbivory pressure index that incorporates all herbivorous species in the region, based on estimates of species-specific consumption rates, providing the first assessment of regime shifts in macroalgal forests driven by cumulative herbivory pressure rather than by individual species. Our results indicate high cumulative herbivore pressure across both regions, despite differences in herbivore community structure. Various positive feedback mechanisms appear to hinder macroalgal recovery across this gradient. Linking cumulative herbivory pressure to forest degradation has conservation and management implications and highlights how climate change and biological invasions are reshaping regime shift dynamics in temperate rocky reefs — insights relevant across the Mediterranean basin and beyond.

---

## Software and Package Requirements

All analyses were conducted in R. Install the following packages before running any script:

```r
install.packages(c(
  "ggplot2", "dplyr", "patchwork", "ggdist", "cowplot", "ggtext",
  "vegan", "changepoint", "strucchange", "mgcv", "gratia", "DHARMa",
  "here", "tidync", "sf", "rnaturalearth", "rnaturalearthdata",
  "robis", "rgbif"
))
```

---

## Repository Structure

The repository is organised into four folders, one per figure group. Scripts must be run in the numbered order within each folder.

```
├── Figure 2 Barplots/
│   ├── 01_barplots.R
│   └── dataset_eastMed_reefs_BARPLOTS.xlsx
│
├── Figure 3 NMDs, regression and changepoint/
│   ├── data/
│   │   └── nmds_fourth_root.rds
│   ├── 02a_nmds_fourth_root.R
│   ├── 02b_nmds_fourth_root_regressions_changepoint.R
│   └── dataset_eastMed_reefs_abundance_density.csv
│
├── Figure 4 CTI plots and estimates_med/
│   ├── data/
│   │   └── dataset_eastMed_reefs_abundance_density.csv
│   ├── 03a_med_fetch_occurrences.R
│   ├── 03b_med_godas_sst.R
│   ├── 03c_med.R
│   └── 03d_raincloud.R
│
└── Figure 5 changepoint and GAMs/
    ├── Changepoint/
    │   └── 04a_breakpoint_analysis.R
    ├── GAMs/
    │   ├── average/
    │   │   └── 04b_bined_gam_avg.R
    │   └── max/
    │       └── 04c_bined_gam_max.R
    ├── ch_data.xlsx
    ├── dataset_eastMed_reefs_BARPLOTS.xlsx
    └── dataset_eastMed_reefs_abundance_density.csv
```

---

## Datasets

### `dataset_eastMed_reefs_BARPLOTS.csv`

**Used by:** `01_barplots.R`

Site-level survey data for Figure 2. Each row is one survey site (n = 60 sites; 30 NAS, 30 SAL).

| Column | Units | Description |
|--------|-------|-------------|
| `id` | — | Site identifier (factor) |
| `region` | — | Region: NAS (North Aegean Sea) or SAL (South Aegean & Levantine Sea) |
| `canopy_avg` | % | Average macroalgal forest canopy cover across the 0–5 m depth zone |
| `pl` | ind/m² | Population density of *Paracentrotus lividus* |
| `al` | ind/m² | Population density of *Arbacia lixula* |
| `ss` | g/m² | Biomass of *Sarpa salpa* |
| `sl` | g/m² | Biomass of *Siganus luridus* |
| `sr` | g/m² | Biomass of *Siganus rivulatus* |
| `sc` | g/m² | Biomass of *Sparisoma cretense* |
| `CHI` | g consumed/day/m² | Cumulative Herbivory Index |

---

### `dataset_eastMed_reefs_abundance_density.csv`

**Used by:** `02a`, `02b`, `03c`, `03d`

Extended site-level dataset including temperature data and computed indices. Each row is one survey site (n = 60).

> **Note:** The three temperature variables (annual average, minimum, and maximum SST) are downloaded from the Copernicus Marine Service product `SST_MED_SST_L4_REP_OBSERVATIONS_010_021` (2016–2021). These data can be freely accessed at [https://data.marine.copernicus.eu](https://data.marine.copernicus.eu).

| Column | Units | Description |
|--------|-------|-------------|
| `id` | — | Site identifier (factor) |
| `region` | — | NAS or SAL (see above) |
| `canopy_avg` | % | Average macroalgal canopy cover (0–5 m depth zone) |
| `pl` | ind/m² | Population density of *Paracentrotus lividus* |
| `al` | ind/m² | Population density of *Arbacia lixula* |
| `ss` | ind/m² | Population density of *Sarpa salpa* |
| `sl` | ind/m² | Population density of *Siganus luridus* |
| `sr` | ind/m² | Population density of *Siganus rivulatus* |
| `sc` | ind/m² | Population density of *Sparisoma cretense* |
| `avg_temp` | °C | Annual average SST (Copernicus, 2016–2021 mean) |
| `min_temp` | °C | Annual minimum SST |
| `max_temp` | °C | Annual maximum SST |
| `CHI` | g consumed/day/m² | Cumulative Herbivory Index |
| `CTI_med` | °C | Community Temperature Index (Mediterranean occurrence data; added by `03c`) |

---

### `ch_data.xlsx`

**Used by:** `04b_bined_gam_avg.R` and `04c_bined_gam_max.R`

Binned data derived from the main survey dataset. Sites are grouped into bins of approximately equal size along the CHI gradient to reduce noise and avoid ties, as required by the breakpoint and GAM analyses. Each row represents one bin.

---

## Script Descriptions

### Figure 2 — `01_barplots.R`
**Purpose:** Produces bar plots of canopy cover, herbivore abundances, and CHI by site (Figure 2).  
**Input:** `dataset_eastMed_reefs_BARPLOTS.csv`

---

### Figure 3 — `02a_nmds_fourth_root.R`
**Purpose:** Runs NMDS ordination on fourth-root transformed herbivore data and saves the NMDS object.  
**Input:** `dataset_eastMed_reefs_abundance_density.csv` (place in `data/` subfolder)

### Figure 3 — `02b_nmds_fourth_root_regressions_changepoint.R`
**Purpose:** Regresses NMDS1 scores against temperature variables and detects changepoints.  
**Input:** `data/nmds_fourth_root.rds` + `dataset_eastMed_reefs_abundance_density.csv`

---

### Figure 4 — `03a_med_fetch_occurrences.R`
**Purpose:** Downloads and cleans species occurrence records from OBIS and GBIF for six herbivore species.  
**Input:** None (downloads from APIs). Requires internet connection and GBIF credentials.

### Figure 4 — `03b_med_godas_sst.R`
**Purpose:** Assigns annual mean sea surface temperature to each occurrence record from GODAS.  
**Input:** `data/med/occurrences_med_raw.rds` + GODAS NetCDF files (downloaded automatically)

### Figure 4 — `03c_med_CTI_corr.R`
**Purpose:** Computes species thermal preference (Ti) and Community Temperature Index (CTI) for all survey sites.  
**Input:** `data/med/occurrences_med_sst.rds` + `data/dataset_eastMed_reefs_abundance_density.csv`

### Figure 4 — `03d_raincloud.R`
**Purpose:** Compares CTI and CHI between NAS and SAL using Wilcoxon tests and raincloud plots.  
**Input:** `data/med/dataset_with_CTI_med.csv` (output of `03c` — not the original CSV)

---

### Figure 5 — `04a_breakpoint_analysis.R`
**Purpose:** Detects structural breakpoints in the canopy cover vs CHI relationship.  
**Input:** `ch_data_avg.csv` and `ch_data_max.csv` (loaded via clipboard or file path)

### Figure 5 — `04b_bined_gam_avg.R` and `04c_bined_gam_max.R`
**Purpose:** Fits GAMs to binned canopy cover data and identifies significant regions of change using SiZer.  
**Input:** `ch_data_avg.csv` (`04b`) or `ch_data_max.csv` (`04c`), loaded via clipboard or file path

---

## Data and Code Availability


All R scripts and datasets required to reproduce the analyses and figures are available in the associated Zenodo repository at: **[DOI: 10.5281/zenodo.20187921]**. 

---

> For questions regarding the data or scripts, please contact the corresponding author via the email address provided in the associated manuscript.
