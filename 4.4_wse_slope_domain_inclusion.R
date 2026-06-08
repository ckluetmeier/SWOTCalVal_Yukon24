# =============================================================================
# SWOT WSE & Slope Domain Inclusion
# -----------------------------------------------------------------------------
# Determines, for each node/reach in the Yukon River domain, whether it
# was observed by Version C, Version D, both, or neither SWOT product, and
# writes the result as shapefiles for mapping. Also produces per-river domain
# statistics and the in situ summary shapefiles + example timeseries
# used in Figures 2 and 3.
# SWOT processing versions:
#   - Version C / PIC0 (SWORD v16, RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# SWOT and in situ measurements are matched in time/space in earlier scripts
# and all data are harmonized to SWORD v17b node / reach IDs before analysis.
# Version inclusion code (per node / reach):
#   -1 = observed only in Version C (RiverSP v16 / PIC0)
#    0 = observed in both versions
#    1 = observed only in Version D (v17b / PGD0)
#    2 = in the YR domain but not observed in either SWOT version
#
# Contains:
#   - Tables: 1
#   - Figures: 2, 3
#   - Shapefiles: node WSE inclusion, reach WSE inclusion, reach slope
#                 inclusion, PT summary, GNSS summary
# =============================================================================

library(sf)
library(dplyr)
library(tidyverse)
library(lubridate)


# =============================================================================
# NODE-LEVEL WSE DOMAIN INCLUSION
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Read in node-level WSE data (time/space matched SWOT vs. in situ)
# -----------------------------------------------------------------------------

# PT — Version C (SWORD v16 / RiverSP) and Version D (SWORD v17b / RiverTile)
node_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  filter(dark_frac < 0.5)

# GNSS — Version C and Version D
node_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 2. Build the in situ node domain (union of PT and GNSS observations)
# -----------------------------------------------------------------------------

# Combine the munged PT files (per-river, per-PT cluster) into one frame.
# SWORD v17b
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_YR"
setwd(wd)

# List all PT node CSVs (toolbox output) and read them in
csv_files <- list.files(wd, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  return(df)
})
combined_PT_df <- bind_rows(data_list)

# Parse the PT timestamp column to POSIXct (UTC)
combined_PT_df$pt_time_UTC <- as.POSIXct(
  combined_PT_df$pt_time_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
)

# GNSS drift node WSEs (SWORD v17b)
GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv"
)

# YR node domain = unique nodes seen by either PT or GNSS
YR_domain <- bind_rows(combined_PT_df, GNSS_df) %>%
  distinct(node_id)


# -----------------------------------------------------------------------------
# 3. Harmonize to SWORD v17b and compute per-node version inclusion
# -----------------------------------------------------------------------------

# Translator: maps v16 node IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv"
)

# Apply translation to Version C node data
node_SWOT_PT_vC <- node_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")
  ) %>%
  rename(node_id = v17_node_id)

node_SWOT_GNSS_vC <- node_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")
  ) %>%
  rename(node_id = v17_node_id)

# Merge all four sources together
node_SWOT_full_insitu <- bind_rows(
  node_SWOT_PT_vC, node_SWOT_PT_vD,
  node_SWOT_GNSS_vC, node_SWOT_GNSS_vD
)

# Compute version inclusion flag per node:
#   -1 = observed only in Version C
#    0 = observed in both versions
#    1 = observed only in Version D
all_nodes <- node_SWOT_full_insitu %>%
  distinct(node_id, source) %>%                  # one row per node_id × source
  group_by(node_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile  ~  0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~  1L,
      TRUE                         ~ NA_integer_
    )
  ) %>%
  select(node_id, version_inclusion)

# Join inclusion flag to the full YR node domain.
# Nodes in the domain but missing from both SWOT versions get version_inclusion = 2.
all_YR_domain_nodes <- YR_domain %>%
  left_join(all_nodes, by = "node_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))


# -----------------------------------------------------------------------------
# 4. Save node inclusion as a shapefile
# -----------------------------------------------------------------------------

# Bring in SWORD shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

# Attach version_inclusion to SWORD geometry
all_YR_domain_nodes_sf <- sword_sf %>%
  left_join(all_YR_domain_nodes, by = "node_id")

# Subset to nodes with a version_inclusion value
all_YR_domain_nodes_sf_subset <- all_YR_domain_nodes_sf %>%
  filter(!is.na(version_inclusion))

# Save shapefile
# NOTE: bad nodes (e.g. no matched overpass) were not removed from the field
# data, so they were manually deleted in QGIS.
st_write(
  all_YR_domain_nodes_sf_subset,
  "/Users/camryn/Desktop/all_YR_domain_nodes_subset.shp",
  delete_layer = TRUE
)


# -----------------------------------------------------------------------------
# 4b. Width version inclusion shapefile
# -----------------------------------------------------------------------------
#
# sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")
# node_SWOT_ortho <- sword_sf %>%
#   left_join(node_SWOT_ortho, by = "node_id")
# node_SWOT_ortho_subset <- node_SWOT_ortho %>%
#   filter(!is.na(version_inclusion))
# st_write(
#   node_SWOT_ortho_subset,
#   "/Users/camryn/Desktop/all_YR_domain_nodes_width_subset.shp",
#   delete_layer = TRUE
# )


# =============================================================================
# REACH-LEVEL WSE DOMAIN INCLUSION
# =============================================================================


# -----------------------------------------------------------------------------
# 5. Read in reach-level WSE data (time/space matched SWOT vs. in situ)
# -----------------------------------------------------------------------------

# PT — Version C and Version D
reach_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile") %>%
  filter(dark_frac < 0.5)

# GNSS — Version C and Version D
reach_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 6. Build the in situ reach domain (union of PT and GNSS reaches)
# -----------------------------------------------------------------------------

# Combine the munged PT reach files
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b"
setwd(wd)

csv_files <- list.files(wd, pattern = "^YR_812.*\\.csv$", full.names = TRUE)
data_list <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  return(df)
})
# Coerce sorted_nodelist to character so bind_rows does not coerce it to NA
data_list <- lapply(data_list, function(df) {
  df %>% mutate(sorted_nodelist = as.character(sorted_nodelist))
})
combined_PT_df <- bind_rows(data_list)

# Parse PT timestamp to POSIXct
combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, tz = "UTC")

# GNSS reach WSE & slope (SWORD v17b)
GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv"
)

# YR reach domain
YR_domain <- bind_rows(combined_PT_df, GNSS_df) %>%
  distinct(reach_id)


# -----------------------------------------------------------------------------
# 7. Harmonize to SWORD v17b and compute per-reach version inclusion
# -----------------------------------------------------------------------------

# Translator: maps v16 reach IDs to v17b
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

# Merge all reach sources together, coalescing in situ WSE from PT or GNSS
reach_SWOT_full_insitu <- bind_rows(
  reach_SWOT_PT_vC, reach_SWOT_PT_vD,
  reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD
) %>%
  mutate(insitu_wse_m        = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m))

# Per-reach version inclusion flag (same coding as nodes)
all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source) %>%
  group_by(reach_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile  ~  0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~  1L,
      TRUE                         ~ NA_integer_
    )
  ) %>%
  select(reach_id, version_inclusion)

# Join to YR reach domain; reaches with no SWOT match get version_inclusion = 2
all_YR_domain_reaches <- YR_domain %>%
  left_join(all_reaches, by = "reach_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))


# -----------------------------------------------------------------------------
# 8. Save reach WSE inclusion as a shapefile
# -----------------------------------------------------------------------------

sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

all_YR_domain_reaches_sf <- sword_sf %>%
  left_join(all_YR_domain_reaches, by = "reach_id")

all_YR_domain_reaches_sf_subset <- all_YR_domain_reaches_sf %>%
  filter(!is.na(version_inclusion))

st_write(
  all_YR_domain_reaches_sf_subset,
  "/Users/camryn/Desktop/all_YR_domain_reaches_subset.shp",
  delete_layer = TRUE
)

# Colors used in the corresponding plot:
#   #e97132  (Version C only)
#   #ececec  (both)
#   #00008b  (Version D only)
#   dashed gray 2 (no SWOT match)


# =============================================================================
# RIVER DOMAIN STATISTICS (TABLE 1)
# =============================================================================


# -----------------------------------------------------------------------------
# 9. Per-river length, median width, median slope (from SWORD attributes)
# -----------------------------------------------------------------------------

sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

# Subset SWORD to the YR domain reaches
sword_subset <- sword_sf %>%
  filter(reach_id %in% YR_domain$reach_id) %>%
  distinct(reach_id, .keep_all = TRUE)

cat("Total km of river:", sum(sword_subset$reach_len) / 1000, "km\n")

# Tag each reach with its river using the SWORD reach_id prefix.
# Some reaches need manual assignment (SJ, BL).
sword_subset <- sword_subset %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                      "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",  # until the Circle bifurcation
      river_code == "812509" ~ "lowerYR",  # past the PR confluence
      river_code == "812705" ~ "upperYR",  # Circle up
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

# Per-river totals from SWORD
river_stats <- sword_subset %>%
  group_by(river) %>%
  summarise(
    total_km         = sum(reach_len, na.rm = TRUE) / 1000,
    median_width     = median(width),
    median_slope     = median(slope) * 100,
    n_distinct(reach_id)
  )


# -----------------------------------------------------------------------------
# 10. GNSS per-river drift count
# -----------------------------------------------------------------------------

# Tag each GNSS reach with its river (same case_when as above)
GNSS_df <- GNSS_df %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                      "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",
      river_code == "812509" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

GNSS_stats <- GNSS_df %>%
  group_by(river) %>%
  summarise(
    num_drift_reaches = sum(!is.na(reach_id))
  )


# =============================================================================
# REACH-LEVEL SLOPE DOMAIN INCLUSION
# =============================================================================


# -----------------------------------------------------------------------------
# 11. Read in reach-level slope data
# -----------------------------------------------------------------------------

# PT — Version C and Version D
reach_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_PT.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# GNSS — Version C and Version D
reach_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 12. Build the in situ slope reach domain (PT + GNSS)
# -----------------------------------------------------------------------------

wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b"
setwd(wd)

csv_files <- list.files(wd, pattern = "^YR_812.*\\.csv$", full.names = TRUE)
data_list <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  return(df)
})
data_list <- lapply(data_list, function(df) {
  df %>% mutate(sorted_nodelist = as.character(sorted_nodelist))
})
combined_PT_df <- bind_rows(data_list)
combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, tz = "UTC")

GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv"
)

YR_domain <- bind_rows(combined_PT_df, GNSS_df) %>%
  distinct(reach_id)


# -----------------------------------------------------------------------------
# 13. Harmonize to SWORD v17b and compute per-reach slope version inclusion
# -----------------------------------------------------------------------------

SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

# Merge all four sources; coalesce in situ slope from PT or GNSS
reach_SWOT_full_insitu <- bind_rows(
  reach_SWOT_PT_vC, reach_SWOT_PT_vD,
  reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD
) %>%
  mutate(insitu_slope_m_m        = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs)) %>%
  mutate(insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias))

all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source) %>%
  group_by(reach_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile  ~  0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~  1L,
      TRUE                         ~ NA_integer_
    )
  ) %>%
  select(reach_id, version_inclusion)

all_YR_domain_reaches <- YR_domain %>%
  left_join(all_reaches, by = "reach_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))


# -----------------------------------------------------------------------------
# 14. Save reach slope inclusion as a shapefile
# -----------------------------------------------------------------------------

sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

all_YR_domain_reaches_sf <- sword_sf %>%
  left_join(all_YR_domain_reaches, by = "reach_id")

all_YR_domain_reaches_sf_subset <- all_YR_domain_reaches_sf %>%
  filter(!is.na(version_inclusion))

st_write(
  all_YR_domain_reaches_sf_subset,
  "/Users/camryn/Desktop/all_YR_domain_reaches_slope_subset.shp",
  delete_layer = TRUE
)

# Colors used in the corresponding plot:
#   #e97132  (Version C only)
#   #ececec  (both)
#   #00008b  (Version D only)
#   dashed gray 2 (no SWOT match)


# =============================================================================
# FIGURE 2 — PT & GNSS SUMMARY SHAPEFILES
# =============================================================================


# -----------------------------------------------------------------------------
# 15. PT summary stats (per PT serial): obs count, deployment days, location
# -----------------------------------------------------------------------------

# Merge all PT node CSVs
base_dir <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b"

csv_files <- list.files(path = base_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

combined_PT_df <- map_dfr(csv_files, read_csv, show_col_types = FALSE)

# Save the merged PT node file
write_csv(combined_PT_df, file.path(base_dir, "flyby_SWOTCalVal_YR_PT_L1_v17b.csv"))

# SWOT node timeseries (RiverSP v17b / PGD0), filtered to the PT node set
SWOT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv"
)

SWOT_df_filtered <- SWOT_df %>%
  semi_join(combined_PT_df, by = c("node_id" = "Node_ID"))

# time_tai is seconds since 2000-01-01, offset 37 seconds from UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# Match SWOT to PT in time (±7.5 min window)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_df_filtered %>%
        filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Match SWOT to PT in space (node level)
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(Node_ID == node_id)

# Drop any duplicated PT observations
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[
  !duplicated(time_space_matched_SWOT_PT[c("pt_wse_m", "pt_time_UTC", "pt_serial")]), ]

# Per-PT summary stats
PT_summary_stats <- time_space_matched_SWOT_PT %>%
  mutate(
    pt_install_UTC   = ymd_hms(pt_install_UTC),
    pt_uninstall_UTC = ymd_hms(pt_uninstall_UTC)) %>%
  group_by(pt_serial, pt_install_UTC, pt_uninstall_UTC) %>%
  summarise(
    n_obs = n(),
    obs_days = as.numeric(
      difftime(first(pt_uninstall_UTC), first(pt_install_UTC), units = "days")),
    avg_pt_lat = mean(pt_lat, na.rm = TRUE),
    avg_pt_lon = mean(pt_lon, na.rm = TRUE),
    Node_ID    = first(Node_ID),
    Reach_ID   = first(Reach_ID),
    .groups = "drop") %>%
  group_by(pt_serial) %>%
  summarise(
    n_obs = sum(n_obs, na.rm = TRUE),
    obs_days = sum(obs_days, na.rm = TRUE),
    pt_install_UTC   = min(pt_install_UTC, na.rm = TRUE),
    pt_uninstall_UTC = max(pt_uninstall_UTC, na.rm = TRUE),
    avg_pt_lat = mean(avg_pt_lat, na.rm = TRUE),
    avg_pt_lon = mean(avg_pt_lon, na.rm = TRUE),
    Node_ID    = first(Node_ID),
    Reach_ID   = first(Reach_ID),
    .groups = "drop")

# Add river name labels to data frame
PT_summary_stats <- PT_summary_stats %>%
  mutate(
    river_code = substr(Reach_ID, 1, 6),
    river = case_when(
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      Reach_ID %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      Reach_ID %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                      "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",
      river_code == "812509" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

river_order <- c("upperYR", "lowerYR", "PR", "CD", "SJ", "CL")

PT_summary_stats <- PT_summary_stats %>%
  mutate(river = factor(river, levels = river_order)) %>%
  arrange(river)

PT_summary_by_river <- PT_summary_stats %>%
  group_by(river) %>%
  summarise(
    earliest_pt_install_UTC   = min(pt_install_UTC, na.rm = TRUE),
    latest_pt_install_UTC     = max(pt_install_UTC, na.rm = TRUE),
    earliest_pt_uninstall_UTC = min(pt_uninstall_UTC, na.rm = TRUE),
    latest_pt_uninstall_UTC   = max(pt_uninstall_UTC, na.rm = TRUE),
    mean_obs_days   = mean(obs_days, na.rm = TRUE),
    median_obs_days = median(obs_days, na.rm = TRUE),
    mean_n_obs      = mean(n_obs, na.rm = TRUE),
    median_n_obs    = median(n_obs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(river)

# NOTE PT 2159244 ERRONEOUSLY HAS 7/31 LISTED AS FIRST INSTALL TIME 
# THIS IS BECAUSE FIRST TIME IS PARSED AS NA
# should be 2024-07-08T00:00:00.000000Z to 2024-07-26T18:50:00.000000Z
# and then 2024-07-31T00:25:00.000000Z to 2024-08-21T15:35:00.000000Z
PT_summary_stats$obs_days[abs(PT_summary_stats$obs_days - 21.63194) < 1e-5] <- 40.416667
PT_summary_stats$obs_days <- as.numeric(PT_summary_stats$obs_days)


# Convert to sf (WGS84) and save shapefile + CSV
PT_summary_sf <- st_as_sf(
  PT_summary_stats,
  coords = c("avg_pt_lon", "avg_pt_lat"),
  crs    = 4326
)


st_write(
  PT_summary_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/PT_summary_stats.shp",
  delete_dsn = TRUE
)

write.csv(
  PT_summary_stats,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/PT_summary_stats.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 16. GNSS summary stats per reach and node
# -----------------------------------------------------------------------------

# REACH
# ---------------------------

GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv")

GNSS_summary_stats <- GNSS_df %>%
  group_by(reach_id) %>%
  summarise(
    n_observations = n(),
    .groups        = "drop")

# Attach to SWORD geometry, keep only reaches with GNSS observations
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

GNSS_summary_stats_sf <- sword_sf %>%
  left_join(GNSS_summary_stats, by = "reach_id") %>%
  filter(!is.na(n_observations))

st_write(
  GNSS_summary_stats_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/GNSS_summary_stats_sf.shp",
  delete_layer = TRUE
)

# NODE
# ---------------------------

GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv")

GNSS_summary_stats <- GNSS_df %>%
  group_by(node_id) %>%
  summarise(
    n_observations = n(),
    .groups        = "drop"
  )

# Attach to SWORD geometry, keep only reaches with GNSS observations
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

GNSS_summary_stats_sf <- sword_sf %>%
  left_join(GNSS_summary_stats, by = "node_id") %>%
  filter(!is.na(n_observations))

st_write(
  GNSS_summary_stats_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/GNSS_summary_stats_node_sf.shp",
  delete_layer = TRUE
)

# =============================================================================
# Table S4 — GNSS METADATA
# =============================================================================

GNSS_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv") %>%
  filter(time_UTC > as.POSIXct("2024-01-01 00:00:00", tz = "UTC"))

# Attach to SWORD geometry, keep only reaches with GNSS observations
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

GNSS_sf <- sword_sf %>%
  right_join(GNSS_df, by = "node_id")

# Add river name labels to data frame
GNSS_sf <- GNSS_sf %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                      "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",
      river_code == "812509" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )


GNSS_summary_stats <- GNSS_sf %>%
  group_by(drift_id) %>%
  summarise(
    river_list       = paste(unique(river), collapse = ", "),
    start_time       = min(time_UTC, na.rm = TRUE),
    end_time         = max(time_UTC, na.rm = TRUE),
    survey_length_km = round(sum(node_len) / 1000, 2),
    reach_list       = paste(unique(reach_id), collapse = ", "),
    .groups          = "drop"
  ) %>%
  mutate(
    survey_time_hours = round(as.numeric(difftime(end_time, start_time, units = "hours")), 2)
  ) %>%
  arrange(start_time) %>%
  filter(survey_time_hours > 0.000000) %>%
  st_drop_geometry() %>%
  select(-drift_id) %>%
  relocate(survey_time_hours, .after = end_time)

# total km of GNSS data collected
sum(GNSS_summary_stats$survey_length_km, na.rm = TRUE)

# mean survey time in hours
mean(GNSS_summary_stats$survey_time_hours, na.rm = TRUE)
max(GNSS_summary_stats$survey_time_hours, na.rm = TRUE)

# unique number of days we have GNSS data from
map2(
  as.Date(GNSS_summary_stats$start_time),
  as.Date(GNSS_summary_stats$end_time),
  ~ seq(.x, .y, by = "day")
) %>%
  unlist() %>%
  as.Date(origin = "1970-01-01") %>%
  n_distinct()

write.csv(
  GNSS_summary_stats,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/GNSS_summary_stats.csv",
  row.names = FALSE
)


# =============================================================================
# FIGURE 3 — EXAMPLE WSE & WIDTH TIMESERIES
# =============================================================================


# -----------------------------------------------------------------------------
# 17. GNSS vs. SWOT longitudinal WSE (Porcupine River, 2024-08-20)
# -----------------------------------------------------------------------------

GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/3_SWOT_examples/data/GNSS_2024-08-20.csv"
)

# Subset to a single PR section between two p_dist_out bookends
GNSS_PR <- GNSS_df %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139)

GNSS_PR_oneday <- GNSS_PR %>%
  filter(as.Date(time_UTC) == as.Date("2024-08-20"))

# Plot: SWOT WSE colored by elevation + GNSS-derived drift WSE with error band
ggplot() +
  # geom_errorbar(...)  # WSE uncertainty bars, currently disabled
  geom_point(data = GNSS_PR_oneday,
             aes(x = p_dist_out / 1000, y = wse, color = wse), size = 4) +
  scale_color_gradient(low = "#2474b7", high = "#d3e3f3") +
  geom_point(data = GNSS_PR_oneday,
             aes(x = p_dist_out / 1000, y = mean_node_drift_wse_no_bias_m),
             color = "#CC79A7", size = 1) +
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000, y = mean_node_drift_wse_no_bias_m, group = 1),
            color = "#CC79A7", linewidth = 0.7) +
  # Upper / lower drift-error envelope
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000,
                y = mean_node_drift_wse_no_bias_m + node_total_error_m, group = 1),
            color = "#CC79A7", alpha = 0.3, linewidth = 2) +
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000,
                y = mean_node_drift_wse_no_bias_m - node_total_error_m, group = 1),
            color = "#CC79A7", alpha = 0.3, linewidth = 2) +
  guides(color = "none") +
  xlab("Distance to river outlet (km)") +
  ylab("WSE (m)") +
  xlim(2138, 2152) +
  theme_minimal(base_size = 30)


# -----------------------------------------------------------------------------
# 18. PT vs. SWOT WSE timeseries (single node)
# -----------------------------------------------------------------------------

# Candidate nodes: 81260300160901, 81260300170011
PT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_PR/flyby_SWOTCalVal_YR_PT_L1_PT230_20240704T120000_20240821T220000_20250714T183224_SWOTCalVal_YR_KEY_20240704_20240826_v17b.csv"
)

node_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  filter(node_id == 81260300160901)

# NOTE: -0.1098754 is a node-specific PT bias correction; SWOT WSE plotted on
# the same vertical reference.
ggplot() +
  geom_errorbar(data = PT_df,
                aes(x = pt_time_UTC,
                    ymin = pt_wse_m - pt_correction_mean_total_error_m - 0.1098754,
                    ymax = pt_wse_m + pt_correction_mean_total_error_m - 0.1098754),
                color = "#99D8C9", linewidth = 5, alpha = 0.4) +
  geom_point(PT_df,
             mapping = aes(y = pt_wse_m - 0.1098754, x = pt_time_UTC),
             color = "#009E73", size = 1.2) +
  geom_point(node_SWOT_PT_vD,
             mapping = aes(y = wse, x = time_utc),
             shape = 21, fill = "#2474b7", color = "black",
             stroke = 1.6, size = 5.5) +
  xlab("Time") + ylab("WSE (m)") +
  theme_minimal(base_size = 30)


# -----------------------------------------------------------------------------
# 19. SWOT vs. ortho width along p_dist_out (Porcupine, 2024-07-10)
# -----------------------------------------------------------------------------

node_SWOT_ortho <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/node_width_SWOT_Ortho.csv") %>%
  filter(river == "PR") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5) %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139) %>%
  filter(as.Date(time_utc) == as.Date("2024-07-10"))

ggplot(node_SWOT_ortho) +
  geom_point(aes(x = p_dist_out / 1000, y = ortho_width_m),
             color = "#E69F00", size = 2.5, shape = 17, alpha = 0.7) +
  geom_point(aes(x = p_dist_out / 1000, y = width),
             color = "#2474b7", size = 2.5) +
  xlab("Distance to river outlet (km)") +
  ylab("Width (m)") +
  theme_minimal(base_size = 30)
# export dimensions: width 8.13 in, height 3.96 in










ortho_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/ortho_summary_stats.csv") %>%
  mutate(
    start_time_UTC = as.POSIXct(start_time_UTC, format = "%m/%d/%y %H:%M", tz = "UTC"),
    end_time_UTC   = as.POSIXct(end_time_UTC,   format = "%m/%d/%y %H:%M", tz = "UTC"),
    SWOT_time_UTC  = as.POSIXct(SWOT_time_UTC,  format = "%m/%d/%y %H:%M", tz = "UTC"))


# Compute midpoint
ortho_df$midpoint_time_UTC <- ortho_df$start_time_UTC + 
  (ortho_df$end_time_UTC - ortho_df$start_time_UTC) / 2

# Compute absolute offset in hours between midpoint and SWOT time
ortho_df$offset_time_UTC <- round(abs(as.numeric(difftime(ortho_df$midpoint_time_UTC, ortho_df$SWOT_time_UTC, units = "hours"))), 2)

ortho_df$survey_length_hours <- round(as.numeric(difftime(ortho_df$end_time_UTC, ortho_df$start_time_UTC, units = "hours")), 2)


# Remove midpoint column and reorder
ortho_df <- ortho_df[, c("River(s)", "start_time_UTC", "end_time_UTC", 
                         "survey_length_hours", "SWOT_pass_id", 
                         "SWOT_time_UTC", "offset_time_UTC")]

# Save to CSV
write.csv(ortho_df, 
          "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/ortho_summary_stats.csv",
          row.names = FALSE)
