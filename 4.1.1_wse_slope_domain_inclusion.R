library(sf)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# SWOT WSE & slope domain inclusion by version
# ---------------------------------------------------------------------------------------------------------------------------

# Final CalVal dataframes of time/space matched SWOT & in situ data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
node_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  rename(old_node_id = node_id)
node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT")

# GNSS
node_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  rename(old_node_id = node_id)
node_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS")

# Original in situ dataframes
# ---------------------------------------------------------------------------------------------------------------------------

# Set working directory to the folder chucked by separate rivers and PT clusters
# SWORD v17b
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_YR"
setwd(wd)

# Get list of all PT CSV files (these are munged PT dataframes created by the toolboxes)
csv_files <- list.files(wd, pattern = "\\.csv$", full.names = TRUE)

# Merge all PT files into a combined dataframe
data_list <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  return(df)
})

combined_PT_df <- bind_rows(data_list)

# Convert time column to POSIXct
combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# SWORD v17b
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv')

YR_domain <- bind_rows(combined_PT_df, GNSS_df) %>%
  distinct(node_id)


# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv")

# translate the vC data to SWORD v17b
node_SWOT_PT_vC <- node_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")
  ) %>%
  rename(node_id = v17_node_id)

node_SWOT_GNSS_vC <- node_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")
  ) %>%
  rename(node_id = v17_node_id)


# ---------------------------------------------------------------------------------------------------------------------------
# Merge all the dataframes together
# ---------------------------------------------------------------------------------------------------------------------------

# merge all dataframes together
node_SWOT_full_insitu <- bind_rows(node_SWOT_PT_vC, node_SWOT_PT_vD, node_SWOT_GNSS_vC, node_SWOT_GNSS_vD)
  


all_nodes <- node_SWOT_full_insitu %>%
  distinct(node_id, source) %>%                 # one row per node_id × source
  group_by(node_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile ~ 0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(node_id, version_inclusion)

all_YR_domain_nodes <- YR_domain %>%
  left_join(all_nodes, by = "node_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))



# ---------------------------------------------------------------------------------------------------------------------------
# Save as a shapefile
# ---------------------------------------------------------------------------------------------------------------------------

# Bring in SWORD shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

# Add version_inclusion to SWORD
all_YR_domain_nodes_sf <- sword_sf %>%
  left_join(all_YR_domain_nodes, by = "node_id")

# Subset to only version_inclusion reaches
all_YR_domain_nodes_sf_subset <- all_YR_domain_nodes_sf %>%
  filter(!is.na(version_inclusion))

# save!
st_write(all_YR_domain_nodes_sf_subset,
         "/Users/camryn/Desktop/all_YR_domain_nodes_subset.shp",
         delete_layer = TRUE
)




















# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# Final CalVal dataframes of time/space matched SWOT & in situ data
# ---------------------------------------------------------------------------------------------------------------------------

# PT version C & D
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP") %>%
  rename(old_reach_id = reach_id)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile")

# GNSS version C & D
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_SWOT_GNSS.csv")


# Original in situ dataframes
# ---------------------------------------------------------------------------------------------------------------------------

# PT with SWORD v17b
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b"
setwd(wd)

# Get list of all PT reach csv files (these are munged PT dataframes created by the toolboxes)
csv_files <- list.files(wd, pattern = "^YR_812.*\\.csv$", full.names = TRUE)

# Merge all PT files into a combined dataframe
data_list <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  return(df)
})
data_list <- lapply(data_list, function(df) {
  df %>% mutate(sorted_nodelist = as.character(sorted_nodelist))
})
combined_PT_df <- bind_rows(data_list)

combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, tz = "UTC")

# GNSS with SWORD v17b
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv')

YR_domain <- bind_rows(combined_PT_df, GNSS_df) %>%
  distinct(reach_id)

# this is the full YR domain not contained to actual sampled reaches
# YR_domain <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv") %>%
#   distinct(v17b_reach_id) %>%
#   rename(reach_id = v17b_reach_id)


# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

# translate the vC data to SWORD v17b
reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)


# ---------------------------------------------------------------------------------------------------------------------------
# Merge all the dataframes together
# ---------------------------------------------------------------------------------------------------------------------------

# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vC, reach_SWOT_PT_vD, reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_wse_m = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m))


all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source) %>%                 # one row per reach_id × source
  group_by(reach_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile ~ 0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(reach_id, version_inclusion)

all_YR_domain_reaches <- YR_domain %>%
  left_join(all_reaches, by = "reach_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))



# ---------------------------------------------------------------------------------------------------------------------------
# Save as a shapefile
# ---------------------------------------------------------------------------------------------------------------------------

# Bring in SWORD shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

# Add version_inclusion to SWORD
all_YR_domain_reaches_sf <- sword_sf %>%
  left_join(all_YR_domain_reaches, by = "reach_id")

# Subset to only version_inclusion reaches
all_YR_domain_reaches_sf_subset <- all_YR_domain_reaches_sf %>%
  filter(!is.na(version_inclusion))

# save!
st_write(all_YR_domain_reaches_sf_subset,
  "/Users/camryn/Desktop/all_YR_domain_reaches_subset.shp",
  delete_layer = TRUE
)

# colors for plot:
# #e97132
# #ececec
# #00008b
# dashed gray 2


