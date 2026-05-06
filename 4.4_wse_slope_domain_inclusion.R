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
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  filter(dark_frac < 0.5)

# GNSS
node_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  filter(dark_frac < 0.5)

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

# didn't remove bad nodes (e.g. no matched overpass) from the field data, so manually deleted these in QGIS















# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# Final CalVal dataframes of time/space matched SWOT & in situ data
# ---------------------------------------------------------------------------------------------------------------------------

# PT version C & D
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile") %>%
  filter(dark_frac < 0.5)

# GNSS version C & D
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)


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











# ---------------------------------------------------------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------------------------------------------------------


# Bring in SWORD shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")


sword_subset <- sword_sf %>%
  filter(reach_id %in% YR_domain$reach_id) %>%
  distinct(reach_id, .keep_all = TRUE)


cat("Total km of river:", sum(sword_subset$reach_len) / 1000, "km\n")



# add river names to df
sword_subset <- sword_subset %>%
  mutate(river_code = substr(reach_id, 1, 6),
         river = case_when(
           reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ", 
           reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141", "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
           river_code == "812701" ~ "lowerYR", # until the Circle bifurcation
           river_code == "812509" ~ "lowerYR", # past the PR confluence
           river_code == "812705" ~ "upperYR", # Circle up
           river_code == "812508" ~ "CD",
           river_code == "812603" ~ "PR",
           river_code == "812605" ~ "PR",
           river_code == "812604" ~ "CL",
           TRUE ~ NA_character_))


river_stats <- sword_subset %>%
  group_by(river) %>%
  summarise(total_km = sum(reach_len, na.rm = TRUE) / 1000,
            median_width = median(width),
            median_slope = median(slope)*100,
            n_distinct(reach_id)
            )



# add river names to df
GNSS_df <- GNSS_df %>%
  mutate(river_code = substr(reach_id, 1, 6),
         river = case_when(
           reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ", 
           reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141", "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
           river_code == "812701" ~ "lowerYR", # until the Circle bifurcation
           river_code == "812509" ~ "lowerYR", # past the PR confluence
           river_code == "812705" ~ "upperYR", # Circle up
           river_code == "812508" ~ "CD",
           river_code == "812603" ~ "PR",
           river_code == "812605" ~ "PR",
           river_code == "812604" ~ "CL",
           TRUE ~ NA_character_))


GNSS_stats <- sword_subset %>%
  group_by(river) %>%
  summarise(num_drift_reaches = sum(!is.na(reach_id)),
  )



# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# Final CalVal dataframes of time/space matched SWOT & in situ data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_PT.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)


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
  mutate(insitu_slope_m_m = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs)) %>%
  mutate(insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias))

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
         "/Users/camryn/Desktop/all_YR_domain_reaches_slope_subset.shp",
         delete_layer = TRUE
)

# colors for plot:
# #e97132
# #ececec
# #00008b
# dashed gray 2









# ---------------------------------------------------------------------------------------------------------------------------
# Figure 3
# ---------------------------------------------------------------------------------------------------------------------------

GNSS_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/3_SWOT_examples/data/GNSS_2024-08-20.csv")


GNSS_PR <- GNSS_df %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139)


GNSS_PR_oneday <- GNSS_PR %>%
  filter(as.Date(time_UTC) == as.Date("2024-08-20"))


# plot SWOT vs GNSS downstream


ggplot() +
  # geom_errorbar(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, ymin = wse - wse_u, ymax = wse + wse_u),
  #   linewidth = 1, width = 0) +
  geom_point(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, y = wse, color = wse), size = 4) +
  scale_color_gradient(low = "#2474b7", high = "#d3e3f3") +
  geom_point(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, y = mean_node_drift_wse_no_bias_m),
    color = "#CC79A7", size = 1) +
  geom_line(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, y = mean_node_drift_wse_no_bias_m, group = 1),
    color = "#CC79A7", linewidth = 0.7) +
  # Upper bound line
  geom_line(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, y = mean_node_drift_wse_no_bias_m + node_total_error_m, group = 1),
    color = "#CC79A7", alpha = 0.3, linewidth = 2) +
   # Lower bound line
  geom_line(data = GNSS_PR_oneday, aes(x = p_dist_out/1000, y = mean_node_drift_wse_no_bias_m - node_total_error_m, group = 1),
    color = "#CC79A7", alpha = 0.3, linewidth = 2) +
  guides(color = "none") +
  # ggtitle("GNSS") +
  xlab("Distance to river outlet (km)") +
  ylab("WSE (m)") +
  xlim(2138, 2152) +
  theme_minimal(base_size = 30)



# PT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_PR/flyby_SWOTCalVal_YR_PT_L1_PT225_20240704T120000_20240821T220000_20241124T014128_SWOTCalVal_YR_KEY_20240704_20240826_v17b.csv')
PT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_PR/flyby_SWOTCalVal_YR_PT_L1_PT230_20240704T120000_20240821T220000_20250714T183224_SWOTCalVal_YR_KEY_20240704_20240826_v17b.csv')

node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  filter(node_id == 81260300160901)

# 81260300160901
# 81260300170011

# plot SWOT vs PT timeseries
ggplot() +
  geom_errorbar(data = PT_df, aes(x = pt_time_UTC,
      ymin = pt_wse_m - pt_correction_mean_total_error_m - 0.1098754,
      ymax = pt_wse_m + pt_correction_mean_total_error_m - 0.1098754),
    color = "#99D8C9", linewidth = 5, alpha=0.4) +
  geom_point(PT_df, mapping=aes(y=pt_wse_m - 0.1098754, x=pt_time_UTC), color="#009E73", size=1.2) +
  # geom_errorbar(data = node_SWOT_PT_vD, mapping = aes(x = time_utc, ymin = wse - wse_u, ymax = wse + wse_u),
  #               linewidth = 1, width = 0) +
  geom_point(node_SWOT_PT_vD, mapping=aes(y=wse, x=time_utc), shape = 21, fill="#2474b7", color = "black", stroke =1.6, size=5.5) +
  # ggtitle("PT") + 
  xlab("Time") + ylab("WSE (m)") +
  theme_minimal(base_size = 30) 



node_SWOT_ortho <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/node_width_SWOT_Ortho.csv') %>%
  filter(river == "PR") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5) %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139) %>%
  filter(as.Date(time_utc) == as.Date("2024-07-10"))

# plot widths along dist_out
ggplot(node_SWOT_ortho) +
  geom_point(aes(x = p_dist_out/1000, y = ortho_width_m), color = "#E69F00", size = 2.5, shape = 17, alpha = 0.7) +
  geom_point(aes(x = p_dist_out/1000, y = width),  color = "#2474b7", size = 2.5) +
  xlab("Distance to river outlet (km)") +
  ylab("Width (m)") +
  theme_minimal(base_size = 30) 


# 8.13, 3.96
