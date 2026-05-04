library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)

# ---------------------------------------------------------------------------------------------------------------------------
# Correct PT & GNSS slopes with bias corrected WSE values
# ---------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------
# Load dfs 

# RiverSP PGD0
# node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv")
# 
# PT_bias <- group_by(node_SWOT_PT_vD, pt_serial) %>% summarise(
#   count = n(),
#   bias = mean(abs(bias), na.rm = TRUE),
#   reach_id = mean(abs(reach_id), na.rm = TRUE),
#   node_id = mean(abs(node_id), na.rm = TRUE),
#   river = first(river)
# )
# 
# YR_PT_reach_slope <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b/YR_PT_reach_slope.csv')


# RiverSP PIC0

node_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv")

PT_bias <- group_by(node_SWOT_PT_vC, pt_serial) %>% summarise(
  count = n(),
  bias = mean(abs(bias), na.rm = TRUE),
  reach_id = mean(abs(reach_id), na.rm = TRUE),
  node_id = mean(abs(node_id), na.rm = TRUE),
  river = first(river)
)

YR_PT_reach_slope <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16/YR_PT_reach_slope.csv')

# ---------------------------------------------------------------------------------------------------------------------------
# Apply the WSE bias correction

# Helper: given a string of alternating pt_serial / node_id extract just pt_serial values
extract_pt_serials <- function(pt_serials_str) {
  tokens <- str_split(pt_serials_str, "_")[[1]]
  # pt_serials are at odd indices (1, 3, 5, ...)
  as.numeric(tokens[seq(1, length(tokens), by = 2)])
}

# Helper: given a vector of pt_serials, look up their biases and return the mean
get_bias_stats <- function(serials, bias_df) {
  biases <- bias_df$bias[match(serials, bias_df$pt_serial)]
  list(
    bias_mean = mean(biases, na.rm = TRUE),
    n_matched = sum(!is.na(biases))
  )
}

# Apply bias correction to each row
YR_PT_reach_slope_corrected <- YR_PT_reach_slope %>%
  rowwise() %>%
  mutate(
    # Extract pt_serials for us and ds
    serials_us = list(extract_pt_serials(pt_serials_us)),
    serials_ds = list(extract_pt_serials(pt_serials_ds)),
    
    # Get bias stats
    bias_stats_us = list(get_bias_stats(serials_us, PT_bias)),
    bias_stats_ds = list(get_bias_stats(serials_ds, PT_bias)),
    
    # Mean bias for output columns
    bias_us = bias_stats_us$bias_mean,
    bias_ds = bias_stats_ds$bias_mean,
    
    # Corrected WSEs: subtract mean bias
    mean_pt_wse_us_boundary_m = mean_pt_wse_us_boundary_m - bias_stats_us$bias_mean,
    mean_pt_wse_ds_boundary_m = mean_pt_wse_ds_boundary_m - bias_stats_ds$bias_mean,
    
    # Recompute slope
    slope_m_m = (mean_pt_wse_us_boundary_m - mean_pt_wse_ds_boundary_m) / reach_length,
    
    # Get river name
    river = PT_bias$river[match(serials_us[[1]][1], PT_bias$pt_serial)]
  ) %>%
  ungroup() %>%
  select(-serials_us, -serials_ds, -bias_stats_us, -bias_stats_ds)



# save to csv
write.csv(YR_PT_reach_slope_corrected, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16/YR_PT_reach_slope_corrected.csv', row.names = FALSE)










# ---------------------------------------------------------------------------------------------------------------------------
# Correct GNSS slopes with bias corrected WSE values
# ---------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------
# Prep reach info

# RiverSP PGD0
# Load YR domain
# YR_domain <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv')

# Load SWORD shapefiles
# sword_reaches <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")
# sword_nodes   <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")


# RiverSP PIC0
# Load YR domain
YR_domain <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v16.csv')

# Load SWORD shapefiles
sword_reaches <- st_read('/Users/camryn/Desktop/SWORD_v16/NA/na_sword_reaches_hb81_v16.shp')
sword_nodes <- st_read('/Users/camryn/Desktop/SWORD_v16/NA/na_sword_nodes_hb81_v16.shp')


# Subset reaches to YR_domain
YR_reaches <- sword_reaches %>%
  filter(reach_id %in% YR_domain$Reach_ID) %>%
  st_drop_geometry() %>%
  select(reach_id, reach_len)

# Subset nodes to YR_domain
YR_nodes <- sword_nodes %>%
  filter(node_id %in% YR_domain$Node_ID) %>%
  st_drop_geometry() %>%
  select(reach_id, node_id)

# For each reach, get upstream (max node_id) and downstream (min node_id)
YR_node_summary <- YR_nodes %>%
  group_by(reach_id) %>%
  summarise(
    us_node_id = max(node_id),
    ds_node_id = min(node_id),
    .groups = "drop"
  )

# Join reach info with node summary
YR_domain_reach_info <- YR_reaches %>%
  left_join(YR_node_summary, by = "reach_id") %>%
  rename(reach_length = reach_len) %>%
  select(reach_id, us_node_id, ds_node_id, reach_length)



# ---------------------------------------------------------------------------------------------------------------------------
# Load GNSS dfs


# RiverSP PGD0
# ---------------------------
# node_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv")

# GNSS_bias <- group_by(node_SWOT_GNSS_vD, drift_id) %>% summarise(
#   count = n(),
#   bias = mean(abs(bias), na.rm = TRUE),
#   river = first(river)
# )
# 
# YR_GNSS_node_wse <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv')


# RiverSP PIC0
# ---------------------------
node_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv")

GNSS_bias <- group_by(node_SWOT_GNSS_vC, drift_id) %>% summarise(
  count = n(),
  bias = mean(abs(bias), na.rm = TRUE),
  river = first(river)
)

YR_GNSS_node_wse <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_node_wses.csv')



# ---------------------------------------------------------------------------------------------------------------------------
# Apply the WSE bias correction


YR_GNSS_node_wse <- YR_GNSS_node_wse %>%
  left_join(GNSS_bias %>% select(drift_id, bias), by = "drift_id") %>%
  mutate(mean_node_drift_nobias_wse_m = mean_node_drift_wse_m - bias) %>%
  select(-bias)



YR_drift_reach_slope_corrected <- YR_domain_reach_info %>%
  # Join upstream node data
  left_join(
    YR_GNSS_node_wse %>%
      select(node_id, drift_id, 
             us_node_drift_nobias_wse_m = mean_node_drift_nobias_wse_m,
             us_node_time_UTC           = time_UTC,
             us_node_total_error_m      = node_total_error_m),
    by = c("us_node_id" = "node_id")
  ) %>%
  # Join downstream node data, matching on the same drift_id
  left_join(
    YR_GNSS_node_wse %>%
      select(node_id, drift_id,
             ds_node_drift_nobias_wse_m = mean_node_drift_nobias_wse_m,
             ds_node_time_UTC           = time_UTC,
             ds_node_total_error_m      = node_total_error_m),
    by = c("ds_node_id" = "node_id", "drift_id")
  ) %>%
  # Calculate slope and mean node error
  mutate(
    slope_m_m              = (us_node_drift_nobias_wse_m - ds_node_drift_nobias_wse_m) / reach_length,
    mean_node_total_error_m = (us_node_total_error_m + ds_node_total_error_m) / 2
  ) %>%
  # Select and order final columns
  select(reach_id, drift_id,
         us_node_drift_nobias_wse_m, us_node_time_UTC,
         ds_node_drift_nobias_wse_m, ds_node_time_UTC,
         reach_length, slope_m_m, mean_node_total_error_m) %>%
  filter(!is.na(slope_m_m))


# save to csv
write.csv(YR_drift_reach_slope_corrected, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_slope_corrected.csv', row.names = FALSE)


