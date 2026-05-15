# =============================================================================
# Bias-Corrected PT and GNSS Reach Slopes
# -----------------------------------------------------------------------------
# Applies WSE bias corrections (derived from SWOT–in situ node comparisons)
# to PT and GNSS reach-level slope estimates. Produces relative slope CSVs.
#
#   Part 1 - PT slopes:   per-sensor (pt_serial) median bias from node WSE
#                         comparisons subtracted from reach boundary WSEs,
#                         then slope recomputed.
#   Part 2 - GNSS slopes: per-drift (drift_id) median bias from node WSE
#                         comparisons subtracted from node WSEs, then slope
#                         computed from upstream/downstream corrected nodes.
#
# Toggle between SWORD v16 (PIC0) and v17b (PGD0).
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)


# =============================================================================
# PART 1 — PT Slope Bias Correction
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

# --- RiverSP PGD0 (SWORD v17b) -----------------------------------------------
# node_SWOT_PT_vD <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv")
# PT_bias <- group_by(node_SWOT_PT_vD, pt_serial) %>%
#   summarise(
#     count    = n(),
#     bias     = mean(abs(bias), na.rm = TRUE),
#     reach_id = mean(abs(reach_id), na.rm = TRUE),
#     node_id  = mean(abs(node_id), na.rm = TRUE),
#     river    = first(river))
#
# YR_PT_reach_slope <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b/YR_PT_reach_slope.csv")

# --- RiverSP PIC0 (SWORD v16) -----------------------------------------------
node_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv")

PT_bias <- group_by(node_SWOT_PT_vC, pt_serial) %>%
  summarise(
    count    = n(),
    bias     = mean(abs(bias), na.rm = TRUE),
    reach_id = mean(abs(reach_id), na.rm = TRUE),
    node_id  = mean(abs(node_id), na.rm = TRUE),
    river    = first(river))

YR_PT_reach_slope <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16/YR_PT_reach_slope.csv")


# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

# Given a combined "pt_serial_node_id_pt_serial_node_id_..." string, extract
# only the pt_serial values (odd-indexed tokens after splitting on "_").
extract_pt_serials <- function(pt_serials_str) {
  tokens <- str_split(pt_serials_str, "_")[[1]]
  as.numeric(tokens[seq(1, length(tokens), by = 2)])
}

# Given a vector of pt_serial values, look up their biases from PT_bias and
# return the mean bias and number of matched sensors.
get_bias_stats <- function(serials, bias_df) {
  biases <- bias_df$bias[match(serials, bias_df$pt_serial)]
  list(
    bias_mean = mean(biases, na.rm = TRUE),
    n_matched = sum(!is.na(biases))
  )
}


# -----------------------------------------------------------------------------
# 3. Apply bias correction to each reach row
# -----------------------------------------------------------------------------

YR_PT_reach_slope_corrected <- YR_PT_reach_slope %>%
  rowwise() %>%
  mutate(
    # Parse upstream and downstream pt_serial lists
    serials_us = list(extract_pt_serials(pt_serials_us)),
    serials_ds = list(extract_pt_serials(pt_serials_ds)),

    # Look up mean bias for each boundary
    bias_stats_us = list(get_bias_stats(serials_us, PT_bias)),
    bias_stats_ds = list(get_bias_stats(serials_ds, PT_bias)),

    # Store boundary bias values
    bias_us = bias_stats_us$bias_mean,
    bias_ds = bias_stats_ds$bias_mean,

    # Subtract mean bias from each boundary WSE
    mean_pt_wse_us_boundary_m = mean_pt_wse_us_boundary_m - bias_stats_us$bias_mean,
    mean_pt_wse_ds_boundary_m = mean_pt_wse_ds_boundary_m - bias_stats_ds$bias_mean,

    # Recompute reach slope from corrected boundary WSEs
    slope_m_m = (mean_pt_wse_us_boundary_m - mean_pt_wse_ds_boundary_m) / reach_length,

    # Look up river name from the first upstream sensor
    river = PT_bias$river[match(serials_us[[1]][1], PT_bias$pt_serial)]
  ) %>%
  ungroup() %>%
  select(-serials_us, -serials_ds, -bias_stats_us, -bias_stats_ds)


# -----------------------------------------------------------------------------
# 4. Save corrected PT reach slopes
# -----------------------------------------------------------------------------

write.csv(
  YR_PT_reach_slope_corrected,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16/YR_PT_reach_slope_corrected.csv",
  row.names = FALSE)


# =============================================================================
# PART 2 — GNSS Reach Slope Bias Correction
# =============================================================================


# -----------------------------------------------------------------------------
# 5. Load domain and SWORD geometry
# -----------------------------------------------------------------------------

# --- RiverSP PGD0 (SWORD v17b) -----------------------------------------------
# YR_domain    <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv")
# 
# sword_reaches <- st_read(
#   "/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")
# 
# sword_nodes   <- st_read(
#   "/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

# --- RiverSP PIC0 (SWORD v16) -----------------------------------------------
YR_domain     <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v16.csv")

sword_reaches <- st_read(
  "/Users/camryn/Desktop/SWORD_v16/NA/na_sword_reaches_hb81_v16.shp")

sword_nodes   <- st_read(
  "/Users/camryn/Desktop/SWORD_v16/NA/na_sword_nodes_hb81_v16.shp")

# Subset reaches and nodes to the Yukon River domain
YR_reaches <- sword_reaches %>%
  filter(reach_id %in% YR_domain$Reach_ID) %>%
  st_drop_geometry() %>%
  select(reach_id, reach_len)

YR_nodes <- sword_nodes %>%
  filter(node_id %in% YR_domain$Node_ID) %>%
  st_drop_geometry() %>%
  select(reach_id, node_id)

# For each reach, identify the upstream (max node_id) and downstream (min node_id) nodes
YR_node_summary <- YR_nodes %>%
  group_by(reach_id) %>%
  summarise(
    us_node_id = max(node_id),
    ds_node_id = min(node_id),
    .groups    = "drop"
  )

# Build reach info table with boundary node IDs and reach length
YR_domain_reach_info <- YR_reaches %>%
  left_join(YR_node_summary, by = "reach_id") %>%
  rename(reach_length = reach_len) %>%
  select(reach_id, us_node_id, ds_node_id, reach_length)


# -----------------------------------------------------------------------------
# 6. Load GNSS node WSE data
# -----------------------------------------------------------------------------

# --- RiverSP PGD0 (SWORD v17b) -----------------------------------------------
# node_SWOT_GNSS_vD <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv")
#
# GNSS_bias <- group_by(node_SWOT_GNSS_vD, drift_id) %>%
#   summarise(
#     count = n(),
#     bias  = mean(abs(bias), na.rm = TRUE),
#     river = first(river))
# 
# YR_GNSS_node_wse <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv")

# --- RiverSP PIC0 (SWORD v16) -----------------------------------------------
node_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv")

GNSS_bias <- group_by(node_SWOT_GNSS_vC, drift_id) %>%
  summarise(
    count = n(),
    bias  = mean(abs(bias), na.rm = TRUE),
    river = first(river))

YR_GNSS_node_wse <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_node_wses.csv")


# -----------------------------------------------------------------------------
# 7. Apply GNSS node WSE bias correction
# -----------------------------------------------------------------------------

YR_GNSS_node_wse <- YR_GNSS_node_wse %>%
  left_join(GNSS_bias %>% select(drift_id, bias), by = "drift_id") %>%
  mutate(mean_node_drift_nobias_wse_m = mean_node_drift_wse_m - bias) %>%
  select(-bias)


# -----------------------------------------------------------------------------
# 8. Compute bias-corrected GNSS reach slopes
# -----------------------------------------------------------------------------

YR_drift_reach_slope_corrected <- YR_domain_reach_info %>%
  # Join corrected upstream node WSE
  left_join(
    YR_GNSS_node_wse %>%
      select(
        node_id,
        drift_id,
        us_node_drift_nobias_wse_m = mean_node_drift_nobias_wse_m,
        us_node_time_UTC           = time_UTC,
        us_node_total_error_m      = node_total_error_m
      ),
    by = c("us_node_id" = "node_id")
  ) %>%
  # Join corrected downstream node WSE, matched on the same drift_id
  left_join(
    YR_GNSS_node_wse %>%
      select(
        node_id,
        drift_id,
        ds_node_drift_nobias_wse_m = mean_node_drift_nobias_wse_m,
        ds_node_time_UTC           = time_UTC,
        ds_node_total_error_m      = node_total_error_m
      ),
    by = c("ds_node_id" = "node_id", "drift_id")
  ) %>%
  mutate(
    # Compute reach slope from corrected boundary WSEs
    slope_m_m               = (us_node_drift_nobias_wse_m - ds_node_drift_nobias_wse_m) / reach_length,
    # Mean total error across both boundary nodes
    mean_node_total_error_m = (us_node_total_error_m + ds_node_total_error_m) / 2
  ) %>%
  select(
    reach_id,
    drift_id,
    us_node_drift_nobias_wse_m,
    us_node_time_UTC,
    ds_node_drift_nobias_wse_m,
    ds_node_time_UTC,
    reach_length,
    slope_m_m,
    mean_node_total_error_m
  ) %>%
  filter(!is.na(slope_m_m))


# -----------------------------------------------------------------------------
# 9. Save corrected GNSS reach slopes
# -----------------------------------------------------------------------------

write.csv(
  YR_drift_reach_slope_corrected,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_slope_corrected.csv",
  row.names = FALSE)
