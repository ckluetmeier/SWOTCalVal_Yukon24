# =============================================================================
# PT Data Consistency Checks
# -----------------------------------------------------------------------------
# Consistency checks to confirm PT observations make physical sense:
#   - WSE decreases downstream
#   - No unexplained jumps between adjacent PT sensors
#
# WARNING: dist_out topology in SWORD v16 is incorrect for the Sheenjek and
# Coleen rivers. Use v17b dist_out (p_dist_out) for those tributaries.
# =============================================================================

library(tidyverse)


# =============================================================================
# 1. Load SWORD node reference data
# =============================================================================

# --- SWORD v16 nodes ---------------------------------------------------------
# SWORD_node_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/SWORD_v16_domain_nodes.csv"
# )

# --- SWORD v17 nodes ---------------------------------------------------------
# SWORD_node_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/SWORD_v17_domain_nodes.csv"
# )

# --- SWORD v17b nodes (from RiverTile) ---------------------------------------
SWORD_node_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_domain_node_timeseries_v17b.csv"
)
# De-duplicate: v17b time series can have repeated node entries; keep one per node
SWORD_node_df <- SWORD_node_df %>%
  distinct(node_id, .keep_all = TRUE)


# =============================================================================
# 2. Load and merge PT data
# =============================================================================

# Working directory containing PT CSVs split by PT cluster
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_YR"
setwd(wd)

# All CSV files in the directory
csv_files <- list.files(wd, pattern = "\\.csv$", full.names = TRUE)

# Read and merge all PT cluster files into one data frame
data_list <- lapply(seq_along(csv_files), function(i) {
  read.csv(csv_files[i])
})
combined_PT_df <- bind_rows(data_list)

# Parse time column to POSIXct (UTC)
# NOTE: Some older toolbox runs format datetime as "%m/%d/%y %H:%M"; swap
# the format string below if the data were created by an older toolbox version.
combined_PT_df$pt_time_UTC <- as.POSIXct(
  combined_PT_df$pt_time_UTC,
  format = "%Y-%m-%d %H:%M:%S",
  tz     = "UTC"
)

# Join PT data to SWORD node attributes
# --- SWORD v16 / v17 join ---------------------------------------------------
# combined_PT_SWORD_df <- combined_PT_df %>%
#   left_join(
#     SWORD_node_df %>% select(Node_ID, dist_out, node_len),
#     by = "Node_ID"
#   )

# --- SWORD v17b join --------------------------------------------------------
combined_PT_SWORD_df <- combined_PT_df %>%
  left_join(
    SWORD_node_df %>% select(node_id, p_dist_out),
    by = c("Node_ID" = "node_id")
  )


# =============================================================================
# 3. Consistency plots
# =============================================================================

# Plot 1: PT WSE time series colored by distance from outlet (km).
ggplot(combined_PT_SWORD_df, aes(x = pt_time_UTC, y = pt_wse_m, color = p_dist_out * 0.001)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    x     = "Time (UTC)",
    y     = "Water Surface Elevation (m)",
    color = "dist_out (km)"
  ) +
  ggtitle("upper_YR — WSE by distance from outlet") +
  theme_minimal(base_size = 15)

# Plot 2: PT WSE time series colored by PT serial number.
# Checking for gaps or jumps
ggplot(combined_PT_SWORD_df, aes(x = pt_time_UTC, y = pt_wse_m, color = factor(pt_serial))) +
  geom_point(size = 0.1) +
  labs(
    x = "Time (UTC)",
    y = "Water Surface Elevation (m)"
  ) +
  ggtitle("upper_YR — WSE by PT serial") +
  theme_minimal(base_size = 15)
