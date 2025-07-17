library(tidyverse)

# ---------------------------------------------------------------------------------------------------------------------------
# Consistency checks to make sure AK PT data make physical sense
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# load SWORD nodes -- make sure to turn on SWORD version that matches PT processing version!

# SWORD v16 nodes
# SWORD_node_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/SWORD_v16_domain_nodes.csv')
# SWORD v17 nodes
# SWORD_node_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/SWORD_v17_domain_nodes.csv')
# SWORD v17b nodes (loaded by pulling in RiverTile SWOT df from JPL as of July, 2025)
SWORD_node_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_domain_node_timeseries_v17b.csv')
SWORD_node_df <- SWORD_node_df %>%
  distinct(node_id, .keep_all = TRUE) #only need this filter with v17b

# ---------------------------------------------------------------------------------------------------------------------------
# Munge PT data & join to SWORD

# Set working directory to a folder with PTs chucked by separate rivers and clusters
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/_node/SWORD_v17b/upper_YR"
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
combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
# watch out for funky datetimes in PT data -- some old toolbox runs vary in how datetime is output
# for example, run this line with old upper_YR PTs
# combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, format = "%m/%d/%y %H:%M", tz = "UTC")

# combine PT df with SWORD nodes
# SWORD v16, v17
# combined_PT_SWORD_df <- combined_PT_df %>%
#   left_join(SWORD_node_df %>% select(Node_ID, dist_out, node_len), by = "Node_ID")
# SWORD v17b
combined_PT_SWORD_df <- combined_PT_df %>%
  left_join(SWORD_node_df %>% select(node_id, p_dist_out), by = c("Node_ID" = "node_id"))

# ---------------------------------------------------------------------------------------------------------------------------
# Plot PT data to make sure PT wse go downstream as expected by install node & look for jumps in PTs

# watch out for topology issues with dist_out in SWORD v16!!
# Sheenjek and Coleen dist_out are incorrect!!!

# Plot PT timeseries for each cluster colored by dist_out
ggplot(combined_PT_SWORD_df, aes(x = pt_time_UTC, y = pt_wse_m, color = p_dist_out*0.001)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(x = "Time (UTC)", y = "Water Surface Elevation (m)", color = "dist_out (km)") +
  theme_minimal(base_size = 15) +
  ggtitle('upper_YR')

# Plot PT timeseries for each cluster colored by pt_serial
ggplot(combined_PT_SWORD_df, aes(x = pt_time_UTC, y = pt_wse_m, color = factor(pt_serial))) +
  geom_point(size = 0.1) +
  labs(x = "Time (UTC)", y = "Water Surface Elevation (m)") +
  theme_minimal(base_size = 15) +
  ggtitle('upper_YR')
