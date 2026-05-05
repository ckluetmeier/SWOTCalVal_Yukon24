library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)

# Load AK SWORD nodes & reaches
sword_nodes <- st_read('/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp')
sword_reaches <- st_read('/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp')
YR_domain <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv')

# add river names to nodes
YR_sword_nodes <- sword_nodes %>%
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
           TRUE ~ NA_character_)) %>% 
  filter(!is.na(river))


# Subset nodes to YR_domain
YR_sword_nodes <- YR_sword_nodes %>%
  filter(node_id %in% YR_domain$Node_ID) %>%
  st_drop_geometry()







# Width
# Trim to 25th-75th percentile per river
trimmed_df <- YR_sword_nodes %>%
  group_by(river) %>%
  mutate(
    q25 = quantile(width, 0.25, na.rm = TRUE),
    q75 = quantile(width, 0.75, na.rm = TRUE),
    q25 = if_else(river == "CL", quantile(width, 0.5, na.rm = TRUE), q25),
    q75 = if_else(river == "CL", quantile(width, 0.95, na.rm = TRUE), q75),
    q25 = if_else(river == "lowerYR", quantile(width, 0.4, na.rm = TRUE), q25),
    q75 = if_else(river == "lowerYR", quantile(width, 0.6, na.rm = TRUE), q75),
    q25 = if_else(river == "upperYR", quantile(width, 0.3, na.rm = TRUE), q25),
    q75 = if_else(river == "upperYR", quantile(width, 0.7, na.rm = TRUE), q75)
  ) %>%
  filter(width >= q25, width <= q75) %>%
  ungroup() %>%
  filter(river != 'BL')


# Compute medians per river (median = 50th percentile)
median_df <- trimmed_df %>%
  group_by(river) %>%
  summarize(median_width = median(width, na.rm = TRUE)) %>%
  ungroup()

color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A",  "#DA627D")

# Density plot with medians as vertical lines; colors/fills match your palette
ggplot(trimmed_df, aes(x = width, fill = river, color = river)) +
  geom_density(alpha = 0.5, size = 0.9, adjust = 1) + 
  # Add median lines (dashed) colored by river
  geom_vline(
    data = median_df,
    aes(xintercept = median_width, color = river),
    linetype = "dashed",
    size = 1 ) +
  labs(x = "Width",y = "Density", fill = "River", color = "River") +
  theme_minimal(base_size = 20) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank())





# Slope

# add river names to reaches
YR_sword_reaches <- sword_reaches %>%
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
           TRUE ~ NA_character_)) %>% 
  filter(!is.na(river))

# Subset reaches to YR_domain
YR_reaches <- YR_sword_reaches %>%
  filter(reach_id %in% YR_domain$Reach_ID) %>%
  st_drop_geometry()







# Trim to 25th-75th percentile per river
trimmed_df <- YR_reaches %>%
  filter(
    river != 'BL',
    !(reach_id > 81250800041 & river == "CD"),
    !(reach_id > 81270500231 & river == "upperYR"),
    (river == "PR" & slope < 0.6) |
      (river == "lowerYR" & slope < 0.7) |
      (!(river %in% c("PR", "lowerYR")) & slope < 1.40)
  )


# group_by(river) %>%
#   mutate(
#     q25 = if_else(river == "", quantile(width, 0.25, na.rm = TRUE), q25),
#     q75 = if_else(river == "", quantile(width, 0.75, na.rm = TRUE), q75)
#   ) %>%
#   filter(slope >= q25, slope <= q75) %>%
#   ungroup() %>%

# Compute medians per river (median = 50th percentile)
median_df <- trimmed_df %>%
  group_by(river) %>%
  summarize(median_slope = median(slope, na.rm = TRUE)) %>%
  ungroup()

color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A",  "#DA627D")

# Density plot with medians as vertical lines; colors/fills match your palette
ggplot(trimmed_df, aes(x = slope*100, fill = river, color = river)) +
  geom_density(alpha = 0.5, size = 0.9, adjust = 1) + 
  # Add median lines (dashed) colored by river
  geom_vline(
    data = median_df,
    aes(xintercept = median_slope*100, color = river),
    linetype = "dashed",
    size = 1 ) +
  labs(x = "Slope",y = "Density", fill = "River", color = "River") +
  theme_minimal(base_size = 20) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()) +
  xlim(0,160)

# 11.05, 2.22

















# ---------------------------------------------------------------------------------------------------------------------------
# SWOT node dark water %
# ---------------------------------------------------------------------------------------------------------------------------

SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_nodes_merged_RiverSP.csv')

# get ride of possible duplicates from hydrocron pull
# this also filters out bad nodes without data (e.g. time = -999999999999, wse = -1.000000e+12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# filter SWOT data by node_q (0=good, 1=suspect, 2=degraded, 3=bad)
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  #filter(wse_u <= 0.5) %>% #this probably won't filter out many nodes beyond what node_q is doing
  filter(node_q <= 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000)
# filter(dark_frac <= 50)

# Compute the average dark_frac for each unique node_id
SWOT_node_dark_frac <- SWOT_df_filtered %>%
  group_by(node_id) %>%
  summarise(
    count = n(),
    dark_frac_sd = sd(dark_frac, na.rm = TRUE),
    dark_frac_median = median(dark_frac, na.rm = TRUE),
    dark_frac_IQR = IQR(dark_frac, na.rm = TRUE),
    dark_frac_min =min(dark_frac, na.rm = TRUE),
    dark_frac_max =max(dark_frac, na.rm = TRUE),
    dark_frac = mean(dark_frac, na.rm = TRUE),
    reach_id = first(reach_id),  
    lat = first(lat),  
    lon = first(lon),
    p_dist_out = first(p_dist_out)
  ) %>%
  ungroup()

# plot SWOT vs PT wse
ggplot(SWOT_node_dark_frac, aes(x = p_dist_out*0.001, y = dark_frac_IQR)) +
  geom_point(size = 4) +
  xlab("distance to outlet (km)") +
  ylab("dark frac IQR") +
  theme_minimal(base_size = 30)

write_csv(SWOT_node_dark_frac,'/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/SWOT_node_dark_frac.csv')

