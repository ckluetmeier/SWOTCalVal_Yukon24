library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare orthomosaic width & SWOT riverSP reach width
# ---------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep ortho data

# RiverSP Version D (SWORD v17b)
# directory with orthos
dir_path <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b"

# Bring in SWORD shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

# reaches that are fully digitized
viable_reaches <- data.frame(reach_id = c(81260300191, 81260300181, 81260400021, 81260400011, 81260500011, 81260300171, 81260300161,
                                          81250800031, 81270100041, 81270100051, 81270100061, 81270500161, 81270500171)) 

# RiverSP Version C (SWORD v16)
# directory with orthos
# dir_path <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16"
# 
# # Bring in SWORD shapefile
# sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v16/NA/na_sword_nodes_hb81_v16.shp")
# 
# # reaches that are fully digitized
# viable_reaches <- data.frame(reach_id = c())

# Get all ortho_ csv files
ortho_files <- list.files(
  path = dir_path,
  pattern = "^ortho_.*\\.csv$",
  full.names = TRUE)

# Function to assign SWOT_date based on filename
assign_swot_date <- function(filename) {
  base <- basename(filename)
  dplyr::case_when(
    grepl("lowerPR_SJ", base)         ~ "2024-07-26",
    grepl("_CD_", base)                 ~ "2024-07-11",
    grepl("lowerYR", base)            ~ "2024-07-16",
    grepl("upperYR", base)            ~ "2024-07-10",
    grepl("upperPR_CL_071024", base)  ~ "2024-07-10",
    grepl("upperPR_CL_071624", base)  ~ "2024-07-16",
    TRUE                              ~ NA_character_)}

# Read each file, add SWOT_date, then bind
ortho_width <- purrr::map_dfr(ortho_files, function(f) {
  readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::mutate(SWOT_date = assign_swot_date(f))})

# Convert to date type
ortho_width <- ortho_width %>%
  dplyr::mutate(SWOT_date = as.Date(SWOT_date))

# Add ortho_width to SWORD
ortho_width <- sword_sf %>%
  right_join(ortho_width, by = "node_id")

# Filter to viable_reaches & add width
ortho_width_filtered <- ortho_width %>%
  filter(reach_id %in% viable_reaches$reach_id) %>%
  mutate(ortho_width_m = water_area_m2/node_len)

# Compute average reach width (with full nodes)
ortho_reach_width <- ortho_width_filtered %>%
  dplyr::group_by(reach_id, SWOT_date) %>%
  dplyr::summarise(
    ortho_width_m = mean(ortho_width_m, na.rm = TRUE),
    n_nodes       = dplyr::n(),
    .groups       = "drop")

# add river names to ortho df
ortho_reach_width <- ortho_reach_width %>%
  mutate(river_code = substr(reach_id, 1, 6),
         river = case_when(
           # putting the reach id first ensures case_when won't overwrite SJ/BL labels
           # SWORD v16: "81260300061", "81260300231", "81260300241", "81260300251"
           # SWORD v17b: 81260300181", "81260300191", "81260300201", "81260300211
           # only SJ reaches need to be adjusted here
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
  mutate(source = "PGD0") ## CHANGE TO CORRECT VERSION!


# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep SWOT data

# RiverSP Version D (SWORD v17b)
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v17b/RiverSP_domain_reach_timeseries_PGD0_v17b.csv')

# RiverSP Version C (SWORD v16)
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# get ride of possible duplicates / empty observations
# (e.g. time = -999999999999, wse = -1.000000e+12)
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE) %>%
  filter(time > 0) %>%
  filter(wse > 0)

# filter SWOT data by reach_q (0=good, 1=suspect, 2=degraded, 3=bad), cross track distance, reach coverage
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) %>%
  filter(dark_frac < 0.8) %>%
  filter(partial_f == 0) # at least 50% node coverage if 0

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset


# ---------------------------------------------------------------------------------------------------------------------------
# match ortho & SWOT observations in time and space


# match ortho and SWOT in time and reach_id
# observations are matched by same day and reach_id
time_matched_SWOT_ortho <- ortho_reach_width %>%
  rowwise() %>%
  mutate(closest_match = list(SWOT_reach_df_filtered %>%
                                filter(as.Date(time_utc) == SWOT_date,
                                       reach_id == .env$reach_id) %>%
                                dplyr::select(-reach_id))) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())


# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: width diff calculation

# Calculate the width diff ortho - SWOT (residuals)
time_matched_SWOT_ortho$residuals = time_matched_SWOT_ortho$ortho_width_m - time_matched_SWOT_ortho$width

# Calculate the width percent diff Ortho/SWOT
# using ortho as truth:
# % diff = ( |swot - ortho| ) / ( ortho ) *100
time_matched_SWOT_ortho <- time_matched_SWOT_ortho %>%
  mutate(percent_diff = ((abs(ortho_width_m - width)) / ortho_width_m) * 100)


time_matched_SWOT_ortho <- time_matched_SWOT_ortho %>%
  filter(dark_frac < 0.5)

# # Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_matched_SWOT_ortho$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_matched_SWOT_ortho$residuals), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.68, na.rm=TRUE)
percentile_50_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

print(paste("68th Percentile Error:", percentile_68_percent))
print(paste("50th Percentile Error:", percentile_50_percent))

# # correlation test
cor_test <- cor.test(time_matched_SWOT_ortho$width, time_matched_SWOT_ortho$ortho_width_m)
# 
# # Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

# Reorder the factor levels for river
time_matched_SWOT_ortho$river <- factor(
  time_matched_SWOT_ortho$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

color_palette <- c("#F2C14E", "#F4845F", "#DA627D", "#9A348E") # full qual filters
# color_palette <- c("#F2C14E", "#8EAD7A", "#F4845F", "#DA627D", "#9A348E") # no dark water filter
# color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E") # no filters
# 
# # plot SWOT vs GNSS width
ggplot(time_matched_SWOT_ortho, aes(x = ortho_width_m, y = width, color = factor(river))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_matched_SWOT_ortho$ortho_width_m, na.rm = TRUE),
           y = max(time_matched_SWOT_ortho$width, na.rm = TRUE),
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(time_matched_SWOT_ortho)),
           hjust = 0, vjust = 1, size = 8) +
  # theme(legend.position = "none") + # comment off to see reach ids
  labs(color = "Reach ID")
# 
# # CDF plot of absolute wse difference
ggplot(time_matched_SWOT_ortho, aes(x = abs(width - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of absolute width residuals") +
  annotate("text", x = 350, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 350, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# CDF plot of percent difference
ggplot(time_matched_SWOT_ortho, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of %diff width residuals") +
  annotate("text", x = 70, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent, 4)), color = "#222222", size = 6) +
  annotate("text", x = 70, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)
