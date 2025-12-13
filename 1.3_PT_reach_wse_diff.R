library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare PT wse/slope & SWOT riverSP reach wse/slope
# ---------------------------------------------------------------------------------------------------------------------------

# WSE
# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# RiverSP
SWOT_reach_df <- read.csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# RiverTile
# SWORD v16
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_domain_reach_timeseries_v16.csv')
# SWORD v17v
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v17b/RiverTile_domain_reach_timeseries_v17b.csv')

# get ride of possible duplicates from hydrocron pull
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE)

# filter SWOT data by reach_q (0=good, 1=suspect, 2=degraded, 3=bad) & cross track distance, >50% reach obs, dark water
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) %>%
  filter(partial_f == 0)  %>%
  filter(dark_frac <= 0.8)

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset



# ---------------------------------------------------------------------------------------------------------------------------
# read in PT data

# Set working directory to the reach df folder from the toolboxes
# SWORD v16
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16"
# SWORD v17b
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b"
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

#change reach_id so there aren't duplicate columns when merging with SWOT data
combined_PT_df <- rename(combined_PT_df, "PT_reach_id" = "reach_id")



# ---------------------------------------------------------------------------------------------------------------------------
# match PT & SWOT observations in time and space

# match PT and SWOT in time
# observation are matched by 7.5min buffer (a SWOT obs should always be within 7.5min of a PT)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(SWOT_reach_df_filtered %>%
                           filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5))
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# match PT and SWOT in space
# reach level
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(PT_reach_id == reach_id)

# get rid of any duplicate PT values
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[!duplicated(time_space_matched_SWOT_PT[c("mean_reach_pt_wse_m","pt_time_UTC","flaglist")]),]

# add river names to df
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      river_code == "812701" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "upperPR",
      river_code == "812605" ~ "upperPR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )


# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (residuals)
time_space_matched_SWOT_PT$residuals = time_space_matched_SWOT_PT$mean_reach_pt_wse_m - time_space_matched_SWOT_PT$wse

# Filter values to sensical residuals (< 10 m diff)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  filter(abs(residuals) < 10) # %>%
# filter(node_total_error_m < 1)

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

# # Calculating the linear regression model 
# model = lm(mean_reach_pt_wse_m~wse, data = time_space_matched_SWOT_PT) 
# # Extracting R-squared parameter from summary 
# summary(model)
# #RMSE
# rmse <- sqrt(mean((time_space_matched_SWOT_PT$mean_reach_pt_wse_m - time_space_matched_SWOT_PT$wse)^2))
# #RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$mean_reach_pt_wse_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

# for rivers
color_palette <- c("#D86A1A", "#6D398B",  "#F8A31B", "#00429D", "#2E7D32",
                   "#00429D",  "#C83232", "#008F7A", "#E3A700", "#124000")

# plot SWOT vs PT wse
ggplot(time_space_matched_SWOT_PT, aes(x = mean_reach_pt_wse_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_PT$mean_reach_pt_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_PT)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "River") 

# for basin codes
color_palette <- c("#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8", "#F8A31B", 
                   "#00429D", "#2E7D32", "#C83232", "#008F7A", "#E3A700", "#124000")

# plot SWOT vs PT wse
# by basin code
ggplot(time_space_matched_SWOT_PT, aes(x = mean_reach_pt_wse_m, y = wse)) +
  geom_point(aes(color = factor(substr(reach_id, 1, 6))), size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", 
           x = min(time_space_matched_SWOT_PT$mean_reach_pt_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), 
                          "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_PT)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "Basin code")


# CDF plot
ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - mean_reach_pt_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
  annotate("text", x = 2, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 
# xlim(0, 1)

# ---------------------------------------------------------------------------------------------------------------------------
# remove bias from PT data

# removed median bias for individual reaches
# need to set a min threshold of obs for us to calc a bias (using 3 currently)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  group_by(reach_id) %>%
  mutate(
    bias = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    pt_wse_nobias_m = if (n() >= 3)
      mean_reach_pt_wse_m - bias
    else NA_real_) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (residuals)
time_space_matched_SWOT_PT$residuals_nobias = time_space_matched_SWOT_PT$pt_wse_nobias_m - time_space_matched_SWOT_PT$wse

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

#csv subset
save_to_csv <- time_space_matched_SWOT_PT %>%
  dplyr::select(time_utc,pt_time_UTC, reach_id, residuals, residuals_nobias, bias, wse, wse_u, 
                mean_reach_pt_wse_m, pt_wse_nobias_m, flaglist, sorted_nodelist, Number_of_nodes,
                slope, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, 
                area_det_u, area_wse, layovr_val, node_dist,
                xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, river, cycle_id, pass_id) #SWOTFileName, p_n_nodes OR #cycle_id, pass_id

save_to_csv <- save_to_csv %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP")

# save to csv
# write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_SWOT_PT.csv', row.names = FALSE)


# correlation test
cor_test_nobias <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_nobias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs PT wse
ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_nobias_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_PT$pt_wse_nobias_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value_nobias, 4), "\np value = ", signif(p_value_nobias, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_PT)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "River") 


# CDF plot
ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_nobias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
  annotate("text", x = 2, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 
#xlim(0, .751)


# ---------------------------------------------------------------------------------------------------------------------------
























# SLOPE
# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# RiverSP
# SWOT_reach_df <- read.csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# RiverTile
# SWORD v16
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_domain_reach_timeseries_v16.csv')
# SWORD v17b
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v17b/RiverTile_domain_reach_timeseries_v17b.csv')

# get ride of possible duplicates from hydrocron pull
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE)

# filter SWOT data by reach_q (0=good, 1=suspect, 2=degraded, 3=bad) & cross track distance
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) %>%
  filter(partial_f == 0) %>%
  filter(dark_frac <= 0.8)

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset

# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep PT reach data

# SWORD v16
# PT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16/YR_PT_reach_slope.csv')
# SWORD v17b
PT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b/YR_PT_reach_slope.csv')


# Convert time column to POSIXct
PT_reach_df$pt_time_UTC <- as.POSIXct(PT_reach_df$pt_time_UTC, format = "%m/%d/%y %H:%M", tz = "UTC")
#change reach_id so there aren't duplicate columns when merging with SWOT data
PT_reach_df <- rename(PT_reach_df, "PT_reach_id" = "reach_id")


# ---------------------------------------------------------------------------------------------------------------------------
# match PT & SWOT observations in time and space

# match PT and SWOT in time
# observation are matched by 7.5min buffer (a SWOT obs should always be within 7.5min of a PT)
time_matched_SWOT_PT_reach <- PT_reach_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(SWOT_reach_df_filtered %>%
                           filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5))) %>%
  unnest(closest_match) %>% 
  dplyr::select(everything())

# match PT and SWOT in space
# node level
time_space_matched_SWOT_PT_reach <- time_matched_SWOT_PT_reach %>%
  filter(PT_reach_id == reach_id)

# get rid of any duplicate PT values
time_space_matched_SWOT_PT_reach <- time_space_matched_SWOT_PT_reach[!duplicated(time_space_matched_SWOT_PT_reach[c("mean_pt_wse_us_boundary_m","pt_time_UTC","pt_serials_us")]),]

# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# fix negative slopes in SWOT & GNSS data
time_space_matched_SWOT_PT_reach$slope_abs <- abs(time_space_matched_SWOT_PT_reach$slope)
time_space_matched_SWOT_PT_reach$slope_m_m_abs <- abs(time_space_matched_SWOT_PT_reach$slope_m_m)


# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (residuals)
time_space_matched_SWOT_PT_reach$residuals = time_space_matched_SWOT_PT_reach$slope_m_m_abs - time_space_matched_SWOT_PT_reach$slope_abs

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT_reach$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT_reach$residuals), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error*100000))
print(paste("50th Percentile Error:", percentile_50_error*100000))


# removed median bias for individual PT
# need to set a min threshold of obs for us to calc a bias (using 3 currently)
time_space_matched_SWOT_PT_reach <- time_space_matched_SWOT_PT_reach %>%
  group_by(reach_id) %>%
  mutate(
    bias = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    mean_reach_PT_slope_no_bias_m_m = if (n() >= 3)
      slope_m_m_abs - bias
    else NA_real_) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_PT_reach$slope_residuals_nobias = time_space_matched_SWOT_PT_reach$mean_reach_PT_slope_no_bias_m_m - time_space_matched_SWOT_PT_reach$slope_abs

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT_reach$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT_reach$slope_residuals_nobias), 0.50, na.rm=TRUE)



#csv subset
save_to_csv <- time_space_matched_SWOT_PT_reach %>%
  select(pt_time_UTC, time_utc, reach_id, slope_m_m, slope_m_m_abs, slope_uncertainty_m_m, slope_residuals_nobias, p_lat, p_lon,
         slope, slope_abs, slope_u, bias, slope_residuals_nobias, mean_reach_PT_slope_no_bias_m_m, residuals, width, width_u, area_total, 
         area_tot_u, layovr_val, node_dist, xtrk_dist, reach_q, reach_q_b, dark_frac, xovr_cal_q, p_dist_out, n_good_nod
         ) #cycle_id, pass_id, OR p_dist_out, n_good_nod

save_to_csv <- save_to_csv %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile") %>%
  rename(slope_residuals = residuals)

# add river names to df
save_to_csv <- save_to_csv %>%
  mutate(river_code = substr(reach_id, 1, 6),
         river = case_when(
           # putting the reach id first ensures case_when won't overwrite SJ/BL labels
           # SWORD v16: "81260300061", "81260300231", "81260300241", "81260300251"
           # SWORD v17b: "81260300181", "81260300191", "81260300201", "81260300211"
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
           TRUE ~ NA_character_))

# save joined_wse_subset to csv
write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_PT.csv', row.names = FALSE)


# Calculating the linear regression model 
model = lm(slope_m_m~slope, data = time_space_matched_SWOT_PT_reach) 

# Extracting R-squared parameter from summary 
summary(model)

#RMSE
rmse <- sqrt(mean((time_space_matched_SWOT_PT_reach$slope_m_m - time_space_matched_SWOT_PT_reach$slope)^2))
#RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_PT_reach$slope, time_space_matched_SWOT_PT_reach$slope_m_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8",  
                   "#E3A700", "#008F7A", "#C83232", "#2E7D32",  
                   "#D81B60", "#00429D", "#A6761D", "#56B4E9", 
                   "#4c64c1","#7ca92f","#9a3c9a", 'orange',
                   'lightyellow', 'yellow', 'pink', 'black')

# plot SWOT vs PT slope
ggplot(time_space_matched_SWOT_PT_reach, aes(x = mean_reach_PT_slope_no_bias_m*100000, y = slope_abs*100000, color = factor(reach_id))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_PT_reach$slope_m_m*100000, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_PT_reach$slope*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_PT_reach)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "Reach ID") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

# CDF plot
ggplot(time_space_matched_SWOT_PT_reach, aes(x = abs(slope_abs - slope_m_m_abs))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (m/m)", y = "Cumulative Probability", title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 0.00018, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error*100000, 8)), color = "#222222", size = 6) +
  annotate("text", x = 0.00018, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error*100000, 8)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) +
  xlim(0, .00025)

# CDF plot no bias
ggplot(time_space_matched_SWOT_PT_reach, aes(x = abs(slope_abs - mean_reach_PT_slope_no_bias_m)*100000)) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (cm/km)", y = "Cumulative Probability", title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 7.4, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias*100000, 8)), color = "#222222", size = 6) +
  annotate("text", x = 7.4, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias*100000, 8)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# plot SWOT vs PT slope
ggplot(time_space_matched_SWOT_PT_reach, aes(x = slope_uncertainty_m_m, y = abs(residuals), color = factor(reach_id))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("pt slope uncertainty") +
  ylab("abs slope difference") +
  theme_minimal(base_size = 30) +
  labs(color = "Reach ID")






# ---------------------------------------------------------------------------------------------------------------------------
# Compare PT & SWOT riverSP/RiverTile reach slope (SWOT Version C vs D)
# ---------------------------------------------------------------------------------------------------------------------------

time_space_matched_riverSP_PT <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/reach_slope_RiverSP_time_space_matched_SWOT_PT.csv")
time_space_matched_riverTile_PT <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/reach_slope_RiverTile_time_space_matched_SWOT_PT.csv")


# ---------------------------------------------------------------------------------------------------------------------------
# Combine SWOT Version C & D results to plot
# ---------------------------------------------------------------------------------------------------------------------------
# label each verion
RiverSP_df <- time_space_matched_riverSP_PT %>%
  mutate(source = "RiverSP")

RiverTile_df <- time_space_matched_riverTile_PT %>%
  mutate(source = "RiverTile")

# Combine both dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_df, RiverTile_df)



# # Compare the same data subset from Version C & D
RiverTile_filtered <- RiverTile_df %>%
  semi_join(
    RiverSP_df %>% select(pt_time_UTC, reach_id, slope_m_m),
    by = c("pt_time_UTC", "reach_id", "slope_m_m")
  )

RiverSP_filtered <- RiverSP_df %>%
  semi_join(
    RiverTile_df %>% select(pt_time_UTC, reach_id, slope_m_m),
    by = c("pt_time_UTC", "reach_id", "slope_m_m")
  )

# Combine both filtered dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_filtered, RiverTile_filtered)


summary <- group_by(RiverSP_filtered, reach_id) %>% summarise(
  count = n(),
  mean = mean(abs(residuals), na.rm = TRUE),
  sd = sd(abs(residuals), na.rm = TRUE),
  median = median(abs(residuals), na.rm = TRUE),
  IQR = IQR(abs(residuals), na.rm= TRUE),
  min =min(abs(residuals), na.rm = TRUE),
  max =max(abs(residuals), na.rm= TRUE)
)


# ---------------------------------------------------------------------------------------------------------------------------
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_df$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_df$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_df$slope_residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_df$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_df$slope_residuals_nobias), 0.50, na.rm=TRUE)




# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_filtered$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_filtered$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_filtered$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_filtered$slope_residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$slope_residuals_nobias), 0.50, na.rm=TRUE)



# correlation test
cor_test <- cor.test(RiverTile_df$SWOT_slope, RiverTile_df$slope_m_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001



# Combo CDF plot no bias
ggplot(SWOT_versionCD_df, aes(x = abs(residuals*100000), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (cm/km)", y = "Cumulative Probability", 
       title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 10, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error*100000, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile*100000, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 10, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error*100000, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile*100000, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))


# Combo CDF plot no bias
ggplot(SWOT_versionCD_df, aes(x = abs(slope_residuals_nobias*100000), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (cm/km)", y = "Cumulative Probability", 
       title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 7.5, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error_nobias*100000, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile_nobias*100000, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 7.5, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error_nobias*100000, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile_nobias*100000, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))


