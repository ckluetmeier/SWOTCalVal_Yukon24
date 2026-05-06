library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare GNSS & SWOT (RiverSP/RiverTile) reach wse & slope
# ---------------------------------------------------------------------------------------------------------------------------

# contents:
# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# read in & prep GNSS data
# match GNSS & SWOT observations in time and space
# WSE
# SLOPE
# save dataframes
# WSE and slope comparisons across SWOT/SWORD versions


# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data

# RiverSP (SWORD v16)
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')
# SWORD v17b
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v17b/RiverSP_domain_reach_timeseries_PGD0_v17b.csv')


# RiverTile
# SWORD v16
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_domain_reach_timeseries_v16.csv')
# SWORD v17b
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v17b/RiverTile_domain_reach_timeseries_v17b.csv')

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
# read in & prep GNSS data

# SWORD v16
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_wse_slope.csv')
GNSS_slope_corr_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_slope_corrected.csv') %>%
  rename(reach_drift_slope_m_m_abs_nobias = slope_m_m)

# SWORD v17b
# GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv')
# GNSS_slope_corr_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_slope_corrected.csv') %>%
#   rename(reach_drift_slope_m_m_abs_nobias = slope_m_m)

GNSS_df <- GNSS_df %>%
  inner_join(GNSS_slope_corr_df %>%
               select(reach_drift_slope_m_m_abs_nobias, drift_id, reach_id),
             by = c(
               "reach_id",
               "drift_id"))


# Convert times to POSIXct
GNSS_df$wse_drift_start_UTC <- as.POSIXct(GNSS_df$wse_drift_start_UTC, tz = "UTC")
GNSS_df$wse_drift_end_UTC <- as.POSIXct(GNSS_df$wse_drift_end_UTC, tz = "UTC")

# Compute midpoint time of the drift for matching to SWOT overpass
GNSS_df$wse_drift_midpoint_UTC <- as.POSIXct(
  (as.numeric(GNSS_df$wse_drift_start_UTC) + as.numeric(GNSS_df$wse_drift_end_UTC)) / 2, origin = "1970-01-01", tz = "UTC")

# Compute total time of drift
GNSS_df$wse_drift_total_time_UTC <- difftime(GNSS_df$wse_drift_end_UTC, GNSS_df$wse_drift_start_UTC, units = "mins")

#change reach_id so there aren't duplicate columns when merging with SWOT data
GNSS_df <- rename(GNSS_df, "GNSS_reach_id" = "reach_id")

# ---------------------------------------------------------------------------------------------------------------------------
# match GNSS & SWOT observations in time and space

# match GNSS and SWOT in time
# observation are matched by 5 hour buffer
time_matched_SWOT_GNSS <- GNSS_df %>%
  rowwise() %>%
  mutate(closest_match = list(SWOT_reach_df_filtered %>%
                           filter(abs(difftime(wse_drift_midpoint_UTC, time_utc, units = "hours")) <= 5))) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# match GNSS and SWOT in space
# reach level
time_space_matched_SWOT_GNSS <- time_matched_SWOT_GNSS %>%
  filter(GNSS_reach_id == reach_id)


# ---------------------------------------------------------------------------------------------------------------------------
# WSE
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$residuals = time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m - time_space_matched_SWOT_GNSS$wse

summary <- group_by(time_space_matched_SWOT_GNSS, drift_id) %>% summarise(
  count = n(),
  mean = mean(abs(residuals), na.rm = TRUE),
  sd = sd(abs(residuals), na.rm = TRUE),
  median = median(abs(residuals), na.rm = TRUE),
  IQR = IQR(abs(residuals), na.rm= TRUE),
  min =min(abs(residuals), na.rm = TRUE),
  max =max(abs(residuals), na.rm= TRUE)
)

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#e73582","#78e964","#5329a8","#e1df28","#8c66f0","#a8dc3e","#7150ce","#48b32e","#cb4bd0",
                   "#5ddd76","#df36ad","#4e227d","#c7e772","#6549a4","#d4c841","#4c64c1","#7ca92f","#9a3c9a",
                   "#409f46","#b67ce0","#95d47e","#bc4797","#df3542","#6389ee","#d99427","#dd80d5","#417a28",
                   "#de6096","#982161","#c87735","#c9405b","#de4b21","#b84b34", "pink", "blue", "lightblue",
                   "darkblue", "darkred")

# plot SWOT vs GNSS wse
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_reach_drift_wse_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")
#labs(color = "Drift ID") 

# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_reach_drift_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 1, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 1, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)
#xlim(0, 3)

# plot residuals vs GNSS reach total error
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(residuals), y = mean_reach_drift_wse_total_error_m, color = factor(substr(drift_id,70,93)))) + #color = factor(substr(reach_id, 1, 6)
  geom_point(size = 1.5) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS - SWOT wse (m)") +
  ylab("GNSS total error (m)") +
  theme_minimal(base_size = 30)  +
  theme(legend.position = "none")
#labs(color = "Drift ID") 

# ---------------------------------------------------------------------------------------------------------------------------
# remove bias from GNSS data

# need to set a min threshold of obs for us to calc a bias (using 3 currently)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  group_by(drift_id) %>%
  mutate(
    bias = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    mean_reach_drift_wse_no_bias_m = if (n() >= 3)
      mean_reach_drift_wse_m - bias
    else
      NA_real_) %>%
  ungroup()

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$residuals_nobias = time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.50, na.rm=TRUE)

print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# correlation test
cor_test_nobias <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs GNSS wse
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_reach_drift_wse_no_bias_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")

# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_reach_drift_wse_no_bias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 0.5, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 0.5, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 


# plot of bias over time (by drift ID)
ggplot(time_space_matched_SWOT_GNSS, aes(x = time_utc, y = bias, color = factor(substr(drift_id,70,93)))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("time") +
  ylab("bias (m)") +
  theme_minimal(base_size = 30)  +
  theme(legend.position = "none") 
  #labs(color = "Drift ID") 

# ---------------------------------------------------------------------------------------------------------------------------
# save dataframes
# ---------------------------------------------------------------------------------------------------------------------------

# add river names to df
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
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
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PIC0") ## CHANGE TO CORRECT VERSION!


# csv subset
# # RiverSP
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  dplyr::select(reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC, wse_drift_midpoint_UTC, wse_drift_total_time_UTC, residuals, residuals_nobias, bias, mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
                mean_reach_drift_wse_no_bias_m, reach_drift_slope_m_m, reach_drift_slope_precision_m, drift_id, wse, wse_u,
                slope, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, area_det_u, area_wse, layovr_val, node_dist,
                xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, cycle_id, pass_id,
                river_code, river, insitu_type, source) # cycle_id, pass_id

# # RiverTile
# save_to_csv <- time_space_matched_SWOT_GNSS %>%
#   dplyr::select(reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC, wse_drift_midpoint_UTC, wse_drift_total_time_UTC, residuals, residuals_nobias, bias, mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
#                 mean_reach_drift_wse_no_bias_m, reach_drift_slope_m_m, reach_drift_slope_precision_m, drift_id, wse, wse_u,
#                 slope, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, area_det_u, area_wse, layovr_val, node_dist,
#                 xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, p_n_nodes, SWOTFileName,
#                 river_code, river, insitu_type, source) # SWOTFileName, p_n_nodes


# save joined_wse_subset to csv
# write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_wse_SWOT_GNSS.csv', row.names = FALSE)














# ---------------------------------------------------------------------------------------------------------------------------
# SLOPE
# ---------------------------------------------------------------------------------------------------------------------------

# fix negative slopes in SWOT & GNSS data
time_space_matched_SWOT_GNSS$slope_abs <- abs(time_space_matched_SWOT_GNSS$slope)
time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs <- abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m)

# remove near 0 GNSS slopes (0.1 cm/km)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  filter(abs(reach_drift_slope_m_m) > 0.000001)

# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT -GNSS(slope_residuals)
time_space_matched_SWOT_GNSS$slope_residuals = time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs - time_space_matched_SWOT_GNSS$slope_abs

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.50, na.rm=TRUE)

#print the result
# *100000 puts m/m into cm/km for slopes
print(paste("68th Percentile Error:", percentile_68_error*100000))
print(paste("50th Percentile Error:", percentile_50_error*100000))

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_GNSS$slope_abs, time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

# color_palette <- c("#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8",
#                    "#E3A700", "#008F7A", "#C83232", "#2E7D32",
#                    "#D81B60", "#00429D", "#A6761D", "#56B4E9")

# plot uncorrected SWOT vs GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_m_m*100000, y = slope*100000)) +
  geom_point(size = 4) +
  xlab("GNSS slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

# plot abs SWOT vs GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(reach_drift_slope_m_m_abs*100000), y = abs(slope_abs*100000))) +
  geom_point(size = 4) +
  xlab("GNSS slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs*100000), na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$slope_abs*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),  "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(slope_abs*100000 - reach_drift_slope_m_m_abs*100000))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - GNSS slope (cm/km)", y = "Cumulative Probability", title = "CDF of SWOT slope - GNSS slope") +
  annotate("text", x = 100, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error*100000, 4)), color = "#222222", size = 6) +
  annotate("text", x = 100, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error*100000, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# plot residuals vs GNSS slope uncertainty
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_precision_m*100000, y = abs(slope_residuals*100000), color = factor(reach_id))) +
  geom_point(size = 4) +
  #scale_color_manual(values = color_palette) +
  xlab("GNSS slope uncertainty") +
  ylab("abs slope difference") +
  theme_minimal(base_size = 30) +
  labs(color = "Reach ID")

# ---------------------------------------------------------------------------------------------------------------------------
# remove bias from GNSS slope data

# # need to set a min threshold of obs for us to calc a bias (using 3 currently)
# time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
#   group_by(drift_id) %>%
#   mutate(
#     bias_slope = if (n() >= 3) median(slope_residuals, na.rm = TRUE) else NA_real_,
#     reach_drift_slope_m_m_abs_nobias = if (n() >= 3)
#       reach_drift_slope_m_m_abs - bias_slope
#     else
#       NA_real_) %>%
#   ungroup()

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$slope_residuals_nobias = time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias - time_space_matched_SWOT_GNSS$slope_abs


# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals_nobias), 0.50, na.rm=TRUE)

print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias*100000))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias*100000))


# correlation test
cor_test_nobias <- cor.test(time_space_matched_SWOT_GNSS$slope_abs, time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_m_m_abs_nobias*100000, y = slope_abs*100000, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias*100000, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$slope_abs*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")

# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(slope_abs*100000 - reach_drift_slope_m_m_abs_nobias*100000))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - GNSS slope (cm/km)", y = "Cumulative Probability", title = "CDF of SWOT slope - GNSS slope") +
  annotate("text", x = 1, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias*100000, 4)), color = "#222222", size = 6) +
  annotate("text", x = 1, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias*100000, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 




# ---------------------------------------------------------------------------------------------------------------------------
# save dataframes
# ---------------------------------------------------------------------------------------------------------------------------

# csv subset
# RiverSP
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  dplyr::select(reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC, wse_drift_midpoint_UTC, wse_drift_total_time_UTC, residuals, residuals_nobias, bias, mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
                mean_reach_drift_wse_no_bias_m, reach_drift_slope_m_m, reach_drift_slope_m_m_abs, reach_drift_slope_precision_m, slope_residuals, slope_residuals_nobias, reach_drift_slope_m_m_abs_nobias, drift_id, wse, wse_u,
                slope, slope_abs, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, area_det_u, area_wse, layovr_val, node_dist,
                xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, cycle_id, pass_id,
                river_code, river, insitu_type, source) # cycle_id, pass_id

# RiverTile
# save_to_csv <- time_space_matched_SWOT_GNSS %>%
#   dplyr::select(reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC, wse_drift_midpoint_UTC, wse_drift_total_time_UTC, residuals, residuals_nobias, bias, mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
#                 mean_reach_drift_wse_no_bias_m, reach_drift_slope_m_m, reach_drift_slope_m_m_abs, reach_drift_slope_precision_m, slope_residuals, slope_residuals_nobias, bias_slope, reach_drift_slope_m_m_abs_nobias, drift_id, wse, wse_u,
#                 slope, slope_abs, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, area_det_u, area_wse, layovr_val, node_dist,
#                 xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, p_n_nodes, SWOTFileName,
#                 river_code, river, insitu_type, source) # SWOTFileName, p_n_nodes


# save joined_wse_subset to csv
write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv', row.names = FALSE)













