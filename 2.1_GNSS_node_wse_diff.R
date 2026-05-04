library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare GNSS wse & SWOT riverSP node wse
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# read in SWOT data

# RiverSP
# Version C
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv')
# Version D
SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv')

# RiverTile
# SWORD v16
# SWOT_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v16/RiverTile_domain_node_timeseries_v16.csv")
# SWORD v17b
# SWOT_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_domain_node_timeseries_v17b.csv")

# get ride of possible duplicates from hydrocron pull
# this also filters out bad nodes without data (e.g. time = -999999999999, wse = -1.000000e+12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# filter SWOT data
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>% # (0=good, 1=suspect, 2=degraded, 3=bad)
  filter(abs(xtrk_dist) >= 10000) %>% # xtrk dist should be 10-60km
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac <= 0.80) # dark water less than 80%

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep GNSS data

# SWORD v16
# GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_node_wses.csv')
# SWORD v17b
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv')

GNSS_df$time_UTC <- as.POSIXct(GNSS_df$time_UTC, tz = "UTC")

#change Node_ID so there aren't duplicate columns when merging with SWOT data
GNSS_df <- rename(GNSS_df, "Node_ID" = "node_id")

# filter GNSS by node total error in m
# GNSS_df <- GNSS_df %>%
#   filter(node_total_error_m < 2.5)

# ---------------------------------------------------------------------------------------------------------------------------
# match GNSS & SWOT observations in time and space

# match GNSS and SWOT in time
# observation are matched by 5 hour buffer
time_matched_SWOT_GNSS <- GNSS_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(SWOT_df_filtered %>%
                           filter(abs(difftime(time_UTC, time_utc, units = "hours")) <= 5))
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# match GNSS and SWOT in space
# node level
time_space_matched_SWOT_GNSS <- time_matched_SWOT_GNSS %>%
  filter(Node_ID == node_id)


# ---------------------------------------------------------------------------------------------------------------------------
# data retention exploration
# 
# # try
# SWOT_df_filtered <- SWOT_df_filtered[, c("time_utc", "node_id", "lat", "lon")]

# # count by node
# matched_all <- group_by(time_space_matched_SWOT_GNSS, node_id) %>% summarise(
#   count = n(),
#   lat = median(lat),
#   lon = median(lon))
# 
# # matched_filtered <- group_by(time_space_matched_SWOT_GNSS, node_id) %>% summarise(
# #   count_filtered = n())
# 
# 
# data_retention <- matched_all %>%
#   left_join(matched_filtered, by = "node_id") %>%
#   mutate(across(everything(), ~ coalesce(.x, 0)))
# 
# data_retention$percent_kept = as.numeric(data_retention$count_filtered / data_retention$count)
# 
# 
# write.csv(data_retention, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/vC_data_retention.csv', row.names = FALSE)


# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$residuals = time_space_matched_SWOT_GNSS$mean_node_drift_wse_m - time_space_matched_SWOT_GNSS$wse

summary <- group_by(time_space_matched_SWOT_GNSS, drift_id) %>% summarise(
  count = n(),
  mean = mean(abs(residuals), na.rm = TRUE),
  sd = sd(abs(residuals), na.rm = TRUE),
  median = median(abs(residuals), na.rm = TRUE),
  IQR = IQR(abs(residuals), na.rm= TRUE),
  min =min(abs(residuals), na.rm = TRUE),
  max =max(abs(residuals), na.rm= TRUE)
)

# Filter values to sensical residuals (< 3 m diff)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  filter(abs(residuals) < 3) # %>%
  # filter(node_total_error_m < 1)

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.50, na.rm=TRUE)

#print the result
# print(paste("68th Percentile Error:", percentile_68_error))
# print(paste("50th Percentile Error:", percentile_50_error))

# # Calculating the linear regression model 
# model = lm(mean_node_drift_wse_m~wse, data = time_space_matched_SWOT_GNSS) 
# # Extracting R-squared parameter from summary 
# summary(model)
# #RMSE
# rmse <- sqrt(mean((time_space_matched_SWOT_GNSS$mean_node_drift_wse_m - time_space_matched_SWOT_GNSS$wse)^2))
# #RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_node_drift_wse_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#e73582","#78e964","#5329a8","#e1df28","#8c66f0","#a8dc3e","#7150ce","#48b32e","#cb4bd0",
                   "#5ddd76","#df36ad","#4e227d","#c7e772","#6549a4","#d4c841","#4c64c1","#7ca92f","#9a3c9a",
                   "#409f46","#b67ce0","#95d47e","#bc4797","#df3542","#6389ee","#d99427","#dd80d5","#417a28",
                   "#de6096","#982161","#c87735","#c9405b","#de4b21","#b84b34", "pink", "blue", "lightblue",
                   "darkblue", "darkred", 'gray')

# plot SWOT vs GNSS wse
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_node_drift_wse_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_GNSS$mean_node_drift_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")
  #labs(color = "Drift ID") 


# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_node_drift_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 2, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)
  #xlim(0, 3)

# plot residuals vs GNSS node total error
# ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(residuals), y = node_total_error_m, color = factor(substr(drift_id,70,93)))) + #color = factor(substr(reach_id, 1, 6)
#   geom_point(size = 1.5) +
#   scale_color_manual(values = color_palette) +
#   xlab("GNSS - SWOT wse (m)") +
#   ylab("wse uncertainty (m)") +
#   theme_minimal(base_size = 30)  +
#   theme(legend.position = "none")
#   #labs(color = "Drift ID") 

# ---------------------------------------------------------------------------------------------------------------------------
# remove bias from GNSS data

# need to set a min threshold of obs for us to calc a bias (using 3 currently)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  group_by(drift_id) %>%
  mutate(
    bias = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    mean_node_drift_wse_no_bias_m = if (n() >= 3)
      mean_node_drift_wse_m - bias
    else
      NA_real_) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$residuals_nobias = time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))


# need to label each version & add river id

# label each verion
node_SWOT_GNSS_vC <- time_space_matched_SWOT_GNSS %>%
  mutate(source = "RiverSP")

node_SWOT_GNSS_vD <- time_space_matched_SWOT_GNSS %>%
  mutate(source = "RiverTile")

# add river names to df
node_SWOT_GNSS_vD <- node_SWOT_GNSS_vD %>%
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
           TRUE ~ NA_character_))
# Check to make sure rivers are labeled as expected
# node_SWOT_GNSS_vC %>%
#   filter(reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141", "81270100151", "81270100161", "81270200011", "81270200021")) %>%
#   select(reach_id, river_code, river)


# csv subset
# SWOT version D
# save_to_csv <- node_SWOT_GNSS_vD %>%
#   dplyr::select(time_utc,time_UTC, residuals, residuals_nobias, bias, wse, wse_u,
#                 mean_node_drift_wse_m, mean_node_drift_wse_no_bias_m, node_total_error_m, drift_id,
#                 width, width_u, node_id, reach_id, p_dist_out,
#                 node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q, lat, lon, source, river_code, river)

# SWOT version C
save_to_csv <- node_SWOT_GNSS_vD %>%
  dplyr::select(time_utc,time_UTC, residuals, residuals_nobias, bias, wse, wse_u,
                mean_node_drift_wse_m, mean_node_drift_wse_no_bias_m, node_total_error_m, drift_id,
                width, width_u, node_id, reach_id, p_dist_out,
                node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q, cycle_id, pass_id, lat, lon, source, river_code, river)

# save joined_wse_subset to csv
write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv', row.names = FALSE)


# Calculating the linear regression model 
# model_nobias = lm(mean_node_drift_wse_no_bias_m~wse, data = time_space_matched_SWOT_GNSS) 
# 
# Extracting R-squared parameter from summary 
# summary(model_nobias)
# 
# RMSE
# rmse_nobias <- sqrt(mean((time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse)^2))
# RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test_nobias <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs GNSS wse
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_node_drift_wse_no_bias_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")


# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_node_drift_wse_no_bias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 2, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 






# ---------------------------------------------------------------------------------------------------------------------------
# Compare GNSS & SWOT riverSP/RiverTile node wse (SWOT Version C vs D)
# ---------------------------------------------------------------------------------------------------------------------------

# OLD
# time_space_matched_riverSP_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/time_space_matched_SWOT_GNSS_5mdiff.csv")
# time_space_matched_riverTile_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v16/time_space_matched_SWOT_GNSS_5mdiff.csv")

time_space_matched_riverSP_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/matched_SWOT_GNSS_3mdiff.csv")
# time_space_matched_riverTile_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v16/matched_SWOT_GNSS_3mdiff.csv")
time_space_matched_riverTile_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/matched_SWOT_GNSS_3mdiff.csv")

# ---------------------------------------------------------------------------------------------------------------------------
# Combine SWOT Version C & D results to plot
# ---------------------------------------------------------------------------------------------------------------------------
# label each verion
RiverSP_df <- time_space_matched_riverSP_GNSS %>%
  mutate(source = "RiverSP")

RiverTile_df <- time_space_matched_riverTile_GNSS %>%
  mutate(source = "RiverTile")

# Combine both dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_df, RiverTile_df)




summary <- group_by(RiverTile_filtered, node_id) %>% summarise(
  count = n(),
  mean = mean(abs(residuals_nobias), na.rm = TRUE),
  sd = sd(abs(residuals_nobias), na.rm = TRUE),
  median = median(abs(residuals_nobias), na.rm = TRUE),
  IQR = IQR(abs(residuals_nobias), na.rm= TRUE),
  min =min(abs(residuals_nobias), na.rm = TRUE),
  max =max(abs(residuals_nobias), na.rm= TRUE),
  lat =mean(lat),
  lon =mean(lon)
)

# summary node dataframe -- export to visualize in QGIS
#write.csv(summary, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/node_summary_time_space_matched_SWOT_GNSS_5mdiff_nobias.csv', row.names = FALSE)


# ---------------------------------------------------------------------------------------------------------------------------
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_df$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_df$residuals_nobias), 0.50, na.rm=TRUE)



# correlation test
cor_test <- cor.test(RiverTile_df$wse, RiverTile_df$mean_node_drift_wse_m) # mean_node_drift_wse_m
# cor_test <- cor.test(RiverSP_df$wse, RiverSP_df$mean_node_drift_wse_no_bias_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# Mean Absolute Error (which is just the mean residual)
MAE <- mean(abs(RiverTile_df$residuals))
MAE <- mean(abs(RiverTile_df$residuals_nobias))


t.test(abs(RiverTile_df$residuals), abs(RiverSP_df$residuals))
t.test(abs(RiverTile_df$residuals_nobias), abs(RiverSP_df$residuals_nobias))



# # SUBSET TO SAME VERSION C/D data points
# # Calculate the 68th & 50th percentile error
# percentile_68_error <- quantile(abs(RiverSP_filtered$residuals), 0.68, na.rm=TRUE)
# percentile_50_error <- quantile(abs(RiverSP_filtered$residuals), 0.50, na.rm=TRUE)
# 
# percentile_68_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.68, na.rm=TRUE)
# percentile_50_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.50, na.rm=TRUE)
# 
# # Calculate the 68th &50th percentile error, no bias
# percentile_68_error_nobias <- quantile(abs(RiverSP_filtered$residuals_nobias), 0.68, na.rm=TRUE)
# percentile_50_error_nobias <- quantile(abs(RiverSP_filtered$residuals_nobias), 0.50, na.rm=TRUE)
# 
# percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$residuals_nobias), 0.68, na.rm=TRUE)
# percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$residuals_nobias), 0.50, na.rm=TRUE)



# ---------------------------------------------------------------------------------------------------------------------------
# Plots

# Combo CDF plot
ggplot(SWOT_versionCD_df, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 3, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 3, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "dashed"))


# Combo CDF plot no bias
ggplot(SWOT_versionCD_df, aes(x = abs(residuals_nobias), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 3, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error_nobias, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile_nobias, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 3, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error_nobias, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile_nobias, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))







# correlation test
cor_test <- cor.test(time_space_matched_riverTile_GNSS$wse, time_space_matched_riverTile_GNSS$mean_node_drift_wse_m)

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
ggplot(time_space_matched_riverTile_GNSS, aes(x = mean_node_drift_wse_no_bias_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_riverTile_GNSS$mean_node_drift_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_riverTile_GNSS$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_riverTile_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")
#labs(color = "Drift ID") 



