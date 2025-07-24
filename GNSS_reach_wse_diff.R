library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare GNSS wse & SWOT RiverTile reach wse/slope
# ---------------------------------------------------------------------------------------------------------------------------

# WSE
# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# RiverSP (SWORD v16)
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# RiverTile
# SWORD v16
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_domain_reach_timeseries_v16.csv')
# SWORD v17b
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v17b/RiverTile_domain_reach_timeseries_v17b.csv')

# get ride of possible duplicates from hydrocron pull
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE) %>%
  filter(time > 0) %>%
  filter(wse > 0)

# filter SWOT data by reach_q (0=good, 1=suspect, 2=degraded, 3=bad) & cross track distance
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) %>%
  filter(partial_f == 0)

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset

# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep GNSS data
# SWORD v16
# GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/SWORD_v16/YR_drift_reach_wse_slope.csv')
# SWORD v17b
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/SWORD_v17b/YR_drift_reach_wse_slope.csv')
  
GNSS_df$wse_drift_start_UTC <- as.POSIXct(GNSS_df$wse_drift_start_UTC, tz = "UTC")
GNSS_df$wse_drift_end_UTC <- as.POSIXct(GNSS_df$wse_drift_end_UTC, tz = "UTC")

#change reach_id so there aren't duplicate columns when merging with SWOT data
GNSS_df <- rename(GNSS_df, "GNSS_reach_id" = "reach_id")

# ---------------------------------------------------------------------------------------------------------------------------
# match GNSS & SWOT observations in time and space

# match GNSS and SWOT in time
# observation are matched by 5 hour buffer
time_matched_SWOT_GNSS <- GNSS_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(SWOT_reach_df_filtered %>%
                           filter(abs(difftime(wse_drift_end_UTC, time_utc, units = "hours")) <= 5))
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# match GNSS and SWOT in space
# reach level
time_space_matched_SWOT_GNSS <- time_matched_SWOT_GNSS %>%
  filter(GNSS_reach_id == reach_id)

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

# Calculating the linear regression model 
model = lm(mean_reach_drift_wse_m~wse, data = time_space_matched_SWOT_GNSS) 

# Extracting R-squared parameter from summary 
summary(model)

#RMSE
rmse <- sqrt(mean((time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m - time_space_matched_SWOT_GNSS$wse)^2))
#RMSE >= MAE, MAE is similar to 50th quantile error

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
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
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

time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  group_by(drift_id) %>%
  mutate(
    bias = median(residuals, na.rm = TRUE),
    mean_reach_drift_wse_no_bias_m = mean_reach_drift_wse_m - bias
  ) %>%
  ungroup()


# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
time_space_matched_SWOT_GNSS$residuals_nobias = time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

#csv subset
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  dplyr::select(reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC, residuals, residuals_nobias, bias, mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m, 
                mean_reach_drift_wse_no_bias_m, reach_drift_slope_m_m, reach_drift_slope_precision_m, drift_id, wse, wse_u, 
                slope, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct, area_det_u, area_wse, layovr_val, node_dist,
                xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q, p_dist_out, p_lat, p_lon, cycle_id, pass_id) # SWOTFileName, p_n_nodes,

# save joined_wse_subset to csv
# write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_time_space_matched_SWOT_GNSS.csv', row.names = FALSE)


# Calculating the linear regression model 
model_nobias = lm(mean_reach_drift_wse_no_bias_m~wse, data = time_space_matched_SWOT_GNSS) 

# Extracting R-squared parameter from summary 
summary(model_nobias)

#RMSE
rmse_nobias <- sqrt(mean((time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse)^2))
#RMSE >= MAE, MAE is similar to 50th quantile error

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
  annotate("text", x = 1, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 1, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
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






library(dplyr)

time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  mutate(
    reach_quality = factor(reach_q,
                     levels = c(0, 1, 2, 3),
                     labels = c("Good", "Suspect", "Degraded", "Bad"))
  )


# Step 1: Sort by absolute error and compute cumulative retention
pareto_df <- time_space_matched_SWOT_GNSS %>%
  arrange(abs(residuals_nobias)) %>%
  mutate(
    abs_error = abs(residuals_nobias),
    retention = seq_along(abs_error) / n() * 100  # % remaining from current threshold onward
  )

pareto_df <- time_space_matched_SWOT_GNSS %>%
  arrange(reach_q) %>%
  mutate(
    abs_error = abs(residuals_nobias),
    retention = seq_along(reach_q) / n() * 100  # % remaining from current threshold onward
  )

# Step 2: Compute Pareto front (unique thresholds and best retention)
pareto_front <- pareto_df %>%
  group_by(abs_error = round(abs_error, 5)) %>%  # round to avoid noise in x-axis
  summarize(retention = max(retention), .groups = "drop") %>%
  arrange(abs_error)




# Step 3: Plot
ggplot(pareto_df, aes(x = abs_error, y = retention, color = reach_quality)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(data = pareto_front, aes(x = abs_error, y = retention), 
            inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = c(
    "Good" = "#1f77b4", 
    "Suspect" = "#17becf", 
    "Degraded" = "#2ca02c", 
    "Bad" = "#d62728"
  )) +
  labs(
    x = "|68%ile| WSE difference (m)",
    y = "Percent of observations remaining (%)",
    title = "Pareto Front?? WSE Difference vs. Data Retention",
    color = "Overall Quality Flag"
  ) +
  theme_minimal(base_size = 18)




library(ggplot2)
d <- data.frame(x = abs(pareto_df$residuals_nobias), y = pareto_df$retention)
D <- d[order(d$x,d$y,decreasing=TRUE),]

pareto_front <- which(!duplicated(cummax(D$y)))
front <- D[pareto_front,]
other <- D[-pareto_front,]

ggplot(mapping = aes(x, y)) +
  geom_point(data = other) +
  geom_line(data = front, colour = 'red') +
  geom_point(data = front, colour = 'red') +
  lims(x = c(0, NA), y = c(0, NA))
















# ---------------------------------------------------------------------------------------------------------------------------
# Compare GNSS wse & SWOT riverSP/RiverTile reach wse
# ---------------------------------------------------------------------------------------------------------------------------

time_space_matched_riverSP_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_time_space_matched_SWOT_GNSS.csv")
time_space_matched_riverTile_GNSS <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_time_space_matched_SWOT_GNSS.csv")


# ---------------------------------------------------------------------------------------------------------------------------
# Combine SWOT Version C & D results to plot
# ---------------------------------------------------------------------------------------------------------------------------
# label each verion
RiverSP_df <- time_space_matched_riverSP_GNSS %>%
  mutate(source = "RiverSP") %>%
  mutate(abs_reach_drift_slope_m_m = abs(reach_drift_slope_m_m)) %>%
  mutate(abs_slope = abs(slope))

RiverTile_df <- time_space_matched_riverTile_GNSS %>%
  mutate(source = "RiverTile") %>%
  mutate(abs_reach_drift_slope_m_m = abs(reach_drift_slope_m_m)) %>%
  mutate(abs_slope = abs(slope))

# Combine both dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_df, RiverTile_df)



# Compare the same data subset from Version C & D
RiverTile_filtered <- RiverTile_df %>%
  semi_join(
    RiverSP_df %>% select(reach_id, wse_drift_start_UTC, mean_reach_drift_wse_m),
    by = c("reach_id", "wse_drift_start_UTC", "mean_reach_drift_wse_m")
  )

RiverSP_filtered <- RiverSP_df %>%
  semi_join(
    RiverTile_df %>% select(reach_id, wse_drift_start_UTC, mean_reach_drift_wse_m),
    by = c("reach_id", "wse_drift_start_UTC", "mean_reach_drift_wse_m")
  )

# Combine both filtered dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_filtered, RiverTile_filtered)




summary <- group_by(RiverTile_filtered, reach_id) %>% summarise(
  count = n(),
  mean = mean(abs(residuals_nobias), na.rm = TRUE),
  sd = sd(abs(residuals_nobias), na.rm = TRUE),
  median = median(abs(residuals_nobias), na.rm = TRUE),
  IQR = IQR(abs(residuals_nobias), na.rm= TRUE),
  min =min(abs(residuals_nobias), na.rm = TRUE),
  max =max(abs(residuals_nobias), na.rm= TRUE),
  lat =mean(p_lat),
  lon =mean(p_lon)
)

#write.csv(summary, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/reach_summary_time_space_matched_SWOT_GNSS_nobias.csv', row.names = FALSE)



# ---------------------------------------------------------------------------------------------------------------------------
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(time_space_matched_riverSP_GNSS$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_riverSP_GNSS$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(time_space_matched_riverTile_GNSS$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(time_space_matched_riverTile_GNSS$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(time_space_matched_riverSP_GNSS$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_riverSP_GNSS$residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(time_space_matched_riverTile_GNSS$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(time_space_matched_riverTile_GNSS$residuals_nobias), 0.50, na.rm=TRUE)




# SUBSET TO SAME VERSION C/D data points
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_filtered$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_filtered$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_filtered$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_filtered$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_filtered$residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$residuals_nobias), 0.50, na.rm=TRUE)



# ---------------------------------------------------------------------------------------------------------------------------
# Plots

# Combo CDF plot
ggplot(SWOT_versionCD_df, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - GNSS WSE (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT WSE - GNSS WSE") +
  annotate("text", x = 0.75, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 0.75, y = 0.53, 
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
  annotate("text", x = 0.8, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error_nobias, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile_nobias, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 0.8, y = 0.53, 
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







# ---------------------------------------------------------------------------------------------------------------------------
# Slope comparisons
# ---------------------------------------------------------------------------------------------------------------------------
# Calculate the wse diff SWOT - PT (slope_residuals)
RiverSP_df$slope_residuals = RiverSP_df$abs_reach_drift_slope_m_m - RiverSP_df$abs_slope


RiverSP_df <- RiverSP_df %>%
  group_by(drift_id) %>%
  mutate(
    bias = median(slope_residuals, na.rm = TRUE),
    mean_abs_reach_drift_slope_no_bias_m = abs_reach_drift_slope_m_m - bias
  ) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
RiverSP_df$slope_residuals_nobias = RiverSP_df$mean_abs_reach_drift_slope_no_bias_m - RiverSP_df$abs_slope


# Calculate the wse diff SWOT - PT (slope_residuals)
RiverTile_df$slope_residuals = RiverTile_df$abs_reach_drift_slope_m_m - RiverTile_df$abs_slope


RiverTile_df <- RiverTile_df %>%
  group_by(drift_id) %>%
  mutate(
    bias = median(slope_residuals, na.rm = TRUE),
    mean_abs_reach_drift_slope_no_bias_m = abs_reach_drift_slope_m_m - bias
  ) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - GNSS (residuals)
RiverTile_df$slope_residuals_nobias = RiverTile_df$mean_abs_reach_drift_slope_no_bias_m - RiverTile_df$abs_slope


# Combine both dataframes
SWOT_versionCD_df <- bind_rows(RiverSP_df, RiverTile_df)

summary <- group_by(RiverSP_df, reach_id) %>% summarise(
  count = n(),
  mean = mean(abs(slope_residuals_nobias), na.rm = TRUE),
  sd = sd(abs(slope_residuals_nobias), na.rm = TRUE),
  median = median(abs(slope_residuals_nobias), na.rm = TRUE),
  IQR = IQR(abs(slope_residuals_nobias), na.rm= TRUE),
  min =min(abs(slope_residuals_nobias), na.rm = TRUE),
  max =max(abs(slope_residuals_nobias), na.rm= TRUE)
)


# ---------------------------------------------------------------------------------------------------------------------------
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_df$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_df$slope_residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$slope_residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_df$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_df$slope_residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_df$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_df$slope_residuals_nobias), 0.50, na.rm=TRUE)




# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_filtered$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_filtered$slope_residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_filtered$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_filtered$slope_residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_filtered$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_filtered$slope_residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_filtered$slope_residuals_nobias), 0.50, na.rm=TRUE)




# correlation test
cor_test <- cor.test(RiverTile_df$slope, RiverTile_df$reach_drift_slope_m_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# Plots

RiverTile_df <- RiverTile_df %>%
  filter(abs(slope_residuals_nobias) < 15/100000)
RiverSP_df <- RiverSP_df %>%
  filter(abs(slope_residuals_nobias) < 15/100000)

SWOT_versionCD_df <- SWOT_versionCD_df %>%
  filter(abs(slope_residuals_nobias) < 15/100000)


# add river names to df
RiverTile_df <- RiverTile_df %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      river_code == "812701" ~ "lower_YR",
      river_code == "812705" ~ "upper_YR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "upper_PR",
      river_code == "812605" ~ "upper_PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

# for rivers
color_palette <- c("#D86A1A", "#F8A31B", "#00429D", "#2E7D32", "#6D398B",
                   "#00429D",  "#C83232", "#008F7A", "#E3A700", "#124000")
# plot SWOT vs PT slope
ggplot(RiverTile_df, aes(x = reach_drift_slope_m_m*100000, y = slope*100000, color = factor(river))) +
  geom_point(size = 4) +
  xlab("PT slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  scale_color_manual(values = color_palette) +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(abs(RiverTile_df$abs_reach_drift_slope_m_m*100000), na.rm = TRUE), 
           y = max(RiverTile_df$abs_slope*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),  "\nn = ", nrow(RiverTile_df)),
           hjust = 0, vjust = 1, size = 8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) 

ggplot(SWOT_versionCD_df, aes(x = abs(slope_residuals*100000), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (cm/km)", y = "Cumulative Probability", 
       title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 7, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error*100000, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile*100000, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 7, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error*100000, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile*100000, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "dashed"))


# Combo CDF plot no bias
ggplot(SWOT_versionCD_df, aes(x = abs(slope_residuals_nobias*100000), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - GNSS slope (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT slope - GNSS slope") +
  annotate("text", x = 7, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error_nobias*100000, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile_nobias*100000, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 7, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error_nobias*100000, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile_nobias*100000, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))

















































# SLOPE
# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (slope_residuals)
time_space_matched_SWOT_GNSS$slope_residuals = time_space_matched_SWOT_GNSS$reach_drift_slope_m_m - time_space_matched_SWOT_GNSS$slope

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error*100000))
print(paste("50th Percentile Error:", percentile_50_error*100000))

#csv subset
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  select(time_utc, slope_residuals, slope, reach_drift_slope_m_m, reach_id)

# save joined_wse_subset to csv
#write.csv(time_space_matched_SWOT_GNSS, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/test.csv', row.names = FALSE)


# Calculating the linear regression model 
model = lm(reach_drift_slope_m_m~slope, data = time_space_matched_SWOT_GNSS) 

# Extracting R-squared parameter from summary 
summary(model)

#RMSE
rmse <- sqrt(mean((time_space_matched_SWOT_GNSS$reach_drift_slope_m_m - time_space_matched_SWOT_GNSS$slope)^2))
#RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_GNSS$slope, time_space_matched_SWOT_GNSS$reach_drift_slope_m_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

# color_palette <- c("#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8",  
#                    "#E3A700", "#008F7A", "#C83232", "#2E7D32",  
#                    "#D81B60", "#00429D", "#A6761D", "#56B4E9")

# plot SWOT vs PT slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_m_m*100000, y = slope*100000)) +
  geom_point(size = 4) +
  xlab("PT slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m*100000), na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$slope*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),  "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)


ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(reach_drift_slope_m_m*100000), y = abs(slope*100000))) +
  geom_point(size = 4) +
  xlab("PT slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m*100000), na.rm = TRUE), 
           y = max(time_space_matched_SWOT_GNSS$slope*100000, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),  "\nn = ", nrow(time_space_matched_SWOT_GNSS)),
           hjust = 0, vjust = 1, size = 8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

# CDF plot
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(slope*100000 - reach_drift_slope_m_m*100000))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT slope - PT slope (cm/km)", y = "Cumulative Probability", title = "CDF of SWOT slope - PT slope") +
  annotate("text", x = 100, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error*100000, 4)), color = "#222222", size = 6) +
  annotate("text", x = 100, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error*100000, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# plot SWOT vs PT slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = slope_uncertainty_m_m, y = abs(residuals), color = factor(SWOT_reach_id))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("pt slope uncertainty") +
  ylab("abs slope difference") +
  theme_minimal(base_size = 30) +
  labs(color = "Reach ID")
