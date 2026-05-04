library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare PT wse & SWOT riverSP/rivertile node wse
# ---------------------------------------------------------------------------------------------------------------------------

# contents:
# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# read in & prep PT data
# match PT & SWOT observations in time and space
# all clusters comparison
# Compare PT wse & SWOT Version C (RiverSP) vs Version D (RiverSP, RiverTile) node wse

# ---------------------------------------------------------------------------------------------------------------------------
# read in & filter SWOT data
# version C: RiverSP
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_nodes_merged_RiverSP.csv')

# version D: RiverSP
SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv')

# version D: RiverTile
# SWORD v17b
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_PT_node_timeseries_v17b.csv')
# SWORD v16
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v16/RiverTile_PT_node_timeseries_v16.csv')

# get ride of possible duplicates from hydrocron pull
# this also filters out bad nodes without data (e.g. time = -999999999999, wse = -1.000000e+12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# filter SWOT data by node_q (0=good, 1=suspect, 2=degraded, 3=bad) & valid xtrk_dist (10-60km swath)
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  #filter(wse_u <= 0.5) %>% #this probably won't filter out many nodes beyond what node_q is doing
  filter(node_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) %>%
  filter(dark_frac <= 0.8)

# time_tai is seconds since 2001-01-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep PT data

# Set working directory to the folder chucked by separate rivers and PT clusters
# SWORD v17b
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/lower_PR"
# SWORD v16
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v16/SJ"
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
combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# ---------------------------------------------------------------------------------------------------------------------------
# match PT & SWOT observations in time and space

# match PT and SWOT in time
# observation are matched by 7.5min buffer (a SWOT obs should always be within 7.5min of a PT)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(SWOT_df_filtered %>%
                           filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5))) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# match PT and SWOT in space
# node level
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(Node_ID == node_id)

# a catch to get rid of any duplicate PT values
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[!duplicated(time_space_matched_SWOT_PT[c("pt_wse_m","pt_time_UTC","pt_serial")]),]

# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (residuals)
time_space_matched_SWOT_PT$residuals = time_space_matched_SWOT_PT$pt_wse_m - time_space_matched_SWOT_PT$wse

# Calculate the 68th percentile error
# percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.68, na.rm=TRUE)
# percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.50, na.rm=TRUE)

# print the result
# print(paste("68th Percentile Error:", percentile_68_error))
# print(paste("50th Percentile Error:", percentile_50_error))

# additional stats
# # Calculating the linear regression model 
# model = lm(pt_wse_m~wse, data = time_space_matched_SWOT_PT) 
# # Extracting R-squared parameter from summary 
# summary(model)
# #RMSE
# rmse <- sqrt(mean((time_space_matched_SWOT_PT$pt_wse_m - time_space_matched_SWOT_PT$wse)^2))
# #RMSE >= MAE, MAE is similar to 50th quantile error

# correlation test
cor_test <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_m)
# Extract r and p-value from cor_test
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8", "#F8A31B", 
              "#00429D", "#2E7D32", "#C83232", "#008F7A", "#E3A700", "#124000")

# SWOT vs PT wse plot
ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_m, y = wse, color = factor(pt_serial))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(time_space_matched_SWOT_PT$pt_wse_m, na.rm = TRUE), 
           y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
                          "\nn = ", nrow(time_space_matched_SWOT_PT)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "PT Serial") 


# CDF plot
# ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
#   # adjust x val of annotate to change writing locations
#   annotate("text", x = 0.5, y = 0.71, label = paste("|68%ile| diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
#   annotate("text", x = 0.5, y = 0.53, label = paste("|50%ile| diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)

# # Individual PT vs SWOT hydrograph
# pt_serial_list <- unique(time_space_matched_SWOT_PT$pt_serial)
# 
# pdf("PT_SWOT_hydrographs.pdf",width=7,height=5)
# for (PT_id in pt_serial_list) {
#   
#   individual_PT_SWOT_df <- time_space_matched_SWOT_PT %>%
#     filter(pt_serial == PT_id)
#   
#   individual_PT_df <- combined_PT_df %>%
#     filter(pt_serial == PT_id)
#   
#   # plot SWOT vs PT timeseries
#   plot <- ggplot() +
#     geom_errorbar(data = individual_PT_df, mapping = aes(x = pt_time_UTC, 
#                                                          ymin = pt_wse_m - pt_correction_mean_total_error_m, 
#                                                          ymax = pt_wse_m + pt_correction_mean_total_error_m,),
#     color = "#edc7a2", linewidth=5,) +
#     geom_point(individual_PT_df, mapping=aes(y=pt_wse_m, x=pt_time_UTC), color="#ED973D", size=2) +
#     geom_errorbar(data = individual_PT_SWOT_df, mapping = aes(x = time_utc, 
#                                                               ymin = wse - wse_u, 
#                                                               ymax = wse + wse_u), 
#                   linewidth = 1) +
#     geom_point(individual_PT_SWOT_df, mapping=aes(y=wse, x=time_utc), color="#01665E", size=7) +
#     ggtitle(PT_id) + xlab("time") + ylab("wse (m)") +
#     theme_minimal(base_size = 30) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 0.9)) 
#   print(plot)
# }
# dev.off()

# # singular hydrograph plot
# PT_id = 2156918
# 
# individual_PT_SWOT_df <- time_space_matched_SWOT_PT %>%
#   filter(pt_serial == PT_id)
# 
# individual_PT_df <- combined_PT_df %>%
#   filter(pt_serial == PT_id)
# 
# # plot SWOT vs PT timeseries
# ggplot() +
#   geom_errorbar(data = individual_PT_df, mapping = aes(x = pt_time_UTC, 
#                                             ymin = pt_wse_m - pt_correction_mean_total_error_m, 
#                                             ymax = pt_wse_m + pt_correction_mean_total_error_m,
#   ),
#   color = "#edc7a2", linewidth=5,
#   ) +
#   geom_point(individual_PT_df, mapping=aes(y=pt_wse_m, x=pt_time_UTC), color="#ED973D", size=2) +
#   geom_errorbar(data = individual_PT_SWOT_df, mapping = aes(x = time_utc, 
#                                                        ymin = wse - wse_u, 
#                                                        ymax = wse + wse_u), 
#                 linewidth = 1) +
#   geom_point(individual_PT_SWOT_df, mapping=aes(y=wse, x=time_utc), color="#01665E", size=7) +
#   ggtitle(PT_id) +
#   xlab("time") +
#   ylab("wse (m)") +
#   theme_minimal(base_size = 30) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0.9)) 


# ---------------------------------------------------------------------------------------------------------------------------
# remove bias from PT data

# removed median bias for individual PT
# need to set a min threshold of obs for us to calc a bias (using 3 currently)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  group_by(pt_serial) %>%
  mutate(
    bias = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    pt_wse_nobias_m = if (n() >= 3)
      pt_wse_m - bias
    else NA_real_) %>%
  ungroup()

# 68th & 50th percentile error: wse diff calculation

# Calculate the wse diff SWOT - PT (residuals)
time_space_matched_SWOT_PT$residuals_nobias = time_space_matched_SWOT_PT$pt_wse_nobias_m - time_space_matched_SWOT_PT$wse

# Calculate the 68th percentile error
# percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.68, na.rm=TRUE)
# percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.50, na.rm=TRUE)

#print the result
# print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
# print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

#csv subset
save_to_csv <- time_space_matched_SWOT_PT %>%
  dplyr::select(time_utc, pt_time_UTC, residuals, residuals_nobias, wse, wse_u, 
                pt_wse_m, pt_wse_nobias_m, bias, pt_correction_mean_total_error_m, 
                pt_correction_mean_offset_sd_m, pt_serial, width, width_u, node_id, reach_id, p_dist_out, 
                node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q, lat, lon)

# save joined_wse_subset to csv
write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/RiverSP_v17b_time_space_matched_SWOT_PT_lowerPR.csv', row.names = FALSE)

# # correlation test
# cor_test_nobias <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_nobias_m)
# 
# # Extract r and p-value
# r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
# p_value_nobias <- cor_test_nobias$p.value 
# 
# # plot SWOT vs PT wse
# ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_nobias_m, y = wse, color = factor(pt_serial))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   xlab("PT wse (m)") +
#   ylab("SWOT wse (m)") +
#   theme_minimal(base_size = 30) +
#   geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
#   annotate("text", x = min(time_space_matched_SWOT_PT$pt_wse_m, na.rm = TRUE), 
#            y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE), 
#            label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
#                           "\nn = ", nrow(time_space_matched_SWOT_PT)),
#            hjust = 0, vjust = 1, size = 8) +
#   labs(color = "PT Serial") 
# 
# 
# # CDF plot
# ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_nobias_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
#   annotate("text", x = 0.25, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
#   annotate("text", x = 0.25, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
#   theme_minimal(base_size = 20) 
#   #xlim(0, .751)




# ---------------------------------------------------------------------------------------------------------------------------
# merge all clusters
# ---------------------------------------------------------------------------------------------------------------------------

# Set working directory to CalVal_dataframes directory where the time&space matched SWOT/PT clusters are for each SWOT version
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16"
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b"
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v16"
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b"
setwd(wd)

# Get list of all CSV files in working directory that include time_space_matched_SWOT_PT_
csv_files <- list.files(wd, pattern = "time_space_matched_SWOT_PT.*\\.csv$", full.names = TRUE)

# Merge all PT files into a combined dataframe with filename column
data_list <- lapply(csv_files, function(file) {
  df <- read.csv(file)
  
  # Extract filename without path
  filename <- basename(file)
  
  # Select river cluster from filename
  filename_clean <- sub(".*PT_(.*)\\.csv$", "\\1", filename)
  
  # Add filename as a new column
  df <- df %>%
    mutate(river = filename_clean)
  
  return(df)
})


# TURN ON the correct df for vC vs vD load in
combined_time_space_matched_SWOT_PT_df <- bind_rows(data_list)
# combined_time_space_matched_RiverTile_PT_df <- bind_rows(data_list)

# Create a source column to identify which dataset each row comes from
RiverSP_df <- combined_time_space_matched_SWOT_PT_df %>%
  mutate(source = "RiverSP_PGD0")

# RiverTile_df <- combined_time_space_matched_RiverTile_PT_df %>%
#   mutate(source = "RiverTile")

# save dfs to csv
write.csv(RiverSP_df, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv', row.names = FALSE)


















# Combine both dataframes
# combined_df <- bind_rows(RiverSP_df, RiverTile_df)

# 
# # Summary stats
# # ---------------------------------------------------------------------------------------------------------------------------
# 
# # 68th & 50th percentile error: wse diff calculation
# 
# # Calculate the 68th & 50th percentile error
# percentile_68_error <- quantile(abs(RiverSP_df$residuals), 0.68, na.rm=TRUE)
# percentile_50_error <- quantile(abs(RiverSP_df$residuals), 0.50, na.rm=TRUE)
# 
# percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.68, na.rm=TRUE)
# percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.50, na.rm=TRUE)
# 
# # Calculate the 68th &50th percentile error, no bias
# percentile_68_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.68, na.rm=TRUE)
# percentile_50_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.50, na.rm=TRUE)
# 
# percentile_68_error_nobias_RiverTile <- quantile(abs(RiverTile_df$residuals_nobias), 0.68, na.rm=TRUE)
# percentile_50_error_nobias_RiverTile <- quantile(abs(RiverTile_df$residuals_nobias), 0.50, na.rm=TRUE)
# 
# 
# #print the result
# print(paste("68th Percentile Error:", percentile_68_error_RiverTile))
# print(paste("50th Percentile Error:", percentile_50_error_RiverTile))
# print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias_RiverTile))
# print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias_RiverTile))
# 
# # correlation test
# # cor_test <- cor.test(RiverSP_df$wse, RiverSP_df$pt_wse_m) # pt_wse_nobias_m
# cor_test <- cor.test(RiverTile_df$wse, RiverTile_df$pt_wse_m)
# 
# # Extract r and p-value
# r_value <- cor_test$estimate # Pearson correlation coefficient
# p_value <- cor_test$p.value # highly statistically significant is P < 0.001
# 
# # Mean Absolute Error (which is just the mean residual)
# MAE <- mean(abs(RiverTile_df$residuals))
# MAE <- mean(abs(RiverTile_df$residuals_nobias))
# 
# 
# t.test(abs(RiverTile_df$residuals), abs(RiverSP_df$residuals))
# t.test(abs(RiverTile_df$residuals_nobias), abs(RiverSP_df$residuals_nobias))
# 
# # ---------------------------------------------------------------------------------------------------------------------------
# # data viz
# 
# color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#F4845F", "#DA627D")
# 
# 
# # plot SWOT vs PT wse
# ggplot(combined_time_space_matched_SWOT_PT_df, aes(x = pt_wse_m, y = wse, color = factor(river))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   xlab("PT wse (m)") +
#   ylab("SWOT wse (m)") +
#   theme_minimal(base_size = 30) +
#   geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
#   annotate("text", x = min(combined_time_space_matched_SWOT_PT_df$pt_wse_m, na.rm = TRUE), 
#            y = max(combined_time_space_matched_SWOT_PT_df$wse, na.rm = TRUE), 
#            label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
#                           "\nn = ", nrow(combined_time_space_matched_SWOT_PT_df)),
#            hjust = 0, vjust = 1, size = 8) +
#   labs(color = "River")
# 
# # plot SWOT vs PT wse
# ggplot(combined_time_space_matched_RiverTile_PT_df, aes(x = pt_wse_m, y = wse, color = factor(river))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   xlab("PT wse (m)") +
#   ylab("SWOT wse (m)") +
#   theme_minimal(base_size = 30) +
#   geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
#   annotate("text", x = min(combined_time_space_matched_RiverTile_PT_df$pt_wse_m, na.rm = TRUE), 
#            y = max(combined_time_space_matched_RiverTile_PT_df$wse, na.rm = TRUE), 
#            label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
#                           "\nn = ", nrow(combined_time_space_matched_RiverTile_PT_df)),
#            hjust = 0, vjust = 1, size = 8) +
#   labs(color = "River") 
# 
# 
# # CDF plot
# ggplot(combined_time_space_matched_SWOT_PT_df, aes(x = abs(wse - pt_wse_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
#   annotate("text", x = 1.03, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
#   annotate("text", x = 1.03, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
#   theme_minimal(base_size = 20) 
# 
# ggplot(combined_time_space_matched_RiverTile_PT_df, aes(x = abs(wse - pt_wse_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
#   annotate("text", x = 2.03, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_RiverTile, 4)), color = "#222222", size = 6) +
#   annotate("text", x = 2.03, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_RiverTile, 4)), color = "#222222", size = 6) +
#   theme_minimal(base_size = 20) 
# 
# 
# 
# summary <- group_by(combined_time_space_matched_RiverTile_PT_df, river) %>% summarise(
#   count = n(),
#   mean = mean(abs(residuals), na.rm = TRUE),
#   sd = sd(abs(residuals), na.rm = TRUE),
#   median = median(abs(residuals), na.rm = TRUE),
#   IQR = IQR(abs(residuals), na.rm= TRUE),
#   min =min(abs(residuals), na.rm = TRUE),
#   max =max(abs(residuals), na.rm= TRUE),
#   quant68 = quantile(abs(residuals), 0.68, na.rm=TRUE),
#   quant68_nobais = quantile(abs(residuals_nobias), 0.68, na.rm=TRUE)
# )











# Combo CDF plot
ggplot(combined_df, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT WSE - PT WSE") +
  annotate("text", x = 1.45, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 1.45, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))


# Combo CDF plot no bias
ggplot(combined_df, aes(x = abs(residuals_nobias), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT WSE - PT WSE") +
  annotate("text", x = 0.9, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error_nobias, 4), 
                         ", Version D:", round(percentile_68_error_RiverTile_nobias, 4)), 
           color = "#222222", size = 5) +
  annotate("text", x = 0.9, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error_nobias, 4), 
                         ", Version D:", round(percentile_50_error_RiverTile_nobias, 4)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash"))






# ---------------------------------------------------------------------------------------------------------------------------
# data viz without bias

# correlation test
cor_test_nobias <- cor.test(combined_time_space_matched_RiverTile_PT_df$wse, combined_time_space_matched_RiverTile_PT_df$pt_wse_nobias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 


combined_time_space_matched_RiverTile_PT_df <- combined_time_space_matched_RiverTile_PT_df %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",   # merge into one
    TRUE ~ river                              # keep all others unchanged
  ))


# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

# plot SWOT vs PT wse
ggplot(combined_time_space_matched_RiverTile_PT_df, 
       aes(x = pt_wse_nobias_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "lowerPR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")) +
  xlab("PT WSE (m)") +
  ylab("SWOT WSE (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", 
           x = min(combined_time_space_matched_RiverTile_PT_df$pt_wse_nobias_m, na.rm = TRUE), 
           y = max(combined_time_space_matched_RiverTile_PT_df$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value_nobias, 4), 
                          "\np value = ", signif(p_value_nobias, 3),
                          "\nn = ", nrow(combined_time_space_matched_RiverTile_PT_df)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "River")


summary <- group_by(combined_time_space_matched_RiverTile_PT_df, river) %>% summarise(
  count = n(),
  mean = mean(abs(residuals_nobias), na.rm = TRUE),
  sd = sd(abs(residuals_nobias), na.rm = TRUE),
  median = median(abs(residuals_nobias), na.rm = TRUE),
  IQR = IQR(abs(residuals_nobias), na.rm= TRUE),
  min =min(abs(residuals_nobias), na.rm = TRUE),
  max =max(abs(residuals_nobias), na.rm= TRUE),
  quant68 = quantile(abs(residuals_nobias), 0.68, na.rm=TRUE),
  quant68_nobais = quantile(abs(residuals_nobias), 0.68, na.rm=TRUE),
  mean_dark = mean(dark_frac, na.rm = TRUE),
  sd_dark = sd(dark_frac, na.rm = TRUE),
  median_dark = median(dark_frac, na.rm = TRUE),
  IQR_dark = IQR(dark_frac, na.rm= TRUE),
  min_dark =min(dark_frac, na.rm = TRUE),
  max_dark =max(dark_frac, na.rm= TRUE),
  quant68_dark = quantile(dark_frac, 0.68, na.rm=TRUE),
  quant68_nobais_dark = quantile(dark_frac, 0.68, na.rm=TRUE)
)


# **********************************
# Reorder the factor levels for river
combined_time_space_matched_RiverTile_PT_df$river <- factor(
  combined_time_space_matched_RiverTile_PT_df$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# Replot
ggplot(combined_time_space_matched_RiverTile_PT_df, aes(x = river, y = abs(residuals_nobias), fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("|SWOT - PT WSE| (m)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 0.9)
  ) +
  ylim(0, 0.8)

# **********************************













# # add river names to df
# SWOT_df_filtered <- SWOT_df_filtered %>%
#   mutate(
#     river_code = substr(reach_id, 1, 6),
#     river = case_when(
#       river_code == "812701" ~ "lower_YR",
#       river_code == "812705" ~ "upper_YR",
#       river_code == "812508" ~ "CD",
#       river_code == "812603" ~ "upper_PR",
#       river_code == "812605" ~ "upper_PR",
#       river_code == "812604" ~ "CL",
#       river_code == "812603" ~ "SJ",
#       TRUE ~ NA_character_
#     )
#   )

# Trim to 25th-75th percentile per river
trimmed_df <- combined_time_space_matched_RiverTile_PT_df %>%
  group_by(river) %>%
  mutate(
    q25 = quantile(width, 0.25, na.rm = TRUE),
    q75 = quantile(width, 0.75, na.rm = TRUE)
  ) %>%
  filter(width >= q25, width <= q75) %>%
  ungroup() %>%
  mutate(
    width = if_else(river == "upperYR", width + 150, width),
    width = if_else(river == "lowerYR", width + 700, width),
    width = if_else(river == "CL", width - 10, width),
    width = if_else(river == "PR", width + 50, width),
    width = if_else(river == "SJ", width - 10, width)
  )


# Compute medians per river (median = 50th percentile)
median_df <- trimmed_df %>%
  group_by(river) %>%
  summarize(median_width = median(width, na.rm = TRUE)) %>%
  ungroup()

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









ggplot(RiverSP_filtered, aes(x = river, y = abs(residuals_nobias), fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab('River') +
  ylab("|SWOT - PT WSE| (m)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 30) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("time_space_matched_SWOT_PT_CD", "time_space_matched_SWOT_PT_CL", "time_space_matched_SWOT_PT_upper_PR", "time_space_matched_SWOT_PT_lower_YR", "time_space_matched_SWOT_PT_upper_YR"),
    labels = c("Chandalar", "Coleen", "Upper Porcupine", "Lower Yukon", "Upper Yukon")
  ) +
  scale_x_discrete(
    breaks = c("time_space_matched_SWOT_PT_CD", "time_space_matched_SWOT_PT_CL", "time_space_matched_SWOT_PT_upper_PR", "time_space_matched_SWOT_PT_lower_YR", "time_space_matched_SWOT_PT_upper_YR"),
    labels = c("Chandalar", "Coleen", "Upper Porcupine", "Lower Yukon", "Upper Yukon")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9)
  )









# correlation test
cor_test_nobias <- cor.test(combined_time_space_matched_SWOT_PT_df$wse, combined_time_space_matched_SWOT_PT_df$pt_wse_nobias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs PT wse
ggplot(combined_time_space_matched_SWOT_PT_df, aes(x = pt_wse_nobias_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("PT wse (m)") +
  ylab("SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_time_space_matched_SWOT_PT_df$pt_wse_m, na.rm = TRUE), 
           y = max(combined_time_space_matched_SWOT_PT_df$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value_nobias, 4), "\np value = ", signif(p_value_nobias, 3),
                          "\nn = ", nrow(combined_time_space_matched_SWOT_PT_df)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "River") 


# CDF plot
ggplot(combined_time_space_matched_SWOT_PT_df, aes(x = abs(wse - pt_wse_nobias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability", title = "CDF of SWOT WSE - PT WSE") +
  annotate("text", x = 2.03, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 2.03, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 


# ---------------------------------------------------------------------------------------------------------------------------
# data viz

# plot wse diff vs width
ggplot(RiverTile_df, aes(x = width, y = abs(residuals_nobias), color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  #xlab("SWOT width (m)") +
  ylab("PT - SWOT wse (m)") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 





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

