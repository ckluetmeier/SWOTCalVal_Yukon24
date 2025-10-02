library(tidyverse)
library(lubridate)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------------
# Compare orthomosaic width & SWOT riverSP node width
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# read in SWOT data
# RiverSP
SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv')

# get ride of possible duplicates from hydrocron pull
# this also filters out bad nodes without data (e.g. time = -999999999999, wse = -1.000000e+12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# filter SWOT data by node_q (0=good, 1=suspect, 2=degraded, 3=bad) & xtrk_dist (10-60km)
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>%
  filter(abs(xtrk_dist) >=10000) %>%
  filter(abs(xtrk_dist) <=60000) #%>%
  #filter(dark_frac <= 0.5)

# time_tai is seconds since 2001-011-01, offset 37 seconds from UTC
tai_epoch <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

# Convert time_tai to UTC
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# ---------------------------------------------------------------------------------------------------------------------------
# read in & prep ortho data

# orthos processed with PIXCVec polygons
# upper YR
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/Orthomosaics/upperYR_240710_summaryStats_node_polygon.csv')
# 7/10 Coleen
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/Orthomosaics/CL_upperPR_240710_summaryStats_node_polygon.csv')
# 7/16 Coleen
ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/Orthomosaics/CL_upperPR_240716_summaryStats_node_polygon.csv')



# ---------------------------------------------------------------------------------------------------------------------------
# join ortho & SWOT data, calculate ortho width with SWORD prior node length

# join ortho and SWOT data by node_id
combined_ortho_SWORD_df <- ortho_df %>%
  left_join(SWOT_df_filtered, by = "node_id")

# filter to same date at ortho collect
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  filter(str_starts(time_str, "2024-07-16")) # CHANGE DATE DEPENDING ON ORTHO COLLECT HERE

# Calculate width from area
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(ortho_width_m = water_area_m2/p_length) #p_length

# # filter to same reaches as in the pixcevec polygons for SWORD polygon comparisons
# # upper YR
# combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
#   filter(reach_id %in% c(81270501181, 81270501191, 81270501201))
# upper PR & CL
# combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
#   filter(reach_id %in% c(81260500011, 81260300221, 81260300211, 81260401181, 81260401011))

# ---------------------------------------------------------------------------------------------------------------------------
# Summary stats

# 68th & 50th percentile error: width diff calculation

# Calculate the width diff ortho - SWOT (residuals)
combined_ortho_SWORD_df$residuals = combined_ortho_SWORD_df$ortho_width_m - combined_ortho_SWORD_df$width

# Calculate the width percent diff Ortho/SWOT
# using ortho as truth:
# % diff = ( |swot - ortho| ) / ( ortho ) *100
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(percent_diff = ((abs(ortho_width_m - width)) / ortho_width_m) * 100)

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.68, na.rm=TRUE)
percentile_50_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

print(paste("68th Percentile Error:", percentile_68_percent))
print(paste("50th Percentile Error:", percentile_50_percent))

# correlation test
cor_test <- cor.test(combined_ortho_SWORD_df$width, combined_ortho_SWORD_df$ortho_width_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#48b32e","#6389ee","#d99427","#6D398B","#C83232")

# plot SWOT vs GNSS width
ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = width, color = factor(reach_id))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_ortho_SWORD_df$ortho_width_m, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$width, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(combined_ortho_SWORD_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") + # comment off to see reach ids
  labs(color = "Reach ID") 

# CDF plot of absolute wse difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(width - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of SWOT Width - Ortho Width") +
  annotate("text", x = 150, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 150, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# CDF plot of percent difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of SWOT-Ortho Width %diff") +
  annotate("text", x = 205, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent, 4)), color = "#222222", size = 6) +
  annotate("text", x = 205, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# for plotting upper PR & CL separately
width_dist_out_df <- combined_ortho_SWORD_df %>%
  filter(reach_id %in% c(81260401181, 81260401011))
# PR reaches
# 81260500011, 81260300221, 81260300211, 
# CL reaches
# 81260401181, 81260401011

# plot widths along dist_out
ggplot(width_dist_out_df) +
  geom_point(aes(x = p_dist_out/1000, y = ortho_width_m), color = "lightblue", size = 2.5, shape = 17) +
  geom_point(aes(x = p_dist_out/1000, y = width),  color = "darkblue", size = 2.5, alpha = 0.7) +
  scale_color_manual(values = color_palette) +
  xlab("Distance to outlet (km)") +
  ylab("width (m)") +
  theme_minimal(base_size = 30) +
  ggtitle('Coleen, 7/16/24')

# large residual investigation
# error_investigation_df <- combined_ortho_SWORD_df %>%
#   filter(width > 400)

# ---------------------------------------------------------------------------------------------------------------------------
# remove bias

# calculate bias by river/sub-basin
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  group_by(substr(reach_id, 0,6)) %>%
  mutate(
    bias = median(residuals, na.rm = TRUE),
    swot_width_nobias_m = width + bias
  ) %>%
  ungroup()


# 68th & 50th percentile error no bias: width diff calculation

# Calculate the width diff SWOT - GNSS (residuals)
combined_ortho_SWORD_df$residuals_nobias = combined_ortho_SWORD_df$swot_width_nobias_m - combined_ortho_SWORD_df$ortho_width_m

# Calculate the width percent diff Ortho/SWOT (residuals)
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(percent_diff_nobias = ((abs(ortho_width_m - swot_width_nobias_m)) / ortho_width_m) * 100)


# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.68, na.rm=TRUE)
percentile_50_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# correlation test
cor_test_nobias <- cor.test(combined_ortho_SWORD_df$ortho_width_m, combined_ortho_SWORD_df$swot_width_nobias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 

# plot SWOT vs GNSS width
ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = swot_width_nobias_m, color = factor(reach_id))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_ortho_SWORD_df$swot_width_nobias_m, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$width, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3),5),
                          "\nn = ", nrow(combined_ortho_SWORD_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")


# CDF plot
ggplot(combined_ortho_SWORD_df, aes(x = abs(ortho_width_m - swot_width_nobias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of SWOT Width - Ortho Width") +
  annotate("text", x = 170, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 170, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 

# CDF plot of percent difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(percent_diff_nobias))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of SWOT-Ortho Width %diff") +
  annotate("text", x = 205, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 205, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# for plotting upper PR & CL separately
width_dist_out_df <- combined_ortho_SWORD_df %>%
  filter(reach_id %in% c(81260500011, 81260300221, 81260300211))
# PR reaches
# 81260500011, 81260300221, 81260300211
# CL reaches
# 81260401181, 81260401011

# plot widths along dist_out
ggplot(width_dist_out_df) +
  geom_point(aes(x = p_dist_out/1000, y = ortho_width_m), color = "lightblue", size = 2.5, shape = 17) +
  geom_point(aes(x = p_dist_out/1000, y = swot_width_nobias_m),  color = "darkblue", size = 2.5, alpha = 0.7) +
  scale_color_manual(values = color_palette) +
  xlab("Distance to outlet (km)") +
  ylab("width (m)") +
  theme_minimal(base_size = 30) +
  ggtitle('Coleen no bias, 7/16/24')


# ---------------------------------------------------------------------------------------------------------------------------

#csv subset
save_to_csv <- combined_ortho_SWORD_df %>%
  dplyr::select(node_id,reach_id,time_utc,residuals,percent_diff,residuals_nobias,percent_diff_nobias,
                number_water_pixels,water_area_m2,ortho_width_m,swot_width_nobias_m,
                lat, lon, wse, wse_u, wse_r_u, width, width_u, area_total, area_tot_u, area_detct,
                area_det_u, area_wse, layovr_val, node_dist, xtrk_dist, node_q, node_q_b, dark_frac,
                n_good_pix, rdr_sig0, xovr_cal_q, cycle_id, pass_id, p_dist_out, p_length)

# save joined_wse_subset to csv
write.csv(save_to_csv, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/Orthomosaics/combined_ortho_SWOT_CL_upperPR_071624.csv', row.names = FALSE)



# ---------------------------------------------------------------------------------------------------------------------------
# all clusters comparison
# ---------------------------------------------------------------------------------------------------------------------------

# Set working directory to the orthomosaics directory with the combined_ortho_SWOT_ csvs
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/Orthomosaics"
setwd(wd)

# Get list of all CSV files in working directory that start with combined_ortho_SWOT
csv_files <- list.files(wd, pattern = "^combined_ortho_SWOT_.*\\.csv$", full.names = TRUE)

# Merge all PT files into a combined dataframe with filename column
data_list <- lapply(csv_files, function(file) {
  df <- read.csv(file)
  return(df)
})

# Combine all dataframes into one
combined_df <- bind_rows(data_list)

# add river names to df
combined_ortho_SWORD_df <- combined_df %>%
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

# optional filter to look at isolated groups (e.g. by river)
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  filter(river == "upper_YR")

# ---------------------------------------------------------------------------------------------------------------------------
# Stats & plots
# ---------------------------------------------------------------------------------------------------------------------------

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.68, na.rm=TRUE)
percentile_50_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_percent))
print(paste("50th Percentile Error:", percentile_50_percent))

# correlation test
cor_test <- cor.test(combined_ortho_SWORD_df$width, combined_ortho_SWORD_df$ortho_width_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# ---------------------------------------------------------------------------------------------------------------------------
# data viz

color_palette <- c("#48b32e","#6389ee","#d99427","#6D398B","#C83232","#d4c841", "#e73582", "#95d47e")

# plot SWOT vs GNSS width
ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = width, color = factor(reach_id))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_ortho_SWORD_df$ortho_width_m, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$width, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(combined_ortho_SWORD_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") + 
  labs(color = "Reach ID") 


# CDF plot of absolute wse difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(width - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of SWOT Width - Ortho Width") +
  annotate("text", x = 350, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 350, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# CDF plot of percent difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of SWOT-Ortho Width %diff") +
  annotate("text", x = 200, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent, 4)), color = "#222222", size = 6) +
  annotate("text", x = 200, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# no bias
# ---------------------------------------------------------------------------------------------------------------------------

# Calculate the 68th percentile error
percentile_68_error_nobias <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.68, na.rm=TRUE)
percentile_50_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.50, na.rm=TRUE)

# plot SWOT vs GNSS width
ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = swot_width_nobias_m, color = factor(reach_id))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_ortho_SWORD_df$swot_width_nobias_m, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$width, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3),5),
                          "\nn = ", nrow(combined_ortho_SWORD_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")


# CDF plot
ggplot(combined_ortho_SWORD_df, aes(x = abs(swot_width_nobias_m - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of SWOT Width - Ortho Width") +
  annotate("text", x = 400, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 400, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) 

# CDF plot of percent difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(percent_diff_nobias))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of SWOT-Ortho Width %diff") +
  annotate("text", x = 140, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent_nobias, 4)), color = "#222222", size = 6) +
  annotate("text", x = 140, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent_nobias, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# ---------------------------------------------------------------------------------------------------------------------------
# Comparison with SWOT parameters plots

# correlation test
cor_test <- cor.test(abs(combined_ortho_SWORD_df$dark_frac), combined_ortho_SWORD_df$percent_diff)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# plot percent diff vs swot df variables
ggplot(combined_ortho_SWORD_df, aes(x = dark_frac, y = percent_diff, color = river)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") +
  ylab("percent_diff (%)") +
  theme_minimal(base_size = 30) +
  annotate("text", x = max(combined_ortho_SWORD_df$dark_frac-.28, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$percent_diff, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3),5)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")

color_palette <- c("#D86A1A", "#6D398B", "#00429D", "#F8A31B", "#2E7D32")
color_palette <- c("#6D398B", "#00429D","#2E7D32")

# violin comparison
# factor(node_q_b), factor(pass_id)
ggplot(combined_ortho_SWORD_df, aes(x = factor(substr(time_utc, 0, 10)), y = percent_diff_nobias, fill=factor(substr(time_utc, 0, 10))), drop = FALSE) + 
  geom_violin(alpha = 0.8, color= NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("time_utc") +
  ylab("Percent Difference (%)") +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")
# violin comparision by river
ggplot(combined_ortho_SWORD_df, aes(x = river, y = percent_diff, fill=river), drop = FALSE) + 
  geom_violin(alpha = 0.8, color= NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("River") +
  ylab("Percent Difference (%)") +
  scale_fill_manual(values=color_palette) +
  theme_minimal(base_size = 30)+
  scale_x_discrete(
    breaks = c("CL", "upper_PR", "upper_YR"),
    labels = c("Coleen", "Upper Porcupine", "Upper Yukon")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9)) +
  ylim(0, 270)

# Combo CDF plot by river
ggplot(combined_ortho_SWORD_df, aes(x = percent_diff, color = river, linetype = river)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho width") +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("CL" = "#6D398B", "upper_PR" = "#00429D", "upper_YR" ="#2E7D32")) +
  theme(legend.position = "none")


# Expand to long df for width type comparisons
long_df <- combined_ortho_SWORD_df %>%
  select(river, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
                             ortho_width_m = "Ortho",
                             width = "SWOT"))
# split violin comparisons
devtools::install_github("psyteachr/introdataviz")
ggplot(long_df, aes(x = river, y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, color= NA) +
  geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_x_discrete(labels = c("Coleen", "Porcupine", "Yukon")) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  ylab("width (m)") +
  ylim(0, 1600) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9))

# Expand to long df for width type comparisons by day
long_df <- combined_ortho_SWORD_df %>%
  select(time_utc, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
                             ortho_width_m = "Ortho",
                             width = "SWOT"))
# split violin comparisons
devtools::install_github("psyteachr/introdataviz")
ggplot(long_df, aes(x = factor(substr(time_utc, 0, 10)), y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, color= NA) +
  geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  ylab("width (m)") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9))

# CDF plot by day (used to look at upper PR differences)
ggplot(combined_ortho_SWORD_df, aes(x = percent_diff, color = factor(substr(time_utc, 0, 10)), linetype = factor(substr(time_utc, 0, 10)))) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho width upper_PR") +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("2024-07-10" = "lightblue", "2024-07-16" = "#00429D")) +
  theme(legend.position = "none")
