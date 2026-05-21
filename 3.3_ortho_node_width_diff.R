# =============================================================================
# Orthomosaic vs SWOT Node Width Comparison
# -----------------------------------------------------------------------------
# Joins manually digitized orthomosaic water area (per SWORD node polygon) to
# SWOT RiverSP/RiverTile node observations, filters to the overpass date
# matching each orthomosaic, derives ortho width from water area and SWORD
# prior node length, computes residuals and percent differences, removes
# per-reach-code median bias, and exports the matched dataset to CSV.
# The second section merges all per-cluster CSVs into one combined dataframe.
#
# Script sections:
#   1.  Read and filter SWOT node data
#   2.  Read orthomosaic width data
#   3.  Join ortho and SWOT data; compute ortho width
#   4.  Raw residuals and percent difference
#   5.  Bias removal (per-reach-code median) and bias-corrected metrics
#   6.  Data visualization
#   7.  Export per-cluster CSV
#   8.  Merge all cluster CSVs into one combined dataframe
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)


# =============================================================================
# 1. Read and filter SWOT node data
# =============================================================================

# --- RiverSP Version C (SWORD v16) -------------------------------------------
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv')

# --- RiverSP Version D (SWORD v17b) ------------------------------------------
SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv')

# --- RiverTile (SWORD v17b) --------------------------------------------------
# SWOT_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_domain_node_timeseries_v17b.csv')
# SWORD_v17b <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GIS/SWORD_v17b_domain/SWORD_YR_domain_v17b.csv')
# # RiverTile lacks p_length; join from SWORD prior
# SWOT_df <- SWOT_df %>%
#   left_join(SWORD_v17b, by = c("node_id", "reach_id"))
# SWOT_df$p_length <- SWOT_df$node_len

# Remove duplicates and fill values (time = -999…, wse = -1e12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# Quality filter: node_q < 2, valid cross-track swath, < 80% dark water
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>%           # 0=good, 1=suspect, 2=degraded, 3=bad
  filter(abs(xtrk_dist) >= 10000) %>%  # cross-track 10–60 km
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac <= 0.8)

# Convert TAI time to UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset


# =============================================================================
# 2. Read orthomosaic width data
# =============================================================================
# Each CSV was produced by script 3.2 (ortho_PIXCVec_node_polygon_widths.ipynb)
# and contains node_id, number_water_pixels, and water_area_m2 for one cluster.

# --- RiverSP Version C (SWORD v16) -------------------------------------------
# chandalar 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_CD_071024.csv')
# sheenjek 7/26
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_lowerPR_SJ_072624.csv')
# lower YR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_lowerYR_071624.csv')
# coleen / upper PR 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_upperPR_CL_071024.csv')
# coleen / upper PR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_upperPR_CL_071624.csv')
# upper YR 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/ortho_upperYR_071024.csv')

# --- RiverSP Version D (SWORD v17b) ------------------------------------------
# chandalar 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_CD_071024.csv')
# sheenjek 7/26
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_lowerPR_SJ_072624.csv')
# lower YR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_lowerYR_071624.csv')
# coleen / upper PR 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_upperPR_CL_071024.csv')
# coleen / upper PR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_upperPR_CL_071624.csv')
# upper YR 7/10
ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/ortho_upperYR_071024.csv')

# --- RiverTile (SWORD v17b) --------------------------------------------------
# chandalar 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_CD_071024.csv')
# sheenjek 7/26
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_lowerPR_SJ_072624.csv')
# lower YR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_lowerYR_071624.csv')
# coleen / upper PR 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_upperPR_CL_071024.csv')
# coleen / upper PR 7/16
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_upperPR_CL_071624.csv')
# upper YR 7/10
# ortho_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/ortho_upperYR_071024.csv')


# =============================================================================
# 3. Join ortho and SWOT data; compute ortho width
# =============================================================================

# Join on node_id (left join keeps all ortho nodes even if no SWOT match)
combined_ortho_SWORD_df <- ortho_df %>%
  left_join(SWOT_df_filtered, by = "node_id")

# Temporal filter: keep only the SWOT overpass matching the orthomosaic date
#   chandalar:        "2024-07-11"
#   sheenjek:         "2024-07-26"
#   lower YR:         "2024-07-16"
#   upper YR:         "2024-07-10"
#   coleen / upper PR: "2024-07-10" or "2024-07-16"
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  filter(str_starts(time_str, "2024-07-10"))  # CHANGE DATE TO MATCH

# Derive ortho width from water area and SWORD prior node length
# width (m) = water_area (m²) / node_length (m)
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(ortho_width_m = water_area_m2 / p_length)

# =============================================================================
# 4. Absolute residuals and percent difference
# =============================================================================

# Residuals: ortho width - SWOT width (positive = SWOT underestimates)
combined_ortho_SWORD_df$residuals <- (
  combined_ortho_SWORD_df$ortho_width_m - combined_ortho_SWORD_df$width)

# Percent difference using ortho as truth:
# % diff = |SWOT - ortho| / ortho × 100
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(percent_diff = ((abs(ortho_width_m - width)) / ortho_width_m) * 100)

# Percentile errors
# percentile_68_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.68, na.rm = TRUE)
# percentile_50_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.50, na.rm = TRUE)
# percentile_68_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.68, na.rm = TRUE)
# percentile_50_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.50, na.rm = TRUE)
# print(paste("68th Percentile Error:", percentile_68_error))
# print(paste("50th Percentile Error:", percentile_50_error))
# print(paste("68th Percentile Error %:", percentile_68_percent))
# print(paste("50th Percentile Error %:", percentile_50_percent))

# Pearson correlation
# cor_test <- cor.test(combined_ortho_SWORD_df$width, combined_ortho_SWORD_df$ortho_width_m)
# r_value  <- cor_test$estimate
# p_value  <- cor_test$p.value


# =============================================================================
# 5. Bias removal — per-reach-code median bias
# =============================================================================

# Group by the first 6 characters of reach_id as a sub-basin proxy
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  group_by(substr(reach_id, 1, 6)) %>%
  mutate(
    bias               = median(residuals, na.rm = TRUE),
    swot_width_nobias_m = width + bias  # shift SWOT width toward ortho
  ) %>%
  ungroup()

# Bias-corrected residuals
combined_ortho_SWORD_df$residuals_nobias <- (
  combined_ortho_SWORD_df$swot_width_nobias_m - combined_ortho_SWORD_df$ortho_width_m
)

# Bias-corrected percent difference
combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
  mutate(percent_diff_nobias = ((abs(ortho_width_m - swot_width_nobias_m)) / ortho_width_m) * 100)

# Bias-corrected percentile errors
# percentile_68_error_nobias   <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias),   0.68, na.rm = TRUE)
# percentile_50_error_nobias   <- quantile(abs(combined_ortho_SWORD_df$residuals_nobias),   0.50, na.rm = TRUE)
# percentile_68_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.68, na.rm = TRUE)
# percentile_50_percent_nobias <- quantile(abs(combined_ortho_SWORD_df$percent_diff_nobias), 0.50, na.rm = TRUE)
# print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
# print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# Bias-corrected correlation
# cor_test_nobias <- cor.test(combined_ortho_SWORD_df$ortho_width_m, combined_ortho_SWORD_df$swot_width_nobias_m)
# r_value_nobias  <- cor_test_nobias$estimate
# p_value_nobias  <- cor_test_nobias$p.value


# =============================================================================
# 6. Data visualization
# =============================================================================

# Scatter: SWOT vs ortho width
# color_palette <- c("#48b32e", "#6389ee", "#d99427", "#6D398B", "#C83232", "gray", "lightyellow", "pink")
# ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = width, color = factor(reach_id))) +
#   geom_point(size = 2.5) +
#   scale_color_manual(values = color_palette) +
#   geom_abline(linetype = "dashed", color = "gray") +
#   xlab("Ortho width (m)") +
#   ylab("SWOT width (m)") +
#   annotate("text",
#     x = min(combined_ortho_SWORD_df$ortho_width_m, na.rm = TRUE),
#     y = max(combined_ortho_SWORD_df$width, na.rm = TRUE),
#     label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 4),
#                    "\nn = ", nrow(combined_ortho_SWORD_df)),
#     hjust = 0, vjust = 1, size = 8) +
#   theme(legend.position = "none") +
#   labs(color = "Reach ID") +
#   theme_minimal(base_size = 30)

# Width along distance to outlet (raw vs bias-corrected)
ggplot(combined_ortho_SWORD_df) +
  geom_point(aes(x = p_dist_out / 1000, y = ortho_width_m),       color = "lightblue", size = 2.5, shape = 17) +
  geom_point(aes(x = p_dist_out / 1000, y = swot_width_nobias_m), color = "darkblue",  size = 2.5, alpha = 0.7) +
  xlab("Distance to outlet (km)") +
  ylab("Width (m)") +
  theme_minimal(base_size = 30)


# -----------------------------------------------------------------------------
# 6a. SWOT parameter comparison plots
# -----------------------------------------------------------------------------

# Pearson correlation: dark_frac vs percent_diff
cor_test <- cor.test(abs(combined_ortho_SWORD_df$dark_frac), combined_ortho_SWORD_df$percent_diff)
r_value  <- cor_test$estimate
p_value  <- cor_test$p.value

color_palette <- c("#6D398B", "#00429D", "#2E7D32")  # 3-river palette

# Scatter: dark_frac vs percent_diff
ggplot(combined_ortho_SWORD_df, aes(x = dark_frac, y = percent_diff, color = river)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Dark water fraction") +
  ylab("Percent difference (%)") +
  annotate("text",
    x     = max(combined_ortho_SWORD_df$dark_frac, na.rm = TRUE) - 0.28,
    y     = max(combined_ortho_SWORD_df$percent_diff, na.rm = TRUE),
    label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 5)),
    hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") +
  theme_minimal(base_size = 30)

# Violin: percent_diff by overpass date (bias-corrected)
ggplot(combined_ortho_SWORD_df,
       aes(x = factor(substr(time_utc, 1, 10)), y = percent_diff_nobias,
           fill = factor(substr(time_utc, 1, 10)))) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("Overpass date") +
  ylab("Percent difference (%)") +
  scale_fill_manual(values = c("lightblue", "darkblue")) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")

# Violin: percent_diff by river (raw)
ggplot(combined_ortho_SWORD_df, aes(x = river, y = percent_diff, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("River") +
  ylab("Percent difference (%)") +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(
    breaks = c("CL", "upper_PR", "upper_YR"),
    labels = c("Coleen", "Upper Porcupine", "Upper Yukon")
  ) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 0.9)) +
  ylim(0, 270)


# =============================================================================
# 7. Export per-cluster CSV
# =============================================================================

save_to_csv <- combined_ortho_SWORD_df %>%
  dplyr::select(
    node_id, reach_id, time_utc,
    residuals, percent_diff, bias, residuals_nobias, percent_diff_nobias,
    number_water_pixels, water_area_m2, ortho_width_m, swot_width_nobias_m,
    lat, lon, wse, wse_u, wse_r_u,
    width, width_u, area_total, area_tot_u, area_detct, area_det_u,
    area_wse, layovr_val, node_dist, xtrk_dist,
    node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q,
    cycle_id, pass_id, p_dist_out, p_length)

write.csv(
  save_to_csv,
  file      = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/combined_ortho_SWOT_upperYR_071024.csv',
  row.names = FALSE)


# =============================================================================
# 8. Merge all cluster CSVs into one combined dataframe
# =============================================================================

# Set working directory containing the combined_ortho_SWOT_*.csv files
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b"
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16"
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b"

setwd(wd)

# Read and bind all cluster CSVs
csv_files  <- list.files(wd, pattern = "^combined_ortho_SWOT_.*\\.csv$", full.names = TRUE)
combined_df <- bind_rows(lapply(csv_files, read.csv))

# Add river name labels
combined_df <- combined_df %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                       "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",
      river_code == "812509" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

write.csv(
  combined_df,
  file      = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/node_width_SWOT_Ortho.csv',
  row.names = FALSE)
