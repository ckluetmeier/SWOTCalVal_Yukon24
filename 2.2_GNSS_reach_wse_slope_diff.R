# =============================================================================
# GNSS vs SWOT Reach WSE and Slope Comparison
# -----------------------------------------------------------------------------
# Matches GNSS drift survey reach WSE and slope observations to SWOT reach data
# in time (±5 hour buffer around drift midpoint) and space (same reach ID),
# computes residuals for both WSE and slope, removes per-drift median bias
# from WSE, and exports the matched dataset to CSV.
#
# Script sections:
#   1.  Read and filter SWOT reach data
#   2.  Read and prepare GNSS data
#   3.  Match GNSS and SWOT in time and space
#   4.  WSE: absolute residuals, summary stats, and plots
#   5.  WSE: bias removal (per-drift median) and plots
#   6.  Save WSE dataframe and add river labels
#   7.  Slope: absolute residuals and plots
#   8.  Slope: bias-corrected relative residuals and plots
#   9.  Export
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)


# =============================================================================
# 1. Read and filter SWOT reach data
# =============================================================================

# --- RiverSP Version C (SWORD v16) -------------------------------------------
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# --- RiverSP Version D (SWORD v17b) -------------------------------------------
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v17b/RiverSP_domain_reach_timeseries_PGD0_v17b.csv')

# --- RiverTile ---------------------------------------------------------------
# SWORD v16
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v16/RiverTile_domain_reach_timeseries_v16.csv')
# SWORD v17b
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverTile_v17b/RiverTile_domain_reach_timeseries_v17b.csv')

# Remove duplicates and sentinel fill values (time = -999…, wse = -1e12)
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE) %>%
  filter(time > 0) %>%
  filter(wse > 0)

# Quality filter: reach_q < 2, valid cross-track swath, < 80% dark water,
# partial_f == 0 means >= 50% of nodes present
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%           # 0=good, 1=suspect, 2=degraded, 3=bad
  filter(abs(xtrk_dist) >= 10000) %>%  # cross-track 10–60 km
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac < 0.8) %>%
  filter(partial_f == 0)

# Convert TAI time to UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset


# =============================================================================
# 2. Read and prepare GNSS data
# =============================================================================

# --- SWORD v16 ---------------------------------------------------------------
GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_wse_slope.csv')
GNSS_slope_corr_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_reach_slope_corrected.csv') %>%
  rename(reach_drift_slope_m_m_abs_nobias = slope_m_m)

# --- SWORD v17b --------------------------------------------------------------
# GNSS_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv')
# GNSS_slope_corr_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_slope_corrected.csv') %>%
#   rename(reach_drift_slope_m_m_abs_nobias = slope_m_m)

# Join bias-corrected slope to main GNSS dataframe
GNSS_df <- GNSS_df %>%
  inner_join(
    GNSS_slope_corr_df %>%
      select(reach_drift_slope_m_m_abs_nobias, drift_id, reach_id),
    by = c("reach_id", "drift_id")
  )

# Convert drift start/end times to POSIXct
GNSS_df$wse_drift_start_UTC <- as.POSIXct(GNSS_df$wse_drift_start_UTC, tz = "UTC")
GNSS_df$wse_drift_end_UTC   <- as.POSIXct(GNSS_df$wse_drift_end_UTC,   tz = "UTC")

# Compute midpoint time of each drift for matching to SWOT overpass
GNSS_df$wse_drift_midpoint_UTC <- as.POSIXct(
  (as.numeric(GNSS_df$wse_drift_start_UTC) + as.numeric(GNSS_df$wse_drift_end_UTC)) / 2,
  origin = "1970-01-01", tz = "UTC"
)

# Compute total duration of each drift
GNSS_df$wse_drift_total_time_UTC <- difftime(
  GNSS_df$wse_drift_end_UTC, GNSS_df$wse_drift_start_UTC, units = "mins"
)

# Rename reach_id to avoid column conflict when merging with SWOT data
GNSS_df <- rename(GNSS_df, "GNSS_reach_id" = "reach_id")


# =============================================================================
# 3. Match GNSS and SWOT in time and space
# =============================================================================

# Temporal match: ±5-hour buffer around drift midpoint
time_matched_SWOT_GNSS <- GNSS_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_reach_df_filtered %>%
        filter(abs(difftime(wse_drift_midpoint_UTC, time_utc, units = "hours")) <= 5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Spatial match: same reach ID
time_space_matched_SWOT_GNSS <- time_matched_SWOT_GNSS %>%
  filter(GNSS_reach_id == reach_id)


# =============================================================================
# 4. WSE - absolute residuals, summary statistics, and plots
# =============================================================================

# Compute raw residuals: GNSS WSE - SWOT WSE
time_space_matched_SWOT_GNSS$residuals <- (
  time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m - time_space_matched_SWOT_GNSS$wse
)

# Per-drift summary
summary_by_drift <- group_by(time_space_matched_SWOT_GNSS, drift_id) %>%
  summarise(
    count  = n(),
    mean   = mean(abs(residuals), na.rm = TRUE),
    sd     = sd(abs(residuals), na.rm = TRUE),
    median = median(abs(residuals), na.rm = TRUE),
    IQR    = IQR(abs(residuals), na.rm = TRUE),
    min    = min(abs(residuals), na.rm = TRUE),
    max    = max(abs(residuals), na.rm = TRUE)
  )

percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

# Pearson correlation
cor_test <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m)
r_value  <- cor_test$estimate  # Pearson r
p_value  <- cor_test$p.value   # p < 0.001 is highly statistically significant


# -----------------------------------------------------------------------------
# 4a. Data visualization — WSE (absolute)
# -----------------------------------------------------------------------------

color_palette <- c(
  "#e73582", "#78e964", "#5329a8", "#e1df28", "#8c66f0", "#a8dc3e", "#7150ce",
  "#48b32e", "#cb4bd0", "#5ddd76", "#df36ad", "#4e227d", "#c7e772", "#6549a4",
  "#d4c841", "#4c64c1", "#7ca92f", "#9a3c9a", "#409f46", "#b67ce0", "#95d47e",
  "#bc4797", "#df3542", "#6389ee", "#d99427", "#dd80d5", "#417a28", "#de6096",
  "#982161", "#c87735", "#c9405b", "#de4b21", "#b84b34", "pink", "blue",
  "lightblue", "darkblue", "darkred"
)

# Scatter: SWOT vs GNSS reach WSE, colored by drift ID
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_reach_drift_wse_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS WSE (m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_GNSS$mean_reach_drift_wse_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_GNSS)
    ),
    hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") +
  # labs(color = "Drift ID")  # show legend
  theme_minimal(base_size = 30)

# CDF: absolute WSE difference
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_reach_drift_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT WSE - GNSS WSE (m)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT WSE - GNSS WSE"
  ) +
  annotate("text", x = 1, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 1, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)
# xlim(0, 3)

# Scatter: residuals vs GNSS reach total error, colored by drift ID
ggplot(time_space_matched_SWOT_GNSS,
       aes(x = abs(residuals), y = mean_reach_drift_wse_total_error_m,
           color = factor(substr(drift_id, 70, 93)))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = color_palette) +
  xlab("GNSS - SWOT WSE (m)") +
  ylab("GNSS total error (m)") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")
# labs(color = "Drift ID")


# =============================================================================
# 5. WSE — bias removal (per-drift median, min 3 observations) and plots
# =============================================================================

time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  group_by(drift_id) %>%
  mutate(
    bias                           = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    mean_reach_drift_wse_no_bias_m = if (n() >= 3) mean_reach_drift_wse_m - bias else NA_real_
  ) %>%
  ungroup()

# Bias-corrected residuals
time_space_matched_SWOT_GNSS$residuals_nobias <- (
  time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse
)

percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# Pearson correlation (bias-corrected)
cor_test_nobias <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m)
r_value_nobias  <- cor_test_nobias$estimate
p_value_nobias  <- cor_test_nobias$p.value

# Scatter: SWOT vs bias-corrected GNSS reach WSE
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_reach_drift_wse_no_bias_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS WSE (bias-corrected, m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_GNSS$mean_reach_drift_wse_no_bias_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value_nobias, 4),
      "\np value = ", signif(p_value_nobias, 3),
      "\nn = ", nrow(time_space_matched_SWOT_GNSS)
    ),
    hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") +
  theme_minimal(base_size = 30)

# CDF: bias-corrected absolute WSE difference
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_reach_drift_wse_no_bias_m))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT WSE - GNSS WSE (m)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT WSE - GNSS WSE (bias-corrected)"
  ) +
  annotate("text", x = 0.5, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 0.5, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# Bias over time (by drift ID)
ggplot(time_space_matched_SWOT_GNSS, aes(x = time_utc, y = bias, color = factor(substr(drift_id, 70, 93)))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("Time (UTC)") +
  ylab("Bias (m)") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")
# labs(color = "Drift ID")


# =============================================================================
# 6. Add river labels and export WSE dataframe
# =============================================================================

time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # Specific reach_id checks take priority to avoid overwriting SJ/BL labels
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                       "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",  # until the Circle bifurcation
      river_code == "812509" ~ "lowerYR",  # past the PR confluence
      river_code == "812705" ~ "upperYR",  # Circle up
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PIC0")  ## CHANGE TO CORRECT VERSION!

# RiverSP column selection for CSV export
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  dplyr::select(
    reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC,
    wse_drift_midpoint_UTC, wse_drift_total_time_UTC,
    residuals, residuals_nobias, bias,
    mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
    mean_reach_drift_wse_no_bias_m,
    reach_drift_slope_m_m, reach_drift_slope_m_m_abs,
    reach_drift_slope_precision_m, slope_residuals, slope_residuals_nobias,
    reach_drift_slope_m_m_abs_nobias,
    drift_id, wse, wse_u,
    slope, slope_abs, slope_u, slope_r_u,
    width, width_u, area_total, area_tot_u, area_detct, area_det_u,
    area_wse, layovr_val, node_dist, xtrk_dist,
    reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q,
    p_dist_out, p_lat, p_lon, cycle_id, pass_id,
    river_code, river, insitu_type, source
  )

# RiverTile column selection (includes p_n_nodes and SWOTFileName):
# save_to_csv <- time_space_matched_SWOT_GNSS %>%
#   dplyr::select(
#     reach_id, time_utc, wse_drift_start_UTC, wse_drift_end_UTC,
#     wse_drift_midpoint_UTC, wse_drift_total_time_UTC,
#     residuals, residuals_nobias, bias,
#     mean_reach_drift_wse_m, mean_reach_drift_wse_total_error_m,
#     mean_reach_drift_wse_no_bias_m,
#     reach_drift_slope_m_m, reach_drift_slope_precision_m, drift_id,
#     wse, wse_u, slope, slope_u, slope_r_u,
#     width, width_u, area_total, area_tot_u, area_detct, area_det_u,
#     area_wse, layovr_val, node_dist, xtrk_dist,
#     reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q,
#     p_dist_out, p_lat, p_lon, p_n_nodes, SWOTFileName,
#     river_code, river, insitu_type, source
#   )

# write.csv(save_to_csv,
#   file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv',
#   row.names = FALSE)


# =============================================================================
# 7. Slope — absolute residuals and plots
# =============================================================================

time_space_matched_SWOT_GNSS$slope_abs              <- abs(time_space_matched_SWOT_GNSS$slope)
time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs <- abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m)

# Remove near-zero GNSS slopes (threshold: 0.1 cm/km = 1e-6 m/m)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  filter(abs(reach_drift_slope_m_m) > 0.000001)

# Compute raw slope residuals: GNSS slope - SWOT slope (both absolute)
time_space_matched_SWOT_GNSS$slope_residuals <- (
  time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs - time_space_matched_SWOT_GNSS$slope_abs
)

# * 100000 converts m/m to cm/km for all slope stats and plots
percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error:", percentile_68_error * 100000))
print(paste("50th Percentile Error:", percentile_50_error * 100000))

# Pearson correlation for slope
cor_test <- cor.test(time_space_matched_SWOT_GNSS$slope_abs, time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs)
r_value  <- cor_test$estimate
p_value  <- cor_test$p.value


# -----------------------------------------------------------------------------
# 7a. Data visualization — slope (absolute)
# -----------------------------------------------------------------------------

# Scatter: absolute (signed) SWOT vs GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_m_m * 100000, y = slope * 100000)) +
  geom_point(size = 4) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 30)

# Scatter: absolute SWOT vs GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(reach_drift_slope_m_m_abs * 100000), y = abs(slope_abs * 100000))) +
  geom_point(size = 4) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS slope (cm/km)") +
  ylab("SWOT slope (cm/km)") +
  annotate("text",
    x     = min(abs(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs * 100000), na.rm = TRUE),
    y     = max(time_space_matched_SWOT_GNSS$slope_abs * 100000, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_GNSS)
    ),
    hjust = 0, vjust = 1, size = 8) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 30)

# CDF: absolute slope difference
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(slope_abs * 100000 - reach_drift_slope_m_m_abs * 100000))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT slope - GNSS slope (cm/km)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT slope - GNSS slope"
  ) +
  annotate("text", x = 100, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error * 100000, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 100, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error * 100000, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# Scatter: residuals vs GNSS slope uncertainty, colored by reach ID
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_precision_m * 100000, y = abs(slope_residuals * 100000), color = factor(reach_id))) +
  geom_point(size = 4) +
  xlab("GNSS slope uncertainty (cm/km)") +
  ylab("Abs slope difference (cm/km)") +
  labs(color = "Reach ID") +
  theme_minimal(base_size = 30)


# =============================================================================
# 8. Slope — bias removal and bias-corrected plots
# =============================================================================

# Bias-corrected slope residuals
time_space_matched_SWOT_GNSS$slope_residuals_nobias <- (
  time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias - time_space_matched_SWOT_GNSS$slope_abs
)

percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$slope_residuals_nobias), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias * 100000))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias * 100000))

# Pearson correlation (bias-corrected slope)
cor_test_nobias <- cor.test(time_space_matched_SWOT_GNSS$slope_abs, time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias)
r_value_nobias  <- cor_test_nobias$estimate
p_value_nobias  <- cor_test_nobias$p.value

# Scatter: SWOT vs bias-corrected GNSS slope
ggplot(time_space_matched_SWOT_GNSS, aes(x = reach_drift_slope_m_m_abs_nobias * 100000, y = slope_abs * 100000, color = factor(drift_id))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS slope (bias-corrected, cm/km)") +
  ylab("SWOT slope (cm/km)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_GNSS$reach_drift_slope_m_m_abs_nobias * 100000, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_GNSS$slope_abs * 100000, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value_nobias, 4),
      "\np value = ", signif(p_value_nobias, 3),
      "\nn = ", nrow(time_space_matched_SWOT_GNSS)
    ),
    hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") +
  theme_minimal(base_size = 30)

# CDF: bias-corrected slope difference
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(slope_abs * 100000 - reach_drift_slope_m_m_abs_nobias * 100000))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT slope - GNSS slope (cm/km)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT slope - GNSS slope (bias-corrected)"
  ) +
  annotate("text", x = 1, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error_nobias * 100000, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 1, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error_nobias * 100000, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# =============================================================================
# 9. Export
# =============================================================================

write.csv(
  save_to_csv,
  file      = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv',
  row.names = FALSE
)
