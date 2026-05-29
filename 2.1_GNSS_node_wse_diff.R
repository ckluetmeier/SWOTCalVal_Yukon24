# =============================================================================
# GNSS vs SWOT Node WSE Comparison
# -----------------------------------------------------------------------------
# Matches GNSS drift survey WSE observations to SWOT node WSE in time
# (±5 hour buffer) and space (same node ID), computes residuals, removes
# per-drift median bias, and exports the matched dataset to CSV.
#
# Script sections:
#   1.  Read and filter SWOT node data
#   2.  Read and prepare GNSS data
#   3.  Match GNSS and SWOT in time and space
#   4.  Raw residuals, summary table, and scatter / CDF plots
#   5.  Bias removal (per-drift median) and export
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)


# =============================================================================
# 1. Read and filter SWOT node data
# =============================================================================

# --- RiverSP Version C -------------------------------------------------------
# SWOT_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv")

# --- RiverSP Version D (SWORD v17b) ------------------------------------------
SWOT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv")

# --- RiverTile ---------------------------------------------------------------
# SWORD v16
# SWOT_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v16/RiverTile_domain_node_timeseries_v16.csv")
# SWORD v17b
# SWOT_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverTile_v17b/RiverTile_domain_node_timeseries_v17b.csv")

# Remove duplicates and sentinel fill values (time = -999…, wse = -1e12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# Quality filter: node_q < 2, valid cross-track swath, < 80% dark water
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>%           # 0=good, 1=suspect, 2=degraded, 3=bad
  filter(abs(xtrk_dist) >= 10000) %>%  # cross-track 10–60 km
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac <= 0.80)

# Convert TAI time to UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset


# =============================================================================
# 2. Read and prepare GNSS data
# =============================================================================

# --- SWORD v16 ---------------------------------------------------------------
# GNSS_df <- read_csv(
#   "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v16/YR_drift_node_wses.csv")

# --- SWORD v17b --------------------------------------------------------------
GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv")

GNSS_df$time_UTC <- as.POSIXct(GNSS_df$time_UTC, tz = "UTC")

# Rename node_id to avoid column conflict with SWOT node_id
GNSS_df <- rename(GNSS_df, Node_ID = node_id)


# =============================================================================
# 3. Match GNSS and SWOT in time and space
# =============================================================================

# Temporal match: ±5-hour buffer
time_matched_SWOT_GNSS <- GNSS_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_df_filtered %>%
        filter(abs(difftime(time_UTC, time_utc, units = "hours")) <= 5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Spatial match: same node ID
time_space_matched_SWOT_GNSS <- time_matched_SWOT_GNSS %>%
  filter(Node_ID == node_id)


# =============================================================================
# 4. Absolute residuals and summary statistics
# =============================================================================

# Compute absolute residuals: GNSS WSE - SWOT WSE
time_space_matched_SWOT_GNSS$residuals <- (
  time_space_matched_SWOT_GNSS$mean_node_drift_wse_m - time_space_matched_SWOT_GNSS$wse
)

# Per-drift summary (all drifts)
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

# Remove implausible residuals (> 3 m absolute difference)
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  filter(abs(residuals) < 3)

percentile_68_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_GNSS$residuals), 0.50, na.rm = TRUE)

# percentile prints:
# print(paste("68th Percentile Error:", percentile_68_error))
# print(paste("50th Percentile Error:", percentile_50_error))

# Pearson correlation
cor_test <- cor.test(time_space_matched_SWOT_GNSS$wse, time_space_matched_SWOT_GNSS$mean_node_drift_wse_m)
r_value  <- cor_test$estimate  # Pearson r
p_value  <- cor_test$p.value   # p < 0.001 is highly statistically significant


# -----------------------------------------------------------------------------
# 4a. Data visualization - absolute WSE
# -----------------------------------------------------------------------------

color_palette <- c(
  "#e73582", "#78e964", "#5329a8", "#e1df28", "#8c66f0", "#a8dc3e", "#7150ce",
  "#48b32e", "#cb4bd0", "#5ddd76", "#df36ad", "#4e227d", "#c7e772", "#6549a4",
  "#d4c841", "#4c64c1", "#7ca92f", "#9a3c9a", "#409f46", "#b67ce0", "#95d47e",
  "#bc4797", "#df3542", "#6389ee", "#d99427", "#dd80d5", "#417a28", "#de6096",
  "#982161", "#c87735", "#c9405b", "#de4b21", "#b84b34", "pink", "blue",
  "lightblue", "darkblue", "darkred", "gray"
)

# Scatter: SWOT vs GNSS WSE, colored by drift ID
ggplot(time_space_matched_SWOT_GNSS, aes(x = mean_node_drift_wse_m, y = wse, color = factor(drift_id))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("GNSS WSE (m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_GNSS$mean_node_drift_wse_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_GNSS)
    ),
    hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") +
  # labs(color = "Drift ID")  # uncomment to show legend
  theme_minimal(base_size = 30)

# CDF: absolute WSE difference
ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_node_drift_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT WSE - GNSS WSE (m)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT WSE - GNSS WSE"
  ) +
  annotate("text", x = 2, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)
# xlim(0, 3)

# residuals vs GNSS node total error
# ggplot(time_space_matched_SWOT_GNSS,
#        aes(x = abs(residuals), y = node_total_error_m,
#            color = factor(substr(drift_id, 70, 93)))) +
#   geom_point(size = 1.5) +
#   scale_color_manual(values = color_palette) +
#   xlab("GNSS - SWOT WSE (m)") + ylab("WSE uncertainty (m)") +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none")


# =============================================================================
# 5. Bias removal — per-drift median bias (min 3 observations) and export
# =============================================================================

time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
  group_by(drift_id) %>%
  mutate(
    bias                          = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    mean_node_drift_wse_no_bias_m = if (n() >= 3) mean_node_drift_wse_m - bias else NA_real_
  ) %>%
  ungroup()

# Bias-corrected residuals
time_space_matched_SWOT_GNSS$residuals_nobias <- (
  time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse
)

percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_GNSS$residuals_nobias), 0.50, na.rm = TRUE)

# print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
# print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))


# Add river name labels to data frame
time_space_matched_SWOT_GNSS <- time_space_matched_SWOT_GNSS %>%
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


# CSV export
save_to_csv <- time_space_matched_SWOT_GNSS %>%
  dplyr::select(
    time_utc, time_UTC, residuals, residuals_nobias, bias, wse, wse_u,
    mean_node_drift_wse_m, mean_node_drift_wse_no_bias_m, node_total_error_m, drift_id,
    width, width_u, node_id, reach_id, p_dist_out,
    node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q,
    cycle_id, pass_id, lat, lon, source, river_code, river
  )

write.csv(
  save_to_csv,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv",
  row.names = FALSE
)

# # Bias-corrected correlation and plots:
# cor_test_nobias <- cor.test(
#   time_space_matched_SWOT_GNSS$wse,
#   time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m
# )
# r_value_nobias <- cor_test_nobias$estimate
# p_value_nobias <- cor_test_nobias$p.value
# 
# # regression / RMSE:
# # model_nobias <- lm(mean_node_drift_wse_no_bias_m ~ wse, data = time_space_matched_SWOT_GNSS)
# # summary(model_nobias)
# # rmse_nobias <- sqrt(mean((time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m - time_space_matched_SWOT_GNSS$wse)^2))
# 
# # Scatter: SWOT vs bias-corrected GNSS WSE
# ggplot(time_space_matched_SWOT_GNSS,
#        aes(x = mean_node_drift_wse_no_bias_m, y = wse, color = factor(drift_id))) +
#   geom_point(size = 2) +
#   scale_color_manual(values = color_palette) +
#   geom_abline(linetype = "dashed", color = "gray") +
#   xlab("GNSS WSE (bias-corrected, m)") +
#   ylab("SWOT WSE (m)") +
#   # NOTE: annotation uses r_value / p_value (raw), not r_value_nobias; update if needed
#   annotate("text",
#     x     = min(time_space_matched_SWOT_GNSS$mean_node_drift_wse_no_bias_m, na.rm = TRUE),
#     y     = max(time_space_matched_SWOT_GNSS$wse, na.rm = TRUE),
#     label = paste0(
#       "r = ", round(r_value, 4),
#       "\np value = ", signif(p_value, 3),
#       "\nn = ", nrow(time_space_matched_SWOT_GNSS)
#     ),
#     hjust = 0, vjust = 1, size = 8) +
#   theme(legend.position = "none") +
#   theme_minimal(base_size = 30)
# 
# # CDF: bias-corrected absolute WSE difference
# ggplot(time_space_matched_SWOT_GNSS, aes(x = abs(wse - mean_node_drift_wse_no_bias_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(
#     x     = "SWOT WSE - GNSS WSE (m)",
#     y     = "Cumulative Probability",
#     title = "CDF of SWOT WSE - GNSS WSE (bias-corrected)"
#   ) +
#   annotate("text", x = 2, y = 0.71,
#     label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 2, y = 0.53,
#     label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)


