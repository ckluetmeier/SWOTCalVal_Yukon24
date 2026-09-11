# =============================================================================
# PT vs SWOT Reach WSE and Slope Comparison
# -----------------------------------------------------------------------------
# Matches PT reach WSE and slope observations to SWOT reach data in time
# (±7.5 min buffer) and space (same reach ID). Computes residuals, removes
# per-reach median bias, and exports matched datasets to CSV.
#
# Script sections:
#   PART A — Reach WSE
#     1.  Read and filter SWOT reach data
#     2.  Read and prepare PT reach WSE data
#     3.  Match PT and SWOT in time and space
#     4.  Raw residuals and scatter / CDF plots
#     5.  Bias removal (per-reach median) and export
#
#   PART B — Reach Slope
#     6.  Read SWOT data
#     7.  Read and prepare PT reach slope data (with bias-corrected slopes)
#     8.  Match, compute slope residuals, and export
#     9.  Slope scatter / CDF plots
# 
# -----------------------------------------------------------------------------
# Script by:
# Camryn Kluetmeier (camryn.kluetmeier@duke.edu)
# 
# Parts of this script were developed with assistance from Claude Code 
# (Anthropic) for debugging, documentation, and related editorial suggestions.
# 
# Last updated: 2026-09-11
# 
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)


# =============================================================================
# Configuration - edit these paths before running
# =============================================================================

# Root of the field-campaign and SWOT data products.
DATA_ROOT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats"


# =============================================================================
# PART A — REACH WSE
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Read and filter SWOT reach data
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Select one SWOT product version.
#   Version C = RiverSP PIC0, SWORD v16
#   Version D = RiverSP PGD0, SWORD v17b
# -----------------------------------------------------------------------------

# --- Version C: RiverSP PIC0 (SWORD v16) -------------------------------------
# SWOT_reach_df <- read.csv(file.path(
#   DATA_ROOT, "SWOT/reach/RiverSP_v16",
#   "RiverSP_domain_reach_timeseries_v16.csv"))

# --- Version D: RiverSP PGD0 (SWORD v17b) ------------------------------------
SWOT_reach_df <- read_csv(file.path(
  DATA_ROOT, "SWOT/reach/RiverSP_v17b",
  "RiverSP_domain_reach_timeseries_PGD0_v17b.csv"))

# Remove duplicate rows / fill-value rows (time = -999…, wse = -1e12)
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE)

# Quality filter: reach_q < 2 (good/suspect), valid cross-track swath,
# >50% node coverage, <80% dark water
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(partial_f == 0) %>%   # partial_f == 0 means >= 50% node coverage
  filter(dark_frac <= 0.8)

# Convert TAI time (seconds since 2000-01-01) to UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset


# -----------------------------------------------------------------------------
# 2. Read and prepare PT reach data
# -----------------------------------------------------------------------------

# Directory of PT reach CSV files (output by the PT toolboxes)
# -----------------------------------------------------------------------------

# --- SWORD v16 ---------------------------------------------------------------
# wd <- file.path(
#   DATA_ROOT,
#   "PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v16")

# --- SWORD v17b --------------------------------------------------------------
wd <- file.path(
  DATA_ROOT,
  "PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b")

setwd(wd)

# Only load files matching the reach-level naming pattern
csv_files <- list.files(wd, pattern = "^YR_812.*\\.csv$", full.names = TRUE)

data_list <- lapply(seq_along(csv_files), function(i) {
  read.csv(csv_files[i])
})
# Ensure sorted_nodelist is character type to prevent bind_rows type conflict
data_list <- lapply(data_list, function(df) {
  df %>% mutate(sorted_nodelist = as.character(sorted_nodelist))
})
combined_PT_df <- bind_rows(data_list)

combined_PT_df$pt_time_UTC <- as.POSIXct(combined_PT_df$pt_time_UTC, tz = "UTC")

# Rename reach_id to avoid column name conflict with SWOT reach_id
combined_PT_df <- rename(combined_PT_df, PT_reach_id = reach_id)


# -----------------------------------------------------------------------------
# 3. Match PT and SWOT in time and space
# -----------------------------------------------------------------------------

# Temporal match: ±7.5-minute buffer (a SWOT overpass should always fall
# within 7.5 min of a PT observation)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_reach_df_filtered %>%
        filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Spatial match: same reach ID
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(PT_reach_id == reach_id)

# Remove duplicate PT values (same WSE + time + flaglist)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[
  !duplicated(time_space_matched_SWOT_PT[c("mean_reach_pt_wse_m", "pt_time_UTC", "flaglist")]),]

# Add river labels and metadata columns
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # Reach IDs are checked first so case_when does not overwrite the
      # SJ/BL labels
      # SJ reach IDs differ between SWORD versions:
      #   SWORD v16:  81260300061, 81260300231, 81260300241, 81260300251
      #   SWORD v17b: 81260300181, 81260300191, 81260300201, 81260300211
      # Must match the selected version above
      reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                       "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",  # downstream to Circle bifurcation
      river_code == "812509" ~ "lowerYR",  # past the Porcupine confluence
      river_code == "812705" ~ "upperYR",  # Circle and above
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0")  # Must match the selected version above


# -----------------------------------------------------------------------------
# 4. Raw residuals and summary statistics
# -----------------------------------------------------------------------------

# Compute raw residuals: PT WSE − SWOT WSE
time_space_matched_SWOT_PT$residuals <- (
  time_space_matched_SWOT_PT$mean_reach_pt_wse_m - time_space_matched_SWOT_PT$wse
)

# Remove physically implausible residuals (> 10 m absolute difference)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  filter(abs(residuals) < 10)

# Percentile error
percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error:", percentile_68_error))
print(paste("50th Percentile Error:", percentile_50_error))

# model <- lm(mean_reach_pt_wse_m ~ wse, data = time_space_matched_SWOT_PT)
# summary(model)
# rmse <- sqrt(mean((time_space_matched_SWOT_PT$mean_reach_pt_wse_m - time_space_matched_SWOT_PT$wse)^2))
# RMSE >= MAE; MAE ≈ 50th-quantile error

# Pearson correlation
cor_test <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$mean_reach_pt_wse_m)
r_value  <- cor_test$estimate  # Pearson r
p_value  <- cor_test$p.value   # p < 0.001 is highly statistically significant


# -----------------------------------------------------------------------------
# 4a. Data visualization - absolute WSE
# -----------------------------------------------------------------------------

# River color palette
color_palette <- c(
  "#D86A1A", "#6D398B", "#F8A31B", "#00429D", "#2E7D32",
  "#C83232", "#008F7A", "#E3A700", "#124000"
)

# Scatter: SWOT vs PT WSE, colored by river
ggplot(time_space_matched_SWOT_PT, aes(x = mean_reach_pt_wse_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("PT WSE (m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_PT$mean_reach_pt_wse_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_PT)
    ),
    hjust = 0, vjust = 1, size = 8) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# Basin-code color palette
color_palette <- c(
  "#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8", "#F8A31B",
  "#00429D", "#2E7D32", "#C83232", "#008F7A", "#E3A700", "#124000"
)

# Scatter: SWOT vs PT WSE, coloured by basin code prefix
ggplot(time_space_matched_SWOT_PT, aes(x = mean_reach_pt_wse_m, y = wse)) +
  geom_point(aes(color = factor(substr(reach_id, 1, 6))), size = 4) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +
  xlab("PT WSE (m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_PT$mean_reach_pt_wse_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_PT)
    ),
    hjust = 0, vjust = 1, size = 8) +
  labs(color = "Basin code") +
  theme_minimal(base_size = 30)

# CDF: absolute WSE difference
ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - mean_reach_pt_wse_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT WSE - PT WSE (m)",
    y     = "Cumulative Probability",
    title = "CDF of SWOT WSE - PT WSE"
  ) +
  annotate("text", x = 2, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 2, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# -----------------------------------------------------------------------------
# 5. Bias removal - per-reach median bias (min 3 observations) and export
# -----------------------------------------------------------------------------

time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  group_by(reach_id) %>%
  mutate(
    bias            = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    pt_wse_nobias_m = if (n() >= 3) mean_reach_pt_wse_m - bias else NA_real_
  ) %>%
  ungroup()

# Bias-corrected residuals
time_space_matched_SWOT_PT$residuals_nobias <- (
  time_space_matched_SWOT_PT$pt_wse_nobias_m - time_space_matched_SWOT_PT$wse
)

percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# Export reach WSE CalVal CSV
save_to_csv <- time_space_matched_SWOT_PT %>%
  dplyr::select(
    time_utc, pt_time_UTC, reach_id, residuals, residuals_nobias, bias, wse, wse_u,
    mean_reach_pt_wse_m, pt_wse_nobias_m, flaglist, sorted_nodelist, Number_of_nodes,
    slope, slope_u, slope_r_u, width, width_u, area_total, area_tot_u, area_detct,
    area_det_u, area_wse, layovr_val, node_dist,
    xtrk_dist, reach_q, reach_q_b, dark_frac, n_good_nod, partial_f, xovr_cal_q,
    p_dist_out, p_lat, p_lon, cycle_id, pass_id,
    source, river_code, river, insitu_type  # cycle_id, pass_id OR p_n_nodes
  )

# Output directory must match the selected version above
write.csv(
  save_to_csv,
  file      = file.path(
    DATA_ROOT, "CalVal_dataframes/wse/reach/RiverSP_v17b",
    "reach_wse_SWOT_PT.csv"),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 5a. Data visualization - bias-corrected WSE
# -----------------------------------------------------------------------------

cor_test_nobias <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_nobias_m)
r_value_nobias  <- cor_test_nobias$estimate
p_value_nobias  <- cor_test_nobias$p.value

ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_nobias_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +
  xlab("PT WSE (bias-corrected, m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_PT$pt_wse_nobias_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value_nobias, 4),
      "\np value = ", signif(p_value_nobias, 3),
      "\nn = ", nrow(time_space_matched_SWOT_PT)
    ),
    hjust = 0, vjust = 1, size = 8) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# Shared river factor levels and color palette used across all
# inter-river plots
river_levels  <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
river_labels  <- c("Coleen", "Sheenjek", "Chandalar", "Porcupine",
                   "Single-channel Yukon", "Braided Yukon")
color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")


ggplot(time_space_matched_SWOT_PT,
       aes(x = river, y = bias * 100, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
            aes(x = river, y = -0.2, label = paste0("n=", n)),
            inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab("SWOT -" ~ italic("in situ") ~ "Slope (cm/km)") +
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25)

# ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_nobias_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(
#     x     = "SWOT WSE - PT WSE (m)",
#     y     = "Cumulative Probability",
#     title = "CDF of SWOT WSE - PT WSE (bias-corrected)"
#   ) +
#   annotate("text", x = 2, y = 0.71,
#     label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 2, y = 0.53,
#     label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)


# =============================================================================
# PART B — REACH SLOPE
# =============================================================================


# -----------------------------------------------------------------------------
# 6. Read and filter SWOT reach data for slope
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Select one SWOT product version
#   Version C = RiverSP PIC0, SWORD v16
#   Version D = RiverSP PGD0, SWORD v17b
# -----------------------------------------------------------------------------

# --- Version C: RiverSP PIC0 (SWORD v16) -------------------------------------
SWOT_reach_df <- read_csv(file.path(
  DATA_ROOT, "SWOT/reach/RiverSP_v16",
  "RiverSP_domain_reach_timeseries_v16.csv"))

# --- Version D: RiverSP PGD0 (SWORD v17b) ------------------------------------
# SWOT_reach_df <- read_csv(file.path(
#   DATA_ROOT, "SWOT/reach/RiverSP_v17b",
#   "RiverSP_domain_reach_timeseries_PGD0_v17b.csv"))

SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE)

SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(partial_f == 0) %>%
  filter(dark_frac <= 0.8)

tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37

SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset


# -----------------------------------------------------------------------------
# 7. Read and prepare PT reach slope data
# -----------------------------------------------------------------------------

# The bias-corrected slope CSV is written by
# 0_fetch_SWOT_data/0.3_PT_GNSS_reach_slope_bias_corr.R
# -----------------------------------------------------------------------------
# Select one PT reach slope pair. Must match the SWORD version selected in
# section 6.
# -----------------------------------------------------------------------------

# --- SWORD v16 ---------------------------------------------------------------
PT_reach_df <- read_csv(file.path(
  DATA_ROOT, "PTs/toolboxes_dataframes/reprocessed_2025_09_02",
  "_reach/SWORD_v16/YR_PT_reach_slope.csv"))
PT_reach_corrected_df <- read_csv(file.path(
  DATA_ROOT, "PTs/toolboxes_dataframes/reprocessed_2025_09_02",
  "_reach/SWORD_v16/YR_PT_reach_slope_corrected.csv")) %>%
  rename(mean_reach_PT_slope_no_bias_m_m = slope_m_m)

# --- SWORD v17b --------------------------------------------------------------
# PT_reach_df <- read_csv(file.path(
#   DATA_ROOT, "PTs/toolboxes_dataframes/reprocessed_2025_09_02",
#   "_reach/SWORD_v17b/YR_PT_reach_slope.csv"))
# PT_reach_corrected_df <- read_csv(file.path(
#   DATA_ROOT, "PTs/toolboxes_dataframes/reprocessed_2025_09_02",
#   "_reach/SWORD_v17b/YR_PT_reach_slope_corrected.csv")) %>%
#   rename(mean_reach_PT_slope_no_bias_m_m = slope_m_m)

# Join corrected slopes back to raw reach data
PT_reach_df <- PT_reach_df %>%
  inner_join(
    PT_reach_corrected_df %>%
      select(
        pt_time_UTC, total_error_pt_wse_us_boundary, pt_serials_us,
        reach_id, mean_reach_PT_slope_no_bias_m_m, bias_us, bias_ds, river
      ),
    by = c("pt_time_UTC", "total_error_pt_wse_us_boundary", "pt_serials_us", "reach_id")
  )

# Parse time to POSIXct (UTC). The format string below must match the
# datetime format in the PT reach slope CSVs, here "%m/%d/%y %H:%M".
PT_reach_df$pt_time_UTC <- as.POSIXct(
  PT_reach_df$pt_time_UTC,
  format = "%m/%d/%y %H:%M",
  tz     = "UTC"
)

# Rename reach_id to avoid conflict with SWOT reach_id
PT_reach_df <- rename(PT_reach_df, PT_reach_id = reach_id)


# -----------------------------------------------------------------------------
# 8. Match, compute slope residuals, and export
# -----------------------------------------------------------------------------

# Temporal match: ±7.5-minute buffer (a SWOT overpass should always fall
# within 7.5 min of a PT observation)
time_matched_SWOT_PT_reach <- PT_reach_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_reach_df_filtered %>%
        filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Spatial match: same reach ID
time_space_matched_SWOT_PT_reach <- time_matched_SWOT_PT_reach %>%
  filter(PT_reach_id == reach_id)

# Remove duplicate PT values (same upstream WSE + time + serials)
time_space_matched_SWOT_PT_reach <- time_space_matched_SWOT_PT_reach[
  !duplicated(time_space_matched_SWOT_PT_reach[c("mean_pt_wse_us_boundary_m", "pt_time_UTC", "pt_serials_us")]),
]

# Force absolute value for slopes
time_space_matched_SWOT_PT_reach$slope_abs     <- abs(time_space_matched_SWOT_PT_reach$slope)
time_space_matched_SWOT_PT_reach$slope_m_m_abs <- abs(time_space_matched_SWOT_PT_reach$slope_m_m)

# Raw slope residuals: PT slope − SWOT slope (m/m)
time_space_matched_SWOT_PT_reach$residuals <- (
  time_space_matched_SWOT_PT_reach$slope_m_m_abs - time_space_matched_SWOT_PT_reach$slope_abs)

percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT_reach$residuals), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT_reach$residuals), 0.50, na.rm = TRUE)

# *100000 converts m/m to cm/km
print(paste("68th Percentile Error:", percentile_68_error * 100000))
print(paste("50th Percentile Error:", percentile_50_error * 100000))

# Bias-corrected slope residuals: bias-corrected PT slope - SWOT slope
time_space_matched_SWOT_PT_reach$slope_residuals_nobias <- (
  time_space_matched_SWOT_PT_reach$mean_reach_PT_slope_no_bias_m_m - time_space_matched_SWOT_PT_reach$slope_abs)

percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT_reach$slope_residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT_reach$slope_residuals_nobias), 0.50, na.rm = TRUE)

# Add river labels and version metadata
save_to_csv <- time_space_matched_SWOT_PT_reach %>%
  select(
    pt_time_UTC, time_utc, reach_id,
    slope_m_m, slope_m_m_abs, slope_uncertainty_m_m, slope_residuals_nobias,
    p_lat, p_lon, slope, slope_abs, slope_u, bias_us, bias_ds,
    slope_residuals_nobias, mean_reach_PT_slope_no_bias_m_m,
    residuals, width, width_u, area_total, area_tot_u, layovr_val, node_dist,
    xtrk_dist, reach_q, reach_q_b, dark_frac, xovr_cal_q, cycle_id, pass_id, river
    # cycle_id, pass_id OR p_dist_out, n_good_nod
  ) %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PIC0") %>%  # Must match the selected version above
  rename(slope_residuals = residuals)

# Add river labels
save_to_csv <- save_to_csv %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # SJ reach IDs differ between SWORD versions:
      #   SWORD v16:  81260300061, 81260300231, 81260300241, 81260300251
      #   SWORD v17b: 81260300181, 81260300191, 81260300201, 81260300211
      # Must match the selected version above
      reach_id %in% c("81260300061", "81260300231", "81260300241", "81260300251") ~ "SJ",
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

# Output directory must match the selected version above
write.csv(
  save_to_csv,
  file      = file.path(
    DATA_ROOT, "CalVal_dataframes/wse/reach/RiverSP_v16",
    "reach_slope_SWOT_PT.csv"),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 9. Slope diagnostic plots
# -----------------------------------------------------------------------------

# Full color palette
# color_palette <- c(
#   "#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8",
#   "#E3A700", "#008F7A", "#C83232", "#2E7D32",
#   "#D81B60", "#00429D", "#A6761D", "#56B4E9",
#   "#4c64c1", "#7ca92f", "#9a3c9a", "orange",
#   "lightyellow", "yellow", "pink", "black"
# )
#
# Regression and correlation stats
# model   <- lm(slope_m_m ~ slope, data = time_space_matched_SWOT_PT_reach)
# summary(model)
# rmse    <- sqrt(mean((time_space_matched_SWOT_PT_reach$slope_m_m - time_space_matched_SWOT_PT_reach$slope)^2))
# RMSE >= MAE; MAE ≈ 50th-quantile error
#
# cor_test <- cor.test(time_space_matched_SWOT_PT_reach$slope, time_space_matched_SWOT_PT_reach$slope_m_m)
# r_value  <- cor_test$estimate
# p_value  <- cor_test$p.value
#
# Scatter: bias-corrected PT slope vs SWOT slope
# ggplot(time_space_matched_SWOT_PT_reach,
#        aes(x = mean_reach_PT_slope_no_bias_m * 100000, y = slope_abs * 100000, color = factor(reach_id))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   geom_abline(linetype = "dashed", color = "gray") +
#   xlab("PT slope (cm/km)") +
#   ylab("SWOT slope (cm/km)") +
#   annotate("text",
#     x     = min(time_space_matched_SWOT_PT_reach$slope_m_m * 100000, na.rm = TRUE),
#     y     = max(time_space_matched_SWOT_PT_reach$slope * 100000, na.rm = TRUE),
#     label = paste0(
#       "r = ", round(r_value, 4),
#       "\np value = ", signif(p_value, 3),
#       "\nn = ", nrow(time_space_matched_SWOT_PT_reach)
#     ),
#     hjust = 0, vjust = 1, size = 8) +
#   labs(color = "Reach ID") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme_minimal(base_size = 30)
#
# CDF: raw slope difference
# ggplot(time_space_matched_SWOT_PT_reach, aes(x = abs(slope_abs - slope_m_m_abs))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(
#     x     = "SWOT slope - PT slope (m/m)",
#     y     = "Cumulative Probability",
#     title = "CDF of SWOT slope - PT slope"
#   ) +
#   annotate("text", x = 0.00018, y = 0.71,
#     label = paste("68% abs diff:", round(percentile_68_error * 100000, 8)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 0.00018, y = 0.53,
#     label = paste("50% abs diff:", round(percentile_50_error * 100000, 8)),
#     color = "#222222", size = 6) +
#   xlim(0, 0.00025) +
#   theme_minimal(base_size = 20)
#
# CDF: bias-corrected slope difference (cm/km)
# ggplot(time_space_matched_SWOT_PT_reach,
#        aes(x = abs(slope_abs - mean_reach_PT_slope_no_bias_m_m) * 100000)) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(
#     x     = "SWOT slope - PT slope (cm/km)",
#     y     = "Cumulative Probability",
#     title = "CDF of SWOT slope - PT slope (bias-corrected)"
#   ) +
#   annotate("text", x = 7.4, y = 0.71,
#     label = paste("68% abs diff:", round(percentile_68_error_nobias * 100000, 8)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 7.4, y = 0.53,
#     label = paste("50% abs diff:", round(percentile_50_error_nobias * 100000, 8)),
#     color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)
#
# Scatter: slope uncertainty vs absolute slope difference
# ggplot(time_space_matched_SWOT_PT_reach,
#        aes(x = slope_uncertainty_m_m, y = abs(residuals), color = factor(reach_id))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   xlab("PT slope uncertainty (m/m)") +
#   ylab("Absolute slope difference (m/m)") +
#   labs(color = "Reach ID") +
#   theme_minimal(base_size = 30)
