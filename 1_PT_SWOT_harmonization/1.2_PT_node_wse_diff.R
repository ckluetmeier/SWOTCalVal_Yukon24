# =============================================================================
# PT vs SWOT Node WSE Comparison
# -----------------------------------------------------------------------------
# Matches pressure transducer (PT) WSE observations to SWOT node WSE in time
# (±7.5 min buffer) and space (same node ID), computes residuals, removes
# per-PT median bias, and exports the matched dataset to CSV.
#
# Script sections:
#   1. Read and filter SWOT node data
#   2. Read and prepare PT data
#   3. Match PT and SWOT in time and space
#   4. Summary statistics (absolute comparisons)
#   5. Scatter plot of SWOT vs PT WSE
#   6. Bias removal (for relative comparisons)
#   7. Export per-cluster CSV
#   8. Merge all clusters and export final CalVal CSV
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
# 1. Read and filter SWOT node data
# =============================================================================

# -----------------------------------------------------------------------------
# Select one SWOT product version. Exactly one block below must be active.
#   Version C = RiverSP PIC0, SWORD v16
#   Version D = RiverSP PGD0, SWORD v17b
# -----------------------------------------------------------------------------

# --- Version C: RiverSP PIC0 (SWORD v16) -------------------------------------
# SWOT_df <- read_csv(file.path(
#   DATA_ROOT, "SWOT/node/hydrocron_timeseries",
#   "YR_nodes_merged_RiverSP.csv"))

# --- Version D: RiverSP PGD0 (SWORD v17b) ------------------------------------
SWOT_df <- read_csv(file.path(
  DATA_ROOT, "SWOT/node/RiverSP_v17b",
  "RiverSP_domain_node_timeseries_PGD0_v17b.csv"))

# Remove duplicates and fill values (time = -999…, wse = -1e12)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# Quality filter: node_q < 2 (good/suspect), valid cross-track, <80% dark
# water
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac <= 0.8)

# Convert TAI time (seconds since 2000-01-01, offset 37 s from UTC) to
# POSIXct UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset


# =============================================================================
# 2. Read and prepare PT data
# =============================================================================

# Directory of PT cluster CSV files (output by the PT toolboxes)
# -----------------------------------------------------------------------------
# Select one PT cluster directory.
# -----------------------------------------------------------------------------

# --- SWORD v17b --------------------------------------------------------------
wd <- file.path(
  DATA_ROOT,
  "PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/lower_PR")

# --- SWORD v16 ---------------------------------------------------------------
# wd <- file.path(
#   DATA_ROOT,
#   "PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v16/SJ")

setwd(wd)

csv_files <- list.files(wd, pattern = "\\.csv$", full.names = TRUE)

data_list <- lapply(seq_along(csv_files), function(i) {
  read.csv(csv_files[i])
})
combined_PT_df <- bind_rows(data_list)

# Time column to POSIXct UTC
combined_PT_df$pt_time_UTC <- as.POSIXct(
  combined_PT_df$pt_time_UTC,
  format = "%Y-%m-%d %H:%M:%S",
  tz     = "UTC")


# =============================================================================
# 3. Match PT and SWOT observations in time and space
# =============================================================================

# Temporal match: ±7.5-minute buffer (a SWOT overpass should always fall
# within 7.5 min of a PT observation)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_df_filtered %>%
        filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Spatial match: same node ID
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(Node_ID == node_id)

# Remove any duplicate PT observations (same WSE + time + serial)
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[
  !duplicated(time_space_matched_SWOT_PT[c("pt_wse_m", "pt_time_UTC", "pt_serial")]),]


# =============================================================================
# 4. Absolute summary statistics
# =============================================================================

# Compute absolute residuals: PT WSE − SWOT WSE
time_space_matched_SWOT_PT$residuals <- (
  time_space_matched_SWOT_PT$pt_wse_m - time_space_matched_SWOT_PT$wse)

# percentile_68_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.68, na.rm = TRUE)
# percentile_50_error <- quantile(abs(time_space_matched_SWOT_PT$residuals), 0.50, na.rm = TRUE)
# print(paste("68th Percentile Error:", percentile_68_error))
# print(paste("50th Percentile Error:", percentile_50_error))

# Pearson correlation: SWOT vs PT WSE
cor_test <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_m)
r_value  <- cor_test$estimate  # Pearson r
p_value  <- cor_test$p.value   # p < 0.001 is highly statistically significant


# =============================================================================
# 5. Data visualization - absolute WSE
# =============================================================================

color_palette <- c(
  "#4A4A4A", "#D86A1A", "#6D398B", "#9EBCD8", "#F8A31B",
  "#00429D", "#2E7D32", "#C83232", "#008F7A", "#E3A700", "#124000"
)

# Scatter: SWOT vs PT WSE, colored by PT
ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_m, y = wse, color = factor(pt_serial))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("PT WSE (m)") +
  ylab("SWOT WSE (m)") +
  annotate("text",
    x     = min(time_space_matched_SWOT_PT$pt_wse_m, na.rm = TRUE),
    y     = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", nrow(time_space_matched_SWOT_PT)
    ),
    hjust = 0, vjust = 1, size = 8) +
  labs(color = "PT Serial") +
  theme_minimal(base_size = 30)

# CDF plot (raw residuals):
# ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability",
#        title = "CDF of SWOT WSE - PT WSE") +
#   annotate("text", x = 0.5, y = 0.71,
#     label = paste("|68%ile| diff:", round(percentile_68_error, 4)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 0.5, y = 0.53,
#     label = paste("|50%ile| diff:", round(percentile_50_error, 4)),
#     color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)

# Individual PT vs SWOT hydrograph loop (PDF export):
# pt_serial_list <- unique(time_space_matched_SWOT_PT$pt_serial)
# pdf("PT_SWOT_hydrographs.pdf", width = 7, height = 5)
# for (PT_id in pt_serial_list) {
#   individual_PT_SWOT_df <- time_space_matched_SWOT_PT %>% filter(pt_serial == PT_id)
#   individual_PT_df      <- combined_PT_df %>% filter(pt_serial == PT_id)
#   plot <- ggplot() +
#     geom_errorbar(data = individual_PT_df,
#       aes(x = pt_time_UTC,
#           ymin = pt_wse_m - pt_correction_mean_total_error_m,
#           ymax = pt_wse_m + pt_correction_mean_total_error_m),
#       color = "#edc7a2", linewidth = 5) +
#     geom_point(individual_PT_df, mapping = aes(y = pt_wse_m, x = pt_time_UTC),
#       color = "#ED973D", size = 2) +
#     geom_errorbar(data = individual_PT_SWOT_df,
#       aes(x = time_utc, ymin = wse - wse_u, ymax = wse + wse_u), linewidth = 1) +
#     geom_point(individual_PT_SWOT_df, mapping = aes(y = wse, x = time_utc),
#       color = "#01665E", size = 7) +
#     ggtitle(PT_id) + xlab("Time") + ylab("WSE (m)") +
#     theme_minimal(base_size = 30) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
#   print(plot)
# }
# dev.off()

# Singular hydrograph (one PT):
# PT_id <- 2156918
# individual_PT_SWOT_df <- time_space_matched_SWOT_PT %>% filter(pt_serial == PT_id)
# individual_PT_df      <- combined_PT_df %>% filter(pt_serial == PT_id)
# ggplot() +
#   geom_errorbar(data = individual_PT_df,
#     aes(x = pt_time_UTC,
#         ymin = pt_wse_m - pt_correction_mean_total_error_m,
#         ymax = pt_wse_m + pt_correction_mean_total_error_m),
#     color = "#edc7a2", linewidth = 5) +
#   geom_point(individual_PT_df, mapping = aes(y = pt_wse_m, x = pt_time_UTC),
#     color = "#ED973D", size = 2) +
#   geom_errorbar(data = individual_PT_SWOT_df,
#     aes(x = time_utc, ymin = wse - wse_u, ymax = wse + wse_u), linewidth = 1) +
#   geom_point(individual_PT_SWOT_df, mapping = aes(y = wse, x = time_utc),
#     color = "#01665E", size = 7) +
#   ggtitle(PT_id) + xlab("Time") + ylab("WSE (m)") +
#   theme_minimal(base_size = 30) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0.9))


# =============================================================================
# 6. Bias removal - per-sensor median bias (min 3 observations)
# =============================================================================

time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT %>%
  group_by(pt_serial) %>%
  mutate(
    bias            = if (n() >= 3) median(residuals, na.rm = TRUE) else NA_real_,
    pt_wse_nobias_m = if (n() >= 3) pt_wse_m - bias else NA_real_
  ) %>%
  ungroup()

# Compute bias-corrected residuals
time_space_matched_SWOT_PT$residuals_nobias <- (
  time_space_matched_SWOT_PT$pt_wse_nobias_m - time_space_matched_SWOT_PT$wse
)

# percentile_68_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.68, na.rm = TRUE)
# percentile_50_error_nobias <- quantile(abs(time_space_matched_SWOT_PT$residuals_nobias), 0.50, na.rm = TRUE)
# print(paste("68th Percentile Error Without Bias:", percentile_68_error_nobias))
# print(paste("50th Percentile Error Without Bias:", percentile_50_error_nobias))

# cor_test_nobias <- cor.test(time_space_matched_SWOT_PT$wse, time_space_matched_SWOT_PT$pt_wse_nobias_m)
# r_value_nobias  <- cor_test_nobias$estimate
# p_value_nobias  <- cor_test_nobias$p.value
#
# ggplot(time_space_matched_SWOT_PT, aes(x = pt_wse_nobias_m, y = wse, color = factor(pt_serial))) +
#   geom_point(size = 4) +
#   scale_color_manual(values = color_palette) +
#   geom_abline(linetype = "dashed", color = "gray") +
#   xlab("PT WSE (bias-corrected, m)") + ylab("SWOT WSE (m)") +
#   annotate("text",
#     x = min(time_space_matched_SWOT_PT$pt_wse_m, na.rm = TRUE),
#     y = max(time_space_matched_SWOT_PT$wse, na.rm = TRUE),
#     label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3),
#                    "\nn = ", nrow(time_space_matched_SWOT_PT)),
#     hjust = 0, vjust = 1, size = 8) +
#   labs(color = "PT Serial") +
#   theme_minimal(base_size = 30)
#
# ggplot(time_space_matched_SWOT_PT, aes(x = abs(wse - pt_wse_nobias_m))) +
#   stat_ecdf(geom = "step", color = "darkblue", size = 1) +
#   geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
#   labs(x = "SWOT WSE - PT WSE (m)", y = "Cumulative Probability",
#        title = "CDF of SWOT WSE - PT WSE (bias-corrected)") +
#   annotate("text", x = 0.25, y = 0.71,
#     label = paste("68% abs diff:", round(percentile_68_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   annotate("text", x = 0.25, y = 0.53,
#     label = paste("50% abs diff:", round(percentile_50_error_nobias, 4)),
#     color = "#222222", size = 6) +
#   theme_minimal(base_size = 20)


# =============================================================================
# 7. Export per-cluster matched dataset
# =============================================================================

save_to_csv <- time_space_matched_SWOT_PT %>%
  dplyr::select(
    time_utc, pt_time_UTC, residuals, residuals_nobias, wse, wse_u,
    pt_wse_m, pt_wse_nobias_m, bias, pt_correction_mean_total_error_m,
    pt_correction_mean_offset_sd_m, pt_serial, width, width_u, node_id, reach_id,
    p_dist_out, node_q, node_q_b, dark_frac, n_good_pix, rdr_sig0, xovr_cal_q,
    lat, lon)

# Output directory and filename must match the selected version above
write.csv(
  save_to_csv,
  file      = file.path(
    DATA_ROOT, "CalVal_dataframes/wse/node/RiverSP_v17b",
    "RiverSP_v17b_time_space_matched_SWOT_PT_lowerPR.csv"),
  row.names = FALSE)


# =============================================================================
# 8. Merge all per-cluster CSVs into one CalVal node dataset
# =============================================================================

# Directory containing all per-cluster matched CSV files
# -----------------------------------------------------------------------------
# Select one directory.
# -----------------------------------------------------------------------------

# --- Version C: RiverSP PIC0 (SWORD v16) -------------------------------------
# wd <- file.path(DATA_ROOT, "CalVal_dataframes/wse/node/RiverSP_v16")

# --- Version D: RiverSP PGD0 (SWORD v17b) ------------------------------------
wd <- file.path(DATA_ROOT, "CalVal_dataframes/wse/node/RiverSP_v17b")

setwd(wd)

csv_files <- list.files(wd, pattern = "time_space_matched_SWOT_PT.*\\.csv$", full.names = TRUE)

# Read each file and tag with the river cluster name derived from the filename
data_list <- lapply(csv_files, function(file) {
  df             <- read.csv(file)
  filename_clean <- sub(".*PT_(.*)\\.csv$", "\\1", basename(file))
  df %>% mutate(river = filename_clean)
})

combined_time_space_matched_SWOT_PT_df <- bind_rows(data_list)

RiverSP_df <- combined_time_space_matched_SWOT_PT_df %>%
  mutate(source = "PGD0")  # Must match the selected version above

# Export merged node CalVal CSV
# Output directory must match the selected version above
write.csv(
  RiverSP_df,
  file      = file.path(
    DATA_ROOT, "CalVal_dataframes/wse/node/RiverSP_v17b",
    "node_SWOT_PT.csv"),
  row.names = FALSE)
