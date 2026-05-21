# =============================================================================
# Orthomosaic vs SWOT Reach Width Comparison
# -----------------------------------------------------------------------------
# Reads per-cluster orthomosaic node-width CSVs (output of script 3.3),
# joins to SWORD nodes to obtain node_len, aggregates to reach-level
# mean widths, matches to SWOT reach observations on the same date, computes
# residuals and percent differences, and produces summary tables and plots.
#
# Script sections:
#   1.  Read and prepare orthomosaic node width data
#   2.  Read and filter SWOT reach data
#   3.  Match ortho and SWOT observations in time and space
#   4.  Raw residuals and percent difference
#   5.  Data visualization
#   6.  Summary table by river
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)


# =============================================================================
# 1. Read and prepare orthomosaic node width data
# =============================================================================

# --- RiverSP Version D (SWORD v17b) ------------------------------------------
dir_path <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b"

# Bring in SWORD node shapefile (provides node_len for width calculation)
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")

# Reaches where orthomosaic coverage is complete enough for reach-level comparison
viable_reaches <- data.frame(reach_id = c(
  81260300191, 81260300181, 81260400021, 81260400011,
  81260500011, 81260300171, 81260300161,
  81250800031, 81270100041, 81270100051, 81270100061,
  81270500161, 81270500171))

# --- RiverSP Version C (SWORD v16) -------------------------------------------
# dir_path  <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16"
# sword_sf  <- st_read("/Users/camryn/Desktop/SWORD_v16/NA/na_sword_nodes_hb81_v16.shp")
# viable_reaches <- data.frame(reach_id = c())  # need v16 reaches

# Read all per-cluster ortho CSVs and assign an overpass date based on filename
ortho_files <- list.files(path = dir_path, pattern = "^ortho_.*\\.csv$", full.names = TRUE)

# Lookup function: map filename substrings to SWOT overpass dates
assign_swot_date <- function(filename) {
  base <- basename(filename)
  dplyr::case_when(
    grepl("lowerPR_SJ",        base) ~ "2024-07-26",
    grepl("_CD_",              base) ~ "2024-07-11",
    grepl("lowerYR",           base) ~ "2024-07-16",
    grepl("upperYR",           base) ~ "2024-07-10",
    grepl("upperPR_CL_071024", base) ~ "2024-07-10",
    grepl("upperPR_CL_071624", base) ~ "2024-07-16",
    TRUE                             ~ NA_character_
  )
}

# Read, date-stamp, and row-bind all cluster CSVs
ortho_width <- purrr::map_dfr(ortho_files, function(f) {
  readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::mutate(SWOT_date = assign_swot_date(f))
}) %>%
  dplyr::mutate(SWOT_date = as.Date(SWOT_date))

# Join SWORD geometry to get node_len for each node
ortho_width <- sword_sf %>%
  right_join(ortho_width, by = "node_id")

# Filter to viable reaches and compute node-level ortho width
ortho_width_filtered <- ortho_width %>%
  filter(reach_id %in% viable_reaches$reach_id) %>%
  mutate(ortho_width_m = water_area_m2 / node_len)

# Aggregate to reach level: mean ortho width across all nodes per reach and date
ortho_reach_width <- ortho_width_filtered %>%
  dplyr::group_by(reach_id, SWOT_date) %>%
  dplyr::summarise(
    ortho_width_m = mean(ortho_width_m, na.rm = TRUE),
    n_nodes       = dplyr::n(),
    .groups       = "drop"
  )

# Add river name labels to ortho reach dataframe
ortho_reach_width <- ortho_reach_width %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      # Specific reach_id checks take priority to avoid overwriting SJ/BL labels
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
  ) %>%
  mutate(source = "PGD0")  ## CHANGE TO CORRECT VERSION!


# =============================================================================
# 2. Read and filter SWOT reach data
# =============================================================================

# --- RiverSP Version D (SWORD v17b) ------------------------------------------
SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v17b/RiverSP_domain_reach_timeseries_PGD0_v17b.csv')

# --- RiverSP Version C (SWORD v16) -------------------------------------------
# SWOT_reach_df <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v16/RiverSP_domain_reach_timeseries_v16.csv')

# Remove duplicates and sentinel fill values (time = -999…, wse = -1e12)
SWOT_reach_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE) %>%
  filter(time > 0) %>%
  filter(wse > 0)

# Quality filter: reach_q < 2, valid cross-track swath, < 80% dark water,
# partial_f == 0 means >= 50% of nodes present
SWOT_reach_df_filtered <- SWOT_reach_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac < 0.8) %>%
  filter(partial_f == 0)

# Convert TAI time to UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds

SWOT_reach_df_filtered$time_utc <- tai_epoch + SWOT_reach_df_filtered$time_tai - tai_utc_offset


# =============================================================================
# 3. Match ortho and SWOT observations in time and space
# =============================================================================

# Match by same calendar date and reach_id
time_matched_SWOT_ortho <- ortho_reach_width %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_reach_df_filtered %>%
        filter(
          as.Date(time_utc) == SWOT_date,
          reach_id == .env$reach_id
        ) %>%
        dplyr::select(-reach_id)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Filtering explorations:
# time_matched_SWOT_ortho <- time_matched_SWOT_ortho %>%
#   filter(reach_q < 2)


# =============================================================================
# 4. Absolute residuals and percent difference
# =============================================================================

# Residuals: ortho width - SWOT width
time_matched_SWOT_ortho$residuals <- (
  time_matched_SWOT_ortho$ortho_width_m - time_matched_SWOT_ortho$width
)

# Percent difference using ortho as truth
time_matched_SWOT_ortho <- time_matched_SWOT_ortho %>%
  mutate(percent_diff = ((abs(ortho_width_m - width)) / ortho_width_m) * 100)

percentile_68_error   <- quantile(abs(time_matched_SWOT_ortho$residuals),    0.68, na.rm = TRUE)
percentile_50_error   <- quantile(abs(time_matched_SWOT_ortho$residuals),    0.50, na.rm = TRUE)
percentile_68_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.68, na.rm = TRUE)
percentile_50_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error (m):",  percentile_68_error))
print(paste("50th Percentile Error (m):",  percentile_50_error))
print(paste("68th Percentile Error (%):", percentile_68_percent))
print(paste("50th Percentile Error (%):", percentile_50_percent))

# Pearson correlation
cor_test <- cor.test(time_matched_SWOT_ortho$width, time_matched_SWOT_ortho$ortho_width_m)
r_value  <- cor_test$estimate
p_value  <- cor_test$p.value


# =============================================================================
# 5. Data visualization
# =============================================================================

# Set river factor levels for consistent plot ordering
time_matched_SWOT_ortho$river <- factor(
  time_matched_SWOT_ortho$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

color_palette <- c("#F2C14E", "#F4845F", "#DA627D", "#9A348E")  # full quality filters
# color_palette <- c("#F2C14E", "#8EAD7A", "#F4845F", "#DA627D", "#9A348E") # no dark water filter
# color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E") # no filters

# Scatter: SWOT vs ortho reach width, colored by river
ggplot(time_matched_SWOT_ortho, aes(x = ortho_width_m, y = width, color = factor(river))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  annotate("text",
    x     = min(time_matched_SWOT_ortho$ortho_width_m, na.rm = TRUE),
    y     = max(time_matched_SWOT_ortho$width, na.rm = TRUE),
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", round(signif(p_value, 3), 4),
      "\nn = ", nrow(time_matched_SWOT_ortho)
    ),
    hjust = 0, vjust = 1, size = 8) +
  # theme(legend.position = "none") +  # comment out to show legend
  labs(color = "River") +
  theme_minimal(base_size = 30)

# CDF: absolute width difference
ggplot(time_matched_SWOT_ortho, aes(x = abs(width - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "SWOT - Ortho Width (m)",
    y     = "Cumulative Probability",
    title = "CDF of absolute width residuals"
  ) +
  annotate("text", x = 350, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_error, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 350, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_error, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)

# CDF: percent difference
ggplot(time_matched_SWOT_ortho, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = "Width Percent Difference (%)",
    y     = "Cumulative Probability",
    title = "CDF of %diff width residuals"
  ) +
  annotate("text", x = 70, y = 0.71,
    label = paste("68% abs diff:", round(percentile_68_percent, 4)),
    color = "#222222", size = 6) +
  annotate("text", x = 70, y = 0.53,
    label = paste("50% abs diff:", round(percentile_50_percent, 4)),
    color = "#222222", size = 6) +
  theme_minimal(base_size = 20)


# =============================================================================
# 6. Summary table by river
# =============================================================================

table_absolute_reach_width <- time_matched_SWOT_ortho %>%
  st_drop_geometry() %>%
  group_by(river) %>%
  summarise(
    med_swot_wd         = round(median(width, na.rm = TRUE), 1),
    med_ortho_wd        = round(median(ortho_width_m, na.rm = TRUE), 1),
    error_abs_68ile     = round(quantile(abs(residuals),    0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals),    0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals), na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff,      0.50, na.rm = TRUE), 2),
    n                   = sum(!is.na(residuals)),
    n_unique_reaches    = n_distinct(reach_id)
  )

# Pearson r and p-value by river
cor_table <- time_matched_SWOT_ortho %>%
  st_drop_geometry() %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  )

# Join correlation columns into summary table
table_absolute_reach_width <- table_absolute_reach_width %>%
  left_join(cor_table, by = "river")
