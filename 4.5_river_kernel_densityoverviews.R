# =============================================================================
# YR Domain — Width / Slope Kernel Densities & SWOT Sampling Overviews
# -----------------------------------------------------------------------------
# Produces per-river kernel density distributions of SWORD width and slope for
# the Yukon River (YR) domain, and characterizes SWOT sampling properties
# (functional repeat time) at the node and reach scale.
# Inputs:
#   - SWORD v17b nodes / reaches (North America HB81)
#   - YR domain definition (node / reach IDs included in the study)
#   - SWOT RiverSP v17b PGD0 timeseries (node and reach)
#   - Hydrocron-merged SWOT node timeseries (for dark fraction)
# Per-river groupings:
#   CL, SJ, CD, PR, upperYR (single-channel), lowerYR (braided), BL (excluded)
#
# Contains:
#   - Figures: width KDE by river, slope KDE by river, dark-fraction IQR plot
#   - Console output: SWOT functional repeat time stats (node & reach)
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)


# =============================================================================
# 1. Load SWORD geometries and YR domain
# =============================================================================

# AK / NA SWORD v17b node + reach shapefiles
sword_nodes   <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp")
sword_reaches <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

# Study domain — list of SWORD node/reach IDs included in the YR analysis
YR_domain <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv"
)


# =============================================================================
# WIDTH KERNEL DENSITY BY RIVER
# =============================================================================


# -----------------------------------------------------------------------------
# 2. Tag SWORD nodes with river and subset to YR domain
# -----------------------------------------------------------------------------

YR_sword_nodes <- sword_nodes %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
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
  filter(!is.na(river))

# Subset to YR domain nodes and drop geometry
YR_sword_nodes <- YR_sword_nodes %>%
  filter(node_id %in% YR_domain$Node_ID) %>%
  st_drop_geometry()


# -----------------------------------------------------------------------------
# 3. Trim per-river width distributions before KDE
# -----------------------------------------------------------------------------

# Default trim is the inter-quartile range; CL / lowerYR / upperYR use custom
# percentile windows to better isolate their typical channel width.
trimmed_df <- YR_sword_nodes %>%
  group_by(river) %>%
  mutate(
    q25 = quantile(width, 0.25, na.rm = TRUE),
    q75 = quantile(width, 0.75, na.rm = TRUE),
    q25 = if_else(river == "CL",      quantile(width, 0.50, na.rm = TRUE), q25),
    q75 = if_else(river == "CL",      quantile(width, 0.95, na.rm = TRUE), q75),
    q25 = if_else(river == "lowerYR", quantile(width, 0.40, na.rm = TRUE), q25),
    q75 = if_else(river == "lowerYR", quantile(width, 0.60, na.rm = TRUE), q75),
    q25 = if_else(river == "upperYR", quantile(width, 0.30, na.rm = TRUE), q25),
    q75 = if_else(river == "upperYR", quantile(width, 0.70, na.rm = TRUE), q75)
  ) %>%
  filter(width >= q25, width <= q75) %>%
  ungroup() %>%
  filter(river != "BL")

# Per-river medians for plot annotation
median_df <- trimmed_df %>%
  group_by(river) %>%
  summarize(median_width = median(width, na.rm = TRUE)) %>%
  ungroup()


# -----------------------------------------------------------------------------
# 4. Width KDE plot
# -----------------------------------------------------------------------------

# Palette order: CD, CL, PR, lowerYR, SJ, upperYR
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

ggplot(trimmed_df, aes(x = width, fill = river, color = river)) +
  geom_density(alpha = 0.5, size = 0.9, adjust = 1) +
  # Median lines (dashed) colored by river
  geom_vline(
    data = median_df,
    aes(xintercept = median_width, color = river),
    linetype = "dashed", size = 1
  ) +
  labs(x = "Width", y = "Density", fill = "River", color = "River") +
  theme_minimal(base_size = 20) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")
  ) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR")
  ) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    panel.grid      = element_blank()
  )
# export dimensions: width 11.05 in, height 2.22 in

# =============================================================================
# SLOPE KERNEL DENSITY BY RIVER
# =============================================================================


# -----------------------------------------------------------------------------
# 5. Tag SWORD reaches with river and subset to YR domain
# -----------------------------------------------------------------------------

YR_sword_reaches <- sword_reaches %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
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
  filter(!is.na(river))

# Subset to YR domain reaches and drop geometry
YR_reaches <- YR_sword_reaches %>%
  filter(reach_id %in% YR_domain$Reach_ID) %>%
  st_drop_geometry()


# -----------------------------------------------------------------------------
# 6. Filter reaches before plotting slope KDE
# -----------------------------------------------------------------------------

# Drop BL, trim outlying upstream segments on CD/upperYR, and apply per-river
# maximum slopes (PR < 0.6, lowerYR < 0.7, all others < 1.40)
trimmed_df <- YR_reaches %>%
  filter(
    river != "BL",
    !(reach_id > 81250800041 & river == "CD"),
    !(reach_id > 81270500231 & river == "upperYR"),
    (river == "PR"      & slope < 0.6) |
    (river == "lowerYR" & slope < 0.7) |
    (!(river %in% c("PR", "lowerYR")) & slope < 1.40)
  )

# Per-river IQR slope trim
# group_by(river) %>%
#   mutate(
#     q25 = if_else(river == "", quantile(width, 0.25, na.rm = TRUE), q25),
#     q75 = if_else(river == "", quantile(width, 0.75, na.rm = TRUE), q75)
#   ) %>%
#   filter(slope >= q25, slope <= q75) %>%
#   ungroup() %>%

# Per-river medians for plot annotation
median_df <- trimmed_df %>%
  group_by(river) %>%
  summarize(median_slope = median(slope, na.rm = TRUE)) %>%
  ungroup()


# -----------------------------------------------------------------------------
# 7. Slope KDE plot (slope shown as cm/km)
# -----------------------------------------------------------------------------

color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

ggplot(trimmed_df, aes(x = slope * 100, fill = river, color = river)) +
  geom_density(alpha = 0.2, size = 1.4, adjust = 1) +
  geom_vline(
    data = median_df,
    aes(xintercept = median_slope * 100, color = river),
    linetype = "dashed", size = 1
  ) +
  labs(x = "Slope", y = "Density", fill = "River", color = "River") +
  theme_minimal(base_size = 20) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")
  ) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "PR", "lowerYR", "SJ", "upperYR")
  ) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 0, hjust = 0.5),
    panel.grid      = element_blank()
  ) +
  xlim(-10, 160)
# export dimensions: width 11.05 in, height 2.22 in


# =============================================================================
# SWOT FUNCTIONAL REPEAT TIMES
# =============================================================================


# -----------------------------------------------------------------------------
# 8. Node-level functional repeat time
# -----------------------------------------------------------------------------

# SWOT node timeseries (RiverSP v17b / PGD0)
SWOT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv"
)

# Drop any rows that are exact duplicates on (node_id, time, wse)
SWOT_df_noduplicates <- SWOT_df %>%
  distinct(node_id, time, wse, .keep_all = TRUE)

# Restrict to good-quality SWOT obs within nominal cross-track distance
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(node_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(dark_frac <= 0.5)

# time_tai is seconds since 2000-01-01, offset 37 seconds from UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# Unique node count
n_unique_nodes <- n_distinct(SWOT_df_filtered$node_id)

# Per-node repeat time (median spacing between consecutive observations).
repeat_time_summary <- SWOT_df_filtered %>%
  arrange(node_id, time_utc) %>%
  group_by(node_id) %>%
  summarise(
    n_obs                  = n(),
    mean_repeat_time_hours = median(diff(time_utc), na.rm = TRUE) / 86400
  ) %>%
  ungroup()

# Median across all nodes
avg_repeat_time_hours <- median(repeat_time_summary$mean_repeat_time_hours, na.rm = TRUE)

cat("Unique node_id:", n_unique_nodes, "\n")
cat("Median repeat time across nodes (days):", round(avg_repeat_time_hours, 2), "\n")
cat("Average observations per node:", round(mean(repeat_time_summary$n_obs), 1), "\n")


# -----------------------------------------------------------------------------
# 9. Reach-level functional repeat time
# -----------------------------------------------------------------------------

SWOT_reach_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/reach/RiverSP_v17b/RiverSP_domain_reach_timeseries_PGD0_v17b.csv"
)

SWOT_df_noduplicates <- SWOT_reach_df %>%
  distinct(reach_id, time, wse, .keep_all = TRUE)

# Reach filter adds partial_f == 0 (full-reach observations only)
SWOT_df_filtered <- SWOT_df_noduplicates %>%
  filter(reach_q < 2) %>%
  filter(abs(xtrk_dist) >= 10000) %>%
  filter(abs(xtrk_dist) <= 60000) %>%
  filter(partial_f == 0) %>%
  filter(dark_frac <= 0.5)

tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

n_unique_reaches <- n_distinct(SWOT_df_filtered$reach_id)

repeat_time_summary <- SWOT_df_filtered %>%
  arrange(reach_id, time_utc) %>%
  group_by(reach_id) %>%
  summarise(
    n_obs                  = n(),
    mean_repeat_time_hours = median(diff(time_utc), na.rm = TRUE) / 86400
  ) %>%
  ungroup()

avg_repeat_time_hours <- median(repeat_time_summary$mean_repeat_time_hours, na.rm = TRUE)

cat("Unique reach_id:", n_unique_reaches, "\n")
cat("Median repeat time across nodes (days):", round(avg_repeat_time_hours, 2), "\n")
cat("Average observations per node:", round(mean(repeat_time_summary$n_obs), 1), "\n")

