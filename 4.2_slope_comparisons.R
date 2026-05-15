# =============================================================================
# SWOT Slope Validation
# -----------------------------------------------------------------------------
# Compares SWOT water surface slope against in situ measurements (PT and GNSS)
# at the reach level.
# SWOT processing versions:
#   - Version C / PIC0 (SWORD v16, RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# SWOT and in situ measurements are matched in time/space in scripts 1.3 & 2.2
# and all data are harmonized to SWORD v17b node/reach IDs before analysis.
# Units: slope residuals are stored in m/m; multiplied by 100,000 to get cm/km.
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)


# =============================================================================
# 1. Read in reach-level slope data
# =============================================================================

# --- PT ----------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
reach_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_PT.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
reach_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# --- GNSS ---------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
reach_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
reach_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)


# =============================================================================
# 2. Harmonize to SWORD v17b reach IDs
# =============================================================================

# Translator: maps v16 reach IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

# Apply translation to Version C PT data
reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% dplyr::select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")) %>%
  rename(reach_id = v17_reach_id)

# Apply translation to Version C GNSS data
reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% dplyr::select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")) %>%
  rename(reach_id = v17_reach_id)

# Merge all four data frames; unify slope and time column names
reach_SWOT_full_insitu <- bind_rows(
  reach_SWOT_PT_vC,
  reach_SWOT_PT_vD,
  reach_SWOT_GNSS_vC,
  reach_SWOT_GNSS_vD
) %>%
  mutate(
    insitu_slope_m_m        = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs),
    insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias),
    insitu_time_utc         = coalesce(pt_time_UTC, wse_drift_midpoint_UTC)
  )

# Compute version inclusion flag per reach:
#   -1 = observed only in Version C
#    0 = observed in both versions
#    1 = observed only in Version D
all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source, insitu_time_utc) %>%
  group_by(reach_id) %>%
  summarise(
    has_PIC0   = any(source == "PIC0"),
    has_PGD0 = any(source == "PGD0"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_PIC0 & has_PGD0  ~  0L,
      has_PIC0 & !has_PGD0 ~ -1L,
      !has_PIC0 & has_PGD0 ~  1L,
      TRUE                         ~ NA_integer_
    )
  ) %>%
  dplyr::select(reach_id, version_inclusion)

# Join version inclusion to full dataset
reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  left_join(all_reaches, by = "reach_id")

# Shared river factor levels and color palette used across all inter-river plots
river_levels  <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
river_labels  <- c("Coleen", "Sheenjek", "Chandalar", "Porcupine",
                   "Single-channel Yukon", "Braided Yukon")
color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")


# =============================================================================
# REACH TABLES — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 3a. Relative reach slope: by SWOT version (C vs D)
# -----------------------------------------------------------------------------

table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals_nobias)),    # count of non-NA residuals
    n_unique_reaches = n_distinct(reach_id)                  # count of unique reaches
  )

# Pearson correlation by source
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Join correlation columns to summary table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "source")

# Relabel and reorder source factor for bar chart
table_relative_reach_slope <- table_relative_reach_slope %>%
  mutate(
    source = factor(source,
      levels = c("PGD0", "PIC0"),
      labels = c("D",  "C")
    )
  )

# Bar chart: observation count by SWOT version
ggplot(table_relative_reach_slope, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(11, 157)) +
  scale_fill_manual(values = c("C" = "#E69F00", "D" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(
    axis.title.x    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    legend.position = "none"
  )
# export dimensions: width 3.16 in, height 6.54 in


# -----------------------------------------------------------------------------
# 3b. Relative reach slope: by in situ type (PT vs GNSS) and version
# -----------------------------------------------------------------------------

table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by source and in situ type
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = c("source", "insitu_type"))


# -----------------------------------------------------------------------------
# 3c. Relative reach slope: by river (D PT or GNSS only)
# -----------------------------------------------------------------------------

table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0", insitu_type == "GNSS") %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",  # merge lower & upper Porcupine
    TRUE ~ river
  )) %>%
  group_by(river) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by river (all sources/types)
cor_table <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0", insitu_type == "GNSS") %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",  # merge lower & upper Porcupine
    TRUE ~ river
  )) %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "river")


# -----------------------------------------------------------------------------
# 3d. Absolute reach slope: by version
# -----------------------------------------------------------------------------

table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation using absolute slope
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 3e. Absolute reach slope: by in situ type and version
# -----------------------------------------------------------------------------

table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by in situ type
# ERROR HERE IF I KEEP THE CORR VALS
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "insitu_type")


# -----------------------------------------------------------------------------
# 3f. Relative reach slope: by version inclusion (reaches unique to vC or vD)
# -----------------------------------------------------------------------------

table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by version inclusion
cor_table <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    n       = sum(complete.cases(slope_abs, insitu_slope_nobias_m_m)),
    r_value = if (n > 1) {
      round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4)
    } else {
      NA_real_},
    p_value = if (n > 1) {
      cor.test(slope_abs, insitu_slope_nobias_m_m, method = "pearson")$p.value
    } else {
      NA_real_},
    .groups = "drop"
  ) %>%
  dplyr::select(-n)

# Join; drop version_inclusion == 0 (reaches present in both versions)
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>%
  mutate(version_inclusion = factor(version_inclusion, labels = c("C", "D")))


# -----------------------------------------------------------------------------
# 3g. Relative reach slope: matched subset (same reaches present in both versions)
# -----------------------------------------------------------------------------

# Retain only reaches present in both vC and vD; drop rows where the paired
# source has a missing slope_residuals_nobias value, then confirm both sources remain.
same_version_subset_reach_SWOT_insitu <- reach_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%
  group_by(reach_id, insitu_time_utc, insitu_type) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  mutate(
    PIC0_resid_na   = any(source == "PIC0" & is.na(slope_residuals_nobias)),
    PGD0_resid_na = any(source == "PGD0" & is.na(slope_residuals_nobias))
  ) %>%
  filter(
    !(source == "PGD0" & PIC0_resid_na),
    !(source == "PIC0" & PGD0_resid_na)
  ) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  dplyr::select(-PIC0_resid_na, -PGD0_resid_na)

# Pearson correlation for matched reach subset
cor_table <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Summary stats for matched reach subset
table_relative_reach_slope <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 2),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 2),
    n              = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id),
    .groups = "drop"
  )

table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "source")


# =============================================================================
# REACH PLOTS — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 4a. CDF: relative reach slope by SWOT version (C vs D)
# -----------------------------------------------------------------------------

# Observation counts for annotation
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    n_unique_reaches = n_distinct(reach_id),
    n                = sum(!is.na(slope_residuals_nobias)),
    .groups = "drop"
  )

ggplot(reach_SWOT_full_insitu,
       aes(x = abs(slope_residuals_nobias) * 100000, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "Slope (cm/km)"),
    y     = "Cumulative Probability",
    title = "By SWOT version"
  ) +
  annotate("text", x = 3, y = 0.71, hjust = 0,
    label = paste("68% C:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
      "cm/km, D:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
      "cm/km"),
    color = "#222222", size = 8) +
  annotate("text", x = 3, y = 0.53, hjust = 0,
    label = paste("50% C:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$slope_residuals_nobias) * 100000, 0.5, na.rm = TRUE), 2),
      "cm/km, D:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$slope_residuals_nobias) * 100000, 0.5, na.rm = TRUE), 2),
      "cm/km"),
    color = "#222222", size = 8) +
  annotate("text", x = Inf, y = 0.1, hjust = 1, vjust = 0,
    label = paste0("Version C: ",
      n_relative_df[n_relative_df$source == "PIC0", ]$n_unique_reaches,
      " unique, ",
      n_relative_df[n_relative_df$source == "PIC0", ]$n, " total"),
    color = "#E69F00", size = 8) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0,
    label = paste0("Version D: ",
      n_relative_df[n_relative_df$source == "PGD0", ]$n_unique_reaches,
      " unique, ",
      n_relative_df[n_relative_df$source == "PGD0", ]$n, " total"),
    color = "#0072B2", size = 8) +
  theme_minimal(base_size = 22) +
  scale_color_manual(values = c("PIC0" = "#E69F00", "PGD0" = "#0072B2")) +
  scale_linetype_manual(values = c("PIC0" = "solid", "PGD0" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 15))
# export dimensions: width 7.17 in, height 6.35 in


# -----------------------------------------------------------------------------
# 4b. CDF: relative reach slope by in situ type (PT vs GNSS), version D
# -----------------------------------------------------------------------------

# Subset to Version D only
reach_SWOT_PGD0_insitu <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0")

# Observation counts for annotation
n_relative_df <- reach_SWOT_PGD0_insitu %>%
  filter(source == "PGD0") %>%
  group_by(insitu_type) %>%
  summarise(
    n_unique_reaches = n_distinct(reach_id[source == "PGD0"]),
    n                = sum(!is.na(slope_residuals_nobias)),
    .groups = "drop"
  )

ggplot(reach_SWOT_PGD0_insitu,
       aes(x = abs(slope_residuals_nobias) * 100000, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("|SWOT -" ~ italic("in situ") ~ "Slope| (cm/km)"),
    y     = "Cumulative Probability",
    title = expression("By" ~ italic("in situ") ~ "measurement type")
  ) +
  annotate("text", x = 3.1, y = 0.71, hjust = 0,
    label = paste("68% PT:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT",  ]$slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
      " GNSS:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 2),
      "cm/km"),
    color = "#222222", size = 7.5) +
  annotate("text", x = 3.1, y = 0.53, hjust = 0,
    label = paste("50% PT:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT",  ]$slope_residuals_nobias) * 100000, 0.5, na.rm = TRUE), 2),
      " GNSS:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias) * 100000, 0.5, na.rm = TRUE), 2),
      "cm/km"),
    color = "#222222", size = 7.5) +
  annotate("text", x = Inf, y = 0.1, hjust = 1, vjust = 0,
    label = paste0("PT: ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_reaches,
      " unique, ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
    color = "#009E73", size = 8) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0,
    label = paste0("GNSS: ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches,
      " unique, ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"),
    color = "#CC79A7", size = 8) +
  theme_minimal(base_size = 22) +
  scale_color_manual(values = c("PT" = "#009E73", "GNSS" = "#CC79A7")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 15))
# export dimensions: width 7.17 in, height 6.35 in


# =============================================================================
# DIAGNOSTIC / EXPLORATORY SECTION
# =============================================================================


# -----------------------------------------------------------------------------
# 5. Porcupine River reach-level diagnostics (D PT only)
# -----------------------------------------------------------------------------

# Subset to PR reaches in D PT for a little lookie
problems <- reach_SWOT_PT_vD %>%
  filter(river == "PR")

# Violin plot of slope residuals per reach
ggplot(problems, aes(x = factor(reach_id), y = slope_residuals_nobias * 100000)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25)

# Per-reach summary stats for problematic PR reaches
problems_table <- problems %>%
  group_by(reach_id) %>%
  summarise(
    error_68ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(slope_residuals_nobias) * 100000, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(slope_residuals_nobias) * 100000, na.rm = TRUE), 1),
    n              = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Other filters to try:
# %>% filter(slope_residuals_nobias * 100000 > 3)
# %>% filter(reach_id != '81270100061')

# Percentile summary
percentile_68_error <- quantile(abs(problems$slope_residuals_nobias), 0.68, na.rm = TRUE)
percentile_50_error <- quantile(abs(problems$slope_residuals_nobias), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error:", percentile_68_error * 100000))
print(paste("50th Percentile Error:", percentile_50_error * 100000))


# -----------------------------------------------------------------------------
# 6. Merge Porcupine sub-reaches in individual version data frames
# -----------------------------------------------------------------------------

reach_SWOT_PT_vD <- reach_SWOT_PT_vD %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  ))

reach_SWOT_GNSS_vD <- reach_SWOT_GNSS_vD %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  ))


# -----------------------------------------------------------------------------
# 7. Covariate scatter plots: slope error vs reach attributes (D PT and GNSS)
# -----------------------------------------------------------------------------

# layovr_val vs slope error
ggplot(reach_SWOT_PT_vD, aes(x = layovr_val, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("layovr_val") + ylab("Slope error (cm/km)") +
  labs(color = "River") +
  theme_minimal(base_size = 30)

ggplot(reach_SWOT_GNSS_vD, aes(x = layovr_val, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("layovr_val") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# abs(xtrk_dist) vs slope error
ggplot(reach_SWOT_PT_vD, aes(x = abs(xtrk_dist), y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("abs(xtrk_dist)") + ylab("Slope error (cm/km)") +
  labs(color = "River") +
  theme_minimal(base_size = 30)

ggplot(reach_SWOT_GNSS_vD, aes(x = abs(xtrk_dist), y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("abs(xtrk_dist)") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# dark_frac vs slope error
ggplot(reach_SWOT_PT_vD, aes(x = dark_frac, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") + ylab("Slope error (cm/km)") +
  labs(color = "River") +
  theme_minimal(base_size = 30)

ggplot(reach_SWOT_GNSS_vD, aes(x = dark_frac, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# In situ slope magnitude vs slope error
ggplot(reach_SWOT_PT_vD, aes(x = slope_m_m_abs, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("In situ slope (m/m)") + ylab("Slope error (cm/km)") +
  labs(color = "River") +
  theme_minimal(base_size = 30)

ggplot(reach_SWOT_GNSS_vD, aes(x = reach_drift_slope_m_m_abs_nobias, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("In situ slope (m/m)") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# Reach width vs slope error
ggplot(reach_SWOT_PT_vD, aes(x = width, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("Width (m)") + ylab("Slope error (cm/km)") +
  labs(color = "River") +
  theme_minimal(base_size = 30)

ggplot(reach_SWOT_GNSS_vD, aes(x = width, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("Width (m)") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)

# n_good_nod vs slope error (GNSS only)
ggplot(reach_SWOT_GNSS_vD, aes(x = n_good_nod, y = abs(slope_residuals_nobias) * 100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("n_good_nod") + ylab("Slope error (cm/km)") +
  ylim(0, 7) +
  labs(color = "River") +
  theme_minimal(base_size = 30)
