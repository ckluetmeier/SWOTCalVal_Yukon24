# =============================================================================
# SWOT WSE Validation
# -----------------------------------------------------------------------------
# Compares SWOT water surface elevation (WSE) against in situ
# measurements (PT and GNSS) at node and reach levels.
# SWOT processing versions:
#   - Version C / PIC0 (SWORD v16, RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# SWOT and in situ measurements are matched in time/space in scripts 1.2 - 2.2
# and all data are harmonized to SWORD v17b node/reach IDs before analysis.
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)


# =============================================================================
# NODE LEVEL
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Read in node-level data
# -----------------------------------------------------------------------------

# --- PT ----------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
node_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
node_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# --- GNSS ---------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
node_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
node_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 2. Harmonize to SWORD v17b node IDs
# -----------------------------------------------------------------------------

# Translator: maps v16 node IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv")

# Apply translation to Version C PT data
node_SWOT_PT_vC <- node_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# Apply translation to Version C GNSS data
node_SWOT_GNSS_vC <- node_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# Merge all four data frames; unify column names for WSE and time
node_SWOT_full_insitu <- bind_rows(
  node_SWOT_PT_vC,
  node_SWOT_PT_vD,
  node_SWOT_GNSS_vC,
  node_SWOT_GNSS_vD) %>%
  mutate(
    insitu_wse_m        = coalesce(pt_wse_m, mean_node_drift_wse_m),
    insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_node_drift_wse_no_bias_m),
    insitu_time_utc     = coalesce(pt_time_UTC, time_UTC))

# Compute version inclusion flag per node:
#   -1 = observed only in Version C
#    0 = observed in both versions
#    1 = observed only in Version D
all_nodes <- node_SWOT_full_insitu %>%
  distinct(node_id, source) %>%
  group_by(node_id) %>%
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
    )) %>%
  select(node_id, version_inclusion)

# Join version inclusion to full dataset
node_SWOT_full_insitu <- node_SWOT_full_insitu %>%
  left_join(all_nodes, by = "node_id")


# =============================================================================
# NODE TABLES — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 2a. Relative node WSE: by SWOT version (C vs D)
# -----------------------------------------------------------------------------

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    # RMSE          = round(sqrt(mean((residuals_nobias * 100)^2, na.rm = TRUE)), 1),
    bias          = round(median(bias, na.rm = TRUE) * 100, 1),
    n             = sum(!is.na(residuals_nobias)),    # count of non-NA residuals
    n_unique_nodes = n_distinct(node_id)              # count of unique nodes
  )

# Pearson correlation and p-value by source
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_),
    .groups = "drop")

# Join correlation columns to summary table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "source")

# Relabel and reorder source for bar chart
table_relative_node_WSE <- table_relative_node_WSE %>%
  mutate(
    source = factor(source,
      levels = c("PGD0", "PIC0"),
      labels = c("D",  "C")))

# Bar chart: observation count by SWOT version
ggplot(table_relative_node_WSE, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(2000, 7150)) +
  scale_fill_manual(values = c("C" = "#E69F00", "D" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(
    axis.title.x  = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.y   = element_blank(),
    legend.position = "none"
  )
# Export dimensions: width 3.16 in, height 6.54 in


# -----------------------------------------------------------------------------
# 2b. Relative node WSE: by version inclusion (nodes unique to C or D)
# -----------------------------------------------------------------------------

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias          = round(median(bias, na.rm = TRUE) * 100, 1),
    n             = sum(!is.na(residuals_nobias)),
    n_unique_nodes = n_distinct(node_id))

# Pearson correlation by version inclusion
cor_table <- node_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_),
    .groups = "drop")

# Join; drop version_inclusion == 0 (nodes present in both versions)
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>%
  mutate(version_inclusion = factor(version_inclusion, labels = c("PIC0", "PGD0")))


# -----------------------------------------------------------------------------
# 2c. Relative node WSE: matched subset (same nodes present in both versions)
# -----------------------------------------------------------------------------

# Retain only nodes present in both C and D; drop rows where the paired
# source has a missing residuals_nobias value, then confirm both sources remain.
same_version_subset_node_SWOT_insitu <- node_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%
  group_by(node_id, insitu_time_utc, insitu_type) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  mutate(
    PIC0_resid_na   = any(source == "PIC0" & is.na(residuals_nobias)),
    PGD0_resid_na = any(source == "PGD0" & is.na(residuals_nobias))
  ) %>%
  filter(
    !(source == "PGD0" & PIC0_resid_na),
    !(source == "PIC0" & PGD0_resid_na)
  ) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  select(-PIC0_resid_na, -PGD0_resid_na)

# Summary stats for matched subset
table_relative_node_WSE <- same_version_subset_node_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias          = round(median(bias, na.rm = TRUE) * 100, 1),
    n             = sum(!is.na(residuals_nobias)),
    n_unique_nodes = n_distinct(node_id)
  )

# Pearson correlation for matched subset
cor_table <- same_version_subset_node_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 2d. Relative node WSE: by in situ type (PT vs GNSS) and version
# -----------------------------------------------------------------------------

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    RMSE          = round(sqrt(mean((residuals_nobias * 100)^2, na.rm = TRUE)), 1),
    bias          = round(median(bias, na.rm = TRUE) * 100, 1),
    n             = sum(!is.na(residuals_nobias)),
    n_unique_nodes = n_distinct(node_id))

# Pearson correlation by in situ type and source
cor_table <- node_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = c("insitu_type", "source"))


# -----------------------------------------------------------------------------
# 2e. Relative node WSE: by river (D PT only)
# -----------------------------------------------------------------------------

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  filter(source == "PGD0", insitu_type == "PT") %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",  # merge lower & upper Porcupine
    TRUE ~ river
  )) %>%
  group_by(river, insitu_type) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias          = round(median(bias, na.rm = TRUE) * 100, 1),
    n             = sum(!is.na(residuals_nobias)),
    n_unique_nodes = n_distinct(node_id)
  )

# Pearson correlation by river
cor_table <- node_SWOT_full_insitu %>%
  filter(source == "PGD0", insitu_type == "PT") %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  )) %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "river")


# -----------------------------------------------------------------------------
# 2f. Absolute node WSE: by version
# -----------------------------------------------------------------------------

table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals) * 100, na.rm = TRUE), 1),
    n             = sum(!is.na(residuals)),
    n_unique_nodes = n_distinct(node_id)
  )

# Pearson correlation using absolute WSE
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 2g. Absolute node WSE: by in situ type and version
# -----------------------------------------------------------------------------

table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    error_68ile   = round(quantile(abs(residuals) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile   = round(quantile(abs(residuals) * 100, 0.50, na.rm = TRUE), 1),
    MAE           = round(mean(abs(residuals) * 100, na.rm = TRUE), 1),
    n             = sum(!is.na(residuals)),
    n_unique_nodes = n_distinct(node_id)
  )

# Pearson correlation by source and in situ type
cor_table <- node_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = c("source", "insitu_type"))


# =============================================================================
# NODE PLOTS — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 3a. CDF: relative node WSE by SWOT version (C vs D)
# -----------------------------------------------------------------------------

# Observation counts for annotation
n_relative_df <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    n_unique_nodes = n_distinct(node_id),
    n              = sum(!is.na(residuals_nobias)),
    .groups = "drop")

ggplot(node_SWOT_full_insitu,
       aes(x = abs(residuals_nobias) * 100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
    y     = "Cumulative Probability",
    title = "By SWOT version"
  ) +
  annotate("text", x = 40, y = 0.72, hjust = 0,
    label = paste("68% C:",
      round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm, D:",
      round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = 40, y = 0.54, hjust = 0,
    label = paste("50% C:",
      round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm, D:",
      round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = Inf, y = 0.1, hjust = 1, vjust = 0,
    label = paste0("Version C: ",
      n_relative_df[n_relative_df$source == "PIC0", ]$n_unique_nodes,
      " unique, ",
      n_relative_df[n_relative_df$source == "PIC0", ]$n, " total"),
    color = "#E69F00", size = 8) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0,
    label = paste0("Version D: ",
      n_relative_df[n_relative_df$source == "PGD0", ]$n_unique_nodes,
      " unique, ",
      n_relative_df[n_relative_df$source == "PGD0", ]$n, " total"),
    color = "#0072B2", size = 8) +
  theme_minimal(base_size = 22) +
  scale_color_manual(values = c("PIC0" = "#E69F00", "PGD0" = "#0072B2")) +
  scale_linetype_manual(values = c("PIC0" = "solid", "PGD0" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# Export dimensions: width 7.17 in, height 6.35 in


# -----------------------------------------------------------------------------
# 3b. CDF: relative node WSE by in situ type (PT vs GNSS), D only
# -----------------------------------------------------------------------------

# Subset to Version D only
node_SWOT_PGD0_insitu <- node_SWOT_full_insitu %>%
  filter(source == "PGD0")

# Observation counts for annotation
n_relative_df <- node_SWOT_PGD0_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    n_unique_nodes = n_distinct(node_id[source == "PGD0"]),
    n              = sum(!is.na(residuals_nobias)),
    .groups = "drop"
  )

ggplot(node_SWOT_PGD0_insitu,
       aes(x = abs(residuals_nobias) * 100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
    y     = "Cumulative Probability",
    title = expression("By" ~ italic("in situ") ~ "measurement type")
  ) +
  annotate("text", x = 31, y = 0.72, hjust = 0,
    label = paste("68% PT:",
      round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "PT",  ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm, GNSS:",
      round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = 31, y = 0.54, hjust = 0,
    label = paste("50% PT:",
      round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "PT",  ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm, GNSS:",
      round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = Inf, y = 0.1, hjust = 1, vjust = 0,
    label = paste0("PT: ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_nodes,
      " unique, ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
    color = "#009E73", size = 8) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0,
    label = paste0("GNSS: ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_nodes,
      " unique, ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"),
    color = "#CC79A7", size = 8) +
  theme_minimal(base_size = 22) +
  scale_color_manual(values = c("PT" = "#009E73", "GNSS" = "#CC79A7")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# Export dimensions: width 7.17 in, height 6.35 in


# =============================================================================
# REACH LEVEL
# =============================================================================


# -----------------------------------------------------------------------------
# 4. Read in reach-level data
# -----------------------------------------------------------------------------

# --- PT ----------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
reach_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv"
) %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PIC0") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
reach_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv"
) %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# --- GNSS ---------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
reach_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  mutate(source = "PIC0") %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
reach_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_GNSS.csv") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 5. Harmonize to SWORD v17b reach IDs
# -----------------------------------------------------------------------------

# Translator: maps v16 reach IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

# Apply translation to Version C PT data
reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

# Apply translation to Version C GNSS data
reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

# Merge all four data frames; unify WSE and time column names
reach_SWOT_full_insitu <- bind_rows(
  reach_SWOT_PT_vC,
  reach_SWOT_PT_vD,
  reach_SWOT_GNSS_vC,
  reach_SWOT_GNSS_vD
) %>%
  mutate(
    insitu_wse_m           = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m),
    insitu_wse_nobias_m    = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m),
    insitu_time_utc        = coalesce(pt_time_UTC, wse_drift_midpoint_UTC)
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
  select(reach_id, version_inclusion)

reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  left_join(all_reaches, by = "reach_id")


# =============================================================================
# REACH TABLES — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 5a. Relative reach WSE: by SWOT version
# -----------------------------------------------------------------------------

table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias           = round(median(bias, na.rm = TRUE) * 100, 1),
    n              = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by source
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "source")

# Relabel and reorder source for bar chart
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  mutate(
    source = factor(source,
      levels = c("PGD0", "PIC0"),
      labels = c("vD0",  "vC0")
    )
  )

# Bar chart: observation count by SWOT version
ggplot(table_relative_reach_WSE, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(11, 185)) +
  scale_fill_manual(values = c("vC0" = "#E69F00", "vD0" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(
    axis.title.x    = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_blank(),
    legend.position = "none"
  )
# Export dimensions: width 3.16 in, height 6.54 in


# -----------------------------------------------------------------------------
# 5b. Relative reach WSE: by version inclusion (reaches unique to vC or vD)
# -----------------------------------------------------------------------------

table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias           = round(median(bias, na.rm = TRUE) * 100, 1),
    n              = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by version inclusion (guards against n <= 1)
cor_table <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    n       = sum(complete.cases(wse, insitu_wse_nobias_m)),
    r_value = if (n > 1) {
      round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4)
    } else {
      NA_real_
    },
    p_value = if (n > 1) {
      cor.test(wse, insitu_wse_nobias_m, method = "pearson")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Join; drop version_inclusion == 0 (reaches present in both versions)
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>%
  mutate(version_inclusion = factor(version_inclusion, labels = c("vC0", "vD0")))


# -----------------------------------------------------------------------------
# 5c. Relative reach WSE: matched subset (same reaches present in both versions)
# -----------------------------------------------------------------------------

# Retain only reaches present in both vC and vD; drop rows where the paired
# source has a missing residuals_nobias value, then confirm both sources remain.
same_version_subset_reach_SWOT_insitu <- reach_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%
  group_by(reach_id, insitu_time_utc, insitu_type) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  mutate(
    RiverSP_resid_na   = any(source == "PIC0" & is.na(residuals_nobias)),
    RiverTile_resid_na = any(source == "PGD0" & is.na(residuals_nobias))
  ) %>%
  filter(
    !(source == "PGD0" & RiverSP_resid_na),
    !(source == "PIC0" & RiverTile_resid_na)
  ) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  select(-RiverSP_resid_na, -RiverTile_resid_na)

# Pearson correlation for matched reach subset
cor_table <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Summary stats for matched reach subset
table_relative_reach_WSE <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias           = round(median(bias, na.rm = TRUE) * 100, 2),
    n              = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id),
    .groups = "drop"
  )

table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 5d. Relative reach WSE: by in situ type (PT vs GNSS) and version
# -----------------------------------------------------------------------------

table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias           = round(median(bias, na.rm = TRUE) * 100, 1),
    n              = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by in situ type and source
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = c("insitu_type", "source"))


# -----------------------------------------------------------------------------
# 5e. Relative reach WSE: by river (vD PT only)
# -----------------------------------------------------------------------------

table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0", insitu_type == "PT") %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  )) %>%
  group_by(river) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals_nobias) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals_nobias) * 100, na.rm = TRUE), 1),
    bias           = round(median(bias, na.rm = TRUE) * 100, 1),
    n              = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by river (all sources/types)
cor_table <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_nobias_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "river")


# -----------------------------------------------------------------------------
# 5f. Absolute reach WSE: by version
# -----------------------------------------------------------------------------

table_absolute_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals) * 100, na.rm = TRUE), 1),
    n              = sum(!is.na(residuals)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation using raw (biased) WSE
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 5g. Absolute reach WSE: by in situ type and version
# -----------------------------------------------------------------------------

table_absolute_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    error_68ile    = round(quantile(abs(residuals) * 100, 0.68, na.rm = TRUE), 1),
    error_50ile    = round(quantile(abs(residuals) * 100, 0.50, na.rm = TRUE), 1),
    MAE            = round(mean(abs(residuals) * 100, na.rm = TRUE), 1),
    n              = sum(!is.na(residuals)),
    n_unique_reaches = n_distinct(reach_id)
  )

# Pearson correlation by source and in situ type (raw WSE)
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(wse, insitu_wse_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
  left_join(cor_table, by = c("source", "insitu_type"))


# =============================================================================
# REACH PLOTS — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 6a. CDF: relative reach WSE by SWOT version (C vs D)
# -----------------------------------------------------------------------------

# Observation counts for annotation
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    n_unique_reaches = n_distinct(reach_id),
    n                = sum(!is.na(residuals_nobias)),
    .groups = "drop"
  )

ggplot(reach_SWOT_full_insitu,
       aes(x = abs(residuals_nobias) * 100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
    y     = "Cumulative Probability",
    title = "By SWOT version"
  ) +
  annotate("text", x = 40, y = 0.72, hjust = 0,
    label = paste("68% C:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm, D:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = 40, y = 0.54, hjust = 0,
    label = paste("50% C:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm, D:",
      round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm"),
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
  coord_cartesian(xlim = c(0, 150))
# export dimensions: width 7.17 in, height 6.35 in


# -----------------------------------------------------------------------------
# 6b. CDF: relative reach WSE by in situ type (PT vs GNSS), D only
# -----------------------------------------------------------------------------

# Subset to Version D only
reach_SWOT_PGD0_insitu <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0")

# Observation counts for annotation
n_relative_df <- reach_SWOT_PGD0_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    n_unique_reaches = n_distinct(reach_id[source == "PGD0"]),
    n                = sum(!is.na(residuals_nobias)),
    .groups = "drop"
  )

ggplot(reach_SWOT_PGD0_insitu,
       aes(x = abs(residuals_nobias) * 100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
    y     = "Cumulative Probability",
    title = expression("By" ~ italic("in situ") ~ "measurement type")
  ) +
  annotate("text", x = 34, y = 0.72, hjust = 0,
    label = paste("68% PT:",
      round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "PT",  ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm, GNSS:",
      round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias) * 100, 0.68, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = 27, y = 0.54, hjust = 0,
    label = paste("50% PT:",
      round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "PT",  ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm, GNSS:",
      round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias) * 100, 0.5, na.rm = TRUE), 1),
      "cm"),
    color = "#222222", size = 8) +
  annotate("text", x = Inf, y = 0.1, hjust = 1, vjust = 0,
    label = paste0("PT: ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_reaches,
      " unique reaches, ",
      n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
    color = "#009E73", size = 8) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0,
    label = paste0("GNSS: ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches,
      " unique reaches, ",
      n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"),
    color = "#CC79A7", size = 8) +
  theme_minimal(base_size = 22) +
  scale_color_manual(values = c("PT" = "#009E73", "GNSS" = "#CC79A7")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export dimensions: width 7.17 in, height 6.35 in


# =============================================================================
# INTER-RIVER COMPARISON PLOTS
# =============================================================================

# Shared river factor levels and color palette used across all inter-river plots
river_levels  <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
river_labels  <- c("Coleen", "Sheenjek", "Chandalar", "Porcupine",
                   "Single-channel Yukon", "Braided Yukon")
color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")


# -----------------------------------------------------------------------------
# 7a. Violin: absolute node WSE residuals by river (D PT)
# -----------------------------------------------------------------------------

# Merge Porcupine River sub-reaches and set factor order
node_SWOT_PT_vD <- node_SWOT_PT_vD %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  ))

node_SWOT_PT_vD$river <- factor(node_SWOT_PT_vD$river, levels = river_levels)

# Per-river sample sizes for plot annotation
counts <- node_SWOT_PT_vD %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()

ggplot(node_SWOT_PT_vD, aes(x = river, y = abs(residuals_nobias) * 100, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
    aes(x = river, y = -1, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab("| SWOT - PT WSE | (cm)") +
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 20, hjust = 0.9),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)  # extra bottom margin for n= labels
  ) +
  coord_cartesian(ylim = c(-5, 60))
# export dimensions: width 9.44 in, height 6.01 in


# -----------------------------------------------------------------------------
# 7b. Violin: node WSE bias by river (D PT)
# -----------------------------------------------------------------------------


ggplot(node_SWOT_PT_vD, aes(x = river, y = bias * 100, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
    aes(x = river, y = -1, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab("SWOT - PT WSE bias (cm)") +  # corrected from absolute-residual label
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 20, hjust = 0.9),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)
  ) +
  coord_cartesian(ylim = c(-25, 50))


# =============================================================================
# INTER-RIVER COMPARISON PLOTS FOR SLOPE VALIDATION
# =============================================================================


# -----------------------------------------------------------------------------
# 8. Read in reach-level slope data (D only)
# -----------------------------------------------------------------------------

reach_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv"
) %>%
  filter(dark_frac < 0.5)

reach_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_PT.csv"
) %>%
  filter(dark_frac < 0.5)

# Merge slope data frames; unify slope and time column names
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vD, reach_SWOT_GNSS_vD) %>%
  mutate(
    insitu_slope_m_m        = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs),
    insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias),
    insitu_time_utc         = coalesce(pt_time_UTC, wse_drift_midpoint_UTC)
  )

# Merge Porcupine River sub-reaches and set factor order
reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",
    TRUE ~ river
  ))

reach_SWOT_full_insitu$river <- factor(reach_SWOT_full_insitu$river, levels = river_levels)

# Per-river sample sizes for plot annotation
counts <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()


# -----------------------------------------------------------------------------
# 8a. Violin: absolute slope residuals by river
# -----------------------------------------------------------------------------

ggplot(reach_SWOT_full_insitu,
       aes(x = river, y = abs(slope_residuals_nobias) * 100000, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
    aes(x = river, y = -0.2, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab("| SWOT -" ~ italic("in situ") ~ "Slope | (cm/km)") +
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 20, hjust = 0.9),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)
  ) +
  coord_cartesian(ylim = c(-0.5, 12))
# export dimensions: width 9.44 in, height 6.01 in
