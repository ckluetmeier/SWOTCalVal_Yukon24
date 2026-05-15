# =============================================================================
# SWOT Width Validation
# -----------------------------------------------------------------------------
# Compares SWOT river width against in situ orthomosaic imagery
# width measurements at the node level.
# SWOT processing versions:
#   - Version C / PIC0 (SWORD v16, RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# SWOT and in situ measurements are matched in time/space in scripts 3.1 - 3.4
# and all data are harmonized to SWORD v17b node IDs before analysis.
# Outliers with |residuals| >= 1500 m are excluded from both data frames.
#
# Contains:
#   - Tables: 7, 8, S5
#   - Figures: 7a, b
# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)


# =============================================================================
# 1. Read in node-level width data
# =============================================================================

# Version C (SWORD v16 / RiverSP PIC0)
node_SWOT_ortho_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/node_width_SWOT_Ortho.csv") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(abs(residuals) < 1500, dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
node_SWOT_ortho_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/node_width_SWOT_Ortho.csv") %>%
  mutate(source = "PGD0") %>%
  filter(abs(residuals) < 1500, dark_frac < 0.5)


# =============================================================================
# 2. Harmonize to SWORD v17b node IDs
# =============================================================================

# Translator: maps v16 node IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv"
)

# Apply translation to Version C data
node_SWOT_ortho_vC <- node_SWOT_ortho_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")
  ) %>%
  rename(node_id = v17_node_id)

# Merge both versions
node_SWOT_ortho <- bind_rows(node_SWOT_ortho_vC, node_SWOT_ortho_vD)

# Compute version inclusion flag per node:
#   -1 = observed only in Version C
#    0 = observed in both versions
#    1 = observed only in Version D
all_nodes <- node_SWOT_ortho %>%
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
    )
  ) %>%
  select(node_id, version_inclusion)

# Join version inclusion to full dataset
node_SWOT_ortho <- node_SWOT_ortho %>%
  left_join(all_nodes, by = "node_id")

# river factor levels and color palette for inter-river plots
river_levels  <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
river_labels  <- c("Coleen", "Sheenjek", "Chandalar", "Porcupine",
                   "Single-channel Yukon", "Braided Yukon")
color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")


# =============================================================================
# NODE TABLES — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 3a. Absolute node width: by SWOT version (C vs D)
# -----------------------------------------------------------------------------

table_absolute_node_width <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    error_abs_68ile     = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals), na.rm = TRUE), 1),
    bias                = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(percent_diff, 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    n                   = sum(!is.na(residuals)),
    n_unique_nodes      = n_distinct(node_id)
  )

# Pearson correlation by source
cor_table <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(width, ortho_width_m)$p.value,
      error = function(e) NA_real_),
    .groups = "drop")

table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "source")


# -----------------------------------------------------------------------------
# 3b. Absolute node width: by river (D only)
# -----------------------------------------------------------------------------

table_absolute_node_width <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    med_swot_wd         = round(median(width), 1),
    med_ortho_wd        = round(median(ortho_width_m), 1),
    error_abs_68ile     = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals), na.rm = TRUE), 1),
    bias                = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    # min_error_percentdiff_68ile = round(min(abs(percent_diff)), 1),
    # max_error_percentdiff_50ile = round(max(abs(percent_diff)), 1),
    n                   = sum(!is.na(residuals_nobias)),
    n_unique_nodes      = n_distinct(node_id)
  )

# Pearson correlation by river (vD only)
cor_table <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(width, ortho_width_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "river")


# -----------------------------------------------------------------------------
# 3c. Absolute node width: by version inclusion (nodes unique to C or D)
# -----------------------------------------------------------------------------

table_absolute_node_width <- node_SWOT_ortho %>%
  group_by(version_inclusion) %>%
  summarise(
    error_abs_68ile     = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias                = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    n                   = sum(!is.na(residuals_nobias)),
    n_unique_nodes      = n_distinct(node_id)
  )

# Pearson correlation by version inclusion
cor_table <- node_SWOT_ortho %>%
  group_by(version_inclusion) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(width, ortho_width_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Join; drop version_inclusion == 0 (nodes present in both versions)
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>%
  mutate(version_inclusion = factor(version_inclusion, labels = c("vC0", "vD0")))


# -----------------------------------------------------------------------------
# 3d. Absolute node width: matched subset (same nodes present in both versions)
# -----------------------------------------------------------------------------

# Retain only nodes present in both vC and vD; drop rows where the paired
# source has a missing residuals_nobias value, then confirm both sources remain.
same_version_subset_node_SWOT_insitu <- node_SWOT_ortho %>%
  filter(version_inclusion == 0) %>%
  group_by(node_id, cycle_id, pass_id) %>%
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

# Summary stats for matched node subset
table_absolute_node_width <- same_version_subset_node_SWOT_insitu %>%
  filter(version_inclusion == 0) %>%
  group_by(source) %>%
  summarise(
    error_abs_68ile     = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias                = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    n                   = sum(!is.na(residuals)),
    n_unique_nodes      = n_distinct(node_id)
  )

# Pearson correlation for matched subset (version_inclusion == 0 only)
cor_table <- node_SWOT_ortho %>%
  filter(version_inclusion == 0) %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(
      cor.test(width, ortho_width_m)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "source")


# =============================================================================
# NODE PLOTS — SUMMARY STATISTICS
# =============================================================================


# -----------------------------------------------------------------------------
# 4a. CDF: absolute node width difference by SWOT version (C vs D)
# -----------------------------------------------------------------------------

# Observation counts for annotation
n_relative_df <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    n_unique_nodes = n_distinct(node_id),
    n              = sum(!is.na(residuals)),
    .groups = "drop"
  )

ggplot(node_SWOT_ortho, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "Width (m)"),
    y     = "Cumulative Probability",
    title = "By absolute difference"
  ) +
  annotate("text", x = 215, y = 0.72, hjust = 0,
    label = paste("68% C:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$residuals), 0.68, na.rm = TRUE), 1),
      "m, D:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$residuals), 0.68, na.rm = TRUE), 1),
      "m"),
    color = "#222222", size = 7) +
  annotate("text", x = 215, y = 0.54, hjust = 0,
    label = paste("50% C:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$residuals), 0.5, na.rm = TRUE), 1),
      "m, D:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$residuals), 0.5, na.rm = TRUE), 1),
      "m"),
    color = "#222222", size = 7) +
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
  coord_cartesian(xlim = c(0, 500))
# export dimensions: width 7.17 in, height 6.35 in


# -----------------------------------------------------------------------------
# 4b. CDF: percent width difference by SWOT version (C vs D)
# -----------------------------------------------------------------------------

ggplot(node_SWOT_ortho, aes(x = abs(percent_diff), color = source, linetype = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(
    x     = expression("SWOT -" ~ italic("in situ") ~ "Width (% difference)"),
    y     = "Cumulative Probability",
    title = "By percent difference"
  ) +
  annotate("text", x = 65, y = 0.72, hjust = 0,
    label = paste("68% C:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$percent_diff), 0.68, na.rm = TRUE), 1),
      "%, D:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$percent_diff), 0.68, na.rm = TRUE), 1),
      "%"),
    color = "#222222", size = 7) +
  annotate("text", x = 65, y = 0.54, hjust = 0,
    label = paste("50% C:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$percent_diff), 0.5, na.rm = TRUE), 1),
      "%, D:",
      round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$percent_diff), 0.5, na.rm = TRUE), 1),
      "%"),
    color = "#222222", size = 7) +
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
# export dimensions: width 7.17 in, height 6.35 in


# =============================================================================
# INTER-RIVER COMPARISON PLOTS (D only)
# =============================================================================

# Set factor order for river and ensure counts match the same ordering
node_SWOT_ortho_vD$river <- factor(node_SWOT_ortho_vD$river, levels = river_levels)

# Per-river sample sizes for plot annotation
counts <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()
counts$river <- factor(counts$river, levels = river_levels)


# -----------------------------------------------------------------------------
# 5a. Violin: percent width difference by river (D)
# -----------------------------------------------------------------------------

ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(percent_diff), fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
    aes(x = river, y = -2, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab(expression(atop("SWOT -" ~ italic("in situ") ~ "Width", "(% difference)"))) +
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 20, hjust = 0.9),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)  # extra bottom margin for n= labels
  ) +
  coord_cartesian(ylim = c(-5, 150))
# export dimensions: width 9.44 in, height 6.01 in


# -----------------------------------------------------------------------------
# 5b. Violin: absolute width difference by river (D)
# -----------------------------------------------------------------------------

ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(residuals), fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = counts,
    aes(x = river, y = -2, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 6) +
  xlab("River") +
  ylab(expression(atop("SWOT -" ~ italic("in situ") ~ "Width (m)"))) +
  scale_fill_manual(values = color_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 20, hjust = 0.9),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)
  ) +
  coord_cartesian(ylim = c(-5, 750))
# export dimensions: width 9.44 in, height 6.01 in


# =============================================================================
# SCATTER PLOTS AND CORRELATION (D)
# =============================================================================

# Re-order river factor with SJ last so Sheenjek nodes are plotted on top
node_SWOT_ortho_vD$river <- factor(
  node_SWOT_ortho_vD$river,
  levels = c("CL", "CD", "PR", "upperYR", "lowerYR", "SJ")
)
# Adjusted palette to match the new level order (SJ moved to end)
color_palette_scatter <- c("#F2C14E", "#3B6064", "#F4845F", "#DA627D", "#9A348E", "#8EAD7A")

# Overall Pearson correlation (SWOT vs ortho width)
cor_test <- cor.test(node_SWOT_ortho_vD$width, node_SWOT_ortho_vD$ortho_width_m)
r_value  <- cor_test$estimate   # Pearson r
p_value  <- cor_test$p.value    # p < 0.001 is highly statistically significant


# -----------------------------------------------------------------------------
# 6a. Scatter: SWOT width vs ortho width (D, colored by river)
# -----------------------------------------------------------------------------

ggplot() +
  # Plot non-SJ rivers first, then SJ on top so Sheenjek points are visible
  geom_point(
    data = subset(node_SWOT_ortho_vD, river != "SJ"),
    aes(x = ortho_width_m, y = width, color = river), size = 2.5
  ) +
  geom_point(
    data = subset(node_SWOT_ortho_vD, river == "SJ"),
    aes(x = ortho_width_m, y = width, color = river), size = 2.5
  ) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_color_manual(values = color_palette_scatter) +
  xlab(expression(italic("In situ") ~ "Width (m)")) +
  ylab("SWOT Width (m)") +
  ylim(0, 2700) +
  xlim(0, 2700) +
  annotate("text",
    x = min(node_SWOT_ortho_vD$ortho_width_m, na.rm = TRUE),
    y = 2700,
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", round(signif(p_value, 3), 4),
      "\nn = ", nrow(node_SWOT_ortho_vD)
    ),
    hjust = 0, vjust = 1, size = 8) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none")
# Suggested export dimensions: width 6.66 in, height 6.01 in






# =============================================================================
# EXPLORATORY SECTION
# =============================================================================


# =============================================================================
# CROSS-TRACK DISTANCE BIAS ANALYSIS (D)
# =============================================================================


# -----------------------------------------------------------------------------
# 7a. Prepare filtered data: trim residuals to 2.5–97.5% per river
# -----------------------------------------------------------------------------

node_SWOT_ortho_vD_filt <- node_SWOT_ortho_vD %>%
  filter(!is.na(residuals), !is.na(xtrk_dist)) %>%
  group_by(river) %>%
  filter(
    residuals >= quantile(residuals, 0.025, na.rm = TRUE),
    residuals <= quantile(residuals, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Correlation stats (residuals vs |xtrk_dist|) by river after filtering
cor_stats <- node_SWOT_ortho_vD_filt %>%
  group_by(river) %>%
  summarise(
    n       = n(),
    r_value = ifelse(n >= 3, cor(residuals, abs(xtrk_dist), method = "pearson"), NA_real_),
    p_value = ifelse(n >= 3, cor.test(residuals, abs(xtrk_dist))$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", round(r_value, 4), "\np value = ", signif(p_value, 3), "\nn = ", n)
  )


# -----------------------------------------------------------------------------
# 7b. Faceted scatter: width residuals vs |cross-track distance|, by river
# -----------------------------------------------------------------------------

ggplot(node_SWOT_ortho_vD_filt, aes(x = abs(xtrk_dist / 1000), y = residuals, color = river)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~ river, scales = "free") +
  scale_color_manual(values = color_palette_scatter) +
  xlab("|Cross-track distance| (km)") +
  ylab("SWOT Width residuals (m)") +
  geom_text(data = cor_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 6) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none")


# -----------------------------------------------------------------------------
# 7c. Overall scatter: width residuals vs |cross-track distance|, all rivers
# -----------------------------------------------------------------------------

# Overall correlation stats across all rivers
cor_stats_all <- node_SWOT_ortho_vD_filt %>%
  summarise(
    n       = n(),
    r_value = ifelse(n >= 3, cor(residuals, abs(xtrk_dist), method = "pearson"), NA_real_),
    p_value = ifelse(n >= 3, cor.test(residuals, abs(xtrk_dist))$p.value, NA_real_)
  ) %>%
  mutate(
    label = paste0("All rivers\nr = ", round(r_value, 4), "\np = ", signif(p_value, 3), "\nn = ", n)
  )

ggplot(node_SWOT_ortho_vD_filt, aes(x = abs(xtrk_dist / 1000), y = residuals, color = river)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 1, color = "black") +
  scale_color_manual(values = color_palette_scatter) +
  xlab("|Cross-track distance| (km)") +
  ylab("SWOT Width residuals (m)") +
  geom_text(data = cor_stats_all,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 6) +
  theme_minimal(base_size = 25)


# =============================================================================
# SPLIT VIOLIN COMPARISONS: SWOT vs ORTHO WIDTH BY RIVER (vD)
# =============================================================================

# NOTE: introdataviz is a GitHub-only package. Install once with:
# devtools::install_github("psyteachr/introdataviz")

# Expand to long format for side-by-side width type comparison
long_df <- node_SWOT_ortho_vD %>%
  select(river, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to  = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
    ortho_width_m = "Ortho",
    width         = "SWOT"
  ))

# Set river factor order for plotting
long_df$river <- factor(long_df$river, levels = river_levels)


# -----------------------------------------------------------------------------
# 8a. Split violin: all rivers combined
# -----------------------------------------------------------------------------

ggplot(long_df, aes(x = river, y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = FALSE, color = NA) +
  geom_boxplot(width = 0.2, alpha = 0.8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange",
               show.legend = FALSE, position = position_dodge(0.175)) +
  scale_x_discrete(labels = river_levels) +
  scale_fill_manual(values = c("lightblue", "darkblue")) +
  ylab("Width (m)") +
  coord_cartesian(ylim = c(0, 3500)) +
  theme_minimal(base_size = 30) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 0.9)
  )


# -----------------------------------------------------------------------------
# 8b. Split violin: one plot per river
# -----------------------------------------------------------------------------

plots <- lapply(river_levels, function(riv) {
  df_riv <- long_df %>% filter(river == riv)

  ggplot(df_riv, aes(x = river, y = width_value, fill = width_type)) +
    introdataviz::geom_split_violin(alpha = 0.5, trim = FALSE, color = NA) +
    geom_boxplot(width = 0.2, alpha = 0.8, fatten = NULL, show.legend = FALSE) +
    stat_summary(fun.data = "mean_se", geom = "pointrange",
                 position = position_dodge(0.175), show.legend = FALSE) +
    scale_fill_manual(values = c("lightblue", "darkblue")) +
    coord_cartesian(ylim = c(0, 300)) +
    labs(title = riv, y = "Width (m)", x = NULL) +
    theme_minimal(base_size = 30) +
    theme(legend.position = "none")
})

names(plots) <- river_levels

# look at the plots like this:
plots$CL




# -----------------------------------------------------------------------------
# 9. Per-river violin: percent width difference with 68th-percentile annotation
# -----------------------------------------------------------------------------

# Per-river n and q68 for annotation
stats <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    n   = n(),
    q68 = quantile(abs(percent_diff), 0.68, na.rm = TRUE)
  ) %>%
  ungroup()

stats$river <- factor(stats$river, levels = river_levels)

ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(percent_diff), fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  geom_text(data = stats,
    aes(x = river, y = -3, label = paste0("q=", round(q68, 1), "%\n", "n=", n)),
    inherit.aes = FALSE, vjust = 1, size = 5, lineheight = 0.9) +
  coord_cartesian(ylim = c(-8, 150)) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none", panel.clip = "off")


# -----------------------------------------------------------------------------
# 10. Export: save node-level means for data visualisation
# -----------------------------------------------------------------------------

node_means <- node_SWOT_ortho_vD %>%
  group_by(node_id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

write.csv(
  node_means,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/node_avg_width_SWOT_Ortho.csv",
  row.names = FALSE
)



# -----------------------------------------------------------------------------
# 11. Scratch: read SWOT tile AOI pass list
# -----------------------------------------------------------------------------

AOI <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/YF24_AOI_SWOT_tiles.csv"
)

length(unique(AOI$pass))
