# =============================================================================
# SWOT WSE Validation
# -----------------------------------------------------------------------------
# Compares SWOT water surface elevation (WSE) against in situ measurements
# (PT and GNSS) at node and reach levels.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# SWOT and in situ measurements are matched in time/space in scripts 1.2 - 2.2
# and harmonized to SWORD v17b node/reach IDs here.
#
# Produces: Tables 2, 3, S6, S7;  Figures 4a-d, 5, 6a
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

# Shared helpers. Run this script with 4_SWOT_validation/ as the wd
source("4.0_comparison_helpers.R")


# =============================================================================
# Configuration - edit these paths before running
# =============================================================================

# Root of the field-campaign and SWOT data products.
DATA_ROOT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats"

# Root of the SWORD distributions (contains SWORD_v16/ and SWORD_v17b/).
SWORD_ROOT <- "/Users/camryn/Desktop"

# Root of the WSE comparison dataframes written by scripts 1.2-2.2.
BASE <- file.path(DATA_ROOT, "CalVal_dataframes/wse")

# Directory holding the SWORD v16 <-> v17b ID translator csvs.
TRANSLATOR_DIR <- file.path(SWORD_ROOT, "SWORD_translation")

DARK_FRAC_MAX <- 0.5

# Reaches shorter than 9 km
SHORT_REACHES_V16  <- c(81260300061, 81270500131, 81270500141)
SHORT_REACHES_V17B <- c(81260300181, 81270500021, 81270500031)
APPLY_SHORT_REACH_EXCLUSION <- TRUE


# -----------------------------------------------------------------------------
# SWORD ID translators
# -----------------------------------------------------------------------------

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)
reach_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_ReachIDs_v17b_vs_v16.csv"),
                             show_col_types = FALSE)

node_lut  <- build_id_lut(node_translator,  "v16_node_id",  "v17_node_id",  "node translator")
reach_lut <- build_id_lut(reach_translator, "v16_reach_id", "v17_reach_id", "reach translator")


# =============================================================================
# Node level
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Read, filter, harmonize
# -----------------------------------------------------------------------------

read_node <- function(path, insitu, version) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D) %>%
    filter(dark_frac < DARK_FRAC_MAX) %>%
    harmonize_ids("node_id", node_lut, version,
                  label = paste("node", insitu, version))
}

node_PT_vC   <- read_node(file.path(BASE, "node/RiverSP_v16/node_SWOT_PT.csv"),          "PT",   "C")
node_PT_vD   <- read_node(file.path(BASE, "node/RiverSP_v17b/node_SWOT_PT.csv"),         "PT",   "D")
node_GNSS_vC <- read_node(file.path(BASE, "node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv"),  "GNSS", "C")
node_GNSS_vD <- read_node(file.path(BASE, "node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv"), "GNSS", "D")


# -----------------------------------------------------------------------------
# 2. Partition each in situ type
# -----------------------------------------------------------------------------

NODE_VALUE <- "residuals_nobias"  # the metric the partition is defined against

node_PT <- bind_rows(node_PT_vC, node_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonized", "pt_serial", "pt_time_UTC"),
                     value_col = NODE_VALUE)

node_GNSS <- bind_rows(node_GNSS_vC, node_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%  # strip version-specific path prefix
  partition_versions(key_cols  = c("id_harmonized", "cycle_id", "pass_id", "drift_file"),
                     value_col = NODE_VALUE)

node_all <- bind_rows(node_PT, node_GNSS) %>%
  mutate(
    insitu_wse_m        = coalesce(pt_wse_m,        mean_node_drift_wse_m),
    insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_node_drift_wse_no_bias_m),
    insitu_time_utc     = coalesce(pt_time_UTC,     time_UTC)  # kept for reference only
  )
# bind_rows() drops attributes, so restore the partition tag used by
# partition_table()
attr(node_all, "partition_value_col") <- NODE_VALUE


# =============================================================================
# Node tables
# =============================================================================

# --- Table 2, node rows: relative node WSE by version -------------------------
table2_node <- node_all %>%
  summarize_errors(NODE_VALUE, by = c("insitu_type", "source"),
                   scale = 100, digits = 1)
print(table2_node)

# --- Table S6, node section: the full partition -------------------------------
tableS6_node <- node_all %>%
  partition_table(NODE_VALUE, by = c("insitu_type", "source"),
                  scale = 100, digits = 1)
print(tableS6_node, n = Inf)

# How many unique nodes are unique to each version
tableS6_node_ids <- node_all %>%
  filter(!is.na(id_bucket)) %>%
  summarize(n_nodes = n_distinct(id_harmonized), .by = c(insitu_type, id_bucket))
print(tableS6_node_ids)

# --- Table S7, node rows: absolute node WSE -----------------------------------
tableS7_node <- node_all %>%
  summarize_errors("residuals", by = c("insitu_type", "source"),
                   scale = 100, digits = 1)
print(tableS7_node)

# --- Table 3, node rows: relative node WSE by river (D0, PT) ------------------
table3_node <- node_all %>%
  filter(source == VERSION_D, insitu_type == "PT") %>%
  merge_porcupine() %>%
  summarize_errors(NODE_VALUE, by = "river", scale = 100, digits = 1)
print(table3_node)

# --- Percent-change stats -----------------------------------------------------
node_change <- node_all %>% version_change(NODE_VALUE)
print(node_change)


# =============================================================================
# Node figures
# =============================================================================

# Annotation counts are taken from table2_node
lab_version <- function(tbl, src, unit = "nodes") {
  r <- tbl %>% filter(source == src)
  sprintf("Version %s: %d unique, %d total",
          if (src == VERSION_C) "C" else "D", r$n_unique, r$n)
}
lab_insitu <- function(tbl, type, unit = "nodes") {
  r <- tbl %>% filter(insitu_type == type)
  sprintf("%s: %d unique, %d total", type, r$n_unique, r$n)
}

# --- Figure 4a: CDF of relative node WSE by version (PT) ----------------------
t2_node_pt <- table2_node %>% filter(insitu_type == "PT")

ggplot(filter(node_all, insitu_type == "PT"),
       aes(x = abs(.data[[NODE_VALUE]]) * 100, color = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = "SWOT - PT WSE (cm)", y = "Cumulative Probability", title = "By SWOT version") +
  annotate("text", x = 40, y = 0.72, hjust = 0, color = "#222222", size = 8,
           label = sprintf("68%% C: %.1f cm, D: %.1f cm",
                           t2_node_pt$error_68ile[t2_node_pt$source == VERSION_C],
                           t2_node_pt$error_68ile[t2_node_pt$source == VERSION_D])) +
  annotate("text", x = 40, y = 0.54, hjust = 0, color = "#222222", size = 8,
           label = sprintf("50%% C: %.1f cm, D: %.1f cm",
                           t2_node_pt$error_50ile[t2_node_pt$source == VERSION_C],
                           t2_node_pt$error_50ile[t2_node_pt$source == VERSION_D])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_C]], label = lab_version(t2_node_pt, VERSION_C)) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_D]], label = lab_version(t2_node_pt, VERSION_D)) +
  scale_color_manual(values = version_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 4b: CDF of relative node WSE by in situ type (D0) -----------------
t2_node_d <- table2_node %>% filter(source == VERSION_D)

ggplot(filter(node_all, source == VERSION_D),
       aes(x = abs(.data[[NODE_VALUE]]) * 100, color = insitu_type)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
       y = "Cumulative Probability",
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 27, y = 0.72, hjust = 0, color = "#222222", size = 8,
           label = sprintf("68%% PT: %.1f cm, GNSS: %.1f cm",
                           t2_node_d$error_68ile[t2_node_d$insitu_type == "PT"],
                           t2_node_d$error_68ile[t2_node_d$insitu_type == "GNSS"])) +
  annotate("text", x = 27, y = 0.54, hjust = 0, color = "#222222", size = 8,
           label = sprintf("50%% PT: %.1f cm, GNSS: %.1f cm",
                           t2_node_d$error_50ile[t2_node_d$insitu_type == "PT"],
                           t2_node_d$error_50ile[t2_node_d$insitu_type == "GNSS"])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           color = insitu_colors[["PT"]],   label = lab_insitu(t2_node_d, "PT")) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           color = insitu_colors[["GNSS"]], label = lab_insitu(t2_node_d, "GNSS")) +
  scale_color_manual(values = insitu_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 5a: node observation count by version -----------------------------
fig4c_data <- table2_node %>%
  summarize(n = sum(n), .by = source) %>%
  mutate(v = factor(source, c(VERSION_D, VERSION_C), c("D0", "C0")))

ggplot(fig4c_data, aes(x = v, y = n, fill = v)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +  # headroom for the label
  scale_fill_manual(values = c(C0 = version_colors[[VERSION_C]],
                               D0 = version_colors[[VERSION_D]])) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")
# export: 3.16 x 6.54 in



# =============================================================================
# Reach level
# =============================================================================

read_reach <- function(path, insitu, version) {
  short <- if (version == "C") SHORT_REACHES_V16 else SHORT_REACHES_V17B
  d <- read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D) %>%
    filter(dark_frac < DARK_FRAC_MAX)
  if (APPLY_SHORT_REACH_EXCLUSION) d <- filter(d, !reach_id %in% short)
  harmonize_ids(d, "reach_id", reach_lut, version,
                label = paste("reach", insitu, version))
}

reach_PT_vC   <- read_reach(file.path(BASE, "reach/RiverSP_v16/reach_wse_SWOT_PT.csv"),    "PT",   "C")
reach_PT_vD   <- read_reach(file.path(BASE, "reach/RiverSP_v17b/reach_wse_SWOT_PT.csv"),   "PT",   "D")
reach_GNSS_vC <- read_reach(file.path(BASE, "reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv"),  "GNSS", "C")
reach_GNSS_vD <- read_reach(file.path(BASE, "reach/RiverSP_v17b/reach_wse_SWOT_GNSS.csv"), "GNSS", "D")

reach_PT <- bind_rows(reach_PT_vC, reach_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonized", "cycle_id", "pass_id", "pt_time_UTC"),
                     value_col = NODE_VALUE)

reach_GNSS <- bind_rows(reach_GNSS_vC, reach_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonized", "cycle_id", "pass_id", "drift_file"),
                     value_col = NODE_VALUE)

reach_all <- bind_rows(reach_PT, reach_GNSS) %>%
  mutate(
    insitu_wse_m        = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m),
    insitu_wse_nobias_m = coalesce(pt_wse_nobias_m,     mean_reach_drift_wse_no_bias_m),
    insitu_time_utc     = coalesce(pt_time_UTC,         wse_drift_midpoint_UTC)
  )
attr(reach_all, "partition_value_col") <- NODE_VALUE


# =============================================================================
# Reach tables
# =============================================================================

table2_reach <- reach_all %>%
  summarize_errors(NODE_VALUE, by = c("insitu_type", "source"),
                   scale = 100, digits = 1)
print(table2_reach)

tableS6_reach <- reach_all %>%
  partition_table(NODE_VALUE, by = c("insitu_type", "source"),
                  scale = 100, digits = 1)
print(tableS6_reach, n = Inf)

tableS6_reach_ids <- reach_all %>%
  filter(!is.na(id_bucket)) %>%
  summarize(n_reaches = n_distinct(id_harmonized), .by = c(insitu_type, id_bucket))
print(tableS6_reach_ids)

tableS7_reach <- reach_all %>%
  summarize_errors("residuals", by = c("insitu_type", "source"),
                   scale = 100, digits = 1)
print(tableS7_reach)

table3_reach <- reach_all %>%
  filter(source == VERSION_D, insitu_type == "PT") %>%
  merge_porcupine() %>%
  summarize_errors(NODE_VALUE, by = "river", scale = 100, digits = 1)
print(table3_reach)

reach_change <- reach_all %>% version_change(NODE_VALUE)
print(reach_change)


# =============================================================================
# Reach figures
# =============================================================================

# --- Figure 4c: CDF of relative reach WSE by version (PT) ---------------------
t2_reach_pt <- table2_reach %>% filter(insitu_type == "PT")

ggplot(filter(reach_all, insitu_type == "PT"),
       aes(x = abs(.data[[NODE_VALUE]]) * 100, color = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = "SWOT - PT WSE (cm)", y = "Cumulative Probability", title = "By SWOT version") +
  annotate("text", x = 40, y = 0.72, hjust = 0, color = "#222222", size = 8,
           label = sprintf("68%% C: %.1f cm, D: %.1f cm",
                           t2_reach_pt$error_68ile[t2_reach_pt$source == VERSION_C],
                           t2_reach_pt$error_68ile[t2_reach_pt$source == VERSION_D])) +
  annotate("text", x = 40, y = 0.54, hjust = 0, color = "#222222", size = 8,
           label = sprintf("50%% C: %.1f cm, D: %.1f cm",
                           t2_reach_pt$error_50ile[t2_reach_pt$source == VERSION_C],
                           t2_reach_pt$error_50ile[t2_reach_pt$source == VERSION_D])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_C]],
           label = sub("nodes", "reaches", lab_version(t2_reach_pt, VERSION_C))) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_D]],
           label = sub("nodes", "reaches", lab_version(t2_reach_pt, VERSION_D))) +
  scale_color_manual(values = version_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 4d: CDF of relative reach WSE by in situ type (D0) ----------------
t2_reach_d <- table2_reach %>% filter(source == VERSION_D)

ggplot(filter(reach_all, source == VERSION_D),
       aes(x = abs(.data[[NODE_VALUE]]) * 100, color = insitu_type)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = expression("SWOT -" ~ italic("in situ") ~ "WSE (cm)"),
       y = "Cumulative Probability",
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 27, y = 0.72, hjust = 0, color = "#222222", size = 8,
           label = sprintf("68%% PT: %.1f cm, GNSS: %.1f cm",
                           t2_reach_d$error_68ile[t2_reach_d$insitu_type == "PT"],
                           t2_reach_d$error_68ile[t2_reach_d$insitu_type == "GNSS"])) +
  annotate("text", x = 27, y = 0.54, hjust = 0, color = "#222222", size = 8,
           label = sprintf("50%% PT: %.1f cm, GNSS: %.1f cm",
                           t2_reach_d$error_50ile[t2_reach_d$insitu_type == "PT"],
                           t2_reach_d$error_50ile[t2_reach_d$insitu_type == "GNSS"])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           color = insitu_colors[["PT"]],
           label = sub("nodes", "reaches", lab_insitu(t2_reach_d, "PT"))) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           color = insitu_colors[["GNSS"]],
           label = sub("nodes", "reaches", lab_insitu(t2_reach_d, "GNSS"))) +
  scale_color_manual(values = insitu_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 5b: reach observation count by version ----------------------------
fig4f_data <- table2_reach %>%
  summarize(n = sum(n), .by = source) %>%
  mutate(v = factor(source, c(VERSION_D, VERSION_C), c("D0", "C0")))

ggplot(fig4f_data, aes(x = v, y = n, fill = v)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c(C0 = version_colors[[VERSION_C]],
                               D0 = version_colors[[VERSION_D]])) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")
# export: 3.16 x 6.54 in


# =============================================================================
# Inter-river figure
# =============================================================================

# --- Figure 6a: node WSE residuals by river (D0, PT) --------------------------
fig6a_data <- node_all %>%
  filter(source == VERSION_D, insitu_type == "PT", !is.na(.data[[NODE_VALUE]])) %>%
  merge_porcupine() %>%
  mutate(river = factor(river, levels = river_levels))

fig6a_counts <- fig6a_data %>% summarize(n = n(), .by = river)

ggplot(fig6a_data, aes(x = river, y = abs(.data[[NODE_VALUE]]) * 100, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, linewidth = 1) +
  geom_text(data = fig6a_counts, aes(x = river, y = -1, label = paste0("n=", n)),
            inherit.aes = FALSE, vjust = 1, size = 6) +
  labs(x = "River", y = "SWOT - PT WSE (cm)") +
  scale_fill_manual(values = river_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.9),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
  coord_cartesian(ylim = c(-5, 60))
# export: 9.44 x 6.01 in
