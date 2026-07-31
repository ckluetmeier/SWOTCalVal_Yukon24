# =============================================================================
# SWOT Width Validation
# -----------------------------------------------------------------------------
# Compares SWOT river width against orthomosaic-derived widths at node scale.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# Outliers with |residuals| >= 1500 m are excluded from both versions.
#
# Produces: Tables 6, 7, S10, S11;  Figures 7a-b, 8a-d
# =============================================================================

library(tidyverse)
source("4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"

DARK_FRAC_MAX   <- 0.5
WIDTH_RESID_MAX <- 1500   # m

# CHECK 1 --------------------------------------------------------------------
# The observation key needs whatever identifies WHICH orthomosaic a node was
# compared against. Most clusters have one water mask, but the upper Porcupine
# and Coleen have two (both acquisition days), so (node, cycle, pass) alone may
# not be unique. Set this to the column that names the mask / acquisition
# (e.g. "ortho_id", "ortho_date", "mask_file"), or to NULL if there is exactly
# one mask per node. partition_versions() warns if the key is not unique.
ORTHO_ID_COL <- "ortho_id"

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)
node_lut <- build_id_lut(node_translator, "v16_node_id", "v17_node_id", "node translator")


# =============================================================================
# 1. Read, filter, harmonise
# =============================================================================

read_width <- function(path, version) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = "Ortho",
           source      = if (version == "C") VERSION_C else VERSION_D) %>%
    filter(abs(residuals) < WIDTH_RESID_MAX, dark_frac < DARK_FRAC_MAX) %>%
    merge_porcupine() %>%
    harmonise_ids("node_id", node_lut, version, label = paste("width", version))
}

width_vC <- read_width(file.path(BASE, "node/RiverSP_v16/node_width_SWOT_Ortho.csv"),  "C")
width_vD <- read_width(file.path(BASE, "node/RiverSP_v17b/node_width_SWOT_Ortho.csv"), "D")

# CHECK 2 --------------------------------------------------------------------
# Confirm the ortho id column exists and see how many masks each node has.
if (!is.null(ORTHO_ID_COL)) {
  if (!ORTHO_ID_COL %in% names(width_vD)) {
    stop("ORTHO_ID_COL '", ORTHO_ID_COL, "' not found. Columns are: ",
         paste(names(width_vD), collapse = ", "))
  }
  width_vD %>%
    summarise(n_masks = n_distinct(.data[[ORTHO_ID_COL]]), .by = id_harmonised) %>%
    count(n_masks) %>% print()
}

width_key <- c("id_harmonised", "cycle_id", "pass_id",
               if (!is.null(ORTHO_ID_COL)) ORTHO_ID_COL)

width_all <- bind_rows(width_vC, width_vD) %>%
  partition_versions(key_cols = width_key, strata = "insitu_type")


# =============================================================================
# 2. Width-specific summarizer
# =============================================================================

summarise_width <- function(df, by) {
  df %>%
    filter(!is.na(residuals)) %>%
    summarise(
      n                   = n(),
      n_unique_nodes      = n_distinct(id_harmonised),
      med_swot_width      = round(median(width, na.rm = TRUE), 1),
      med_ortho_width     = round(median(ortho_width_m, na.rm = TRUE), 1),
      error_abs_68ile     = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
      error_abs_50ile     = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
      MAE                 = round(mean(abs(residuals), na.rm = TRUE), 1),
      RMSE                = round(sqrt(mean(residuals^2, na.rm = TRUE)), 1),
      bias                = round(median(bias, na.rm = TRUE), 1),
      error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
      error_perdiff_50ile = round(quantile(abs(percent_diff), 0.50, na.rm = TRUE), 2),
      .by = all_of(by)
    ) %>%
    arrange(across(all_of(by)))
}


# =============================================================================
# 3. TABLES
# =============================================================================

# --- Table 6: absolute node width by version ----------------------------------
table6 <- width_all %>% summarise_width(by = "source")
print(table6)

# --- Table S10: exhaustive partition ------------------------------------------
tableS10 <- bind_rows(
  width_all %>% summarise_width(by = "source")               %>% mutate(bucket = "total",  .before = 1),
  width_all %>% summarise_width(by = c("source", "obs_bucket")) %>% rename(bucket = obs_bucket)
) %>%
  mutate(bucket = factor(bucket, levels = c("total", "same", "C0_only", "D0_only"))) %>%
  arrange(source, bucket)
print(tableS10, n = Inf)

# CHECK 3 --------------------------------------------------------------------
# Hard reconciliation. If this stops, the observation key is still wrong --
# almost certainly ORTHO_ID_COL.
recon <- tableS10 %>%
  summarise(total = sum(n[bucket == "total"]),
            parts = sum(n[bucket != "total"]), .by = source)
print(recon)
stopifnot(all(recon$total == recon$parts))
message("[4.3] reconciliation OK: same + C0_only + D0_only == total")

# Node-level membership, for the "N% fewer unique nodes in D0" statement
tableS10_ids <- width_all %>%
  filter(!is.na(residuals)) %>%
  summarise(n_nodes = n_distinct(id_harmonised), .by = id_bucket)
print(tableS10_ids)

# --- Table 7: absolute node width by river (D0) -------------------------------
table7 <- width_all %>%
  filter(source == VERSION_D) %>%
  summarise_width(by = "river")
print(table7)

# Section 3.3 attributes a 47.7% |68%ile| percent difference to the braided
# Yukon, but that number is the all-rivers C0-only value from Table S10. This
# gives the river breakdown that claim actually needs.
tableS10_by_river <- width_all %>%
  filter(obs_bucket != "same") %>%
  summarise_width(by = c("source", "obs_bucket", "river"))
print(tableS10_by_river, n = Inf)

# --- Percent-change claim in section 3.3 --------------------------------------
# Manuscript: "5.1% fewer nodes and 7.0% fewer unique nodes in version D0".
# From the published Table 6 those are 5.0% and 6.8%.
width_change <- width_all %>% version_change("residuals")
print(width_change)


# =============================================================================
# 4. FIGURES
# =============================================================================

# --- Figure 7a: CDF of absolute width difference by version -------------------
ggplot(width_all, aes(x = abs(residuals), colour = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", colour = "grey") +
  labs(x = "SWOT - Orthomosaic Width (m)", y = "Cumulative Probability",
       title = "By absolute difference") +
  annotate("text", x = 215, y = 0.72, hjust = 0, colour = "#222222", size = 7,
           label = sprintf("68%% C: %.1f m, D: %.1f m",
                           table6$error_abs_68ile[table6$source == VERSION_C],
                           table6$error_abs_68ile[table6$source == VERSION_D])) +
  annotate("text", x = 215, y = 0.54, hjust = 0, colour = "#222222", size = 7,
           label = sprintf("50%% C: %.1f m, D: %.1f m",
                           table6$error_abs_50ile[table6$source == VERSION_C],
                           table6$error_abs_50ile[table6$source == VERSION_D])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           colour = version_colours[[VERSION_C]],
           label = sprintf("Version C: %d unique, %d total",
                           table6$n_unique_nodes[table6$source == VERSION_C],
                           table6$n[table6$source == VERSION_C])) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           colour = version_colours[[VERSION_D]],
           label = sprintf("Version D: %d unique, %d total",
                           table6$n_unique_nodes[table6$source == VERSION_D],
                           table6$n[table6$source == VERSION_D])) +
  scale_colour_manual(values = version_colours) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 500))
# export: 7.17 x 6.35 in

# --- Figure 7b: CDF of percent width difference by version --------------------
ggplot(width_all, aes(x = abs(percent_diff), colour = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", colour = "grey") +
  labs(x = "SWOT - Orthomosaic Width (% difference)", y = "Cumulative Probability",
       title = "By percent difference") +
  annotate("text", x = 65, y = 0.72, hjust = 0, colour = "#222222", size = 7,
           label = sprintf("68%% C: %.1f%%, D: %.1f%%",
                           table6$error_perdiff_68ile[table6$source == VERSION_C],
                           table6$error_perdiff_68ile[table6$source == VERSION_D])) +
  annotate("text", x = 65, y = 0.54, hjust = 0, colour = "#222222", size = 7,
           label = sprintf("50%% C: %.1f%%, D: %.1f%%",
                           table6$error_perdiff_50ile[table6$source == VERSION_C],
                           table6$error_perdiff_50ile[table6$source == VERSION_D])) +
  scale_colour_manual(values = version_colours) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 8a/8b: inter-river width difference (D0) --------------------------
fig8_data <- width_all %>%
  filter(source == VERSION_D, !is.na(residuals)) %>%
  mutate(river = factor(river, levels = river_levels))

fig8_counts <- fig8_data %>% summarise(n = n(), .by = river)

violin_by_river <- function(dat, counts, yvar, ylab, ylim, ylab_y) {
  ggplot(dat, aes(x = river, y = {{ yvar }}, fill = river)) +
    geom_violin(alpha = 0.8, colour = NA) +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, linewidth = 1) +
    geom_text(data = counts, aes(x = river, y = ylab_y, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = 1, size = 6) +
    labs(x = "River", y = ylab) +
    scale_fill_manual(values = river_palette, breaks = river_levels, labels = river_labels) +
    scale_x_discrete(breaks = river_levels, labels = river_labels) +
    theme_minimal(base_size = 25) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 0.9),
          plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
    coord_cartesian(ylim = ylim)
}

print(violin_by_river(fig8_data, fig8_counts, abs(percent_diff),
                      expression(atop("SWOT - Orthomosaic", " Width (% difference)")),
                      c(-5, 150), -2))
print(violin_by_river(fig8_data, fig8_counts, abs(residuals),
                      "SWOT - Orthomosaic Width (m)", c(-20, 900), -10))
# export: 9.44 x 6.01 in each


# =============================================================================
# 5. EXPORT
# =============================================================================
OUT <- file.path(BASE, "tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_csv(table6,            file.path(OUT, "table6_absolute_width.csv"))
write_csv(table7,            file.path(OUT, "table7_width_by_river.csv"))
write_csv(tableS10,          file.path(OUT, "tableS10_partition_width.csv"))
write_csv(tableS10_by_river, file.path(OUT, "tableS10_partition_width_by_river.csv"))
