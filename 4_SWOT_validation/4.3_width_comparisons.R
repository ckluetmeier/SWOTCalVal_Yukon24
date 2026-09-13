# =============================================================================
# SWOT Width Validation
# -----------------------------------------------------------------------------
# Compares SWOT river width against orthomosaic-derived widths at node scale.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
#
# Produces: Tables 6, 7, S10, S11;  Figures 7a-b, 8a-d
#
# -----------------------------------------------------------------------------
# Script by:
# Camryn Kluetmeier (camryn.kluetmeier@duke.edu)
# 
# Parts of this script were developed with assistance from Claude Code 
# (Anthropic) for debugging, documentation, and related editorial suggestions.
# 
# Last updated: 2026-09-13
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

# Root of the width comparison dataframes written by script 3.2.
BASE <- file.path(DATA_ROOT, "CalVal_dataframes/width")

# Directory holding the SWORD v16 <-> v17b ID translator CSVs.
TRANSLATOR_DIR <- file.path(SWORD_ROOT, "SWORD_translation")

# Combined SWOT/orthomosaic node width table, both SWORD versions in one file.
WIDTH_CSV <- file.path(BASE, "node/node_width_SWOT_Ortho.csv")

SWORD_VERSIONS <- c(C = "v16", D = "v17b")

DARK_FRAC_MAX   <- 0.5
WIDTH_RESID_MAX <- 1500  # m
WIDTH_VALUE     <- "residuals"

ORTHO_ID_COL <- "survey"


# -----------------------------------------------------------------------------
# Local helper and SWORD ID translator
# -----------------------------------------------------------------------------

as_id_chr <- function(x) {
  if (is.numeric(x)) ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  else               ifelse(is.na(x), NA_character_, trimws(as.character(x)))
}

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)

node_lut <- build_id_lut(node_translator, "v16_node_id", "v17_node_id", "node translator")
node_ambiguous <- attr(node_lut, "ambiguous_from_ids")
node_lut <- node_lut %>% mutate(from_id = as_id_chr(from_id),
                                to_id   = as_id_chr(to_id))
attr(node_lut, "ambiguous_from_ids") <- as_id_chr(node_ambiguous)


# =============================================================================
# 1. Read, filter, harmonize
# =============================================================================

width_raw <- read_csv(WIDTH_CSV,
                      col_types = cols(node_id  = col_character(),
                                       reach_id = col_character()),
                      show_col_types = FALSE)

NEEDED <- c("sword_version", "node_id", "river", "residuals", "percent_diff",
            "width", "ortho_width_m", "bias", "dark_frac", "cycle_id",
            "pass_id", ORTHO_ID_COL)
missing_cols <- setdiff(NEEDED, names(width_raw))
if (length(missing_cols)) {
  stop(basename(WIDTH_CSV), " is missing column(s): ",
       paste(missing_cols, collapse = ", "),
       ".\nIt should be the combined output of 3.2_ortho_node_width_diff.R. ",
       "Columns present: ", paste(names(width_raw), collapse = ", "))
}
missing_ver <- setdiff(unname(SWORD_VERSIONS), unique(width_raw$sword_version))
if (length(missing_ver)) {
  stop("sword_version value(s) not in the data: ",
       paste(missing_ver, collapse = ", "),
       ". Present: ", paste(sort(unique(width_raw$sword_version)), collapse = ", "),
       ".\nRe-run 3.3 with both versions in VERSIONS_TO_RUN, or drop the ",
       "missing one from SWORD_VERSIONS here.")
}
message(sprintf("[4.3] read %d rows: %s",
                nrow(width_raw),
                paste(sprintf("%s=%d", names(table(width_raw$sword_version)),
                              as.integer(table(width_raw$sword_version))),
                      collapse = ", ")))

read_width <- function(df, version) {
  out <- df %>%
    filter(sword_version == SWORD_VERSIONS[[version]]) %>%
    mutate(insitu_type = "Ortho",
           source      = if (version == "C") VERSION_C else VERSION_D) %>%
    filter(abs(residuals) < WIDTH_RESID_MAX, dark_frac < DARK_FRAC_MAX) %>%
    merge_porcupine() %>%
    harmonize_ids("node_id", node_lut, version, label = paste("width", version))

  if (version == "C" && nrow(out) > 0 && all(out$xlate_missing)) {
    stop("no node_id in the version-C data matched the translator, so every ",
         "row would be bucketed 'unmappable'.\n",
         "  data node_id  e.g. ", paste(utils::head(out$node_id, 2), collapse = ", "),
         "\n  translator    e.g. ", paste(utils::head(node_lut$from_id, 2), collapse = ", "),
         "\nCheck that both are the same kind of SWORD v16 node id.")
  }
  out
}

width_vC <- read_width(width_raw, "C")
width_vD <- read_width(width_raw, "D")

# -----------------------------------------------------------------------------
# 1a. Orthomosaics per node
# -----------------------------------------------------------------------------
if (!is.null(ORTHO_ID_COL)) {
  message("[4.3] orthomosaics per node (version D):")
  width_vD %>%
    summarize(n_masks = n_distinct(.data[[ORTHO_ID_COL]]), .by = id_harmonized) %>%
    count(n_masks, name = "n_nodes") %>% print()
  message("[4.3] surveys present: ",
          paste(sort(unique(width_vD[[ORTHO_ID_COL]])), collapse = ", "))
}

# -----------------------------------------------------------------------------
# 1b. River labels not listed in river_levels
# -----------------------------------------------------------------------------
unlisted_rivers <- setdiff(unique(c(width_vC$river, width_vD$river)), river_levels)
if (length(unlisted_rivers)) {
  warning("river label(s) not in river_levels: ",
          paste(unlisted_rivers, collapse = ", "),
          ". They appear in the TABLES but are dropped from figure 8. Add them ",
          "to river_levels / river_labels / river_palette in ",
          "4.0_comparison_helpers.R, or fold them into an existing river in 3.3.")
  print(bind_rows(width_vC, width_vD) %>%
          filter(river %in% unlisted_rivers) %>%
          summarize(n = n(), .by = c(source, river)))
}
if (any(is.na(c(width_vC$river, width_vD$river)))) {
  warning(sum(is.na(c(width_vC$river, width_vD$river))),
          " row(s) have river = NA. Check the river_code mapping in 3.3.")
}

# id_harmonized leads the key so that both members of a matched pair refer to
# the same node and n_unique_nodes is symmetric between versions.
width_key <- c("id_harmonized", "cycle_id", "pass_id",
               if (!is.null(ORTHO_ID_COL)) ORTHO_ID_COL)

width_all <- bind_rows(width_vC, width_vD) %>%
  partition_versions(key_cols  = width_key,
                     value_col = WIDTH_VALUE,
                     strata    = "insitu_type")


# =============================================================================
# 2. Width-specific summarize
# =============================================================================

summarize_width <- function(df, by) {
  df %>%
    filter(!is.na(.data[[WIDTH_VALUE]])) %>%
    summarize(
      n                   = n(),
      n_unique_nodes      = n_distinct(id_harmonized),
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
# 3. Tables
# =============================================================================

# --- Table 6: absolute node width by version ----------------------------------
table6 <- width_all %>% summarize_width(by = "source")
print(table6)

# --- Table S10: full partition ------------------------------------------------
tableS10 <- bind_rows(
  width_all %>% summarize_width(by = "source") %>% mutate(bucket = "total", .before = 1),
  width_all %>% filter(!is.na(obs_bucket)) %>%
    summarize_width(by = c("source", "obs_bucket")) %>% rename(bucket = obs_bucket)
) %>%
  mutate(bucket = factor(bucket, levels = BUCKET_LEVELS)) %>%
  arrange(source, bucket)
stopifnot(!any(is.na(tableS10$bucket)))
print(tableS10, n = Inf)

# -----------------------------------------------------------------------------
# 3a. Bucket reconciliation and symmetry
# -----------------------------------------------------------------------------
recon <- tableS10 %>%
  summarize(total = sum(n[bucket == "total"]),
            parts = sum(n[bucket != "total"]), .by = source)
print(recon)
stopifnot(all(recon$total == recon$parts))

same_rows <- tableS10 %>% filter(bucket == "same")
if (n_distinct(same_rows$n) > 1 || n_distinct(same_rows$n_unique_nodes) > 1) {
  print(same_rows)
  stop("The 'same' bucket is not symmetric between versions. Check width_key ",
       "(id_harmonized must be included) and ORTHO_ID_COL.")
}
message("[4.3] OK: buckets reconcile, and the matched subset is symmetric")

# Node-level membership, for the "N% fewer unique nodes in D0" statement
tableS10_ids <- width_all %>%
  filter(!is.na(id_bucket)) %>%
  summarize(n_nodes = n_distinct(id_harmonized), .by = id_bucket)
print(tableS10_ids)

# --- Table 7: absolute node width by river (D0) -------------------------------
table7 <- width_all %>%
  filter(source == VERSION_D) %>%
  summarize_width(by = "river")
print(table7)

# Percent and absolute width statistics for the version-unique buckets, broken
# down by river.
tableS10_by_river <- width_all %>%
  filter(!is.na(obs_bucket), obs_bucket != "same") %>%
  summarize_width(by = c("source", "obs_bucket", "river"))
print(tableS10_by_river, n = Inf)

# Untranslatable v16 nodes in the width data, reported separately from the
# C0_only bucket.
width_unmappable <- width_all %>%
  filter(obs_bucket == "unmappable") %>%
  summarize(n_obs = n(), n_nodes = n_distinct(id_harmonized), .by = river)
if (nrow(width_unmappable)) print(width_unmappable) else
  message("[4.3] no untranslatable v16 nodes in the width data")

# --- Percent-change stats ----------------------------------------------------
width_change <- width_all %>% version_change(WIDTH_VALUE)
print(width_change)


# =============================================================================
# 4. Figures
# =============================================================================

# --- Figure 7a: CDF of absolute width difference by version -------------------
ggplot(width_all, aes(x = abs(residuals), color = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Orthomosaic Width (m)", y = "Cumulative Probability",
       title = "By absolute difference") +
  annotate("text", x = 215, y = 0.72, hjust = 0, color = "#222222", size = 7,
           label = sprintf("68%% C: %.1f m, D: %.1f m",
                           table6$error_abs_68ile[table6$source == VERSION_C],
                           table6$error_abs_68ile[table6$source == VERSION_D])) +
  annotate("text", x = 215, y = 0.54, hjust = 0, color = "#222222", size = 7,
           label = sprintf("50%% C: %.1f m, D: %.1f m",
                           table6$error_abs_50ile[table6$source == VERSION_C],
                           table6$error_abs_50ile[table6$source == VERSION_D])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_C]],
           label = sprintf("Version C: %d unique, %d total",
                           table6$n_unique_nodes[table6$source == VERSION_C],
                           table6$n[table6$source == VERSION_C])) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           color = version_colors[[VERSION_D]],
           label = sprintf("Version D: %d unique, %d total",
                           table6$n_unique_nodes[table6$source == VERSION_D],
                           table6$n[table6$source == VERSION_D])) +
  scale_color_manual(values = version_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 500))
# export: 7.17 x 6.35 in

# --- Figure 7b: CDF of percent width difference by version --------------------
ggplot(width_all, aes(x = abs(percent_diff), color = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Orthomosaic Width (% difference)", y = "Cumulative Probability",
       title = "By percent difference") +
  annotate("text", x = 65, y = 0.72, hjust = 0, color = "#222222", size = 7,
           label = sprintf("68%% C: %.1f%%, D: %.1f%%",
                           table6$error_perdiff_68ile[table6$source == VERSION_C],
                           table6$error_perdiff_68ile[table6$source == VERSION_D])) +
  annotate("text", x = 65, y = 0.54, hjust = 0, color = "#222222", size = 7,
           label = sprintf("50%% C: %.1f%%, D: %.1f%%",
                           table6$error_perdiff_50ile[table6$source == VERSION_C],
                           table6$error_perdiff_50ile[table6$source == VERSION_D])) +
  scale_color_manual(values = version_colors) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# export: 7.17 x 6.35 in

# --- Figure 8a/8b: inter-river width difference (D0) --------------------------
fig8_data <- width_all %>%
  filter(source == VERSION_D, !is.na(residuals)) %>%
  mutate(river = factor(river, levels = river_levels))

fig8_counts <- fig8_data %>% summarize(n = n(), .by = river)

violin_by_river <- function(dat, counts, yvar, ylab, ylim, ylab_y) {
  ggplot(dat, aes(x = river, y = {{ yvar }}, fill = river)) +
    geom_violin(alpha = 0.8, color = NA) +
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
                      expression(atop("SWOT - Orthomosaic", " Width (m)")),
                      c(-20, 800), -10))
# export: 9.44 x 6.01 in each

# =============================================================================
# 5. Scatter plots and correlation (D0)
# =============================================================================

# Re-order river factor with SJ last so Sheenjek nodes are plotted on top
fig8_data$river <- factor(
  fig8_data$river,
  levels = c("CL", "CD", "PR", "upperYR", "lowerYR", "SJ")
)
# Palette reordered to match the level order above (SJ last)
color_palette_scatter <- c("#F2C14E", "#3B6064", "#F4845F", "#DA627D", "#9A348E", "#8EAD7A")

# Overall Pearson correlation (SWOT vs ortho width)
cor_test <- cor.test(fig8_data$width, fig8_data$ortho_width_m)
r_value  <- cor_test$estimate  # Pearson r
p_value  <- cor_test$p.value  # p < 0.001 is highly statistically significant


# -----------------------------------------------------------------------------
# 5a. Scatter: SWOT width vs ortho width (D0, colored by river)
# -----------------------------------------------------------------------------

ggplot() +
  # Plot non-SJ rivers first, then SJ on top so Sheenjek points are visible
  geom_point(
    data = subset(fig8_data, river != "SJ"),
    aes(x = ortho_width_m, y = width, color = river), size = 2.5
  ) +
  geom_point(
    data = subset(fig8_data, river == "SJ"),
    aes(x = ortho_width_m, y = width, color = river), size = 2.5
  ) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_color_manual(values = color_palette_scatter) +
  xlab(expression("Orthomosaic Width (m)")) +
  ylab("SWOT Width (m)") +
  ylim(0, 2700) +
  xlim(0, 2700) +
  annotate("text",
           x = min(fig8_data$ortho_width_m, na.rm = TRUE),
           y = 2700,
           label = paste0(
             "r = ", round(r_value, 4),
             "\np value = ", round(signif(p_value, 3), 4),
             "\nn = ", nrow(fig8_data)
           ),
           hjust = 0, vjust = 1, size = 8) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none")
# export dimensions: width 6.66 in, height 6.01 in
