# =============================================================================
# SWOT Width Validation
# -----------------------------------------------------------------------------
# Compares SWOT river width against orthomosaic-derived widths at node scale.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# Outliers with |residuals| >= 1500 m are excluded from both versions.
#
# Produces: Tables 6, 7, S10, S11;  Figures 7a-b, 8a-d
# (The original header numbered these 7, 8, S5 -- renumbered to match
#  YR_CalVal_D7.)
#
# !! NOT YET VERIFIED AGAINST DATA !!
# Unlike 4.1 and 4.2 the numbers here have not been reproduced against the
# published tables. The CHECK blocks below announce the things to confirm.
#
# UPDATED FOR THE RIVEROBS PIPELINE (script 3.3)
#   * 3.3 now writes ONE node_width_SWOT_Ortho.csv covering every survey and
#     both prior-database versions, with a `sword_version` column, instead of
#     one file per version in RiverSP_v16/ and RiverSP_v17b/. This script reads
#     that file and splits it.
#   * CHECK 1 is answered: the column identifying WHICH orthomosaic a node was
#     compared against is `survey`. That is what disambiguates the two Coleen
#     acquisitions, so ORTHO_ID_COL is set to it rather than left NULL.
#   * A preflight check names any missing column, or a version present in
#     SWORD_VERSIONS but absent from the data, before any analysis runs.
#
# WHAT CHANGED IN THIS REWRITE
#   1. Exhaustive version partitioning. Table 6 vs Table S10 currently loses 23
#      D0 and 9 C0 observations. Note that the width UNIQUE NODE counts already
#      reconcile exactly (655 + 34 = 689; 655 + 84 = 739) because 4.3 was the
#      only script keying on (node_id, cycle_id, pass_id) rather than on the in
#      situ timestamp -- that was the right instinct, it just needed the
#      orthomosaic identity added and the leftover observations bucketed.
#   2. The matched subset is now symmetric by construction: id_harmonised leads
#      the key, and partition_versions() buckets only rows that already carry a
#      usable residual, so a node-overpass where one version is NA becomes
#      version-unique rather than "same".
#   3. Translator join no longer fans out rows.
#   4. Metric columns are now internally consistent. In the old 3c/3d blocks
#      the MAE was computed from residuals_nobias while the 68%ile, 50%ile and
#      RMSE in the same row were computed from residuals, so Table S10's MAE
#      column was not comparable to Table 6's.
#   5. `n` is now defined the same way in every block (the old 3a used
#      residuals, 3c used residuals_nobias, 3d used residuals).
# =============================================================================

library(tidyverse)
source("/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"

# One combined file from 3.3, both versions inside it.
WIDTH_CSV <- file.path(BASE, "node/node_width_SWOT_Ortho.csv")

# Which sword_version value is which manuscript version. Version C is the
# SWORD v16 / RiverSP PIC0 product, version D is SWORD v17b / PGD0.
SWORD_VERSIONS <- c(C = "v16", D = "v17b")

DARK_FRAC_MAX   <- 0.5
WIDTH_RESID_MAX <- 1500      # m
WIDTH_VALUE     <- "residuals"   # metric the partition is defined against

# CHECK 1 -- answered by the 3.3 rewrite --------------------------------------
# The observation key needs whatever identifies WHICH orthomosaic a node was
# compared against. Most clusters have one water mask, but the upper Porcupine
# and Coleen have two (both acquisition days), so (node, cycle, pass) alone is
# NOT unique for those. 3.3 carries `survey` -- the acquisition name, e.g.
# upperPR_CL_071024 vs upperPR_CL_071624 -- so that is the column.
# Set to NULL only if you are certain there is one mask per node.
ORTHO_ID_COL <- "survey"

# Node ids are handled as TEXT from here on, because a 14-digit id read as a
# double is one digit away from as.character() returning scientific notation and
# silently breaking every join. sprintf("%.0f") is exact for integers up to
# 2^53, which as.character() is not guaranteed to be.
as_id_chr <- function(x) {
  if (is.numeric(x)) ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  else               ifelse(is.na(x), NA_character_, trimws(as.character(x)))
}

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)

# IMPORTANT: build the lookup from the translator's NATIVE types, then convert.
# build_id_lut() resolves a v16 id that split into several v17 ids by taking the
# lowest v17 id, via arrange(). On text that ordering is lexicographic, which
# differs from numeric ordering the moment two ids have different digit counts.
# Converting first would make 4.3 pick a different child than 4.1 and 4.2 do,
# for the same v16 node.
node_lut <- build_id_lut(node_translator, "v16_node_id", "v17_node_id", "node translator")
node_ambiguous <- attr(node_lut, "ambiguous_from_ids")
node_lut <- node_lut %>% mutate(from_id = as_id_chr(from_id),
                                to_id   = as_id_chr(to_id))
attr(node_lut, "ambiguous_from_ids") <- as_id_chr(node_ambiguous)


# =============================================================================
# 1. Read, filter, harmonise
# =============================================================================

# Identifiers are read as text so a 14-digit node_id cannot be turned into a
# double and then into scientific notation, which would silently break the
# translator join.
width_raw <- read_csv(WIDTH_CSV,
                      col_types = cols(node_id  = col_character(),
                                       reach_id = col_character()),
                      show_col_types = FALSE)

# --- preflight: fail here, naming the problem, not ten lines down ------------
NEEDED <- c("sword_version", "node_id", "river", "residuals", "percent_diff",
            "width", "ortho_width_m", "bias", "dark_frac", "cycle_id",
            "pass_id", ORTHO_ID_COL)
missing_cols <- setdiff(NEEDED, names(width_raw))
if (length(missing_cols)) {
  stop(basename(WIDTH_CSV), " is missing column(s): ",
       paste(missing_cols, collapse = ", "),
       ".\nIt should be the combined output of 3.3_ortho_node_width_diff.R. ",
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
    harmonise_ids("node_id", node_lut, version, label = paste("width", version))

  # harmonise_ids() does not return NA for an unmatched id -- it builds a
  # surrogate "v16_<id>" so the node still counts as its own feature -- so a
  # total translator failure looks like every row being xlate_missing, not like
  # a column of NAs. Every such row lands in the "unmappable" bucket, which
  # would quietly empty C0_only and every version comparison built on it.
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

# CHECK 2 --------------------------------------------------------------------
# How many orthomosaics does each node appear in? A count above 1 is expected
# for the upper Porcupine and Coleen nodes covered on both acquisition days,
# and is exactly why ORTHO_ID_COL belongs in the observation key.
if (!is.null(ORTHO_ID_COL)) {
  message("[4.3] orthomosaics per node (version D):")
  width_vD %>%
    summarise(n_masks = n_distinct(.data[[ORTHO_ID_COL]]), .by = id_harmonised) %>%
    count(n_masks, name = "n_nodes") %>% print()
  message("[4.3] surveys present: ",
          paste(sort(unique(width_vD[[ORTHO_ID_COL]])), collapse = ", "))
}

# CHECK 4 --------------------------------------------------------------------
# 3.3 labels rivers from the reach_id prefix and can emit labels that
# river_levels in 4.0_comparison_helpers.R does not list -- "BL" (Black River)
# is one. Figure 8 factors `river` on river_levels, so an unlisted label becomes
# NA and those rows vanish from the figure while still counting in the tables.
unlisted_rivers <- setdiff(unique(c(width_vC$river, width_vD$river)), river_levels)
if (length(unlisted_rivers)) {
  warning("river label(s) not in river_levels: ",
          paste(unlisted_rivers, collapse = ", "),
          ". They appear in the TABLES but are dropped from figure 8. Add them ",
          "to river_levels / river_labels / river_palette in ",
          "4.0_comparison_helpers.R, or fold them into an existing river in 3.3.")
  print(bind_rows(width_vC, width_vD) %>%
          filter(river %in% unlisted_rivers) %>%
          summarise(n = n(), .by = c(source, river)))
}
if (any(is.na(c(width_vC$river, width_vD$river)))) {
  warning(sum(is.na(c(width_vC$river, width_vD$river))),
          " row(s) have river = NA. Check the river_code mapping in 3.3.")
}

# id_harmonised leads the key so that both members of a matched pair refer to
# the same node and n_unique_nodes is symmetric between versions.
width_key <- c("id_harmonised", "cycle_id", "pass_id",
               if (!is.null(ORTHO_ID_COL)) ORTHO_ID_COL)

width_all <- bind_rows(width_vC, width_vD) %>%
  partition_versions(key_cols  = width_key,
                     value_col = WIDTH_VALUE,
                     strata    = "insitu_type")


# =============================================================================
# 2. Width-specific summariser
# =============================================================================
# Width tables report both metric (m) and percent columns, so summarise_errors()
# from the helper file is extended here rather than reused directly. Every
# statistic in a row is computed from the same set of rows: those with a
# non-missing `residuals`.
#
# NOTE: the old code reported error_perdiff_50ile from the SIGNED percent_diff
# while error_perdiff_68ile used abs(percent_diff), so the two percentiles in
# the same row were not on the same scale. Both are absolute here. If the
# signed median was intentional, it belongs in the `bias` column, not next to
# an absolute 68th percentile.

summarise_width <- function(df, by) {
  df %>%
    filter(!is.na(.data[[WIDTH_VALUE]])) %>%
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
# BUCKET_LEVELS comes from the helper and includes "unmappable" -- version-C
# nodes with no v17b counterpart, which cannot pair with D0 for reasons
# unrelated to quality filtering. None of the node WSE PT or reach data has
# any; whether the width data does will show up here.
tableS10 <- bind_rows(
  width_all %>% summarise_width(by = "source") %>% mutate(bucket = "total", .before = 1),
  width_all %>% filter(!is.na(obs_bucket)) %>%
    summarise_width(by = c("source", "obs_bucket")) %>% rename(bucket = obs_bucket)
) %>%
  mutate(bucket = factor(bucket, levels = BUCKET_LEVELS)) %>%
  arrange(source, bucket)
stopifnot(!any(is.na(tableS10$bucket)))   # guards against a new bucket label
print(tableS10, n = Inf)

# CHECK 3 --------------------------------------------------------------------
# The same two assertions partition_table() makes in 4.1 and 4.2, applied here
# to the width-specific summary. If either stops, the observation key is wrong
# -- almost certainly ORTHO_ID_COL.
recon <- tableS10 %>%
  summarise(total = sum(n[bucket == "total"]),
            parts = sum(n[bucket != "total"]), .by = source)
print(recon)
stopifnot(all(recon$total == recon$parts))

same_rows <- tableS10 %>% filter(bucket == "same")
if (n_distinct(same_rows$n) > 1 || n_distinct(same_rows$n_unique_nodes) > 1) {
  print(same_rows)
  stop("The 'same' bucket is not symmetric between versions. Check width_key ",
       "(id_harmonised must be included) and ORTHO_ID_COL.")
}
message("[4.3] OK: buckets reconcile, and the matched subset is symmetric")

# Node-level membership, for the "N% fewer unique nodes in D0" statement
tableS10_ids <- width_all %>%
  filter(!is.na(id_bucket)) %>%
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
  filter(!is.na(obs_bucket), obs_bucket != "same") %>%
  summarise_width(by = c("source", "obs_bucket", "river"))
print(tableS10_by_river, n = Inf)

# Untranslatable v16 nodes, if the width data has any (the node WSE GNSS file
# has 39). Report them, do not let them sit inside the C0_only claim.
width_unmappable <- width_all %>%
  filter(obs_bucket == "unmappable") %>%
  summarise(n_obs = n(), n_nodes = n_distinct(id_harmonised), .by = river)
if (nrow(width_unmappable)) print(width_unmappable) else
  message("[4.3] no untranslatable v16 nodes in the width data")

# --- Percent-change claim in section 3.3 --------------------------------------
# Manuscript: "5.1% fewer nodes and 7.0% fewer unique nodes in version D0".
# From the published Table 6 those are 5.0% and 6.8%.
width_change <- width_all %>% version_change(WIDTH_VALUE)
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
