# =============================================================================
# SWOT Slope Validation
# -----------------------------------------------------------------------------
# Compares SWOT water surface slope against in situ PT and GNSS measurements at
# the reach level.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
# Slope residuals are stored in m/m; multiplied by 1e5 for cm/km.
#
# Produces: Tables 4, 5, S8, S9;  Figures 5a-c, 6b
# (The original header numbered these 5, 6, S3, S4 -- that numbering no longer
#  matches the manuscript. Renumbered to match YR_CalVal_D7.)
#
# WHAT CHANGED IN THIS REWRITE (see 4.0_comparison_helpers.R)
#   1. Exhaustive version partitioning. The old Table S8 did not reconcile with
#      Table 4: GNSS D0 37 same + 10 unique = 47 against a total of 56, and the
#      PT rows had no "unique" row at all (49 same against a total of 89).
#   2. The matched subset is now symmetric. partition_versions() takes a
#      required `value_col` and buckets only rows that already carry a usable
#      slope residual, so a key where one version is NA becomes version-unique
#      rather than "same". Without this the PT "same" row came out as 54/17 for
#      PGD0 against 49/14 for PIC0. It is now 49/14 on both sides -- which is
#      exactly what the published Table S8 reports.
#   3. version_inclusion computed within in situ type, not pooled.
#   4. Translator join no longer fans out (no reach was affected in this
#      dataset, but the join is now guarded so it cannot silently start to).
#   5. Metric definitions come from one shared function.
#
# NOTE ON REACH SETS: the <9 km exclusion applied here is the same one the
# published slope tables used, and it IS the correct behaviour. Be aware the
# published reach WSE tables (Table 2/S6/S7) were built WITHOUT it, so the
# slope and WSE results currently describe slightly different reach sets.
# 4.1 now applies it in both places.
# =============================================================================

library(tidyverse)
source("4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"

DARK_FRAC_MAX <- 0.5
SHORT_REACHES_V16  <- c(81260300061, 81270500131, 81270500141)
SHORT_REACHES_V17B <- c(81260300181, 81270500021, 81270500031)

reach_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_ReachIDs_v17b_vs_v16.csv"),
                             show_col_types = FALSE)
reach_lut <- build_id_lut(reach_translator, "v16_reach_id", "v17_reach_id", "reach translator")

SCALE_SLOPE <- 100000   # m/m -> cm/km
DIGITS      <- 2
SLOPE_VALUE <- "slope_residuals_nobias"   # metric the partition is defined against


# =============================================================================
# 1. Read, filter, harmonise
# =============================================================================

read_slope <- function(path, insitu, version) {
  short <- if (version == "C") SHORT_REACHES_V16 else SHORT_REACHES_V17B
  read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D) %>%
    filter(dark_frac < DARK_FRAC_MAX) %>%
    filter(!reach_id %in% short) %>%
    merge_porcupine() %>%
    harmonise_ids("reach_id", reach_lut, version,
                  label = paste("slope", insitu, version))
}

slope_PT_vC   <- read_slope(file.path(BASE, "reach/RiverSP_v16/reach_slope_SWOT_PT.csv"),    "PT",   "C")
slope_PT_vD   <- read_slope(file.path(BASE, "reach/RiverSP_v17b/reach_slope_SWOT_PT.csv"),   "PT",   "D")
slope_GNSS_vC <- read_slope(file.path(BASE, "reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv"),  "GNSS", "C")
slope_GNSS_vD <- read_slope(file.path(BASE, "reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv"), "GNSS", "D")


# =============================================================================
# 2. Partition on version-independent keys
# =============================================================================
# Both keys lead with id_harmonised so that the two members of a matched pair
# always refer to the same reach and n_unique is symmetric between versions.

slope_PT <- bind_rows(slope_PT_vC, slope_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "pt_time_UTC"),
                     value_col = SLOPE_VALUE)

slope_GNSS <- bind_rows(slope_GNSS_vC, slope_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "drift_file"),
                     value_col = SLOPE_VALUE)

slope_all <- bind_rows(slope_PT, slope_GNSS) %>%
  mutate(
    insitu_slope_m_m        = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs),
    insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m,
                                       reach_drift_slope_m_m_abs_nobias),
    insitu_time_utc         = coalesce(pt_time_UTC, wse_drift_midpoint_UTC)  # reference only
  )
# bind_rows() drops attributes, so restore the partition tag for partition_table()
attr(slope_all, "partition_value_col") <- SLOPE_VALUE


# =============================================================================
# 3. TABLES
# =============================================================================
# bias_col = NULL: the slope products carry no single reach-level slope bias
# (the PT file has bias_us / bias_ds; the GNSS file's `bias` column is the WSE
# bias, not a slope bias). Do not report it as if it were one.

# --- Table 4: relative reach slope by version ---------------------------------
table4 <- slope_all %>%
  summarise_errors(SLOPE_VALUE, by = c("insitu_type", "source"),
                   scale = SCALE_SLOPE, digits = DIGITS, bias_col = NULL)
print(table4)

# --- Table S8: exhaustive partition -------------------------------------------
# partition_table() stops if the buckets do not sum to the total, or if the
# "same" row is not identical between the two versions.
tableS8 <- slope_all %>%
  partition_table(SLOPE_VALUE, by = c("insitu_type", "source"),
                  scale = SCALE_SLOPE, digits = DIGITS, bias_col = NULL)
print(tableS8, n = Inf)

tableS8_ids <- slope_all %>%
  filter(!is.na(id_bucket)) %>%
  summarise(n_reaches = n_distinct(id_harmonised), .by = c(insitu_type, id_bucket))
print(tableS8_ids)

# --- Table S9: absolute reach slope by version --------------------------------
# Not partitioned: `slope_residuals` has a different missingness pattern from
# `slope_residuals_nobias`, so the buckets above do not apply to it.
tableS9 <- slope_all %>%
  summarise_errors("slope_residuals", by = c("insitu_type", "source"),
                   scale = SCALE_SLOPE, digits = DIGITS, bias_col = NULL)
print(tableS9)

# --- Table 5: relative reach slope by river (D0) ------------------------------
table5 <- slope_all %>%
  filter(source == VERSION_D) %>%
  summarise_errors(SLOPE_VALUE, by = c("insitu_type", "river"),
                   scale = SCALE_SLOPE, digits = DIGITS, bias_col = NULL)
print(table5)

# --- Percent-change claim in section 3.2 --------------------------------------
# The manuscript's "30% more observations / 24% more unique reaches, aggregated
# across PT and GNSS" is actually the GNSS-only figure. This computes it
# aggregated (n_distinct pooled across in situ types) AND per in situ type, so
# whichever is quoted can be labelled correctly.
slope_change_pooled <- slope_all %>% version_change(SLOPE_VALUE)
print(slope_change_pooled)

slope_change_by_type <- slope_all %>%
  group_split(insitu_type) %>%
  map_dfr(~ version_change(.x, SLOPE_VALUE) %>%
            mutate(insitu_type = unique(.x$insitu_type), .before = 1))
print(slope_change_by_type)


# =============================================================================
# 4. FIGURES
# =============================================================================

# --- Figure 5a: CDF of relative reach slope by version (GNSS) -----------------
t4_gnss <- table4 %>% filter(insitu_type == "GNSS")

ggplot(filter(slope_all, insitu_type == "GNSS"),
       aes(x = abs(.data[[SLOPE_VALUE]]) * SCALE_SLOPE, colour = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", colour = "grey") +
  labs(x = expression("SWOT - GNSS Slope (cm km"^{-1}*")"),
       y = "Cumulative Probability", title = "By SWOT version") +
  annotate("text", x = 3, y = 0.71, hjust = 0, colour = "#222222", size = 8,
           label = sprintf("68%% C: %.2f D: %.2f cm/km",
                           t4_gnss$error_68ile[t4_gnss$source == VERSION_C],
                           t4_gnss$error_68ile[t4_gnss$source == VERSION_D])) +
  annotate("text", x = 3, y = 0.53, hjust = 0, colour = "#222222", size = 8,
           label = sprintf("50%% C: %.2f D: %.2f cm/km",
                           t4_gnss$error_50ile[t4_gnss$source == VERSION_C],
                           t4_gnss$error_50ile[t4_gnss$source == VERSION_D])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           colour = version_colours[[VERSION_C]],
           label = sprintf("Version C: %d unique, %d total",
                           t4_gnss$n_unique[t4_gnss$source == VERSION_C],
                           t4_gnss$n[t4_gnss$source == VERSION_C])) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           colour = version_colours[[VERSION_D]],
           label = sprintf("Version D: %d unique, %d total",
                           t4_gnss$n_unique[t4_gnss$source == VERSION_D],
                           t4_gnss$n[t4_gnss$source == VERSION_D])) +
  scale_colour_manual(values = version_colours) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 15))
# export: 7.17 x 6.35 in

# --- Figure 5b: CDF of relative reach slope by in situ type (D0) --------------
t4_d <- table4 %>% filter(source == VERSION_D)

ggplot(filter(slope_all, source == VERSION_D),
       aes(x = abs(.data[[SLOPE_VALUE]]) * SCALE_SLOPE, colour = insitu_type)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = c(0.50, 0.68), linetype = "dashed", colour = "grey") +
  labs(x = expression("SWOT -" ~ italic("in situ") ~ "Slope (cm km"^{-1}*")"),
       y = "Cumulative Probability",
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 3.1, y = 0.71, hjust = 0, colour = "#222222", size = 7.5,
           label = sprintf("68%% PT: %.2f GNSS: %.2f cm/km",
                           t4_d$error_68ile[t4_d$insitu_type == "PT"],
                           t4_d$error_68ile[t4_d$insitu_type == "GNSS"])) +
  annotate("text", x = 3.1, y = 0.53, hjust = 0, colour = "#222222", size = 7.5,
           label = sprintf("50%% PT: %.2f GNSS: %.2f cm/km",
                           t4_d$error_50ile[t4_d$insitu_type == "PT"],
                           t4_d$error_50ile[t4_d$insitu_type == "GNSS"])) +
  annotate("text", x = Inf, y = 0.10, hjust = 1, vjust = 0, size = 8,
           colour = insitu_colours[["PT"]],
           label = sprintf("PT: %d unique, %d total",
                           t4_d$n_unique[t4_d$insitu_type == "PT"],
                           t4_d$n[t4_d$insitu_type == "PT"])) +
  annotate("text", x = Inf, y = 0.01, hjust = 1, vjust = 0, size = 8,
           colour = insitu_colours[["GNSS"]],
           label = sprintf("GNSS: %d unique, %d total",
                           t4_d$n_unique[t4_d$insitu_type == "GNSS"],
                           t4_d$n[t4_d$insitu_type == "GNSS"])) +
  scale_colour_manual(values = insitu_colours) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 15))
# export: 7.17 x 6.35 in

# --- Figure 5c: observation count by version ----------------------------------
ggplot(table4 %>% mutate(v = factor(source, c(VERSION_D, VERSION_C), c("D", "C"))),
       aes(x = v, y = n, fill = v)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  scale_fill_manual(values = c(C = version_colours[[VERSION_C]], D = version_colours[[VERSION_D]])) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "none")
# export: 3.16 x 6.54 in

# --- Figure 6b: inter-river slope residuals (D0, PT) --------------------------
fig6b_data <- slope_all %>%
  filter(source == VERSION_D, insitu_type == "PT", !is.na(.data[[SLOPE_VALUE]])) %>%
  mutate(river = factor(river, levels = river_levels))

fig6b_counts <- fig6b_data %>% summarise(n = n(), .by = river)

ggplot(fig6b_data, aes(x = river, y = abs(.data[[SLOPE_VALUE]]) * SCALE_SLOPE, fill = river)) +
  geom_violin(alpha = 0.8, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, linewidth = 1) +
  geom_text(data = fig6b_counts, aes(x = river, y = -0.2, label = paste0("n=", n)),
            inherit.aes = FALSE, vjust = 1, size = 6) +
  labs(x = "River", y = expression("SWOT - PT Slope (cm km"^{-1}*")")) +
  scale_fill_manual(values = river_palette, breaks = river_levels, labels = river_labels) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.9),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
  coord_cartesian(ylim = c(-0.5, 16))
# export: 9.44 x 6.01 in


# =============================================================================
# 5. DIAGNOSTICS — covariate scatter plots (D0)
# =============================================================================
# Section 3.2 states slope error is uncorrelated with dark_frac, cross-track
# distance, layover, node count, slope and width. These plots support that
# claim; add a printed correlation so the statement is quantified rather than
# eyeballed.

slope_d <- slope_all %>% filter(source == VERSION_D, !is.na(.data[[SLOPE_VALUE]]))

covariates <- c("layovr_val", "xtrk_dist", "dark_frac", "width", "n_good_nod",
                "insitu_slope_m_m")

slope_correlations <- map_dfr(covariates, function(v) {
  if (!v %in% names(slope_d)) return(NULL)
  slope_d %>%
    filter(!is.na(.data[[v]])) %>%
    summarise(
      covariate = v,
      n         = n(),
      spearman  = round(cor(abs(.data[[v]]), abs(.data[[SLOPE_VALUE]]),
                            method = "spearman", use = "complete.obs"), 3),
      .by = insitu_type
    )
})
print(slope_correlations)

walk(covariates, function(v) {
  if (!v %in% names(slope_d)) return(invisible(NULL))
  p <- ggplot(slope_d, aes(x = abs(.data[[v]]),
                           y = abs(.data[[SLOPE_VALUE]]) * SCALE_SLOPE,
                           colour = river)) +
    geom_point(size = 3) +
    facet_wrap(~ insitu_type, scales = "free_y") +
    scale_colour_manual(values = river_palette, breaks = river_levels, labels = river_labels) +
    labs(x = v, y = expression("Slope error (cm km"^{-1}*")"), colour = "River") +
    theme_minimal(base_size = 20)
  print(p)
})


# =============================================================================
# 6. EXPORT
# =============================================================================
OUT <- file.path(BASE, "tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_csv(table4,  file.path(OUT, "table4_relative_slope.csv"))
write_csv(tableS8, file.path(OUT, "tableS8_partition_slope.csv"))
write_csv(tableS9, file.path(OUT, "tableS9_absolute_slope.csv"))
write_csv(table5,  file.path(OUT, "table5_slope_by_river.csv"))
