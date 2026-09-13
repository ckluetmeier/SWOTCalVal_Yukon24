# =============================================================================
# Orthomosaic vs SWOT Reach Width Comparison — all surveys in one pass
# -----------------------------------------------------------------------------
# Reads ortho_riverobs_reaches_all.csv, the REACH-level output of the RiverObs
# pipeline (3.1.2_run_RiverObs.py), keeps the reaches listed in viable_reaches,
# matches each survey to its own SWOT overpass, and computes residuals, percent
# differences, summary tables and plots.
#
# Reach width is read from the RiverObs reach output and is a LENGTH-WEIGHTED
# mean,
#
#     width = sum(node area over observed nodes)
#             / sum(node p_length over observed nodes)
#
# which is how the operational RiverSP reach product defines it.
#
# viable_reaches is the filter that decides which reaches have enough
# orthomosaic coverage to compare. It is keyed by prior-database version,
# because a v16 reach_id does not identify the same reach in v17b.
#
# Caution - ortho_area_total_m2:
#   RiverObs writes reach area as width x the FULL prior reach length, i.e. it
#   extrapolates across unobserved nodes. It is NOT the measured water area of
#   the surveyed part of the reach. Use ortho_width_m; if you need a measured
#   area, sum the node-level ortho_area_total_m2 from
#   ortho_riverobs_nodes_all.csv instead. (Verified in SWOTRiverEstimator.py:
#   `width = masked_area / reach_area_length`, then
#   `reach_stats['area'] = width * reach_stats['length']`.)
#
# Script sections:
#   0.  Configuration - paths, filters, and reach lists
#   1.  Read orthomosaic reach widths
#   2.  Read and filter SWOT reach data
#   3.  Match ortho and SWOT observations in time and space
#   4.  Residuals and percent difference
#   5.  Data visualization
#   6.  Summary table by river
#   7.  Export
# =============================================================================

library(tidyverse)
library(lubridate)


# =============================================================================
# 0. Configuration — edit these paths before running
# =============================================================================

# Root of the field-campaign and SWOT data products.
DATA_ROOT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats"

REACH_CSV <- file.path(
  DATA_ROOT, "CalVal_dataframes/width/reach",
  "ortho_riverobs_reaches_all.csv")
OUT_CSV   <- file.path(
  DATA_ROOT, "CalVal_dataframes/width/reach",
  "reach_width_SWOT_Ortho.csv")

# SWOT reach timeseries, one per prior-database version. The key must equal the
# sword_version value in the ortho CSV.
SWOT_SOURCES <- c(
  v16  = file.path(DATA_ROOT, "SWOT/reach/RiverSP_v16",
                   "RiverSP_domain_reach_timeseries_v16.csv"),
  v17b = file.path(DATA_ROOT, "SWOT/reach/RiverSP_v17b",
                   "RiverSP_domain_reach_timeseries_PGD0_v17b.csv")
)

# Label used in the `source` column, matching VERSION_C / VERSION_D in
# 4.0_comparison_helpers.R.
SOURCE_LABEL <- c(v16 = "PIC0", v17b = "PGD0")

# --- viable reaches ----------------------------------------------------------
# Reaches where orthomosaic coverage is complete enough for a reach-level
# comparison. This is the filter for the ortho data and is applied before any
# matching. Keyed by prior-database version: a v16 reach_id does not identify
# the same stretch of river in v17b, so one list cannot serve both.
VIABLE_REACHES <- list(
  v17b = c("81260300191", "81260300181", "81260400021", "81260400011",
           "81260500011", "81260300171", "81260300161",
           "81250800031", "81270100041", "81270100051", "81270100061",
           "81270500161", "81270500171"),
  v16  = character(0)  # v16 equivalents of the reaches above are not defined
)

# Versions to process: every version whose viable_reaches list above is
# non-empty. Filling in the v16 list is therefore all that is needed to
# enable v16; assign a character vector here to override.
VERSIONS_TO_RUN <- names(VIABLE_REACHES)[lengths(VIABLE_REACHES) > 0]

# --- SWOT quality filters ----------------------------------------------------
REACH_Q_MAX   <- 2  # keep reach_q < this
XTRK_MIN      <- 10000  # cross-track distance limits, m
XTRK_MAX      <- 60000
PARTIAL_F_MAX <- 0  # 0 means >= 50% of the reach's nodes were observed
DARK_FRAC_MAX <- NA  # set to e.g. 0.8 to filter; NA disables

# --- overpass dates ----------------------------------------------------------
# Overrides the SWOT_date column in the ortho CSV, which is only as good as the
# SURVEYS list in 3.1.2_run_RiverObs.py - a date written there as "07-10-24"
# parses in R as the year 7, and the match then fails silently. Set to NULL to
# trust the CSV column; either way the resolved dates are validated below.
#
# Chandalar is the one that is not the flight date: flown 7/10, overpass 7/11.
SURVEY_DATES <- c(
  CD_071024         = "2024-07-11",
  upperPR_CL_071024 = "2024-07-10",
  upperPR_CL_071624 = "2024-07-16",
  lowerPR_SJ_072624 = "2024-07-26",
  lowerYR_071624    = "2024-07-16",
  upperYR_071024    = "2024-07-10"
)

MATCH_TOLERANCE_DAYS <- 0  # days either side of SWOT_date that still match

# --- reach ids that need a name the river_code prefix cannot give ------------
SJ_REACHES <- list(
  v16  = c("81260300061", "81260300231", "81260300241", "81260300251"),
  v17b = c("81260300181", "81260300191", "81260300201", "81260300211")
)
BL_REACHES <- list(
  v16  = c("81270100111", "81270100121", "81270100131", "81270100141",
           "81270100151", "81270100161", "81270200011", "81270200021"),
  v17b = c("81270100111", "81270100121", "81270100131", "81270100141",
           "81270100151", "81270100161", "81270200011", "81270200021")
)
# The BL list is the v16 list duplicated for v17b; only one set is defined.
# Verify the v17b ids before trusting a "BL" label there.

river_levels <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")


# =============================================================================
# 1. Read orthomosaic reach widths
# =============================================================================

# Identifiers as TEXT: an 11-digit reach_id read as a double is one step from
# scientific notation, and the joins below would then match nothing.
norm_id <- function(x) sub("\\.0+$", "", trimws(as.character(x)))

ortho_all <- read_csv(REACH_CSV, col_types = cols(reach_id = col_character()),
                      show_col_types = FALSE) %>%
  mutate(reach_id     = norm_id(reach_id),
         SWOT_date_in = as.character(SWOT_date))

REQ <- c("survey", "SWOT_date", "sword_version", "reach_id", "ortho_width_m",
         "n_good_nod", "obs_frac_n", "partial_f", "p_length")
miss <- setdiff(REQ, names(ortho_all))
if (length(miss)) {
  stop(basename(REACH_CSV), " is missing column(s): ",
       paste(miss, collapse = ", "),
       ".\nIt should be the combined reach CSV from 3.1.2_run_RiverObs.py.")
}

message(sprintf("[3.3] read %d reach rows | survey(s): %d | version(s): %s",
                nrow(ortho_all), n_distinct(ortho_all$survey),
                paste(sort(unique(ortho_all$sword_version)), collapse = ", ")))

# --- resolve and validate the overpass date ----------------------------------
if (!is.null(SURVEY_DATES)) {
  unknown <- setdiff(unique(ortho_all$survey), names(SURVEY_DATES))
  if (length(unknown)) {
    stop("SURVEY_DATES has no entry for survey(s): ",
         paste(unknown, collapse = ", "),
         ". Add them, or set SURVEY_DATES <- NULL to use the CSV column.")
  }
  ortho_all$SWOT_date <- as.Date(unname(SURVEY_DATES[ortho_all$survey]))
} else {
  ortho_all$SWOT_date <- suppressWarnings(as.Date(ortho_all$SWOT_date_in))
}
date_problem <- ortho_all %>%
  distinct(survey, SWOT_date_in, SWOT_date) %>%
  filter(is.na(SWOT_date) |
           as.integer(format(SWOT_date, "%Y")) < 2022 |
           as.integer(format(SWOT_date, "%Y")) > 2100)
if (nrow(date_problem) > 0) {
  print(as.data.frame(date_problem))
  stop("The overpass date did not parse to a plausible date for the survey(s) ",
       "above. A value like \"07-10-24\" parses as the year 7. Fix the dates ",
       "in SURVEY_DATES here, and in the SURVEYS list in ",
       "3.1.2_run_RiverObs.py, using ISO format: 2024-07-10.")
}

# --- unobserved reaches ------------------------------------------------------
# The RiverTile carries every reach in the prior database. Ones the mask never
# reached have no width.
n_all <- nrow(ortho_all)
ortho_all <- ortho_all %>% filter(!is.na(ortho_width_m), ortho_width_m > 0)
message(sprintf("[3.3] %d unobserved prior reach(es) dropped -> %d",
                n_all - nrow(ortho_all), nrow(ortho_all)))

# --- viable-reach filter -----------------------------------------------------
versions <- Reduce(intersect, list(VERSIONS_TO_RUN,
                                   names(SWOT_SOURCES),
                                   unique(ortho_all$sword_version)))
if (length(versions) == 0) {
  stop("No sword_version to process. In ", basename(REACH_CSV), ": ",
       paste(sort(unique(ortho_all$sword_version)), collapse = ", "),
       " | in SWOT_SOURCES: ", paste(names(SWOT_SOURCES), collapse = ", "),
       " | in VERSIONS_TO_RUN: ", paste(VERSIONS_TO_RUN, collapse = ", "))
}
for (v in versions) {
  if (length(VIABLE_REACHES[[v]]) == 0) {
    stop("VIABLE_REACHES has no reaches for version '", v,
         "'. Fill in the list to process this version.")
  }
}

ortho_reach <- bind_rows(lapply(versions, function(v) {
  keep_ids <- norm_id(VIABLE_REACHES[[v]])
  sub <- ortho_all %>% filter(sword_version == v)

  # A viable reach with no row here was never observed by the orthomosaic,
  # worth naming because the usual cause is a typo in the id list.
  absent <- setdiff(keep_ids, unique(sub$reach_id))
  if (length(absent)) {
    warning("version ", v, ": viable reach(es) with no orthomosaic row: ",
            paste(absent, collapse = ", "),
            ". Check the ids against ", basename(REACH_CSV), ".")
  }
  out <- sub %>% filter(reach_id %in% keep_ids)
  message(sprintf("[3.3] %s: %d of %d reach-survey row(s) kept by viable_reaches (%d distinct reaches)",
                  v, nrow(out), nrow(sub), n_distinct(out$reach_id)))
  out
}))

# Coverage of what survived. RiverObs computes reach width over the OBSERVED
# nodes only, so a partly covered reach is not biased narrow - but it does rest
# on fewer nodes, and viable_reaches exists to exclude exactly that. If anything
# here is far from 1.0, the list needs another look.
cover <- ortho_reach %>%
  select(survey, sword_version, reach_id, obs_frac_n, n_good_nod, partial_f) %>%
  arrange(obs_frac_n)
message("[3.3] lowest orthomosaic coverage among viable reaches:")
print(as.data.frame(head(cover, 5)))
if (any(ortho_reach$partial_f == 1, na.rm = TRUE)) {
  warning(sum(ortho_reach$partial_f == 1, na.rm = TRUE),
          " viable reach-survey row(s) have partial_f == 1 (under half the ",
          "reach's nodes observed in the orthomosaic).")
}

# --- river labels and version label ------------------------------------------
ortho_reach <- ortho_reach %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      sword_version == "v16"  & reach_id %in% SJ_REACHES$v16  ~ "SJ",
      sword_version == "v17b" & reach_id %in% SJ_REACHES$v17b ~ "SJ",
      sword_version == "v16"  & reach_id %in% BL_REACHES$v16  ~ "BL",
      sword_version == "v17b" & reach_id %in% BL_REACHES$v17b ~ "BL",
      river_code == "812701" ~ "lowerYR",
      river_code == "812509" ~ "lowerYR",
      river_code == "812705" ~ "upperYR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_),
    source = unname(SOURCE_LABEL[sword_version])
  )

unlabelled <- ortho_reach %>% filter(is.na(river)) %>% count(river_code)
if (nrow(unlabelled)) {
  message("[3.3] river_code with no river label:"); print(as.data.frame(unlabelled))
}
outside <- setdiff(na.omit(unique(ortho_reach$river)), river_levels)
if (length(outside)) {
  warning("river label(s) not in river_levels: ", paste(outside, collapse = ", "),
          ". They appear in the tables but are dropped from the plots.")
}

# Rename the prior-database columns the SWOT table also carries, so the join
# produces no .x / .y ambiguity.
ortho_reach <- ortho_reach %>%
  rename_with(~ paste0(.x, "_ortho"),
              any_of(c("p_length", "p_width", "obs_frac_n", "partial_f",
                       "n_good_nod", "area_total", "area_detct", "width_u")))


# =============================================================================
# 2. Read and filter SWOT reach data
# =============================================================================

tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI - UTC, seconds

# Columns this script owns; dropped from the SWOT side so they cannot collide.
PROTECTED_COLS <- c("survey", "sword_version", "SWOT_date", "SWOT_date_in",
                    "river", "river_code", "source", "ortho_width_m")

read_swot_reach <- function(path) {

  swot <- read_csv(path, col_types = cols(reach_id = col_character()),
                   show_col_types = FALSE) %>%
    mutate(reach_id = norm_id(reach_id))

  clash <- intersect(names(swot), PROTECTED_COLS)
  if (length(clash)) {
    message("  dropped ", paste(clash, collapse = ", "),
            " from the SWOT table (this script owns those names)")
    swot <- swot %>% dplyr::select(-all_of(clash))
  }

  out <- swot %>%
    # Remove duplicates and sentinel fill values (time = -999..., wse = -1e12)
    distinct(reach_id, time, wse, .keep_all = TRUE) %>%
    filter(time > 0, wse > 0) %>%
    # Quality filter. partial_f == 0 means >= 50% of the reach's nodes present.
    filter(reach_q < REACH_Q_MAX,
           abs(xtrk_dist) >= XTRK_MIN,
           abs(xtrk_dist) <= XTRK_MAX,
           partial_f <= PARTIAL_F_MAX) %>%
    mutate(time_utc = tai_epoch + time_tai - tai_utc_offset,
           obs_date = suppressWarnings(as.Date(substr(time_str, 1, 10)))) %>%
    filter(!is.na(obs_date))

  if (!is.na(DARK_FRAC_MAX)) out <- out %>% filter(dark_frac < DARK_FRAC_MAX)
  out
}


# =============================================================================
# 3. Match ortho and SWOT observations in time and space
# =============================================================================
# The join produces one row per (viable reach, matching overpass) and keeps
# the match key explicit.

matched_list <- lapply(versions, function(v) {

  message(sprintf("\n--- SWORD %s ---", v))
  swot <- read_swot_reach(SWOT_SOURCES[[v]])
  message(sprintf("  SWOT reach rows after quality filter: %d", nrow(swot)))

  ortho_v <- ortho_reach %>% filter(sword_version == v)

  out <- ortho_v %>%
    left_join(swot, by = "reach_id", suffix = c("", "_swot"),
              relationship = "many-to-many") %>%
    filter(!is.na(obs_date),
           abs(as.numeric(obs_date - SWOT_date)) <= MATCH_TOLERANCE_DAYS)

  if (nrow(out) == 0) {
    o_ids <- unique(ortho_v$reach_id); s_ids <- unique(swot$reach_id)
    message(sprintf("\n  ***** SWORD %s matched nothing -- diagnosing *****", v))
    message(sprintf("  ortho reach_ids: %d, e.g. %s", length(o_ids),
                    paste(head(o_ids, 3), collapse = ", ")))
    message(sprintf("  SWOT  reach_ids: %d, e.g. %s", length(s_ids),
                    paste(head(s_ids, 3), collapse = ", ")))
    message(sprintf("  in common: %d", length(intersect(o_ids, s_ids))))
    if (length(intersect(o_ids, s_ids)) == 0) {
      message("  CAUSE: the id sets are disjoint -- wrong SWOT file for this ",
              "prior-database version, or viable_reaches lists v16 ids here.")
    } else {
      have <- sort(unique(swot$obs_date[swot$reach_id %in% o_ids]))
      message("  CAUSE: ids match, the DATE filter empties it.")
      message("  wanted: ", paste(sort(unique(ortho_v$SWOT_date)), collapse = ", "))
      message("  available: ", paste(head(have, 10), collapse = ", "))
    }
    return(out)
  }

  dup <- out %>% count(survey, reach_id, name = "n_rows") %>% filter(n_rows > 1)
  if (nrow(dup)) {
    message(sprintf("  NOTE: %d reach(es) matched more than one SWOT row on ",
                    nrow(dup)),
            "their survey date (up to ", max(dup$n_rows), "). Two passes on the ",
            "same day, or duplicates in the SWOT file -- they double-weight ",
            "that reach in every statistic below.")
  }

  message(sprintf("  matched reach-overpass rows: %d across %d survey(s)",
                  nrow(out), n_distinct(out$survey)))
  print(as.data.frame(out %>% count(survey, name = "n_matched")))
  out
})

time_matched_SWOT_ortho <- bind_rows(matched_list)
if (nrow(time_matched_SWOT_ortho) == 0) {
  stop("Nothing matched. See the diagnostic above.")
}


# =============================================================================
# 4. Residuals and percent difference
# =============================================================================

time_matched_SWOT_ortho <- time_matched_SWOT_ortho %>%
  mutate(
    # Positive residual = SWOT narrower than the orthomosaic
    residuals    = ortho_width_m - width,
    percent_diff = (abs(ortho_width_m - width) / ortho_width_m) * 100
  )

percentile_68_error   <- quantile(abs(time_matched_SWOT_ortho$residuals),    0.68, na.rm = TRUE)
percentile_50_error   <- quantile(abs(time_matched_SWOT_ortho$residuals),    0.50, na.rm = TRUE)
percentile_68_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.68, na.rm = TRUE)
percentile_50_percent <- quantile(abs(time_matched_SWOT_ortho$percent_diff), 0.50, na.rm = TRUE)

print(paste("68th Percentile Error (m):", round(percentile_68_error, 2)))
print(paste("50th Percentile Error (m):", round(percentile_50_error, 2)))
print(paste("68th Percentile Error (%):", round(percentile_68_percent, 2)))
print(paste("50th Percentile Error (%):", round(percentile_50_percent, 2)))

cor_test <- cor.test(time_matched_SWOT_ortho$width,
                     time_matched_SWOT_ortho$ortho_width_m)
r_value  <- cor_test$estimate
p_value  <- cor_test$p.value


# =============================================================================
# 5. Data visualization
# =============================================================================

plot_df <- time_matched_SWOT_ortho %>%
  mutate(river = factor(river, levels = river_levels))

# Named so a missing river cannot silently shift every colour by one.
river_palette <- c(CL = "#F2C14E", SJ = "#8EAD7A", CD = "#3B6064",
                   PR = "#F4845F", upperYR = "#DA627D", lowerYR = "#9A348E")

# Scatter: SWOT vs ortho reach width, coloured by river
ggplot(plot_df, aes(x = ortho_width_m, y = width, color = river)) +
  geom_point(size = 2.5) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_color_manual(values = river_palette, drop = FALSE) +
  xlab("Orthomosaic Width (m)") +
  ylab("SWOT Width (m)") +
  annotate("text",
           x = min(plot_df$ortho_width_m, na.rm = TRUE),
           y = max(plot_df$width, na.rm = TRUE),
           label = paste0("r = ", round(r_value, 4),
                          "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(plot_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none")
# Export dimensions: width 6.66 in, height 6.01 in

# CDF: absolute width difference
ggplot(plot_df, aes(x = abs(residuals))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = c(0.68, 0.50), linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability",
       title = "CDF of absolute width residuals") +
  annotate("text", x = 350, y = 0.71, color = "#222222", size = 6,
           label = paste("68% abs diff:", round(percentile_68_error, 1))) +
  annotate("text", x = 350, y = 0.53, color = "#222222", size = 6,
           label = paste("50% abs diff:", round(percentile_50_error, 1))) +
  theme_minimal(base_size = 20)

# CDF: percent difference
ggplot(plot_df, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = c(0.68, 0.50), linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability",
       title = "CDF of %diff width residuals") +
  annotate("text", x = 70, y = 0.71, color = "#222222", size = 6,
           label = paste("68% abs diff:", round(percentile_68_percent, 1))) +
  annotate("text", x = 70, y = 0.53, color = "#222222", size = 6,
           label = paste("50% abs diff:", round(percentile_50_percent, 1))) +
  theme_minimal(base_size = 20)


# =============================================================================
# 6. Summary table by river
# =============================================================================

table_absolute_reach_width <- time_matched_SWOT_ortho %>%
  group_by(source, river) %>%
  summarise(
    n                   = sum(!is.na(residuals)),
    n_unique_reaches    = n_distinct(reach_id),
    med_swot_wd         = round(median(width, na.rm = TRUE), 1),
    med_ortho_wd        = round(median(ortho_width_m, na.rm = TRUE), 1),
    error_abs_68ile     = round(quantile(abs(residuals),    0.68, na.rm = TRUE), 1),
    error_abs_50ile     = round(quantile(abs(residuals),    0.50, na.rm = TRUE), 1),
    MAE                 = round(mean(abs(residuals), na.rm = TRUE), 1),
    RMSE                = round(sqrt(mean(residuals^2, na.rm = TRUE)), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(abs(percent_diff), 0.50, na.rm = TRUE), 2),
    MAE_percent         = round(mean(abs(percent_diff), na.rm = TRUE), 1),
    .groups = "drop"
  )
print(as.data.frame(table_absolute_reach_width))


# =============================================================================
# 7. Export
# =============================================================================

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write_csv(time_matched_SWOT_ortho, OUT_CSV)
message(sprintf("\nwrote %d rows to %s", nrow(time_matched_SWOT_ortho), OUT_CSV))
