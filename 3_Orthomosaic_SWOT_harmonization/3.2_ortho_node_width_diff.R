# =============================================================================
# Orthomosaic vs SWOT Node Width Comparison
# -----------------------------------------------------------------------------
# Reads the outputs from the RiverObs pipeline:
#
#   ortho_riverobs_nodes_all.csv   every survey x SWORD version, one row per
#                                  node, with RiverObs width
#   node_qc_all.csv                the QGIS exclusion list, one row
#                                  per survey x version x node
#
# joins them to the SWOT node timeseries, matches each survey to its own
# overpass date, computes residuals and percent differences, removes per-group
# median bias, and writes one matched CSV covering every survey.
#
# Script sections:
#   0.  Configuration - paths, filters, and overpass dates
#   1.  Read orthomosaic node widths (RiverObs output)
#   2.  Apply the manual node QC
#   3.  Read and filter SWOT node data
#   4.  Join, and match each survey to its own overpass date
#   5.  Residuals and percent difference
#   6.  Bias removal and bias-corrected metrics
#   7.  River labels
#   8.  Export ONE matched CSV
#   9.  Data visualization
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
library(lubridate)


# =============================================================================
# 0. Configuration - edit these paths before running
# =============================================================================

# Root of the field-campaign and SWOT data products.
DATA_ROOT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats"

ORTHO_NODES_CSV <- file.path(
  DATA_ROOT, "CalVal_dataframes/width/node", "ortho_riverobs_nodes_all.csv")
NODE_QC_CSV     <- file.path(
  DATA_ROOT, "CalVal_dataframes/width/node", "node_qc_all.csv")

SWOT_SOURCES <- c(
  v16  = file.path(DATA_ROOT, "SWOT/node/hydrocron_timeseries",
                   "YR_domain_nodes_merged_RiverSP.csv"),
  v17b = file.path(DATA_ROOT, "SWOT/node/RiverSP_v17b",
                   "RiverSP_domain_node_timeseries_PGD0_v17b.csv")
)

VERSIONS_TO_RUN <- c("v16", "v17b")

OUT_CSV <- file.path(
  DATA_ROOT, "CalVal_dataframes/width/node", "node_width_SWOT_Ortho.csv")

# --- filters -----------------------------------------------------------------
APPLY_NODE_QC <- TRUE
MIN_GOOD_PIX  <- 1  # minimum ortho water pixels for a node to be usable

NODE_Q_MAX    <- 2  # node_q  (0 good, 1 suspect, 2 degraded, 3 bad)
XTRK_MIN      <- 10000  # cross-track distance limits, m
XTRK_MAX      <- 60000
DARK_FRAC_MAX <- 0.8

# --- overpass dates ----------------------------------------------------------
# The SWOT overpass date for each survey.
#
# (Chandalar flown 7/10, overpass 7/11)
SURVEY_DATES <- c(
  CD_071024         = "2024-07-11",
  upperPR_CL_071024 = "2024-07-10",
  upperPR_CL_071624 = "2024-07-16",
  lowerPR_SJ_072624 = "2024-07-26",
  lowerYR_071624    = "2024-07-16",
  upperYR_071024    = "2024-07-10"
)

MATCH_TOLERANCE_DAYS <- 0

# --- bias grouping -----------------------------------------------------------
# One median bias per survey, SWORD version and reach code.
BIAS_GROUP <- c("survey", "sword_version", "river_code")

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


# =============================================================================
# 1. Read orthomosaic node widths (RiverObs output)
# =============================================================================

norm_id <- function(x) sub("\\.0+$", "", trimws(as.character(x)))

ID_COLS <- cols(node_id = col_character(), reach_id = col_character())

ortho_all <- read_csv(ORTHO_NODES_CSV, col_types = ID_COLS,
                      show_col_types = FALSE) %>%
  mutate(
    node_id      = norm_id(node_id),
    reach_id     = norm_id(reach_id),
    SWOT_date_in = as.character(SWOT_date)
  )

ORTHO_REQUIRED <- c("survey", "SWOT_date", "sword_version", "reach_id",
                    "node_id", "lat", "lon", "ortho_width_m",
                    "ortho_area_total_m2", "n_good_pix", "p_length",
                    "p_width", "p_dist_out")
missing_ortho <- setdiff(ORTHO_REQUIRED, names(ortho_all))
if (length(missing_ortho)) {
  stop(basename(ORTHO_NODES_CSV), " is missing column(s): ",
       paste(missing_ortho, collapse = ", "),
       ".\nIt should be the combined CSV written by 3.1.2_run_RiverObs.py.")
}

# --- confirm the overpass date --------------------------------------------
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
message(sprintf("overpass dates in use: %s",
                paste(sprintf("%s=%s",
                              sort(unique(ortho_all$survey)),
                              ortho_all$SWOT_date[match(
                                sort(unique(ortho_all$survey)),
                                ortho_all$survey)]),
                      collapse = "  ")))

message(sprintf("ortho nodes read: %d rows, %d survey(s), version(s): %s",
                nrow(ortho_all), n_distinct(ortho_all$survey),
                paste(sort(unique(ortho_all$sword_version)), collapse = ", ")))

ortho <- ortho_all %>%
  filter(!is.na(ortho_width_m), n_good_pix >= MIN_GOOD_PIX)

message(sprintf("  %d unobserved prior nodes dropped -> %d measured nodes",
                nrow(ortho_all) - nrow(ortho), nrow(ortho)))

# RiverObs writes width = area_total / p_length.
chk <- ortho %>%
  mutate(implied = ortho_area_total_m2 / p_length,
         rel_err = abs(implied - ortho_width_m) / pmax(ortho_width_m, 1e-9))
if (max(chk$rel_err, na.rm = TRUE) > 1e-6) {
  warning(sprintf(
    "ortho_width_m does not equal area_total / p_length (max rel. err %.2e). ",
    max(chk$rel_err, na.rm = TRUE)),
    "Check that the CSV came from 3.1.2_run_RiverObs.py.")
}


# =============================================================================
# 2. Apply the manual node QC
# =============================================================================
# node_qc_all.csv is keyed on survey + sword_version + node_id

if (APPLY_NODE_QC && file.exists(NODE_QC_CSV)) {

  node_qc <- read_csv(NODE_QC_CSV, col_types = cols(node_id = col_character()),
                      show_col_types = FALSE) %>%
    mutate(node_id = norm_id(node_id)) %>%
    distinct(survey, sword_version, node_id, .keep_all = TRUE)

  ortho <- ortho %>%
    left_join(node_qc %>% dplyr::select(survey, sword_version, node_id,
                                        keep, any_of("excl_frac")),
              by = c("survey", "sword_version", "node_id"))

  n_unreviewed <- sum(is.na(ortho$keep))
  if (n_unreviewed > 0) {
    warning(sprintf(
      "%d node(s) have no row in node_qc_all.csv and were kept. Re-run 3.1.4_manual_node_qc.py ",
      n_unreviewed),
      "if the QC should cover them.")
  }

  n_before <- nrow(ortho)
  ortho <- ortho %>% filter(is.na(keep) | keep == 1)
  message(sprintf("  node QC: %d node(s) excluded -> %d remain",
                  n_before - nrow(ortho), nrow(ortho)))

} else if (APPLY_NODE_QC) {
  warning("APPLY_NODE_QC is TRUE but ", NODE_QC_CSV, " does not exist. ",
          "No nodes were excluded.")
  ortho$keep <- NA_integer_
} else {
  message("  node QC skipped (APPLY_NODE_QC = FALSE)")
  ortho$keep <- NA_integer_
}

# Rename the prior-database columns the SWOT table also has
ortho <- ortho %>%
  rename(reach_id_ortho   = reach_id,
         lat_ortho         = lat,
         lon_ortho         = lon,
         n_good_pix_ortho  = n_good_pix,
         p_length_ortho    = p_length,
         p_width_ortho     = p_width,
         p_dist_out_ortho  = p_dist_out)


# =============================================================================
# 3. Read and filter SWOT node data
# =============================================================================

tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI - UTC, seconds

PROTECTED_COLS <- c("survey", "sword_version", "SWOT_date", "SWOT_date_in",
                    "keep", "excl_frac")

read_swot <- function(path) {

  swot <- read_csv(path, col_types = ID_COLS, show_col_types = FALSE) %>%
    mutate(node_id  = norm_id(node_id),
           reach_id = norm_id(reach_id))

  clash <- intersect(names(swot), PROTECTED_COLS)
  if (length(clash)) {
    message(sprintf("  dropped %s from the SWOT table (this script owns %s)",
                    paste(clash, collapse = ", "),
                    if (length(clash) == 1) "that name" else "those names"))
    swot <- swot %>% dplyr::select(-all_of(clash))
  }

  swot %>%
    # Remove duplicates and fill values (time = -999..., wse = -1e12)
    distinct(node_id, time, wse, .keep_all = TRUE) %>%
    # Quality filter: node_q < 2, valid cross-track swath, < 80% dark water
    filter(node_q < NODE_Q_MAX,
           abs(xtrk_dist) >= XTRK_MIN,
           abs(xtrk_dist) <= XTRK_MAX,
           dark_frac <= DARK_FRAC_MAX) %>%
    mutate(
      time_utc = tai_epoch + time_tai - tai_utc_offset,
      # The overpass date used for matching comes from time_str, the string
      # SWOT wrote.
      obs_date = suppressWarnings(as.Date(substr(time_str, 1, 10)))
    ) %>%
    filter(!is.na(obs_date))
}

# -----------------------------------------------------------------------------
# Bonk check
# -----------------------------------------------------------------------------

diagnose_no_match <- function(ortho_v, swot, v) {

  o_ids  <- unique(ortho_v$node_id)
  s_ids  <- unique(swot$node_id)
  common <- intersect(o_ids, s_ids)

  message(sprintf("\n  ***** SWORD %s matched nothing -- diagnosing *****", v))
  message(sprintf("  ortho node_ids: %d unique, e.g. %s  (%s, %d chars)",
                  length(o_ids), paste(utils::head(o_ids, 3), collapse = ", "),
                  class(ortho_v$node_id)[1], nchar(o_ids[1])))
  message(sprintf("  SWOT  node_ids: %d unique, e.g. %s  (%s, %d chars)",
                  length(s_ids), paste(utils::head(s_ids, 3), collapse = ", "),
                  class(swot$node_id)[1], nchar(s_ids[1])))
  message(sprintf("  node_ids in common: %d", length(common)))

  if (length(common) == 0) {
    message("\n  CAUSE: the id sets are disjoint. The join can never produce a row.")
    if (nchar(o_ids[1]) != nchar(s_ids[1])) {
      message("  The two are different lengths, so they are not the same kind of")
      message("  identifier -- check that both files use SWORD node_id and that")
      message("  neither was mangled into scientific notation on read.")
    } else {
      message("  Same length but no overlap: almost certainly the SWOT file for")
      message(sprintf("  '%s' is built against a different prior database.", v))
      message("  Check SWOT_SOURCES -- a v16 SWOT file cannot match v17b nodes.")
    }
    return(invisible(NULL))
  }

  # ids are fine, so it is the dates
  sub  <- swot %>% filter(node_id %in% common)
  have <- sort(unique(sub$obs_date))
  want <- sort(unique(ortho_v$SWOT_date))
  message("\n  CAUSE: node_ids match, so the DATE filter is what empties it.")
  message(sprintf("  SWOT_date wanted (from the ortho CSV): %s",
                  paste(want, collapse = ", ")))
  message(sprintf("  overpass dates present for those nodes: %d distinct, range %s to %s",
                  length(have), min(have), max(have)))
  message(sprintf("  nearest available dates: %s",
                  paste(utils::head(have[order(abs(as.numeric(
                    have - want[1])))], 8), collapse = ", ")))
  message("  Fix the SWOT_date values in the SURVEYS list in")
  message("  3.1.2_run_RiverObs.py, or set MATCH_TOLERANCE_DAYS below.")
  invisible(NULL)
}


# =============================================================================
# 4. Join, and match each survey to its own overpass date
# =============================================================================
#   chandalar (CD_071024)         2024-07-11
#   coleen / upper PR 7/10        2024-07-10
#   coleen / upper PR 7/16        2024-07-16
#   sheenjek (lowerPR_SJ_072624)  2024-07-26
#   lower YR (lowerYR_071624)     2024-07-16
#   upper YR (upperYR_071024)     2024-07-10

versions <- Reduce(intersect, list(VERSIONS_TO_RUN,
                                   names(SWOT_SOURCES),
                                   unique(ortho$sword_version)))
if (length(versions) == 0) {
  stop("No sword_version to process. In ", basename(ORTHO_NODES_CSV), ": ",
       paste(sort(unique(ortho$sword_version)), collapse = ", "),
       " | in SWOT_SOURCES: ", paste(names(SWOT_SOURCES), collapse = ", "),
       " | in VERSIONS_TO_RUN: ", paste(VERSIONS_TO_RUN, collapse = ", "))
}

matched_list <- lapply(versions, function(v) {

  message(sprintf("\n--- SWORD %s ---", v))
  swot <- read_swot(SWOT_SOURCES[[v]])
  message(sprintf("  SWOT rows after quality filter: %d", nrow(swot)))

  ortho_v <- ortho %>% filter(sword_version == v)

  out <- ortho_v %>%
    left_join(swot, by = "node_id",
              suffix = c("", "_swot"),
              relationship = "many-to-many") %>%
    # Keep only the overpass matching each survey's own date.
    # MATCH_TOLERANCE_DAYS is 0 by default, i.e. the same calendar day.
    filter(!is.na(obs_date),
           abs(as.numeric(obs_date - SWOT_date)) <= MATCH_TOLERANCE_DAYS)

  if (nrow(out) == 0) {
    diagnose_no_match(ortho_v, swot, v)
    return(out)
  }

  # After the date filter each survey/version/node should appear at most once.
  dup <- out %>%
    count(survey, sword_version, node_id, name = "n_rows") %>%
    filter(n_rows > 1)
  if (nrow(dup) > 0) {
    message(sprintf(
      "  NOTE: %d node(s) matched more than one SWOT row on their survey date ",
      nrow(dup)))
    message("  (up to ", max(dup$n_rows), " rows). Two passes on the same day, ",
            "or duplicate rows in the")
    message("  SWOT file. Inspect with: out %>% count(survey, node_id) %>% ",
            "filter(n > 1)")
  }

  bad <- sum(out$reach_id_ortho != out$reach_id, na.rm = TRUE)
  if (bad > 0) {
    warning(sprintf(
      "SWORD %s: %d node(s) have different reach_id in the ortho and SWOT ",
      v, bad),
      "tables. The SWOT file may not be the ", v, " product.")
  }

  message(sprintf("  matched node-overpass rows: %d across %d survey(s)",
                  nrow(out), n_distinct(out$survey)))
  if (nrow(out) > 0) {
    per_survey <- out %>% count(survey, name = "n_matched")
    for (i in seq_len(nrow(per_survey))) {
      message(sprintf("     %-20s %5d", per_survey$survey[i],
                      per_survey$n_matched[i]))
    }
  }
  # Surveys that matched nothing
  missing <- setdiff(unique(ortho_v$survey), unique(out$survey))
  if (length(missing)) {
    warning("SWORD ", v, ": no SWOT match for survey(s): ",
            paste(missing, collapse = ", "),
            ". Check SWOT_date against the overpass dates in the SWOT file.")
  }
  out
})

# Two SWOT products are two files written by two pipelines
harmonise <- function(frames) {
  frames <- Filter(function(f) nrow(f) > 0, frames)
  if (length(frames) < 2) return(frames)
  classes <- lapply(frames, function(f) vapply(f, function(x) class(x)[1], ""))
  all_cols <- unique(unlist(lapply(classes, names)))
  conflict <- all_cols[vapply(all_cols, function(cl) {
    seen <- unlist(lapply(classes, function(k) {
      if (cl %in% names(k)) unname(k[[cl]]) else NULL
    }))
    length(unique(seen)) > 1
  }, logical(1))]
  if (length(conflict)) {
    message("\ncolumn(s) with different types across SWORD versions, coerced ",
            "to text: ", paste(conflict, collapse = ", "))
    frames <- lapply(frames, function(f) {
      f %>% mutate(across(any_of(conflict), as.character))
    })
  }
  frames
}

combined <- bind_rows(harmonise(matched_list))

if (nrow(combined) == 0) {
  stop("Nothing matched. Check SWOT_SOURCES paths and the SWOT_date values.")
}


# =============================================================================
# 5. Residuals and percent difference
# =============================================================================

combined <- combined %>%
  mutate(
    # Positive residual = SWOT narrower than the orthomosaic
    residuals    = ortho_width_m - width,
    percent_diff = (abs(ortho_width_m - width) / ortho_width_m) * 100
  )


# =============================================================================
# 6. Bias removal and bias-corrected metrics
# =============================================================================
# river_code is the first 6 digits of reach_id

combined <- combined %>%
  mutate(river_code = substr(reach_id, 1, 6)) %>%
  group_by(across(all_of(BIAS_GROUP))) %>%
  mutate(
    n_in_bias_group     = n(),
    bias                = median(residuals, na.rm = TRUE),
    swot_width_nobias_m = width + bias  # shift SWOT width toward ortho
  ) %>%
  ungroup() %>%
  mutate(
    residuals_nobias    = swot_width_nobias_m - ortho_width_m,
    percent_diff_nobias = (abs(ortho_width_m - swot_width_nobias_m) /
                             ortho_width_m) * 100
  )

thin <- combined %>%
  distinct(across(all_of(BIAS_GROUP)), n_in_bias_group) %>%
  filter(n_in_bias_group < 5)
if (nrow(thin) > 0) {
  message(sprintf("\n%d bias group(s) have fewer than 5 nodes:", nrow(thin)))
  print(as.data.frame(thin))
}


# =============================================================================
# 7. River labels
# =============================================================================

combined <- combined %>%
  mutate(
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
      TRUE ~ NA_character_
    )
  )

unlabelled <- combined %>% filter(is.na(river)) %>% count(river_code, name = "n")
if (nrow(unlabelled) > 0) {
  message("\nriver_code values with no river label (river = NA):")
  print(as.data.frame(unlabelled))
}


# =============================================================================
# 8. Export ONE matched CSV
# =============================================================================

REQUIRED_COLS <- c(
  "survey", "sword_version", "SWOT_date", "SWOT_date_in", "obs_date",
  "node_id", "reach_id", "river", "river_code",
  "residuals", "percent_diff", "bias", "n_in_bias_group",
  "residuals_nobias", "percent_diff_nobias",
  "ortho_width_m", "ortho_area_total_m2", "n_good_pix_ortho", "p_length_ortho",
  "width"
)

# Columns from SWOT versions
OPTIONAL_COLS <- c(
  "time_utc", "ortho_area_detct_m2", "p_width_ortho", "p_dist_out_ortho",
  "keep", "excl_frac",
  "lat", "lon", "wse", "wse_u", "wse_r_u",
  "width_u", "area_total", "area_tot_u", "area_detct", "area_det_u",
  "area_wse", "layovr_val", "node_dist", "xtrk_dist",
  "node_q", "node_q_b", "dark_frac", "n_good_pix", "rdr_sig0", "xovr_cal_q",
  "cycle_id", "pass_id", "p_dist_out", "p_length"
)

missing_req <- setdiff(REQUIRED_COLS, names(combined))
if (length(missing_req)) {
  stop("these columns should exist by now and do not: ",
       paste(missing_req, collapse = ", "),
       ". Something upstream in this script did not run as expected.")
}
missing_opt <- setdiff(OPTIONAL_COLS, names(combined))
if (length(missing_opt)) {
  message("\nnot in the SWOT file(s), omitted from the output: ",
          paste(missing_opt, collapse = ", "))
}

save_to_csv <- combined %>%
  dplyr::select(all_of(REQUIRED_COLS), any_of(OPTIONAL_COLS))

write.csv(save_to_csv, file = OUT_CSV, row.names = FALSE)

message(sprintf("\nwrote %d rows to %s", nrow(save_to_csv), OUT_CSV))
message("summary by survey and version:")
print(as.data.frame(
  combined %>%
    group_by(survey, sword_version) %>%
    summarise(n            = n(),
              med_pct_diff = round(median(percent_diff, na.rm = TRUE), 2),
              med_pct_nb   = round(median(percent_diff_nobias, na.rm = TRUE), 2),
              bias_m       = round(median(bias, na.rm = TRUE), 2),
              .groups = "drop")
))


# =============================================================================
# 9. Data visualization
# =============================================================================

# --- SWOT vs ortho width, colored by river ----------------------------------
cor_test <- cor.test(combined$width, combined$ortho_width_m)
r_value  <- cor_test$estimate
p_value  <- cor_test$p.value

ggplot(combined, aes(x = ortho_width_m, y = width, color = river)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(linetype = "dashed", color = "gray") +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  annotate("text",
           x = min(combined$ortho_width_m, na.rm = TRUE),
           y = max(combined$width, na.rm = TRUE),
           label = paste0("r = ", round(r_value, 4),
                          "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(combined)),
           hjust = 0, vjust = 1, size = 6) +
  labs(color = "River") +
  theme_minimal(base_size = 20)

# --- width along distance to outlet, one panel per survey --------------------
ggplot(combined) +
  geom_point(aes(x = p_dist_out / 1000, y = ortho_width_m),
             color = "lightblue", size = 2.5, shape = 17) +
  geom_point(aes(x = p_dist_out / 1000, y = swot_width_nobias_m),
             color = "darkblue", size = 2.5, alpha = 0.7) +
  facet_wrap(~ survey, scales = "free_x") +
  xlab("Distance to outlet (km)") +
  ylab("Width (m)") +
  theme_minimal(base_size = 16)

# --- dark water fraction vs percent difference -------------------------------
cor_dark <- cor.test(abs(combined$dark_frac), combined$percent_diff)

ggplot(combined, aes(x = dark_frac, y = percent_diff, color = river)) +
  geom_point(size = 2.5, alpha = 0.8) +
  xlab("Dark water fraction") +
  ylab("Percent difference (%)") +
  annotate("text",
           x = max(combined$dark_frac, na.rm = TRUE) - 0.28,
           y = max(combined$percent_diff, na.rm = TRUE),
           label = paste0("r = ", round(cor_dark$estimate, 4),
                          "\np value = ",
                          round(signif(cor_dark$p.value, 3), 5)),
           hjust = 0, vjust = 1, size = 6) +
  theme_minimal(base_size = 20)

# --- percent difference by overpass date (bias-corrected) --------------------
ggplot(combined, aes(x = factor(obs_date), y = percent_diff_nobias,
                     fill = factor(obs_date))) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 2, lwd = 0.8) +
  xlab("Overpass date") +
  ylab("Percent difference (%)") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.9))

# --- percent difference by river (raw) ---------------------------------------
ggplot(combined %>% filter(!is.na(river)),
       aes(x = river, y = percent_diff, fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 2, lwd = 0.8) +
  xlab("River") +
  ylab("Percent difference (%)") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.9))
