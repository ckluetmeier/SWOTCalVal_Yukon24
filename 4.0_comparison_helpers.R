# =============================================================================
# Shared helpers for SWOT version C0 / D0 comparison scripts (4.1, 4.2, 4.3)
# -----------------------------------------------------------------------------
# Written for R >= 4.4 and dplyr >= 1.1 (uses .by=, relationship=, unmatched=).
# =============================================================================

library(tidyverse)

VERSION_C <- "PIC0"   # SWOT version C0, SWORD v16
VERSION_D <- "PGD0"   # SWOT version D0, SWORD v17b


# -----------------------------------------------------------------------------
# build_id_lut() — collapse a SWORD translator into a strict 1:1 v16 -> v17 map
# -----------------------------------------------------------------------------
# The translator files (NA_NodeIDs_v17b_vs_v16.csv, NA_ReachIDs_v17b_vs_v16.csv)
# carry ONE ROW PER v17 ID. A v16 id that was split into several v17 ids
# therefore appears on several rows. Joining that table directly duplicates the
# version-C observation once per v17, which inflates the C0 counts.
#
# In the Yukon domain this affected 25 v16 nodes (each splitting into exactly 2
# v17 nodes) and no reaches. The translator provides no field that can pick a
# winner (boundary_percent is 0 for all of them), so we take the lowest v17 id
# deterministically and flag the affected rows.
#
# Returns a two-column tibble: from_id, to_id  (guaranteed unique on from_id).
# -----------------------------------------------------------------------------
build_id_lut <- function(translator, from_col, to_col, label = "translator") {

  lut_all <- translator %>%
    select(from_id = all_of(from_col), to_id = all_of(to_col)) %>%
    filter(!is.na(from_id), !is.na(to_id)) %>%
    distinct()

  ambiguous <- lut_all %>%
    summarise(n_children = n_distinct(to_id), .by = from_id) %>%
    filter(n_children > 1)

  lut <- lut_all %>%
    arrange(from_id, to_id) %>%
    slice_head(n = 1, by = from_id)

  message(sprintf(
    "[%s] %s unique v16 ids; %s split into >1 v17 id (lowest v17 id used, rows flagged)",
    label, format(nrow(lut), big.mark = ","), format(nrow(ambiguous), big.mark = ",")))

  attr(lut, "ambiguous_from_ids") <- ambiguous$from_id
  lut
}


# -----------------------------------------------------------------------------
# harmonise_ids() — apply the 1:1 map without ever changing the row count
# -----------------------------------------------------------------------------
# Adds three columns:
#   id_harmonised    v17b id as character; for v16 ids with no v17 counterpart,
#                    "v16_<original id>" so each unmatched id still counts as
#                    its own distinct feature instead of collapsing into a
#                    single NA pseudo-node (which is what n_distinct() did
#                    before -- 39 untranslatable v16 nodes / 53 GNSS node
#                    observations were being counted as ONE node).
#   xlate_ambiguous  TRUE if the v16 id split into >1 v17 id
#   xlate_missing    TRUE if the v16 id has no v17 counterpart at all
#
# `version` = "C" applies the map; "D" is already v17b and is passed through.
# -----------------------------------------------------------------------------
harmonise_ids <- function(df, id_col, lut, version, label = "") {

  n_before <- nrow(df)

  if (version == "D") {
    out <- df %>%
      mutate(
        id_harmonised   = as.character(.data[[id_col]]),
        xlate_ambiguous = FALSE,
        xlate_missing   = FALSE
      )
  } else {
    amb <- attr(lut, "ambiguous_from_ids")
    out <- df %>%
      left_join(lut, by = setNames("from_id", id_col),
                relationship = "many-to-one") %>%   # errors if the map is not 1:1
      mutate(
        xlate_missing   = is.na(to_id),
        xlate_ambiguous = .data[[id_col]] %in% amb,
        id_harmonised   = if_else(xlate_missing,
                                  paste0("v16_", .data[[id_col]]),
                                  as.character(to_id))
      ) %>%
      select(-to_id)

    message(sprintf(
      "[%s] %d rows | %d untranslatable v16 ids (%d rows, kept with surrogate ids) | %d rows on ambiguous splits",
      label, nrow(out), n_distinct(out[[id_col]][out$xlate_missing]),
      sum(out$xlate_missing), sum(out$xlate_ambiguous)))
  }

  stopifnot(nrow(out) == n_before)   # harmonisation must never fan out
  out
}


# -----------------------------------------------------------------------------
# add_obs_key() — a version-INDEPENDENT identifier for one paired observation
# -----------------------------------------------------------------------------
# Do NOT key on the in situ timestamp. For GNSS the timestamp is the drift
# midpoint computed INSIDE the node/reach footprint, and those footprints move
# between SWORD v16 and v17b, so the same physical observation gets a different
# timestamp in each version and the match silently fails. Measured on the node
# GNSS data: of 3,064 rows sharing node+cycle+pass across versions, only 2,207
# (72%) had an identical timestamp.
#
# Use instead, per product:
#   node GNSS   : id_harmonised + cycle_id + pass_id + basename(drift_id)
#   node PT     : pt_serial + pt_time_UTC
#                 (the node-level PT file carries no cycle_id/pass_id; the PT
#                  serial is a physical sensor and is version-independent, so
#                  it is a better key than the node id anyway)
#   reach GNSS  : id_harmonised + cycle_id + pass_id + basename(drift_id)
#   reach PT    : id_harmonised + cycle_id + pass_id + pt_time_UTC
#   node width  : id_harmonised + cycle_id + pass_id + ortho id
#
# basename() strips the version-specific directory prefix that drift_id carries
# ("Munged drifts/reprocessed_2025_09_02/..." vs ".../v17b_reprocessed_...").
# -----------------------------------------------------------------------------
add_obs_key <- function(df, key_cols) {
  missing_cols <- setdiff(key_cols, names(df))
  if (length(missing_cols)) {
    stop("add_obs_key(): missing key column(s): ", paste(missing_cols, collapse = ", "))
  }
  df %>% tidyr::unite("obs_key", all_of(key_cols), sep = "|", remove = FALSE)
}


# -----------------------------------------------------------------------------
# partition_versions() — an EXHAUSTIVE three-way split
# -----------------------------------------------------------------------------
# Adds:
#   obs_bucket  "same" / "C0_only" / "D0_only", at OBSERVATION level
#   id_bucket   "same" / "C0_only" / "D0_only", at NODE/REACH level
#
# These are computed separately because they answer different questions, and
# conflating them is what broke the original tables. They are also computed
# WITHIN each in situ type (`strata`), because every published table is
# stratified by in situ type while the old version_inclusion flag pooled PT and
# GNSS together.
#
# "same" means both versions passed quality filtering at the same overpass.
# By construction: n(same) + n(C0_only) + n(D0_only) == n(total), per version.
# -----------------------------------------------------------------------------
partition_versions <- function(df, key_cols, strata = "insitu_type",
                               value_col = NULL) {

  df <- add_obs_key(df, key_cols)

  # the key must identify at most one row per version, or "same" is ambiguous
  dup <- df %>%
    summarise(n = n(), .by = c(all_of(strata), source, obs_key)) %>%
    filter(n > 1)
  if (nrow(dup)) {
    warning(sprintf(
      "partition_versions(): observation key is not unique -- %d duplicated key(s), max %d rows. Widen key_cols.",
      nrow(dup), max(dup$n)))
  }

  classify <- function(has_c, has_d) {
    case_when(has_c &  has_d ~ "same",
              has_c & !has_d ~ "C0_only",
              !has_c & has_d ~ "D0_only")
  }

  obs_map <- df %>%
    distinct(across(all_of(strata)), obs_key, source) %>%
    summarise(has_c = any(source == VERSION_C),
              has_d = any(source == VERSION_D),
              .by = c(all_of(strata), obs_key)) %>%
    mutate(obs_bucket = classify(has_c, has_d)) %>%
    select(all_of(strata), obs_key, obs_bucket)

  id_map <- df %>%
    distinct(across(all_of(strata)), id_harmonised, source) %>%
    summarise(has_c = any(source == VERSION_C),
              has_d = any(source == VERSION_D),
              .by = c(all_of(strata), id_harmonised)) %>%
    mutate(id_bucket = classify(has_c, has_d)) %>%
    select(all_of(strata), id_harmonised, id_bucket)

  df %>%
    left_join(obs_map, by = c(strata, "obs_key"),       relationship = "many-to-one") %>%
    left_join(id_map,  by = c(strata, "id_harmonised"), relationship = "many-to-one")
}


# -----------------------------------------------------------------------------
# summarise_errors() — the single canonical metric definition
# -----------------------------------------------------------------------------
# One function for every table in 4.1/4.2/4.3, so that a percentile and an MAE
# in the same row can never be computed from different columns again. (In the
# original 4.3, Table S10's MAE used residuals_nobias while the percentiles and
# RMSE in the same row used residuals.)
#
#   scale   100 for m -> cm (WSE), 1e5 for m/m -> cm/km (slope), 1 for width (m)
#   n            number of non-missing values of `value_col`
#   n_unique     number of distinct features among THOSE SAME ROWS
#                (the original counted distinct ids over all rows, including
#                rows whose residual was NA, so n and n_unique described
#                different row sets)
# -----------------------------------------------------------------------------
summarise_errors <- function(df, value_col, by, scale = 100, digits = 1,
                             bias_col = "bias", id_col = "id_harmonised") {

  has_bias <- !is.null(bias_col) && bias_col %in% names(df)

  df %>%
    filter(!is.na(.data[[value_col]])) %>%
    summarise(
      n           = n(),
      n_unique    = n_distinct(.data[[id_col]]),
      error_68ile = round(quantile(abs(.data[[value_col]]) * scale, 0.68, na.rm = TRUE), digits),
      error_50ile = round(quantile(abs(.data[[value_col]]) * scale, 0.50, na.rm = TRUE), digits),
      MAE         = round(mean(abs(.data[[value_col]]) * scale, na.rm = TRUE), digits),
      RMSE        = round(sqrt(mean((.data[[value_col]] * scale)^2, na.rm = TRUE)), digits),
      bias        = if (has_bias) round(median(.data[[bias_col]], na.rm = TRUE) * scale, digits) else NA_real_,
      .by = all_of(by)
    ) %>%
    arrange(across(all_of(by)))
}


# -----------------------------------------------------------------------------
# partition_table() — total + the three buckets, in one table, with a check
# -----------------------------------------------------------------------------
# This is what should be published in place of the old "Same" / "Unique"
# supplementary tables. It is self-reconciling: the TOTAL row always equals the
# sum of the three bucket rows, and the function stops if it does not.
# -----------------------------------------------------------------------------
partition_table <- function(df, value_col, by = c("insitu_type", "source"),
                            scale = 100, digits = 1, bias_col = "bias") {

  totals <- df %>%
    summarise_errors(value_col, by = by, scale = scale, digits = digits, bias_col = bias_col) %>%
    mutate(bucket = "total", .before = 1)

  buckets <- df %>%
    summarise_errors(value_col, by = c(by, "obs_bucket"),
                     scale = scale, digits = digits, bias_col = bias_col) %>%
    rename(bucket = obs_bucket)

  out <- bind_rows(totals, buckets) %>%
    mutate(bucket = factor(bucket, levels = c("total", "same", "C0_only", "D0_only"))) %>%
    arrange(across(all_of(by)), bucket)

  check <- out %>%
    summarise(total = sum(n[bucket == "total"]),
              parts = sum(n[bucket != "total"]),
              .by = all_of(by))

  if (!all(check$total == check$parts)) {
    print(check)
    stop("partition_table(): buckets do not sum to the total. The observation key is wrong.")
  }
  message("[partition_table] reconciliation OK: same + C0_only + D0_only == total for all strata")
  out
}


# -----------------------------------------------------------------------------
# version_change() — the "X% more data in D0" claims, computed one way
# -----------------------------------------------------------------------------
# The manuscript currently computes these three different ways: the reach WSE
# claim sums PT and GNSS, the slope claim uses GNSS only despite saying
# "aggregated", and the node claim matches neither. This returns both the
# observation count and the pooled n_distinct feature count, always the same way.
# -----------------------------------------------------------------------------
version_change <- function(df, value_col, id_col = "id_harmonised") {
  df %>%
    filter(!is.na(.data[[value_col]])) %>%
    summarise(n_obs = n(), n_unique = n_distinct(.data[[id_col]]), .by = source) %>%
    tidyr::pivot_wider(names_from = source, values_from = c(n_obs, n_unique)) %>%
    mutate(
      pct_change_obs    = round(100 * (.data[[paste0("n_obs_",    VERSION_D)]] - .data[[paste0("n_obs_",    VERSION_C)]]) /
                                        .data[[paste0("n_obs_",    VERSION_C)]], 1),
      pct_change_unique = round(100 * (.data[[paste0("n_unique_", VERSION_D)]] - .data[[paste0("n_unique_", VERSION_C)]]) /
                                        .data[[paste0("n_unique_", VERSION_C)]], 1)
    )
}


# -----------------------------------------------------------------------------
# Plot constants, shared by all three scripts
# -----------------------------------------------------------------------------
river_levels  <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
river_labels  <- c("Coleen", "Sheenjek", "Chandalar", "Porcupine",
                   "Single-channel Yukon", "Braided Yukon")
river_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

version_colours <- c(PIC0 = "#E69F00", PGD0 = "#0072B2")
insitu_colours  <- c(PT   = "#009E73", GNSS = "#CC79A7")

merge_porcupine <- function(df) {
  df %>% mutate(river = if_else(river %in% c("lowerPR", "upperPR"), "PR", river))
}
