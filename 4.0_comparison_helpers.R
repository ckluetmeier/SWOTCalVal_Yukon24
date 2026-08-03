# =============================================================================
# Shared helpers for SWOT version C0 / D0 comparison scripts (4.1, 4.2, 4.3)
# -----------------------------------------------------------------------------
# Written for R >= 4.4 and dplyr >= 1.1 (uses .by=, relationship=, unmatched=).
# Verified against dplyr 1.2.1 / R 4.6.0 idiom.
#
# WHY THIS FILE EXISTS
# --------------------
# The original 4.1/4.2/4.3 each re-implemented the same three operations
# (SWORD v16 -> v17b harmonisation, version partitioning, error summaries)
# with small differences, and two of those operations were wrong:
#
#   1. left_join() to the SWORD translator FANNED OUT rows, because the
#      translator has one row per v17 id and some v16 ids split into several
#      v17 ids. This silently inflated every version-C count
#      (+36 GNSS node obs, +9 PT node obs).
#
#   2. The "same subset" was keyed on (id, insitu_time_utc) while the
#      version_inclusion flag was keyed on id alone, so the two were not
#      complementary and Same + Unique != Total. The in situ timestamp is
#      NOT stable across SWORD versions (node/reach footprints move, so the
#      GNSS drift midpoint inside them moves too).
#
# REVISION: the "same" bucket must be SYMMETRIC
# ---------------------------------------------
# A first version of this file bucketed every row and then let
# summarise_errors() drop rows with a missing residual. That made the "same"
# bucket asymmetric between versions, which is wrong by definition -- a matched
# subset must contain the same number of observations on both sides. Two causes,
# both of which showed up in the node PT table (269/57 vs 263/52):
#
#   a) 6 pairs had a usable residual in one version and NA in the other. The
#      pair was bucketed as "same" but only contributed to one version's n.
#      FIX: bucket only rows that already carry a usable residual, so a key
#      where one side is NA becomes C0_only / D0_only. That is also the correct
#      reading -- one version produced a usable comparison and the other did not.
#      partition_versions() now takes `value_col` for this reason.
#
#   b) 10 of 270 PT pairs had the two versions disagreeing on the harmonised
#      node id, because the node-PT key was (pt_serial, pt_time_UTC) and did
#      not mention the node. n matched but n_unique did not.
#      FIX: include id_harmonised in every observation key, node PT included.
#      A pair that straddles two different nodes is not the same node
#      observation and is now correctly bucketed as version-unique.
#
# partition_table() asserts symmetry and stops if it is ever violated again.
#
# REVISION 3: the "unmappable" bucket
# -----------------------------------
# 39 v16 nodes (53 GNSS node observations) have no v17b counterpart in the
# translator at all. They cannot pair with D0 by construction, so they were
# falling into C0_only -- which is the row the manuscript uses to argue that
# nodes dropped by D0's quality filter were poor quality. That inference does
# not hold for these 39: SWORD v17b simply re-noded those reaches, and for 32
# of the 39 the physically-nearest v17b node IS present in the D0 dataset. D0
# observed that stretch of river, under a different node id.
#
# They also carry anomalously high error (|68%ile| 56.9 cm against 16.2 cm for
# C0 overall), so folding them into C0_only moved that statistic from 33.0 to
# 33.7 cm -- a small effect, but in the direction that flatters the argument.
#
# They are now reported as their own bucket: still inside the C0 total (so
# Table 2 reconciles unchanged at 6,048), but visible and excluded from the
# C0_only claim. Rows are identified by xlate_missing == TRUE.
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
# version-C observation once per v17 child, which inflates the C0 counts.
#
# In the Yukon domain this affects 25 v16 nodes (each splitting into exactly 2
# v17 nodes) and no reaches. The translator provides no field that can pick a
# winner (boundary_percent is 0 for all of them), so we take the lowest v17 id
# deterministically and flag the affected rows so they stay auditable.
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
# ALWAYS include id_harmonised, so that both members of a matched pair refer to
# the same node/reach and n_unique is symmetric by construction.
#
# Use, per product:
#   node GNSS   : id_harmonised + cycle_id + pass_id + basename(drift_id)
#   node PT     : id_harmonised + pt_serial + pt_time_UTC
#                 (the node-level PT file carries no cycle_id/pass_id; pt_serial
#                  is a physical sensor and so is version-independent)
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
  if (!"id_harmonised" %in% key_cols) {
    warning("add_obs_key(): id_harmonised is not in key_cols. n_unique will not be ",
            "symmetric between versions in the 'same' bucket.")
  }
  df %>% tidyr::unite("obs_key", all_of(key_cols), sep = "|", remove = FALSE)
}


# -----------------------------------------------------------------------------
# partition_versions() — an EXHAUSTIVE, SYMMETRIC three-way split
# -----------------------------------------------------------------------------
# Adds:
#   obs_bucket  "same" / "C0_only" / "D0_only" / "unmappable", OBSERVATION level
#   id_bucket   "same" / "C0_only" / "D0_only" / "unmappable", NODE/REACH level
# Rows with a missing `value_col` get NA in both and are excluded from all
# tables (summarise_errors() drops them anyway).
#
# "unmappable" takes precedence over everything else and marks version-C rows
# whose v16 id has no v17b counterpart (xlate_missing). They cannot pair with
# D0 for a reason that has nothing to do with quality filtering, so they must
# not be counted as evidence in the C0_only row. They remain inside the C0
# total, so reconciliation is unaffected.
#
# `value_col` is REQUIRED: the partition is defined with respect to one metric.
# Bucketing is done on rows that already carry a usable value, so a key where
# one version is NA is correctly labelled version-unique rather than "same".
# Without this, the "same" bucket is asymmetric -- see the header note.
#
# obs_bucket and id_bucket are computed separately because they answer different
# questions, and conflating them is what broke the original tables. Both are
# computed WITHIN each in situ type (`strata`), because every published table is
# stratified by in situ type while the old version_inclusion flag pooled PT and
# GNSS together.
#
# "same" means both versions passed quality filtering, produced a usable value,
# and did so at the same overpass on the same feature.
# By construction: n(same) + n(C0_only) + n(D0_only) == n(total), per version,
# and n(same) is identical for the two versions.
# -----------------------------------------------------------------------------
partition_versions <- function(df, key_cols, value_col, strata = "insitu_type") {

  if (missing(value_col)) {
    stop("partition_versions(): `value_col` is required -- the partition is ",
         "defined with respect to one metric (e.g. \"residuals_nobias\").")
  }

  df <- add_obs_key(df, key_cols)

  # Bucket on rows that carry a usable value; NA rows cannot be part of a pair.
  usable <- df %>% filter(!is.na(.data[[value_col]]))

  dup <- usable %>%
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

  obs_map <- usable %>%
    distinct(across(all_of(strata)), obs_key, source) %>%
    summarise(has_c = any(source == VERSION_C),
              has_d = any(source == VERSION_D),
              .by = c(all_of(strata), obs_key)) %>%
    mutate(obs_bucket = classify(has_c, has_d)) %>%
    select(all_of(strata), obs_key, obs_bucket)

  id_map <- usable %>%
    distinct(across(all_of(strata)), id_harmonised, source) %>%
    summarise(has_c = any(source == VERSION_C),
              has_d = any(source == VERSION_D),
              .by = c(all_of(strata), id_harmonised)) %>%
    mutate(id_bucket = classify(has_c, has_d)) %>%
    select(all_of(strata), id_harmonised, id_bucket)

  out <- df %>%
    left_join(obs_map, by = c(strata, "obs_key"),       relationship = "many-to-one") %>%
    left_join(id_map,  by = c(strata, "id_harmonised"), relationship = "many-to-one") %>%
    # unmappable overrides: no v17b counterpart, so pairing is impossible for a
    # reason unrelated to quality filtering
    mutate(
      obs_bucket = if_else(xlate_missing, "unmappable", obs_bucket),
      id_bucket  = if_else(xlate_missing, "unmappable", id_bucket)
    ) %>%
    mutate(
      obs_bucket = if_else(is.na(.data[[value_col]]), NA_character_, obs_bucket),
      id_bucket  = if_else(is.na(.data[[value_col]]), NA_character_, id_bucket)
    )

  n_unmap <- sum(out$obs_bucket == "unmappable", na.rm = TRUE)
  if (n_unmap > 0) {
    message(sprintf(
      "[partition_versions] %d observations on %d features have no v17b counterpart -> bucket 'unmappable' (reported separately from C0_only)",
      n_unmap, n_distinct(out$id_harmonised[which(out$obs_bucket == "unmappable")])))
  }

  attr(out, "partition_value_col") <- value_col
  out
}


# Bucket display order, shared by partition_table() and 4.3's width tables
BUCKET_LEVELS <- c("total", "same", "C0_only", "D0_only", "unmappable")


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
# partition_table() — total + the buckets, in one table, with two checks
# -----------------------------------------------------------------------------
# This is what should be published in place of the old "Same" / "Unique"
# supplementary tables. It enforces:
#   RECONCILIATION  the TOTAL row equals the sum of same + C0_only + D0_only
#                   + unmappable
#   SYMMETRY        the "same" bucket has identical n and n_unique in both
#                   versions -- the property that was broken for PT
# and stops if either fails.
# -----------------------------------------------------------------------------
partition_table <- function(df, value_col, by = c("insitu_type", "source"),
                            scale = 100, digits = 1, bias_col = "bias") {

  pv <- attr(df, "partition_value_col")
  if (!is.null(pv) && !identical(pv, value_col)) {
    stop("partition_table(): this data was partitioned on '", pv,
         "' but you are tabulating '", value_col,
         "'. Re-run partition_versions() with value_col = \"", value_col, "\".")
  }

  totals <- df %>%
    summarise_errors(value_col, by = by, scale = scale, digits = digits, bias_col = bias_col) %>%
    mutate(bucket = "total", .before = 1)

  buckets <- df %>%
    filter(!is.na(obs_bucket)) %>%
    summarise_errors(value_col, by = c(by, "obs_bucket"),
                     scale = scale, digits = digits, bias_col = bias_col) %>%
    rename(bucket = obs_bucket)

  out <- bind_rows(totals, buckets) %>%
    mutate(bucket = factor(bucket, levels = BUCKET_LEVELS)) %>%
    arrange(across(all_of(by)), bucket)

  if (any(is.na(out$bucket))) {
    stop("partition_table(): unrecognised bucket label -- BUCKET_LEVELS is out of date.")
  }

  # --- check 1: buckets sum to the total ---------------------------------
  recon <- out %>%
    summarise(total = sum(n[bucket == "total"]),
              parts = sum(n[bucket != "total"]),
              .by = all_of(by))
  if (!all(recon$total == recon$parts)) {
    print(recon)
    stop("partition_table(): buckets do not sum to the total. The observation key is wrong.")
  }

  # --- check 2: the matched subset is symmetric between versions ---------
  strata_only <- setdiff(by, "source")
  sym <- out %>%
    filter(bucket == "same") %>%
    summarise(n_versions        = n(),
              distinct_n        = n_distinct(n),
              distinct_n_unique = n_distinct(n_unique),
              .by = all_of(strata_only))
  if (any(sym$distinct_n > 1) || any(sym$distinct_n_unique > 1)) {
    print(out %>% filter(bucket == "same"))
    stop("partition_table(): the 'same' bucket is not symmetric between versions. ",
         "Either id_harmonised is missing from key_cols, or partition_versions() ",
         "was given the wrong value_col.")
  }

  message("[partition_table] OK: buckets reconcile, and the matched subset is symmetric")
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
