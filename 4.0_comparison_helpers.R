# =============================================================================
# Shared helpers for SWOT version C0 / D0 comparison scripts (4.1, 4.2, 4.3)
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
# deterministically and flag the affected rows.
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

  stopifnot(nrow(out) == n_before)   # harmonization must never fan out
  out
}


# -----------------------------------------------------------------------------
# add_obs_key() — a version-INDEPENDENT identifier for one paired observation
# -----------------------------------------------------------------------------
# Use, per product:
#   node GNSS   : id_harmonised + cycle_id + pass_id + basename(drift_id)
#   node PT     : id_harmonised + pt_serial + pt_time_UTC
#                 (the node-level PT file carries no cycle_id/pass_id; pt_serial
#                  is a physical sensor and so is version-independent)
#   reach GNSS  : id_harmonised + cycle_id + pass_id + basename(drift_id)
#   reach PT    : id_harmonised + cycle_id + pass_id + pt_time_UTC
#   node width  : id_harmonised + cycle_id + pass_id + ortho id
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
# tables.
#
# "unmappable" takes precedence over everything else and marks version-C rows
# whose v16 id has no v17b counterpart (xlate_missing).
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
# summarise_errors() — metric definition
# -----------------------------------------------------------------------------
#
#   scale   100 for m -> cm (WSE), 1e5 for m/m -> cm/km (slope), 1 for width (m)
#   n            number of non-missing values of `value_col`
#   n_unique     number of distinct features among THOSE SAME ROWS
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
# partition_table() — total + the buckets, in one table
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
# Plotting constants
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
