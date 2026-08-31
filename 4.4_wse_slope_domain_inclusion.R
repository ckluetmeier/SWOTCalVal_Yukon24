# =============================================================================
# SWOT WSE & Slope Domain Inclusion
# -----------------------------------------------------------------------------
# Determines, for each node/reach in the Yukon River domain, whether it was
# observed by Version C, Version D, both, or neither SWOT product, and writes
# the result as shapefiles for mapping. Also produces per-river domain
# statistics and the in situ summary shapefiles + example timeseries used in
# Figures 2 and 3.
#
# SWOT processing versions:
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
#
# Version inclusion code (per node / reach, in SWORD v17b space):
#   -1 = observed only in Version C (RiverSP v16 / PIC0)
#    0 = observed in both versions
#    1 = observed only in Version D (v17b / PGD0)
#    2 = in the in situ YR domain but not observed in either SWOT version
#
# Contains:
#   - Tables: 1, S4
#   - Figures: 2, 3
#   - Shapefiles: node WSE inclusion, reach WSE inclusion, reach slope
#                 inclusion, PT summary, GNSS summary
#
# =============================================================================

library(sf)
library(tidyverse)
library(lubridate)

source("/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

WSE_BASE       <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"
SWORD_NODES    <- "/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_nodes_hb81_v17b.shp"
SWORD_REACHES  <- "/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp"
INCLUSION_OUT  <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/inclusion_maps"

PT_NODE_DIR    <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b"
PT_REACH_DIR   <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_reach/SWORD_v17b"
GNSS_NODE_CSV  <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_node_wses.csv"
GNSS_REACH_CSV <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GNSS/_processed_data/reprocessed_2025_09_02/SWORD_v17b/YR_drift_reach_wse_slope.csv"

# Section 6 of this script (Figure 2) writes a merged PT node file back into
# PT_NODE_DIR. Exclude it from the domain glob so a second run does not read
# its own output.
PT_MERGED_NODE_FILE <- "flyby_SWOTCalVal_YR_PT_L1_v17b.csv"

# --- filters: must match 4.1 and 4.2 -----------------------------------------
DARK_FRAC_MAX <- 0.5   # upstream scripts 1.2-2.2 only filter at <= 0.8,
                       # so this filter is NOT redundant and does remove rows

# Reaches shorter than 9 km, listed by ID in each SWORD version because
# p_length is not carried through. (Manually checked in QGIS.) Character, so
# the comparison happens on the same type the ids are carried in here.
SHORT_REACHES_V16  <- c("81260300061", "81270500131", "81270500141")
SHORT_REACHES_V17B <- c("81260300181", "81270500021", "81270500031")
APPLY_SHORT_REACH_EXCLUSION <- TRUE

# --- metrics the partitions are defined against: must match 4.1 and 4.2 ------
NODE_VALUE  <- "residuals_nobias"
REACH_VALUE <- "residuals_nobias"
SLOPE_VALUE <- "slope_residuals_nobias"

# --- inclusion codes ----------------------------------------------------------
INCL_C_ONLY      <- -1L
INCL_BOTH        <-  0L
INCL_D_ONLY      <-  1L
INCL_DOMAIN_ONLY <-  2L


# =============================================================================
# 0b. LOCAL HELPERS
# =============================================================================

# --- as_id_chr(): one canonical string form for a SWORD id --------------------
# Same helper as 4.3. A 14-digit node id read as a double would render as
# "8.126e+13" under as.character(), which matches nothing. sprintf("%.0f") is
# exact for ids well under 2^53.
as_id_chr <- function(x) {
  if (is.numeric(x)) ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  else               ifelse(is.na(x), NA_character_, trimws(as.character(x)))
}

# --- chr_lut(): put a build_id_lut() map into character space -----------------
# build_id_lut() returns whatever type read_csv() gave the translator (numeric).
# The data ids are coerced to character before harmonise_ids(), so the map has
# to move with them -- including the ambiguous-id attribute.
chr_lut <- function(lut) {
  amb <- attr(lut, "ambiguous_from_ids")
  out <- lut %>% mutate(from_id = as_id_chr(from_id), to_id = as_id_chr(to_id))
  attr(out, "ambiguous_from_ids") <- as_id_chr(amb)
  out
}

# --- collect_domain_ids(): the in situ domain, robust to column casing --------
# PT node toolbox files use Node_ID / Reach_ID; PT reach files and the GNSS
# files use node_id / reach_id. Reading the whole frame invites type conflicts
# across clusters, so only the id column is read, and it is read as character.
collect_domain_ids <- function(paths, candidates, label = "") {
  paths <- unique(paths)
  if (!length(paths)) stop("collect_domain_ids(): no files given for ", label)

  out <- map_dfr(paths, function(f) {
    hdr <- suppressWarnings(
      read_csv(f, n_max = 0, show_col_types = FALSE, progress = FALSE))
    hit <- intersect(candidates, names(hdr))
    if (!length(hit)) {
      warning("[domain ", label, "] no id column (", paste(candidates, collapse = "/"),
              ") in ", basename(f), " -- file contributes nothing to the domain")
      return(tibble(id = character(), src_file = character(), src_col = character()))
    }
    vals <- suppressWarnings(
      read_csv(f, col_select = all_of(hit[1]), col_types = cols(.default = col_character()),
               progress = FALSE))
    tibble(id = as_id_chr(vals[[1]]), src_file = basename(f), src_col = hit[1])
  })

  out <- out %>% filter(!is.na(id), id != "")

  message(sprintf("[domain %s] %d file(s), %d row(s), %d distinct id(s); id column(s) used: %s",
                  label, length(paths), nrow(out), n_distinct(out$id),
                  paste(sort(unique(out$src_col)), collapse = ", ")))
  out
}

# --- report_unmappable(): what the v17b maps cannot show ----------------------
report_unmappable <- function(df, value_col, raw_id_col, label) {
  u <- df %>% filter(!is.na(.data[[value_col]]), xlate_missing)
  if (!nrow(u)) {
    message(sprintf("[unmappable %s] none -- every version-C id has a v17b counterpart", label))
    return(invisible(NULL))
  }
  ids <- sort(unique(as_id_chr(u[[raw_id_col]])))
  message(sprintf(
    "[unmappable %s] %d observation(s) on %d version-C feature(s) have no v17b counterpart.",
    label, nrow(u), length(ids)))
  message(sprintf(
    "[unmappable %s]   These are EXCLUDED from the shapefile: there is no v17b geometry to draw them on.",
    label))
  message(sprintf("[unmappable %s]   v16 ids: %s", label, paste(ids, collapse = ", ")))
  invisible(tibble(v16_id = ids))
}

# --- inclusion_from_partition(): feature-level buckets, pooled over in situ ----
inclusion_from_partition <- function(df, value_col, label = "") {
  usable <- df %>% filter(!is.na(.data[[value_col]]), !xlate_missing)

  out <- usable %>%
    summarise(n_obs_C = sum(source == VERSION_C),
              n_obs_D = sum(source == VERSION_D),
              .by = id_harmonised) %>%
    mutate(version_inclusion = case_when(
      n_obs_C >  0 & n_obs_D >  0 ~ INCL_BOTH,
      n_obs_C >  0 & n_obs_D == 0 ~ INCL_C_ONLY,
      n_obs_C == 0 & n_obs_D >  0 ~ INCL_D_ONLY
    )) %>%
    rename(id = id_harmonised)

  if (any(is.na(out$version_inclusion))) {
    stop("inclusion_from_partition(): a feature has neither a C nor a D ",
         "observation, which cannot happen -- check `source` values against ",
         "VERSION_C / VERSION_D.")
  }
  message(sprintf("[inclusion %s] %d observed feature(s): %d both, %d C-only, %d D-only",
                  label, nrow(out),
                  sum(out$version_inclusion == INCL_BOTH),
                  sum(out$version_inclusion == INCL_C_ONLY),
                  sum(out$version_inclusion == INCL_D_ONLY)))
  out
}

# --- build_inclusion_table(): observed features + the unobserved domain -------
build_inclusion_table <- function(observed, domain_ids, id_name, label = "") {
  domain_ids <- unique(domain_ids[!is.na(domain_ids)])

  observed_not_in_domain <- setdiff(observed$id, domain_ids)
  domain_not_observed    <- setdiff(domain_ids, observed$id)

  if (length(observed_not_in_domain)) {
    message(sprintf(
      "[domain check %s] %d observed feature(s) are NOT in the in situ domain (kept in the map, coded normally): %s",
      label, length(observed_not_in_domain),
      paste(utils::head(observed_not_in_domain, 10), collapse = ", ")))
  }

  out <- tibble(id = union(domain_ids, observed$id)) %>%
    left_join(observed, by = "id") %>%
    mutate(
      n_obs_C = replace_na(n_obs_C, 0L),
      n_obs_D = replace_na(n_obs_D, 0L),
      version_inclusion = replace_na(version_inclusion, INCL_DOMAIN_ONLY)
    ) %>%
    rename(!!id_name := id)

  message(sprintf("[domain check %s] %d domain feature(s) with no SWOT observation -> code 2",
                  label, length(domain_not_observed)))
  message(sprintf("[inclusion %s] final table: %d feature(s)", label, nrow(out)))
  print(count(out, version_inclusion))
  out
}

# --- SHP_CODE_FIELD: the on-disk name of the inclusion code -------------------
SHP_CODE_FIELD <- "vrsn_nc"

# --- write_inclusion_shapefile(): join to SWORD geometry and write ------------
write_inclusion_shapefile <- function(tbl, sword_path, id_name, out_path, label = "") {
  geom <- st_read(sword_path, quiet = TRUE)

  if (!id_name %in% names(geom)) {
    stop("write_inclusion_shapefile(): '", id_name, "' is not a field in ",
         basename(sword_path), ". Fields: ", paste(names(geom), collapse = ", "))
  }

  geom  <- geom %>% mutate(.join_id = as_id_chr(.data[[id_name]]))
  tbl_j <- tbl  %>% mutate(.join_id = as_id_chr(.data[[id_name]])) %>% select(-all_of(id_name))

  matched <- sum(tbl_j$.join_id %in% geom$.join_id)
  if (matched == 0) {
    stop("write_inclusion_shapefile(): no ", label, " id matched ",
         basename(sword_path), ".\n",
         "  table e.g.     ", paste(utils::head(tbl_j$.join_id, 2), collapse = ", "), "\n",
         "  shapefile e.g. ", paste(utils::head(geom$.join_id, 2), collapse = ", "),
         "\nCheck that both are SWORD v17b ids.")
  }
  if (matched < nrow(tbl_j)) {
    unmatched <- setdiff(tbl_j$.join_id, geom$.join_id)
    warning(sprintf(
      "[shapefile %s] %d of %d feature(s) have no geometry in %s and are NOT in the shapefile: %s",
      label, length(unmatched), nrow(tbl_j), basename(sword_path),
      paste(utils::head(unmatched, 10), collapse = ", ")))
  }

  # many-to-one, not one-to-one: the SWORD reach layer can carry duplicate
  # reach_id rows (which is why Table 1 below still needs distinct()), and each
  # of those geometries should receive the same code.
  out <- geom %>%
    left_join(tbl_j, by = ".join_id", relationship = "many-to-one") %>%
    filter(!is.na(version_inclusion)) %>%
    # Overwrite SWORD's numeric id with the canonical string. st_read() brings
    # an 11/14-digit id in as a double and st_write() then puts it back as a
    # Real with 12 decimals ("81250800021.000000000000"), which is unusable as
    # a join key in QGIS.
    mutate(!!id_name := .join_id) %>%
    rename(!!SHP_CODE_FIELD := version_inclusion) %>%
    select(-any_of(".join_id"))

  # --- guard: nothing may exceed the shapefile field-name limit --------------
  fields  <- setdiff(names(out), attr(out, "sf_column"))
  too_long <- fields[nchar(fields) > 10]
  if (length(too_long)) {
    stop("write_inclusion_shapefile(): field name(s) longer than 10 characters: ",
         paste(too_long, collapse = ", "),
         ".\nsf would abbreviate EVERY field in the layer, not just these. ",
         "Shorten them before writing.")
  }

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  st_write(out, out_path, delete_layer = TRUE)

  # --- read back and verify what actually landed on disk ---------------------
  chk <- st_read(out_path, quiet = TRUE)
  renamed <- setdiff(fields, names(chk))
  if (length(renamed)) {
    stop("[shapefile ", label, "] the driver renamed field(s) on write: ",
         paste(renamed, collapse = ", "), ". On-disk names: ",
         paste(names(chk), collapse = ", "))
  }
  if (nrow(chk) != nrow(out)) {
    stop("[shapefile ", label, "] wrote ", nrow(out), " feature(s) but read back ",
         nrow(chk), ".")
  }
  codes <- sort(unique(chk[[SHP_CODE_FIELD]]))
  if (!all(codes %in% c(INCL_C_ONLY, INCL_BOTH, INCL_D_ONLY, INCL_DOMAIN_ONLY))) {
    print(table(chk[[SHP_CODE_FIELD]], useNA = "ifany"))
    stop("[shapefile ", label, "] ", SHP_CODE_FIELD,
         " on disk contains value(s) outside {-1, 0, 1, 2}: ",
         paste(codes, collapse = ", "))
  }

  # --- CSV twin, with full unabbreviated names, as the source of truth -------
  csv_path <- sub("\\.shp$", "_attributes.csv", out_path)
  write_csv(st_drop_geometry(out) %>% rename(version_inclusion = all_of(SHP_CODE_FIELD)),
            csv_path)

  message(sprintf("[shapefile %s] wrote %d feature(s) to %s; %s verified in {-1,0,1,2}: %s",
                  label, nrow(chk), basename(out_path), SHP_CODE_FIELD,
                  paste(codes, collapse = ", ")))
  message(sprintf("[shapefile %s] field names on disk: %s",
                  label, paste(names(chk), collapse = ", ")))
  invisible(out)
}


# =============================================================================
# 0c. TRANSLATORS
# =============================================================================

node_translator  <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                             show_col_types = FALSE)
reach_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_ReachIDs_v17b_vs_v16.csv"),
                             show_col_types = FALSE)

node_lut  <- chr_lut(build_id_lut(node_translator,  "v16_node_id",  "v17_node_id",  "node translator"))
reach_lut <- chr_lut(build_id_lut(reach_translator, "v16_reach_id", "v17_reach_id", "reach translator"))


# =============================================================================
# 1. THE IN SITU YR DOMAIN
# =============================================================================
# A node/reach is "in the domain" if any in situ instrument observed it

pt_node_files <- list.files(PT_NODE_DIR, pattern = "\\.csv$",
                            full.names = TRUE, recursive = TRUE)
pt_node_files <- pt_node_files[basename(pt_node_files) != PT_MERGED_NODE_FILE]

pt_reach_files <- list.files(PT_REACH_DIR, pattern = "^YR_812.*\\.csv$",
                             full.names = TRUE)

domain_nodes <- bind_rows(
  collect_domain_ids(pt_node_files, c("Node_ID", "node_id"), "PT nodes"),
  collect_domain_ids(GNSS_NODE_CSV, c("node_id", "Node_ID"), "GNSS nodes")
) %>% distinct(id) %>% pull(id)

domain_reaches <- bind_rows(
  collect_domain_ids(pt_reach_files, c("reach_id", "Reach_ID"), "PT reaches"),
  collect_domain_ids(GNSS_REACH_CSV, c("reach_id", "Reach_ID"), "GNSS reaches")
) %>% distinct(id) %>% pull(id)

message(sprintf("[domain] %d node(s), %d reach(es) with in situ data",
                length(domain_nodes), length(domain_reaches)))

# The short reaches excluded from 4.1/4.2 are also removed from the reach
# domain, so they cannot come back as code 2.
if (APPLY_SHORT_REACH_EXCLUSION) {
  domain_reaches <- setdiff(domain_reaches, SHORT_REACHES_V17B)
  message(sprintf("[domain] %d reach(es) after removing sub-9 km reaches",
                  length(domain_reaches)))
}


# =============================================================================
# 2. NODE-LEVEL WSE DOMAIN INCLUSION
# =============================================================================

# --- 2a. Read, filter, harmonise (mirrors read_node() in 4.1) ----------------
read_node_wse <- function(path, insitu, version) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D,
           node_id     = as_id_chr(node_id)) %>%
    filter(dark_frac < DARK_FRAC_MAX) %>%
    harmonise_ids("node_id", node_lut, version,
                  label = paste("node", insitu, version))
}

node_PT_vC   <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v16/node_SWOT_PT.csv"),           "PT",   "C")
node_PT_vD   <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v17b/node_SWOT_PT.csv"),          "PT",   "D")
node_GNSS_vC <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv"),  "GNSS", "C")
node_GNSS_vD <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv"), "GNSS", "D")

# --- 2b. Partition on version-independent keys (identical to 4.1) ------------
node_PT <- bind_rows(node_PT_vC, node_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonised", "pt_serial", "pt_time_UTC"),
                     value_col = NODE_VALUE)

node_GNSS <- bind_rows(node_GNSS_vC, node_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "drift_file"),
                     value_col = NODE_VALUE)

node_all <- bind_rows(node_PT, node_GNSS)
attr(node_all, "partition_value_col") <- NODE_VALUE

# Cross-check against tableS6_node_ids in 4.1 -- these numbers should match.
message("[check] node id_bucket by in situ type (compare with tableS6_node_ids in 4.1):")
node_all %>%
  filter(!is.na(id_bucket)) %>%
  summarise(n_nodes = n_distinct(id_harmonised), .by = c(insitu_type, id_bucket)) %>%
  arrange(insitu_type, id_bucket) %>%
  print()

report_unmappable(node_all, NODE_VALUE, "node_id", "node WSE")

# --- 2c. Inclusion table and shapefile ---------------------------------------
node_inclusion <- inclusion_from_partition(node_all, NODE_VALUE, "node WSE") %>%
  build_inclusion_table(domain_nodes, "node_id", "node WSE")

write_inclusion_shapefile(
  node_inclusion, SWORD_NODES, "node_id",
  file.path(INCLUSION_OUT, "all_YR_domain_nodes_subset.shp"), "node WSE")

# =============================================================================
# 3. REACH-LEVEL WSE DOMAIN INCLUSION
# =============================================================================

# --- 3a. Read, filter, harmonise (mirrors read_reach() in 4.1) ---------------
read_reach_wse <- function(path, insitu, version) {
  short <- if (version == "C") SHORT_REACHES_V16 else SHORT_REACHES_V17B
  d <- read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D,
           reach_id    = as_id_chr(reach_id)) %>%
    filter(dark_frac < DARK_FRAC_MAX)
  if (APPLY_SHORT_REACH_EXCLUSION) d <- filter(d, !reach_id %in% short)
  harmonise_ids(d, "reach_id", reach_lut, version,
                label = paste("reach", insitu, version))
}

reach_PT_vC   <- read_reach_wse(file.path(WSE_BASE, "reach/RiverSP_v16/reach_wse_SWOT_PT.csv"),    "PT",   "C")
reach_PT_vD   <- read_reach_wse(file.path(WSE_BASE, "reach/RiverSP_v17b/reach_wse_SWOT_PT.csv"),   "PT",   "D")
reach_GNSS_vC <- read_reach_wse(file.path(WSE_BASE, "reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv"),  "GNSS", "C")
reach_GNSS_vD <- read_reach_wse(file.path(WSE_BASE, "reach/RiverSP_v17b/reach_wse_SWOT_GNSS.csv"), "GNSS", "D")

# --- 3b. Partition (identical to 4.1) ----------------------------------------
reach_PT <- bind_rows(reach_PT_vC, reach_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "pt_time_UTC"),
                     value_col = REACH_VALUE)

reach_GNSS <- bind_rows(reach_GNSS_vC, reach_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "drift_file"),
                     value_col = REACH_VALUE)

reach_all <- bind_rows(reach_PT, reach_GNSS)
attr(reach_all, "partition_value_col") <- REACH_VALUE

message("[check] reach WSE id_bucket by in situ type (compare with tableS6_reach_ids in 4.1):")
reach_all %>%
  filter(!is.na(id_bucket)) %>%
  summarise(n_reaches = n_distinct(id_harmonised), .by = c(insitu_type, id_bucket)) %>%
  arrange(insitu_type, id_bucket) %>%
  print()

report_unmappable(reach_all, REACH_VALUE, "reach_id", "reach WSE")

# --- 3c. Inclusion table and shapefile ---------------------------------------
reach_wse_inclusion <- inclusion_from_partition(reach_all, REACH_VALUE, "reach WSE") %>%
  build_inclusion_table(domain_reaches, "reach_id", "reach WSE")

write_inclusion_shapefile(
  reach_wse_inclusion, SWORD_REACHES, "reach_id",
  file.path(INCLUSION_OUT, "all_YR_domain_reaches_subset.shp"), "reach WSE")

# Colors used in the corresponding plot:
#   #E69F00  (Version C only, -1)
#   #ececec  (both, 0)
#   #0072B2  (Version D only, 1)
#   dashed gray 2 (no SWOT match, 2)


# =============================================================================
# 4. REACH-LEVEL SLOPE DOMAIN INCLUSION
# =============================================================================

# --- 4a. Read, filter, harmonise (mirrors read_slope() in 4.2) ---------------
# 4.2 always excludes the short reaches; the flag here is honoured so the two
# reach maps are built on the same reach set.
read_reach_slope <- function(path, insitu, version) {
  short <- if (version == "C") SHORT_REACHES_V16 else SHORT_REACHES_V17B
  d <- read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D,
           reach_id    = as_id_chr(reach_id)) %>%
    filter(dark_frac < DARK_FRAC_MAX)
  if (APPLY_SHORT_REACH_EXCLUSION) d <- filter(d, !reach_id %in% short)
  d %>%
    merge_porcupine() %>%
    harmonise_ids("reach_id", reach_lut, version,
                  label = paste("slope", insitu, version))
}

slope_PT_vC   <- read_reach_slope(file.path(WSE_BASE, "reach/RiverSP_v16/reach_slope_SWOT_PT.csv"),    "PT",   "C")
slope_PT_vD   <- read_reach_slope(file.path(WSE_BASE, "reach/RiverSP_v17b/reach_slope_SWOT_PT.csv"),   "PT",   "D")
slope_GNSS_vC <- read_reach_slope(file.path(WSE_BASE, "reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv"),  "GNSS", "C")
slope_GNSS_vD <- read_reach_slope(file.path(WSE_BASE, "reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv"), "GNSS", "D")

# --- 4b. Partition (identical to 4.2) ----------------------------------------
slope_PT <- bind_rows(slope_PT_vC, slope_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "pt_time_UTC"),
                     value_col = SLOPE_VALUE)

slope_GNSS <- bind_rows(slope_GNSS_vC, slope_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "drift_file"),
                     value_col = SLOPE_VALUE)

slope_all <- bind_rows(slope_PT, slope_GNSS)
attr(slope_all, "partition_value_col") <- SLOPE_VALUE

message("[check] slope id_bucket by in situ type (compare with tableS8_ids in 4.2):")
slope_all %>%
  filter(!is.na(id_bucket)) %>%
  summarise(n_reaches = n_distinct(id_harmonised), .by = c(insitu_type, id_bucket)) %>%
  arrange(insitu_type, id_bucket) %>%
  print()

report_unmappable(slope_all, SLOPE_VALUE, "reach_id", "reach slope")

# --- 4c. Inclusion table and shapefile ---------------------------------------
slope_inclusion <- inclusion_from_partition(slope_all, SLOPE_VALUE, "reach slope") %>%
  build_inclusion_table(domain_reaches, "reach_id", "reach slope")

write_inclusion_shapefile(
  slope_inclusion, SWORD_REACHES, "reach_id",
  file.path(INCLUSION_OUT, "all_YR_domain_reaches_slope_subset.shp"), "reach slope")

# Colors used in the corresponding plot:
#   #e97132  (Version C only, -1)
#   #ececec  (both, 0)
#   #00008b  (Version D only, 1)
#   dashed gray 2 (no SWOT match, 2)


# =============================================================================
# 5. INCLUSION SUMMARY (all three products side by side)
# =============================================================================

inclusion_summary <- bind_rows(
  node_inclusion       %>% count(version_inclusion) %>% mutate(product = "node WSE",    .before = 1),
  reach_wse_inclusion  %>% count(version_inclusion) %>% mutate(product = "reach WSE",   .before = 1),
  slope_inclusion      %>% count(version_inclusion) %>% mutate(product = "reach slope", .before = 1)
) %>%
  mutate(label = case_when(
    version_inclusion == INCL_C_ONLY      ~ "C0 only",
    version_inclusion == INCL_BOTH        ~ "both",
    version_inclusion == INCL_D_ONLY      ~ "D0 only",
    version_inclusion == INCL_DOMAIN_ONLY ~ "in domain, no SWOT"
  )) %>%
  pivot_wider(id_cols = product, names_from = label, values_from = n, values_fill = 0L)
print(inclusion_summary)

# =============================================================================
# 6. RIVER DOMAIN STATISTICS (TABLE 1)
# =============================================================================

YR_domain <- tibble(reach_id = domain_reaches)

sword_sf <- st_read(SWORD_REACHES)

# Subset SWORD to the YR domain reaches
sword_subset <- sword_sf %>%
  mutate(reach_id_chr = as_id_chr(reach_id)) %>%
  filter(reach_id_chr %in% YR_domain$reach_id) %>%
  distinct(reach_id_chr, .keep_all = TRUE)

if (nrow(sword_subset) == 0) {
  stop("Table 1: no SWORD reach matched the in situ domain. Check as_id_chr() ",
       "against the reach_id field type in ", basename(SWORD_REACHES), ".")
}

cat("Total km of river:", sum(sword_subset$reach_len) / 1000, "km\n")

# Name each river
sword_subset <- sword_subset %>%
  mutate(
    river_code = substr(reach_id_chr, 1, 6),
    river = case_when(
      reach_id_chr %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id_chr %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                          "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
      river_code == "812701" ~ "lowerYR",  # until the Circle bifurcation
      river_code == "812509" ~ "lowerYR",  # past the PR confluence
      river_code == "812705" ~ "upperYR",  # Circle up
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "PR",
      river_code == "812605" ~ "PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

# Per-river totals from SWORD
river_stats <- sword_subset %>%
  st_drop_geometry() %>%
  group_by(river) %>%
  summarise(
    total_km     = sum(reach_len, na.rm = TRUE) / 1000,
    median_width = median(width, na.rm = TRUE),
    median_slope = median(slope, na.rm = TRUE) * 100,
    n_reaches    = n_distinct(reach_id_chr),
    .groups      = "drop"
  )
print(river_stats)


# -----------------------------------------------------------------------------
# 6b. GNSS per-river drift count
# -----------------------------------------------------------------------------

GNSS_df <- read_csv(GNSS_REACH_CSV, show_col_types = FALSE)

# To get node counts
# GNSS_df <- GNSS_df %>%
#   mutate(node_id = as.character(node_id),
#          reach_id = paste0(substr(node_id, 1, 10), substr(node_id, 14, 14)))

# Tag each GNSS reach with its river (same case_when as above)
GNSS_df <- GNSS_df %>%
  mutate(
    reach_id_chr = as_id_chr(reach_id),
    river_code   = substr(reach_id_chr, 1, 6),
    river = case_when(
      reach_id_chr %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id_chr %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                          "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
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

GNSS_stats <- GNSS_df %>%
  group_by(river) %>%
  summarise(num_drift_reaches = sum(!is.na(reach_id)), .groups = "drop")
print(GNSS_stats)


# =============================================================================
# 7. FIGURE 2 — PT & GNSS SUMMARY SHAPEFILES
# =============================================================================
# Unchanged from the previous version of this script.


# -----------------------------------------------------------------------------
# 7a. PT summary stats (per PT serial): obs count, deployment days, location
# -----------------------------------------------------------------------------

# Merge all PT node CSVs
base_dir <- PT_NODE_DIR

csv_files <- list.files(path = base_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
csv_files <- csv_files[basename(csv_files) != PT_MERGED_NODE_FILE]

combined_PT_df <- map_dfr(csv_files, read_csv, show_col_types = FALSE)

# Save the merged PT node file
write_csv(combined_PT_df, file.path(base_dir, PT_MERGED_NODE_FILE))

# SWOT node timeseries (RiverSP v17b / PGD0), filtered to the PT node set
SWOT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node/RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv"
)

SWOT_df_filtered <- SWOT_df %>%
  semi_join(combined_PT_df, by = c("node_id" = "Node_ID"))

# time_tai is seconds since 2000-01-01, offset 37 seconds from UTC
tai_epoch      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
tai_utc_offset <- 37  # TAI-UTC offset in seconds
SWOT_df_filtered$time_utc <- tai_epoch + SWOT_df_filtered$time_tai - tai_utc_offset

# Match SWOT to PT in time (±7.5 min window)
time_matched_SWOT_PT <- combined_PT_df %>%
  rowwise() %>%
  mutate(
    closest_match = list(
      SWOT_df_filtered %>%
        filter(abs(difftime(pt_time_UTC, time_utc, units = "mins")) <= 7.5)
    )
  ) %>%
  unnest(closest_match) %>%
  dplyr::select(everything())

# Match SWOT to PT in space (node level)
time_space_matched_SWOT_PT <- time_matched_SWOT_PT %>%
  filter(Node_ID == node_id)

# Drop any duplicated PT observations
time_space_matched_SWOT_PT <- time_space_matched_SWOT_PT[
  !duplicated(time_space_matched_SWOT_PT[c("pt_wse_m", "pt_time_UTC", "pt_serial")]), ]

# Per-PT summary stats
PT_summary_stats <- time_space_matched_SWOT_PT %>%
  mutate(
    pt_install_UTC   = ymd_hms(pt_install_UTC),
    pt_uninstall_UTC = ymd_hms(pt_uninstall_UTC)) %>%
  group_by(pt_serial, pt_install_UTC, pt_uninstall_UTC) %>%
  summarise(
    n_obs = n(),
    obs_days = as.numeric(
      difftime(first(pt_uninstall_UTC), first(pt_install_UTC), units = "days")),
    avg_pt_lat = mean(pt_lat, na.rm = TRUE),
    avg_pt_lon = mean(pt_lon, na.rm = TRUE),
    Node_ID    = first(Node_ID),
    Reach_ID   = first(Reach_ID),
    .groups = "drop") %>%
  group_by(pt_serial) %>%
  summarise(
    n_obs = sum(n_obs, na.rm = TRUE),
    obs_days = sum(obs_days, na.rm = TRUE),
    pt_install_UTC   = min(pt_install_UTC, na.rm = TRUE),
    pt_uninstall_UTC = max(pt_uninstall_UTC, na.rm = TRUE),
    avg_pt_lat = mean(avg_pt_lat, na.rm = TRUE),
    avg_pt_lon = mean(avg_pt_lon, na.rm = TRUE),
    Node_ID    = first(Node_ID),
    Reach_ID   = first(Reach_ID),
    .groups = "drop")

# Add river name labels to data frame
PT_summary_stats <- PT_summary_stats %>%
  mutate(
    Reach_ID_chr = as_id_chr(Reach_ID),
    river_code   = substr(Reach_ID_chr, 1, 6),
    river = case_when(
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      Reach_ID_chr %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      Reach_ID_chr %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                          "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
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

river_order <- c("upperYR", "lowerYR", "PR", "CD", "SJ", "CL")

PT_summary_stats <- PT_summary_stats %>%
  mutate(river = factor(river, levels = river_order)) %>%
  arrange(river)

PT_summary_by_river <- PT_summary_stats %>%
  group_by(river) %>%
  summarise(
    earliest_pt_install_UTC   = min(pt_install_UTC, na.rm = TRUE),
    latest_pt_install_UTC     = max(pt_install_UTC, na.rm = TRUE),
    earliest_pt_uninstall_UTC = min(pt_uninstall_UTC, na.rm = TRUE),
    latest_pt_uninstall_UTC   = max(pt_uninstall_UTC, na.rm = TRUE),
    mean_obs_days   = mean(obs_days, na.rm = TRUE),
    median_obs_days = median(obs_days, na.rm = TRUE),
    mean_n_obs      = mean(n_obs, na.rm = TRUE),
    median_n_obs    = median(n_obs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(river)

# NOTE PT 2159244 ERRONEOUSLY HAS 7/31 LISTED AS FIRST INSTALL TIME
# THIS IS BECAUSE FIRST TIME IS PARSED AS NA
# should be 2024-07-08T00:00:00.000000Z to 2024-07-26T18:50:00.000000Z
# and then 2024-07-31T00:25:00.000000Z to 2024-08-21T15:35:00.000000Z
PT_summary_stats$obs_days[abs(PT_summary_stats$obs_days - 21.63194) < 1e-5] <- 40.416667
PT_summary_stats$obs_days <- as.numeric(PT_summary_stats$obs_days)


# Convert to sf (WGS84) and save shapefile + CSV
PT_summary_sf <- st_as_sf(
  PT_summary_stats,
  coords = c("avg_pt_lon", "avg_pt_lat"),
  crs    = 4326
)

st_write(
  PT_summary_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/PT_summary_stats.shp",
  delete_dsn = TRUE
)

write.csv(
  PT_summary_stats,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/PT_summary_stats.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 7b. GNSS summary stats per reach and node
# -----------------------------------------------------------------------------
# NOT FIXED (left alone deliberately -- see note (h) in the header):
# "n_observations" is 14 characters, so these two layers still get their WHOLE
# attribute table abbreviated by sf, exactly as the inclusion shapefiles did.
# Renaming it to something <= 10 characters (e.g. n_obs) fixes it, but it also
# changes the on-disk field name, which the Figure 2 QGIS styling is keyed to.
# The same applies to PT_summary_sf below (pt_install_UTC, pt_uninstall_UTC,
# avg_pt_lat, avg_pt_lon are all over the limit).

# REACH
# ---------------------------

GNSS_df <- read_csv(GNSS_REACH_CSV, show_col_types = FALSE)

GNSS_summary_stats <- GNSS_df %>%
  group_by(reach_id) %>%
  summarise(n_observations = n(), .groups = "drop")

# Attach to SWORD geometry, keep only reaches with GNSS observations
sword_sf <- st_read(SWORD_REACHES)

GNSS_summary_stats_sf <- sword_sf %>%
  mutate(.join_id = as_id_chr(reach_id)) %>%
  left_join(GNSS_summary_stats %>% mutate(.join_id = as_id_chr(reach_id)) %>% select(-reach_id),
            by = ".join_id", relationship = "many-to-one") %>%
  filter(!is.na(n_observations)) %>%
  select(-any_of(".join_id"))

st_write(
  GNSS_summary_stats_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/GNSS_summary_stats_sf.shp",
  delete_layer = TRUE
)

# NODE
# ---------------------------

GNSS_df <- read_csv(GNSS_NODE_CSV, show_col_types = FALSE)

GNSS_summary_stats <- GNSS_df %>%
  group_by(node_id) %>%
  summarise(n_observations = n(), .groups = "drop")

# Attach to SWORD geometry, keep only nodes with GNSS observations
sword_sf <- st_read(SWORD_NODES)

GNSS_summary_stats_sf <- sword_sf %>%
  mutate(.join_id = as_id_chr(node_id)) %>%
  left_join(GNSS_summary_stats %>% mutate(.join_id = as_id_chr(node_id)) %>% select(-node_id),
            by = ".join_id", relationship = "many-to-one") %>%
  filter(!is.na(n_observations)) %>%
  select(-any_of(".join_id"))

st_write(
  GNSS_summary_stats_sf,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/2_domain_map/data/GNSS_summary_stats_node_sf.shp",
  delete_layer = TRUE
)


# =============================================================================
# 8. Table S4 — GNSS METADATA
# =============================================================================
# Unchanged from the previous version of this script.

GNSS_df <- read_csv(GNSS_NODE_CSV, show_col_types = FALSE) %>%
  filter(time_UTC > as.POSIXct("2024-01-01 00:00:00", tz = "UTC"))

sword_sf <- st_read(SWORD_NODES)

GNSS_sf <- sword_sf %>%
  mutate(.join_id = as_id_chr(node_id)) %>%
  right_join(GNSS_df %>% mutate(.join_id = as_id_chr(node_id)) %>% select(-node_id),
             by = ".join_id")

# Add river name labels to data frame
GNSS_sf <- GNSS_sf %>%
  mutate(
    reach_id_chr = as_id_chr(reach_id),
    river_code   = substr(reach_id_chr, 1, 6),
    river = case_when(
      # SWORD v16 SJ reaches: "81260300061", "81260300231", "81260300241", "81260300251"
      # SWORD v17b SJ reaches: "81260300181", "81260300191", "81260300201", "81260300211"
      reach_id_chr %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ",
      reach_id_chr %in% c("81270100111", "81270100121", "81270100131", "81270100141",
                          "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
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

GNSS_summary_stats <- GNSS_sf %>%
  group_by(drift_id) %>%
  summarise(
    river_list       = paste(unique(river), collapse = ", "),
    start_time       = min(time_UTC, na.rm = TRUE),
    end_time         = max(time_UTC, na.rm = TRUE),
    survey_length_km = round(sum(node_len) / 1000, 2),
    reach_list       = paste(unique(reach_id_chr), collapse = ", "),
    .groups          = "drop"
  ) %>%
  mutate(
    survey_time_hours = round(as.numeric(difftime(end_time, start_time, units = "hours")), 2)
  ) %>%
  arrange(start_time) %>%
  filter(survey_time_hours > 0.000000) %>%
  st_drop_geometry() %>%
  select(-drift_id) %>%
  relocate(survey_time_hours, .after = end_time)

# total km of GNSS data collected
sum(GNSS_summary_stats$survey_length_km, na.rm = TRUE)

# mean survey time in hours
mean(GNSS_summary_stats$survey_time_hours, na.rm = TRUE)
max(GNSS_summary_stats$survey_time_hours, na.rm = TRUE)

# unique number of days we have GNSS data from
map2(
  as.Date(GNSS_summary_stats$start_time),
  as.Date(GNSS_summary_stats$end_time),
  ~ seq(.x, .y, by = "day")
) %>%
  unlist() %>%
  as.Date(origin = "1970-01-01") %>%
  n_distinct()

write.csv(
  GNSS_summary_stats,
  file      = "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/GNSS_summary_stats.csv",
  row.names = FALSE
)


# =============================================================================
# 9. FIGURE 3 — EXAMPLE WSE & WIDTH TIMESERIES
# =============================================================================
# Unchanged from the previous version of this script.


# -----------------------------------------------------------------------------
# 9a. GNSS vs. SWOT longitudinal WSE (Porcupine River, 2024-08-20)
# -----------------------------------------------------------------------------

GNSS_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/3_SWOT_examples/data/GNSS_2024-08-20.csv"
)

# Subset to a single PR section between two p_dist_out bookends
GNSS_PR <- GNSS_df %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139)

GNSS_PR_oneday <- GNSS_PR %>%
  filter(as.Date(time_UTC) == as.Date("2024-08-20"))

# Plot: SWOT WSE colored by elevation + GNSS-derived drift WSE with error band
ggplot() +
  # geom_errorbar(...)  # WSE uncertainty bars, currently disabled
  geom_point(data = GNSS_PR_oneday,
             aes(x = p_dist_out / 1000, y = wse, color = wse), size = 4) +
  scale_color_gradient(low = "#2474b7", high = "#d3e3f3") +
  geom_point(data = GNSS_PR_oneday,
             aes(x = p_dist_out / 1000, y = mean_node_drift_wse_no_bias_m),
             color = "#CC79A7", size = 1) +
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000, y = mean_node_drift_wse_no_bias_m, group = 1),
            color = "#CC79A7", linewidth = 0.7) +
  # Upper / lower drift-error envelope
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000,
                y = mean_node_drift_wse_no_bias_m + node_total_error_m, group = 1),
            color = "#CC79A7", alpha = 0.3, linewidth = 2) +
  geom_line(data = GNSS_PR_oneday,
            aes(x = p_dist_out / 1000,
                y = mean_node_drift_wse_no_bias_m - node_total_error_m, group = 1),
            color = "#CC79A7", alpha = 0.3, linewidth = 2) +
  guides(color = "none") +
  xlab("Distance to river outlet (km)") +
  ylab("WSE (m)") +
  xlim(2138, 2152) +
  theme_minimal(base_size = 30)


# -----------------------------------------------------------------------------
# 9b. PT vs. SWOT WSE timeseries (single node)
# -----------------------------------------------------------------------------

# Candidate nodes: 81260300160901, 81260300170011
PT_df <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/PTs/toolboxes_dataframes/reprocessed_2025_09_02/_node/SWORD_v17b/upper_PR/flyby_SWOTCalVal_YR_PT_L1_PT230_20240704T120000_20240821T220000_20250714T183224_SWOTCalVal_YR_KEY_20240704_20240826_v17b.csv"
)

example_node_SWOT_PT <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  filter(node_id == 81260300160901)

# NOTE: -0.1098754 is a node-specific PT bias correction; SWOT WSE plotted on
# the same vertical reference.
ggplot() +
  geom_errorbar(data = PT_df,
                aes(x = pt_time_UTC,
                    ymin = pt_wse_m - pt_correction_mean_total_error_m - 0.1098754,
                    ymax = pt_wse_m + pt_correction_mean_total_error_m - 0.1098754),
                color = "#99D8C9", linewidth = 5, alpha = 0.4) +
  geom_point(PT_df,
             mapping = aes(y = pt_wse_m - 0.1098754, x = pt_time_UTC),
             color = "#009E73", size = 1.2) +
  geom_point(example_node_SWOT_PT,
             mapping = aes(y = wse, x = time_utc),
             shape = 21, fill = "#2474b7", color = "black",
             stroke = 1.6, size = 5.5) +
  xlab("Time") + ylab("WSE (m)") +
  theme_minimal(base_size = 30)


# -----------------------------------------------------------------------------
# 9c. SWOT vs. ortho width along p_dist_out (Porcupine, 2024-07-10)
# -----------------------------------------------------------------------------

node_SWOT_ortho <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/node_width_SWOT_Ortho.csv") %>%
  filter(river == "PR") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5) %>%
  filter(p_dist_out > 2138400) %>%
  filter(p_dist_out < 2152139) %>%
  filter(as.Date(time_utc) == as.Date("2024-07-10"))

ggplot(node_SWOT_ortho) +
  geom_point(aes(x = p_dist_out / 1000, y = ortho_width_m),
             color = "#E69F00", size = 2.5, shape = 17, alpha = 0.7) +
  geom_point(aes(x = p_dist_out / 1000, y = width),
             color = "#2474b7", size = 2.5) +
  xlab("Distance to river outlet (km)") +
  ylab("Width (m)") +
  theme_minimal(base_size = 30)
# export dimensions: width 8.13 in, height 3.96 in


# =============================================================================
# 10. ORTHO SURVEY METADATA TABLE
# =============================================================================
# Unchanged from the previous version of this script.
# NOTE: this block reads and overwrites the same file, so a second run parses
# already-reformatted timestamps and re-derives the columns from them.

ortho_df <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/ortho_summary_stats.csv") %>%
  mutate(
    start_time_UTC = as.POSIXct(start_time_UTC, format = "%m/%d/%y %H:%M", tz = "UTC"),
    end_time_UTC   = as.POSIXct(end_time_UTC,   format = "%m/%d/%y %H:%M", tz = "UTC"),
    SWOT_time_UTC  = as.POSIXct(SWOT_time_UTC,  format = "%m/%d/%y %H:%M", tz = "UTC"))

# Compute midpoint
ortho_df$midpoint_time_UTC <- ortho_df$start_time_UTC +
  (ortho_df$end_time_UTC - ortho_df$start_time_UTC) / 2

# Compute absolute offset in hours between midpoint and SWOT time
ortho_df$offset_time_UTC <- round(abs(as.numeric(difftime(ortho_df$midpoint_time_UTC, ortho_df$SWOT_time_UTC, units = "hours"))), 2)

ortho_df$survey_length_hours <- round(as.numeric(difftime(ortho_df$end_time_UTC, ortho_df$start_time_UTC, units = "hours")), 2)

# Remove midpoint column and reorder
ortho_df <- ortho_df[, c("River(s)", "start_time_UTC", "end_time_UTC",
                         "survey_length_hours", "SWOT_pass_id",
                         "SWOT_time_UTC", "offset_time_UTC")]

# Save to CSV
write.csv(ortho_df,
          "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/_figures/Tables/ortho_summary_stats.csv",
          row.names = FALSE)
