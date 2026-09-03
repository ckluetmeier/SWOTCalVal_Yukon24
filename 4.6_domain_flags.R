# =============================================================================
# SWOT Version-Unique Node Observation Diagnostics
# -----------------------------------------------------------------------------
# Asks WHY a paired observation appears in one SWOT version's validation set but
# not the other's. For every version-unique observation, this script goes back
# to the unfiltered node timeseries of the OTHER version, finds that
# same node AT THAT SAME OVERPASS, and works out what happened to it there.
#   - Version C / PIC0 (SWORD v16,  RiverSP)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
#
# -----------------------------------------------------------------------------
# "Unique to C0" does not mean "failed D0's quality filtering". There are five
# ways an observation can end up version-unique:
#
#   no_sword_counterpart    the node does not exist in the other SWORD version
#   node_absent             the node exists but the other version's timeseries
#                           has no row for it at all
#   overpass_absent         the node is in the other timeseries, but not on this
#                           overpass -- SWOT did not deliver that node/pass pair
#   failed_qc               the row exists at this overpass and fails the
#                           quality cascade
#   passed_qc_no_pairing    the row exists at this overpass and passes every
#                           filter, so the SWOT side was fine and the pairing
#                           failed on the IN SITU side
#
# The last category is the interesting one at observation level: it means SWOT
# had good data at that node and pass in both versions, and only the in situ
# assignment differs. Section 6 tests whether that is SWORD renumbering moving
# the in situ point to a neighboring node.
#
# -----------------------------------------------------------------------------
# THE QUALITY CASCADE
# -----------------------------------------------------------------------------
# Reproduced from 1.2_PT_node_wse_diff.R and 2.1_GNSS_node_wse_diff.R (identical
# in both), plus the tighter dark_frac used from 4.1 onwards:
#
#   distinct(node_id, time, wse)   drop repeated rows
#   node_q < 2                     0=good, 1=suspect, 2=degraded, 3=bad
#   abs(xtrk_dist) >= 10 km        inner swath edge
#   abs(xtrk_dist) <= 60 km        outer swath edge
#   dark_frac <= 0.80              matching stage (1.2 / 2.1)
#   dark_frac <  0.50              comparison stage (4.1 / 4.2 / 4.4)
#
#
# -----------------------------------------------------------------------------
# LOCATING AN OBSERVATION IN THE OTHER VERSION'S TIMESERIES
# -----------------------------------------------------------------------------
#   GNSS  the 4.1 key already carries cycle_id and pass_id, so the lookup is a
#         direct join on (node, cycle, pass)
#   PT    the node-level PT file carries pt_serial and pt_time_UTC and no
#         cycle/pass, so the overpass is found by time, reusing 1.2's +-7.5 min
#         window. "The overpass exists but failed QC" therefore means exactly
#         what it means in 1.2
#
# Contains:
#   - Tables: Table S6 rebuild and check, observation-level decomposition,
#             node_q_b bit prevalence, attribute comparison, per-overpass
#             adjacency test and the node pairs behind it
#   - Figures: bit prevalence, node_q class mix, attribute distributions,
#              decomposition by river, adjacency against a permutation null
#
# =============================================================================

library(tidyverse)

source("/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

WSE_BASE       <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse"
SWOT_NODE_DIR  <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"

RAW_C0 <- file.path(SWOT_NODE_DIR, "hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv")
RAW_D0 <- file.path(SWOT_NODE_DIR, "RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv")

OUT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/domain_flags"

# --- filters: must match 4.1 --------------------------------------------------
DARK_FRAC_CMP   <- 0.50    # 4.1 / 4.2 / 4.4 comparison stage
DARK_FRAC_MATCH <- 0.80    # 1.2 / 2.1 matching stage
NODE_Q_MAX      <- 2       # keep node_q < 2
XTRK_MIN        <- 10000   # m
XTRK_MAX        <- 60000   # m
FILL_WSE        <- -1e10   # anything below this is a product fill value

NODE_VALUE   <- "residuals_nobias"   # the metric 4.1 partitions against
PT_MATCH_MIN <- 7.5                  # minutes; 1.2's PT-to-overpass window

# --- TAI -> UTC ---------------------------------------------------------------
TAI_EPOCH      <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
TAI_UTC_OFFSET <- 37   # seconds

# --- published Table S6, node half --------------------------------------------
# The rebuild in section 1 must reproduce these exactly
TABLE_S6_NODE <- tribble(
  ~insitu_type, ~version, ~bucket,    ~n,     ~n_unique,
  "PT",         "D0",     "same",      253L,   51L,
  "PT",         "C0",     "same",      253L,   51L,
  "GNSS",       "D0",     "same",     5120L, 2946L,
  "GNSS",       "C0",     "same",     5120L, 2946L,
  "PT",         "D0",     "D0_only",    90L,   46L,
  "PT",         "C0",     "C0_only",    46L,   23L,
  "GNSS",       "D0",     "D0_only",  1476L, 1057L,
  "GNSS",       "C0",     "C0_only",   875L,  626L
)
STOP_ON_TABLE_S6_MISMATCH <- TRUE   # set FALSE if Table S6 is knowingly stale

# --- adjacency test -----------------------------------------------------------
ADJ_N_PERM    <- 499L      # permutation replicates
ADJ_IMMEDIATE <- 1L        # "adjacent" means within this many node positions
ADJ_SEED      <- 20260831L

# --- node_q_b bit labels -----------------------------------------------------
# node_q_b is the bitwise node quality field. 390 distinct values appear in C0
# and 538 in D0, so the raw integers are not interpretable on their own; the
# bits are.
#
# Table 12. Measurement Quality Flag Bit Definitions in the
# SWOT RiverSP product documentation (node_q_b column).
#
# Bits 5, 6, 8, 12, 16, 17, 20 and 21 are undefined in the table, and bit 15
# (partially_observed) is defined for reach_q_b only, not node_q_b. Those are
# reported as "bit_NN" wherever they appear.
NODE_Q_B_LABELS <- tribble(
  ~bit, ~decimal,     ~label,
   0L,          1,    "sig0_qual_suspect",
   1L,          2,    "classification_qual_suspect",
   2L,          4,    "geolocation_qual_suspect",
   3L,          8,    "water_fraction_suspect",
   4L,         16,    "blocking_width_suspect",
   7L,        128,    "bright_land",
   9L,        512,    "few_sig0_observations",
  10L,       1024,    "few_area_observations",
  11L,       2048,    "few_wse_observations",
  13L,       8192,    "far_range_suspect",
  14L,      16384,    "near_range_suspect",
  18L,     262144,    "classification_qual_degraded",
  19L,     524288,    "geolocation_qual_degraded",
  22L,    4194304,    "lake_flagged",
  23L,    8388608,    "wse_outlier",
  24L,   16777216,    "wse_bad",
  25L,   33554432,    "no_sig0_observations",
  26L,   67108864,    "no_area_observations",
  27L,  134217728,    "no_wse_observations",
  28L,  268435456,    "no_pixels"
)

# The reach flags differ from the node flags at bits 0, 4, 9, 15, 23, 24 and 25
REACH_Q_B_LABELS <- tribble(
  ~bit, ~decimal,     ~label,
   1L,          2,    "classification_qual_suspect",
   2L,          4,    "geolocation_qual_suspect",
   3L,          8,    "water_fraction_suspect",
   7L,        128,    "bright_land",
  10L,       1024,    "few_area_observations",
  11L,       2048,    "few_wse_observations",
  13L,       8192,    "far_range_suspect",
  14L,      16384,    "near_range_suspect",
  15L,      32768,    "partially_observed",
  18L,     262144,    "classification_qual_degraded",
  19L,     524288,    "geolocation_qual_degraded",
  22L,    4194304,    "lake_flagged",
  25L,   33554432,    "below_min_fit_points",
  26L,   67108864,    "no_area_observations",
  27L,  134217728,    "no_wse_observations",
  28L,  268435456,    "no_pixels"
)

# Transcription check: every documented decimal must equal 2^bit. This catches
# a mistyped row in either table before any of it reaches a figure.
walk2(list(NODE_Q_B_LABELS, REACH_Q_B_LABELS), c("node_q_b", "reach_q_b"),
      function(tbl, nm) {
        bad <- tbl %>% filter(decimal != 2^bit)
        if (nrow(bad)) {
          print(bad)
          stop("bit/decimal mismatch in ", nm,
               " labels -- the flag table was transcribed wrong.")
        }
      })

N_QB_BITS <- 32L   # node_q_b is a 32-bit field; the largest value observed in
                   # either file (520093696) occupies 29 bits

# --- columns carried out of the raw files -------------------------------------
# The two products do not share a schema exactly (the D0 file adds crid and
# collection_version), and sword_version is "16" in C0 but "7b" in D0 -- readr
# types the first as a double and the second as character, so binding the whole
# frames fails. Only the columns actually used are carried through, and
# sword_version is forced to character.
ID_COLS        <- c("node_id", "reach_id", "sword_version")
OVERPASS_COLS  <- c("time", "time_tai", "cycle_id", "pass_id")
QC_COLS        <- c("wse", "node_q", "node_q_b", "xtrk_dist", "dark_frac")
ATTR_QC_COLS   <- c("node_q", "dark_frac", "xtrk_dist", "n_good_pix",
                    "layovr_val", "xovr_cal_q")
ATTR_GEOM_COLS <- c("wse_u", "wse_r_u", "width", "width_u", "area_total",
                    "area_detct", "rdr_sig0", "node_dist", "p_length", "p_dist_out")
KEEP_COLS <- unique(c(ID_COLS, OVERPASS_COLS, QC_COLS,
                      ATTR_QC_COLS, ATTR_GEOM_COLS))


# =============================================================================
# 0b. LOCAL HELPERS
# =============================================================================

# --- as_id_chr(): one canonical string form for a SWORD id --------------------
# Same helper as 4.3 / 4.4. A 14-digit node id read as a double renders as
# "8.126e+13" under as.character(), which matches nothing.
as_id_chr <- function(x) {
  if (is.numeric(x)) ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  else               ifelse(is.na(x), NA_character_, trimws(as.character(x)))
}

# --- chr_lut(): put a build_id_lut() map into character space -----------------
chr_lut <- function(lut) {
  amb <- attr(lut, "ambiguous_from_ids")
  out <- lut %>% mutate(from_id = as_id_chr(from_id), to_id = as_id_chr(to_id))
  attr(out, "ambiguous_from_ids") <- as_id_chr(amb)
  out
}

# --- reach_id_from_node(): SWORD node id -> its reach id ----------------------
reach_id_from_node <- function(node_id) {
  n <- as_id_chr(node_id)
  if_else(nchar(n) == 14L,
          paste0(substr(n, 1, 10), substr(n, 14, 14)),
          NA_character_)
}

node_number <- function(node_id) as.integer(substr(as_id_chr(node_id), 11, 13))

# --- label_river(): the river naming used in 4.4 and 4.5 ----------------------
# Takes SWORD v17b reach ids as character.
label_river <- function(reach_id_chr) {
  river_code <- substr(reach_id_chr, 1, 6)
  case_when(
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
}

# --- qb_bit(): is bit k set in a bitwise field? -------------------------------
qb_bit <- function(v, k) {
  out <- (floor(v / 2^k) %% 2) == 1
  out[is.na(v)] <- NA
  out
}

bit_label <- function(bit) {
  lab <- NODE_Q_B_LABELS$label[match(bit, NODE_Q_B_LABELS$bit)]
  if_else(is.na(lab), sprintf("bit_%02d", bit), lab)
}

# --- decode_q_b(): a whole node_q_b value as its list of set flags -------------
# 520093696, the most common value in both products, decodes to
# wse_bad + no_sig0_observations + no_area_observations + no_wse_observations +
# no_pixels -- i.e. "nothing was measured at this node on this pass", which is
# why every fill-WSE row carries it.
decode_q_b <- function(v) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    set <- which(vapply(0:(N_QB_BITS - 1L),
                        function(k) isTRUE(qb_bit(x, k)), logical(1))) - 1L
    if (!length(set)) "none" else paste(bit_label(set), collapse = " + ")
  }, character(1), USE.NAMES = FALSE)
}


# =============================================================================
# 0c. TRANSLATOR
# =============================================================================

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)
node_lut <- chr_lut(build_id_lut(node_translator, "v16_node_id", "v17_node_id",
                                 "node translator"))

# v17b ids that have a v16 match A D0-only observation on a node outside
# this set is on a node that is new in v17b.
v17b_with_v16_ancestor <- unique(node_lut$to_id)


# =============================================================================
# 1. REBUILD THE 4.1 NODE PARTITION, THEN CHECK IT AGAINST TABLE S6
# =============================================================================

# --- 1a. read, filter, harmonize (mirrors read_node() in 4.1) ----------------
# KEEP IN SYNC WITH 4.1!!
read_node_wse <- function(path, insitu, version) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(insitu_type = insitu,
           source      = if (version == "C") VERSION_C else VERSION_D,
           node_id     = as_id_chr(node_id)) %>%
    filter(dark_frac < DARK_FRAC_CMP) %>%
    harmonise_ids("node_id", node_lut, version,
                  label = paste("node", insitu, version))
}

node_PT_vC   <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v16/node_SWOT_PT.csv"),           "PT",   "C")
node_PT_vD   <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v17b/node_SWOT_PT.csv"),          "PT",   "D")
node_GNSS_vC <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv"),  "GNSS", "C")
node_GNSS_vD <- read_node_wse(file.path(WSE_BASE, "node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv"), "GNSS", "D")

# --- 1b. partition on the 4.1 keys -------------------------------------------
node_PT <- bind_rows(node_PT_vC, node_PT_vD) %>%
  partition_versions(key_cols  = c("id_harmonised", "pt_serial", "pt_time_UTC"),
                     value_col = NODE_VALUE)

node_GNSS <- bind_rows(node_GNSS_vC, node_GNSS_vD) %>%
  mutate(drift_file = basename(drift_id)) %>%
  partition_versions(key_cols  = c("id_harmonised", "cycle_id", "pass_id", "drift_file"),
                     value_col = NODE_VALUE)

node_all <- bind_rows(node_PT, node_GNSS)
attr(node_all, "partition_value_col") <- NODE_VALUE

# --- 1c. reproduce Table S6 and assert ---------------------------------------
tableS6_rebuilt <- node_all %>%
  partition_table(NODE_VALUE, by = c("insitu_type", "source"),
                  scale = 100, digits = 1)
print(tableS6_rebuilt, n = Inf)

stopifnot(all(TABLE_S6_NODE$version %in% c("C0", "D0")))
s6_check <- TABLE_S6_NODE %>%
  mutate(source = if_else(version == "D0", VERSION_D, VERSION_C)) %>%
  left_join(tableS6_rebuilt %>%
              mutate(bucket = as.character(bucket)) %>%
              select(insitu_type, source, bucket, n_got = n, n_unique_got = n_unique),
            by = c("insitu_type", "source", "bucket")) %>%
  mutate(n_ok = !is.na(n_got) & n_got == n,
         n_unique_ok = !is.na(n_unique_got) & n_unique_got == n_unique)

if (!all(s6_check$n_ok & s6_check$n_unique_ok)) {
  print(s6_check %>% filter(!n_ok | !n_unique_ok), n = Inf, width = Inf)
  msg <- paste0(
    "the rebuilt partition does not reproduce the published Table S6.\n",
    "  Either the comparison CSVs have changed since the table was made, or\n",
    "  read_node_wse() here has drifted from read_node() in 4.1.\n",
    "  Reconcile before trusting anything below; set ",
    "STOP_ON_TABLE_S6_MISMATCH <- FALSE to proceed anyway.")
  if (STOP_ON_TABLE_S6_MISMATCH) stop(msg) else warning(msg)
} else {
  message("[4.6] OK: rebuild reproduces every published Table S6 node count")
}

# --- 1d. the observation-level sets ------------------------------------------
# One row per version-unique paired observation, carrying everything needed to
# find it again in the other version's timeseries.
obs_cols <- c("obs_key", "id_harmonised", "insitu_type", "source", "obs_bucket",
              "cycle_id", "pass_id", "pt_time_UTC", "drift_file")

observations <- node_all %>%
  filter(!is.na(.data[[NODE_VALUE]]), !is.na(obs_bucket)) %>%
  select(any_of(obs_cols)) %>%
  mutate(reach_id_v17b = reach_id_from_node(id_harmonised),
         river         = label_river(reach_id_v17b))

message("[4.6] observations by bucket (the Table S6 `Count` column):")
observations %>% count(insitu_type, source, obs_bucket) %>% print(n = Inf)

unique_obs <- observations %>% filter(obs_bucket %in% c("C0_only", "D0_only"))
same_obs   <- observations %>% filter(obs_bucket == "same")

# Version-C observations whose v16 node has no v17b counterpart cannot be looked
# up by v17b id at all. Reported, not silently dropped.
unmappable_obs <- node_all %>%
  filter(!is.na(.data[[NODE_VALUE]]), obs_bucket == "unmappable")
message(sprintf("[4.6] %d observation(s) on %d node(s) are 'unmappable' (no v17b counterpart) and are excluded from the lookup",
                nrow(unmappable_obs), n_distinct(unmappable_obs$id_harmonised)))


# =============================================================================
# 2. THE ORIGINAL TIMESERIES, HARMONISED TO SWORD v17b
# =============================================================================

read_raw <- function(path, version) {
  raw <- read_csv(path, show_col_types = FALSE, guess_max = 100000)

  needed <- c("node_id", "reach_id", "time", "time_tai", "wse", "node_q",
              "node_q_b", "xtrk_dist", "dark_frac", "cycle_id", "pass_id")
  missing_cols <- setdiff(needed, names(raw))
  if (length(missing_cols)) {
    stop(basename(path), " is missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  absent <- setdiff(c(ATTR_QC_COLS, ATTR_GEOM_COLS), names(raw))
  if (length(absent)) {
    message("[4.6] ", basename(path), " has no ", paste(absent, collapse = ", "),
            " -- those attributes are dropped from the comparison")
  }

  raw %>%
    select(any_of(KEEP_COLS)) %>%
    mutate(source        = if (version == "C") VERSION_C else VERSION_D,
           node_id       = as_id_chr(node_id),
           reach_id      = as_id_chr(reach_id),
           sword_version = as.character(sword_version)) %>%
    distinct(node_id, time, wse, .keep_all = TRUE) %>%   # 1.2 / 2.1 start here
    harmonise_ids("node_id", node_lut, version, label = paste("raw", version)) %>%
    mutate(time_utc = TAI_EPOCH + time_tai - TAI_UTC_OFFSET)
}

add_qc_flags <- function(df) {
  df %>%
    mutate(
      is_fill        = wse < FILL_WSE,
      pass_node_q    = node_q < NODE_Q_MAX,
      pass_xtrk_near = abs(xtrk_dist) >= XTRK_MIN,
      pass_xtrk_far  = abs(xtrk_dist) <= XTRK_MAX,
      pass_dark_080  = dark_frac <= DARK_FRAC_MATCH,
      pass_dark_050  = dark_frac <  DARK_FRAC_CMP,
      pass_match     = pass_node_q & pass_xtrk_near & pass_xtrk_far & pass_dark_080,
      pass_all       = pass_match & pass_dark_050,
      first_fail = case_when(
        !pass_node_q    ~ "node_q >= 2",
        !pass_xtrk_near ~ "xtrk < 10 km",
        !pass_xtrk_far  ~ "xtrk > 60 km",
        !pass_dark_080  ~ "dark_frac > 0.80",
        !pass_dark_050  ~ "dark_frac >= 0.50",
        TRUE            ~ NA_character_
      ),
      reach_id_v17b = reach_id_from_node(id_harmonised),
      river         = label_river(reach_id_v17b)
    )
}

raw_C0 <- add_qc_flags(read_raw(RAW_C0, "C"))
raw_D0 <- add_qc_flags(read_raw(RAW_D0, "D"))

message("[4.6] cascade survival, whole files:")
bind_rows(raw_C0, raw_D0) %>%
  summarise(rows = n(), fill_rows = sum(is_fill, na.rm = TRUE),
            pass_match = sum(pass_match, na.rm = TRUE),
            pass_all = sum(pass_all, na.rm = TRUE), .by = source) %>%
  print()

fill_survivors <- bind_rows(raw_C0, raw_D0) %>% filter(is_fill, pass_all)
if (nrow(fill_survivors)) {
  print(count(fill_survivors, source, node_q))
  stop("fill-WSE rows survived the quality cascade. 1.2 / 2.1 assume node_q ",
       "removes them; that assumption no longer holds.")
}
message("[4.6] OK: no fill-WSE row survives the cascade in either version")


# =============================================================================
# 3. WHAT HAPPENED TO EACH OBSERVATION IN THE OTHER VERSION
# =============================================================================

# --- 3a. locate one set of observations in one raw file ----------------------
# GNSS joins on (node, cycle, pass). PT has no cycle/pass at node level, so the
# overpass is found by time within 1.2's +-7.5 min window. Returns one row per
# (observation x matching raw row); an observation with no match appears once
# with NA raw columns.
RAW_KEEP <- c("id_harmonised", "cycle_id", "pass_id", "time_utc", "is_fill",
              "pass_all", "pass_match", "first_fail", "node_q", "node_q_b",
              ATTR_QC_COLS, ATTR_GEOM_COLS)

locate_in_raw <- function(obs, raw, label = "") {
  raw_small <- raw %>% select(any_of(RAW_KEEP))

  gnss <- obs %>%
    filter(insitu_type == "GNSS") %>%
    left_join(raw_small, by = c("id_harmonised", "cycle_id", "pass_id"),
              relationship = "many-to-many")

  pt <- obs %>% filter(insitu_type == "PT")
  if (nrow(pt)) {
    pt_hits <- pt %>%
      select(obs_key, id_harmonised, pt_time_UTC) %>%
      inner_join(raw_small, by = "id_harmonised", relationship = "many-to-many") %>%
      filter(abs(as.numeric(difftime(pt_time_UTC, time_utc, units = "mins")))
             <= PT_MATCH_MIN) %>%
      select(-id_harmonised, -pt_time_UTC)
    pt <- pt %>%
      select(-any_of(c("cycle_id", "pass_id"))) %>%
      left_join(pt_hits, by = "obs_key", relationship = "many-to-many")
  }

  out <- bind_rows(gnss, pt)
  message(sprintf("[locate %s] %d observation(s) -> %d row(s); %d observation(s) found no row at their overpass",
                  label, nrow(obs), nrow(out),
                  sum(is.na(out$pass_all[!duplicated(out$obs_key)]))))
  out
}

# C0-only observations are looked for in the D0 timeseries, and vice versa.
unique_C_rows <- locate_in_raw(unique_obs %>% filter(obs_bucket == "C0_only"),
                               raw_D0, "C0_only in D0")
unique_D_rows <- locate_in_raw(unique_obs %>% filter(obs_bucket == "D0_only"),
                               raw_C0, "D0_only in C0")

# The control: observations that paired in BOTH versions, located in the same
# file by the same route. This is a properly matched reference -- same nodes,
# same overpasses, same lookup -- not "all rows of every node".
same_in_D_rows <- locate_in_raw(same_obs %>% filter(source == VERSION_D),
                                raw_D0, "same in D0")
same_in_C_rows <- locate_in_raw(same_obs %>% filter(source == VERSION_C),
                                raw_C0, "same in C0")

# --- 3b. classify each version-unique observation ----------------------------
nodes_in_raw_D <- unique(raw_D0$id_harmonised)
nodes_in_raw_C <- unique(raw_C0$id_harmonised)

classify_obs <- function(rows, nodes_present, check_ancestor, label) {
  rows %>%
    summarise(
      insitu_type   = first(insitu_type),
      obs_bucket    = first(obs_bucket),
      id_harmonised = first(id_harmonised),
      river         = first(river),
      reach_id_v17b = first(reach_id_v17b),
      n_raw_rows    = sum(!is.na(pass_all)),
      n_pass_all    = sum(pass_all,   na.rm = TRUE),
      n_pass_match  = sum(pass_match, na.rm = TRUE),
      n_fill        = sum(is_fill,    na.rm = TRUE),
      .by = obs_key
    ) %>%
    mutate(
      has_ancestor = if (check_ancestor) id_harmonised %in% v17b_with_v16_ancestor else TRUE,
      node_present = id_harmonised %in% nodes_present,
      reason = case_when(
        !has_ancestor    ~ "no_sword_counterpart",
        !node_present    ~ "node_absent",
        n_raw_rows == 0  ~ "overpass_absent",
        n_pass_all == 0  ~ "failed_qc",
        TRUE             ~ "passed_qc_no_pairing"
      ),
      looked_in = label
    )
}

why_C_only <- classify_obs(unique_C_rows, nodes_in_raw_D,
                           check_ancestor = FALSE, label = VERSION_D)
why_D_only <- classify_obs(unique_D_rows, nodes_in_raw_C,
                           check_ancestor = TRUE,  label = VERSION_C)

why_all <- bind_rows(why_C_only, why_D_only) %>%
  mutate(reason = factor(reason, levels = c("no_sword_counterpart", "node_absent",
                                            "overpass_absent", "failed_qc",
                                            "passed_qc_no_pairing")))

message("[4.6] why each version-unique OBSERVATION is unique:")
table3_reasons <- why_all %>%
  count(obs_bucket, insitu_type, reason) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(obs_bucket, insitu_type))
print(table3_reasons, n = Inf)

message("[4.6] ... and the distinct nodes those observations fall on:")
table3_reason_nodes <- why_all %>%
  summarise(n_obs = n(), n_nodes = n_distinct(id_harmonised),
            .by = c(obs_bucket, insitu_type, reason))
print(table3_reason_nodes, n = Inf)

message("[4.6] ... broken down by river:")
table3_reasons_river <- why_all %>% count(obs_bucket, reason, river)
print(table3_reasons_river, n = Inf)

# --- 3c. within failed_qc, which filter did the removing? --------------------
qc_rows <- bind_rows(
  unique_C_rows %>% semi_join(why_all %>% filter(reason == "failed_qc",
                                                 obs_bucket == "C0_only"), by = "obs_key"),
  unique_D_rows %>% semi_join(why_all %>% filter(reason == "failed_qc",
                                                 obs_bucket == "D0_only"), by = "obs_key")
) %>% filter(!is.na(pass_all))

message("[4.6] failed_qc rows, first binding filter:")
table4_binding <- qc_rows %>%
  count(obs_bucket, insitu_type, first_fail) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(obs_bucket, insitu_type))
print(table4_binding, n = Inf)

message("[4.6] failed_qc rows, each filter failed at all (not mutually exclusive):")
table4_anyfail <- qc_rows %>%
  summarise(
    n_rows              = n(),
    `node_q >= 2`       = sum(node_q >= NODE_Q_MAX, na.rm = TRUE),
    `xtrk out of range` = sum(abs(xtrk_dist) < XTRK_MIN | abs(xtrk_dist) > XTRK_MAX, na.rm = TRUE),
    `dark_frac > 0.80`  = sum(dark_frac >  DARK_FRAC_MATCH, na.rm = TRUE),
    `dark_frac >= 0.50` = sum(dark_frac >= DARK_FRAC_CMP,   na.rm = TRUE),
    fill_wse            = sum(is_fill, na.rm = TRUE),
    .by = c(obs_bucket, insitu_type)) %>%
  pivot_longer(-c(obs_bucket, insitu_type, n_rows),
               names_to = "filter", values_to = "n_failing") %>%
  mutate(pct_of_rows = round(100 * n_failing / n_rows, 1))
print(table4_anyfail, n = Inf)


# =============================================================================
# 4. node_q_b BIT DECOMPOSITION
# =============================================================================
# failed_qc rows against the matched control: the rows the SAME lookup returns
# for observations that paired in both versions.

ref_rows <- bind_rows(
  same_in_D_rows %>% filter(!is.na(pass_all)),
  same_in_C_rows %>% filter(!is.na(pass_all))
)

cmp_rows <- bind_rows(
  qc_rows  %>% mutate(group = "failed_qc"),
  ref_rows %>% mutate(group = "both_versions")
) %>%
  mutate(group = factor(group, levels = c("failed_qc", "both_versions")),
         looked_in = if_else(source == VERSION_C, VERSION_D, VERSION_C))

message(sprintf("[4.6] bit analysis on %d failed_qc rows vs %d matched reference rows",
                sum(cmp_rows$group == "failed_qc"),
                sum(cmp_rows$group == "both_versions")))

bit_prevalence <- map_dfr(0:(N_QB_BITS - 1L), function(k) {
  cmp_rows %>%
    summarise(n_rows = n(), n_set = sum(qb_bit(node_q_b, k), na.rm = TRUE),
              .by = c(insitu_type, group)) %>%
    mutate(bit = k, .before = 1)
}) %>%
  mutate(prevalence = n_set / n_rows, flag = bit_label(bit))

bits_used <- bit_prevalence %>%
  summarise(any_set = sum(n_set) > 0, .by = bit) %>%
  filter(any_set) %>% pull(bit)
bit_prevalence <- bit_prevalence %>% filter(bit %in% bits_used)
message(sprintf("[4.6] %d of %d bits are set at least once: %s",
                length(bits_used), N_QB_BITS, paste(bits_used, collapse = ", ")))

table5_bits <- bit_prevalence %>%
  select(insitu_type, bit, flag, group, prevalence) %>%
  pivot_wider(names_from = group, values_from = prevalence)
for (col in c("failed_qc", "both_versions")) {
  if (!col %in% names(table5_bits)) {
    warning("no '", col, "' rows to compare -- the bit contrast is undefined.")
    table5_bits[[col]] <- NA_real_
  }
}
table5_bits <- table5_bits %>%
  mutate(diff  = round(failed_qc - both_versions, 3),
         ratio = round(failed_qc / both_versions, 2),
         failed_qc = round(failed_qc, 3),
         both_versions = round(both_versions, 3)) %>%
  arrange(insitu_type, desc(abs(diff)))
print(table5_bits, n = Inf)

table5_qb_values <- cmp_rows %>%
  count(insitu_type, group, node_q_b) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(insitu_type, group)) %>%
  slice_max(order_by = n, n = 10, by = c(insitu_type, group)) %>%
  arrange(insitu_type, group, desc(n)) %>%
  mutate(flags = decode_q_b(node_q_b))
print(table5_qb_values, n = Inf, width = Inf)

table5_node_q <- cmp_rows %>%
  count(insitu_type, group, node_q) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(insitu_type, group))
print(table5_node_q, n = Inf)


# =============================================================================
# 5. ATTRIBUTE COMPARISON
# =============================================================================
# Continuous attributes on NON-FILL rows only: a fill row carries -1e12 in wse
# and matching fill elsewhere, and mixing those into a median describes the fill
# convention rather than the river. The fill fraction is its own diagnostic.

attr_cols <- intersect(c(ATTR_QC_COLS, ATTR_GEOM_COLS), names(cmp_rows))

table6_attrs <- cmp_rows %>%
  filter(!is_fill) %>%
  pivot_longer(all_of(attr_cols), names_to = "attribute", values_to = "value") %>%
  filter(!is.na(value)) %>%
  summarise(n = n(),
            p25 = round(quantile(value, 0.25), 3),
            median = round(median(value), 3),
            p75 = round(quantile(value, 0.75), 3),
            .by = c(insitu_type, group, attribute)) %>%
  pivot_wider(names_from = group, values_from = c(n, p25, median, p75)) %>%
  arrange(insitu_type, attribute)
print(table6_attrs, n = Inf)

table6_fill <- cmp_rows %>%
  summarise(n_rows = n(), n_fill = sum(is_fill, na.rm = TRUE),
            pct_fill = round(100 * mean(is_fill, na.rm = TRUE), 1),
            .by = c(insitu_type, group))
print(table6_fill)

table7_reach <- why_all %>%
  filter(reason %in% c("failed_qc", "passed_qc_no_pairing")) %>%
  count(obs_bucket, reason, river, reach_id_v17b, sort = TRUE)
print(table7_reach %>% slice_head(n = 30), n = Inf)



# =============================================================================
# 6. FIGURES
# =============================================================================

group_colours <- c(failed_qc = "#D55E00", both_versions = "#999999")

# --- Figure A: decomposition of version-unique observations -------------------
ggplot(table3_reasons, aes(x = reason, y = n, fill = insitu_type)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = n), position = position_dodge(width = 0.85),
            vjust = -0.3, size = 4.5) +
  facet_wrap(~ obs_bucket, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = insitu_colours) +
  labs(x = NULL, y = "Version-unique observations", fill = NULL,
       title = "Why each version-unique observation is unique") +
  theme_minimal(base_size = 17) +
  theme(axis.text.x = element_text(angle = 20, hjust = 0.9),
        legend.position = "top")
# export: 10 x 8 in

# --- Figure B: node_q_b bit contrast ------------------------------------------
figB_data <- table5_bits %>% filter(!is.na(diff)) %>%
  mutate(flag = fct_reorder(flag, diff))

ggplot(figB_data, aes(x = diff, y = flag, fill = diff > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ insitu_type, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(x = "Prevalence in failed_qc minus prevalence in matched control",
       y = NULL, title = "Which quality bits separate the version-unique observations") +
  theme_minimal(base_size = 17)
# export: 11 x 7 in

# --- Figure C: node_q class mix ----------------------------------------------
ggplot(table5_node_q, aes(x = factor(node_q), y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ insitu_type, ncol = 2) +
  scale_fill_manual(values = group_colours) +
  labs(x = "node_q  (0 good, 1 suspect, 2 degraded, 3 bad)",
       y = "Percent of rows", fill = NULL) +
  theme_minimal(base_size = 17) +
  theme(legend.position = "top")
# export: 9 x 5 in

# --- Figure D: quality-driver attributes --------------------------------------
figD_data <- cmp_rows %>%
  filter(!is_fill) %>%
  select(insitu_type, group, all_of(intersect(ATTR_QC_COLS, attr_cols))) %>%
  pivot_longer(-c(insitu_type, group), names_to = "attribute", values_to = "value") %>%
  filter(!is.na(value))

ggplot(figD_data, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.8, colour = NA) +
  geom_boxplot(width = 0.18, fill = "white", outlier.size = 1) +
  facet_grid(attribute ~ insitu_type, scales = "free_y") +
  scale_fill_manual(values = group_colours, guide = "none") +
  labs(x = NULL, y = NULL, title = "Quality-driver attributes, non-fill rows") +
  theme_minimal(base_size = 15)
# export: 9 x 12 in

# --- Figure E: version-unique observations by river ---------------------------
figE_data <- why_all %>%
  filter(!is.na(river)) %>%
  count(obs_bucket, reason, river) %>%
  mutate(river = factor(river, levels = river_levels))

ggplot(figE_data, aes(x = river, y = n, fill = reason)) +
  geom_col(width = 0.8) +
  facet_wrap(~ obs_bucket, ncol = 2) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  labs(x = "River", y = "Version-unique observations", fill = NULL) +
  theme_minimal(base_size = 17) +
  theme(axis.text.x = element_text(angle = 25, hjust = 0.9),
        legend.position = "top")
# export: 11 x 6 in


# =============================================================================
# 8. EXPORT
# =============================================================================

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

write_csv(tableS6_rebuilt,       file.path(OUT, "tableS6_node_rebuilt.csv"))
write_csv(why_all,               file.path(OUT, "version_unique_obs_reasons.csv"))
write_csv(table3_reasons,        file.path(OUT, "reason_counts.csv"))
write_csv(table3_reason_nodes,   file.path(OUT, "reason_counts_with_nodes.csv"))
write_csv(table3_reasons_river,  file.path(OUT, "reason_counts_by_river.csv"))
write_csv(table4_binding,        file.path(OUT, "failed_qc_binding_filter.csv"))
write_csv(table4_anyfail,        file.path(OUT, "failed_qc_filter_failures.csv"))
write_csv(table5_bits,           file.path(OUT, "node_q_b_bit_prevalence.csv"))
write_csv(table5_qb_values,      file.path(OUT, "node_q_b_top_values.csv"))
write_csv(table5_node_q,         file.path(OUT, "node_q_class_mix.csv"))
write_csv(table6_attrs,          file.path(OUT, "attribute_comparison.csv"))
write_csv(table6_fill,           file.path(OUT, "fill_fraction.csv"))
write_csv(table7_reach,          file.path(OUT, "unique_obs_by_reach.csv"))
write_csv(table8_adjacency,      file.path(OUT, "adjacency_test.csv"))
write_csv(table8_pairs,          file.path(OUT, "adjacency_pairs.csv"))
write_csv(table8_pair_summary,   file.path(OUT, "adjacency_pair_summary.csv"))
write_csv(table8_offsets,        file.path(OUT, "adjacency_reach_offsets.csv"))
write_csv(adj_nodes,             file.path(OUT, "adjacency_node_overpass_events.csv"))

message("[4.6] wrote ", length(list.files(OUT, pattern = "\\.csv$")),
        " table(s) to ", OUT)
