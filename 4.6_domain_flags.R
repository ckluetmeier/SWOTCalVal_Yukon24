# =============================================================================
# SWOT Version-Unique Node Diagnostics
# -----------------------------------------------------------------------------
# Asks WHY a node appears in one SWOT version's validation set but not the
# other's. For each version-unique node, this script goes back to the ORIGINAL
# (unfiltered) node timeseries of the OTHER version and works out what happened
# to it there.
#   - Version C / PIC0 (SWORD v16,  RiverSP, hydrocron timeseries)
#   - Version D / PGD0 (SWORD v17b, RiverSP)
#
# The version-unique node sets come from 4.4_wse_slope_domain_inclusion.R, so
# "unique" here means exactly what it means in Tables 2/S6 and the inclusion
# maps: the node produced a usable relative WSE residual in one version and not
# the other, pooled over PT and GNSS.
#
# -----------------------------------------------------------------------------
# WHY THE FOUR-WAY DECOMPOSITION COMES FIRST
# -----------------------------------------------------------------------------
# It is tempting to read "unique to C0" as "failed D0's quality filtering", but
# there are four ways a node can end up version-unique, and only one of them is
# a quality-flag story:
#
#   no_sword_counterpart       the node does not exist in the other SWORD
#                              version at all, so no amount of processing could
#                              have produced an observation
#   absent_from_timeseries     the node exists in that SWORD version but the
#                              other version's timeseries has no row for it
#   failed_qc                  rows exist, but none survive the quality cascade
#                              -- THIS is the group the flag analysis is about
#   passed_qc_no_insitu_match  rows exist and survive the cascade, but never
#                              matched an in situ observation in time/space
#                              (+-7.5 min for PT in 1.2, +-5 hr for GNSS in 2.1)
#                              so they never reached the comparison file
#
# Pooling these together would mix nodes that were never observable with nodes
# that were flagged, and any pattern in the flags would be diluted by nodes that
# have no flags to speak of. Section 3 splits them; sections 4-6 analyse
# `failed_qc` against a reference group of nodes that passed in BOTH versions.
#
# ASYMMETRY WORTH KNOWING: a C0-only node reached its v17b id through the
# translator in 4.1, and 4.4 drops the untranslatable ones, so
# `no_sword_counterpart` cannot occur when looking for C0-only nodes in the D0
# timeseries. It CAN occur in the other direction: a D0-only node may be new in
# v17b with no v16 ancestor. Section 3 checks and reports this.
#
# -----------------------------------------------------------------------------
# THE QUALITY CASCADE
# -----------------------------------------------------------------------------
# Reproduced from 1.2_PT_node_wse_diff.R and 2.1_GNSS_node_wse_diff.R (both use
# an identical cascade), plus the tighter dark_frac used from 4.1 onwards:
#
#   distinct(node_id, time, wse)   drop repeated rows
#   node_q < 2                     0=good, 1=suspect, 2=degraded, 3=bad
#   abs(xtrk_dist) >= 10 km        inner swath edge
#   abs(xtrk_dist) <= 60 km        outer swath edge
#   dark_frac <= 0.80              matching stage (1.2 / 2.1)
#   dark_frac <  0.50              comparison stage (4.1 / 4.2 / 4.4)
#
# VERIFIED, not assumed: the comment in 1.2 says fill values are removed, but
# the code only calls distinct(). It turns out not to matter -- every row with
# a fill WSE (< -1e10) carries node_q = 3 in both files (36,913 of 120,146 rows
# in C0; 30,351 of 122,876 in D0), so `node_q < 2` removes all of them and none
# reach the comparison. They ARE kept here, flagged as `is_fill`, because "the
# node only ever returned fill values" is a diagnosis, not noise.
#
# Section 6 then tests the leading explanation for the dominant category. If
# version-uniqueness comes from SWORD renumbering moving the IN SITU assignment
# between neighbouring nodes, C0-only and D0-only nodes should sit next to each
# other far more often than chance allows. That is tested against a permutation
# null, with nodes usable in both versions as a control.
#
# Contains:
#   - Tables: version-unique decomposition, node_q_b bit prevalence,
#             attribute comparison, spatial/temporal clustering,
#             adjacency test and the node pairs behind it
#   - Figures: bit prevalence, node_q class mix, attribute distributions,
#              unique-node counts by river and overpass, adjacency against the
#              null, and the along-river layout of every node
#
# DEPENDS ON: 4.4 having been run (it writes the inclusion attribute CSV).
# =============================================================================

library(tidyverse)

source("/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/4.0_comparison_helpers.R")


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

SWOT_NODE_DIR  <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/node"
TRANSLATOR_DIR <- "/Users/camryn/Desktop/SWORD_translation"

RAW_C0 <- file.path(SWOT_NODE_DIR, "hydrocron_timeseries/YR_domain_nodes_merged_RiverSP.csv")
RAW_D0 <- file.path(SWOT_NODE_DIR, "RiverSP_v17b/RiverSP_domain_node_timeseries_PGD0_v17b.csv")

# Written by 4.4_wse_slope_domain_inclusion.R
INCLUSION_CSV <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/inclusion_maps/all_YR_domain_nodes_subset_attributes.csv"

OUT <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/domain_flags"

# --- inclusion codes, as written by 4.4 --------------------------------------
INCL_C_ONLY      <- -1L
INCL_BOTH        <-  0L
INCL_D_ONLY      <-  1L
INCL_DOMAIN_ONLY <-  2L

# --- quality cascade ---------------------------------------------------------
NODE_Q_MAX      <- 2       # keep node_q < 2
XTRK_MIN        <- 10000   # m
XTRK_MAX        <- 60000   # m
DARK_FRAC_MATCH <- 0.80    # 1.2 / 2.1
DARK_FRAC_CMP   <- 0.50    # 4.1 / 4.2 / 4.4
FILL_WSE        <- -1e10   # anything below this is a product fill value

# --- node_q_b bit labels -----------------------------------------------------
# node_q_b is the bitwise node quality field. 390 distinct values appear in C0
# and 538 in D0, so the raw integers are not interpretable on their own; the
# bits are.
#
# Transcribed from "Table 12. Measurement Quality Flag Bit Definitions" in the
# SWOT product documentation (node_q_b column). `decimal` is the value given in
# that table and is carried here only so the transcription can be checked
# against 2^bit -- see the assertion below.
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

# Kept for reference only -- this script is node-level. The reach flags differ
# from the node flags at bits 0, 4, 9, 15, 23, 24 and 25, so they are NOT
# interchangeable if this analysis is ever extended to reach_q_b in 4.2.
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
ID_COLS      <- c("node_id", "reach_id", "sword_version")
OVERPASS_COLS<- c("time", "cycle_id", "pass_id")
QC_COLS      <- c("wse", "node_q", "node_q_b", "xtrk_dist", "dark_frac")
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
# A SWORD node id is <10-digit reach base><3-digit node number><1-digit type>,
# and the reach id is <10-digit reach base><1-digit type>. This is the same
# reconstruction that sits commented out in section 6b of 4.4.
reach_id_from_node <- function(node_id) {
  n <- as_id_chr(node_id)
  if_else(nchar(n) == 14L,
          paste0(substr(n, 1, 10), substr(n, 14, 14)),
          NA_character_)
}

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
# Deliberately NOT bitwAnd(): node_q_b is read as a double, and as.integer()
# returns NA above .Machine$integer.max. Double arithmetic is exact well beyond
# the 29 bits this field actually uses.
qb_bit <- function(v, k) {
  out <- (floor(v / 2^k) %% 2) == 1
  out[is.na(v)] <- NA
  out
}

# --- bit_label(): apply NODE_Q_B_LABELS, falling back to bit_NN ---------------
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
# 1. THE VERSION-UNIQUE NODE SETS (from 4.4)
# =============================================================================

if (!file.exists(INCLUSION_CSV)) {
  stop("4.6 needs the node inclusion table written by 4.4:\n  ", INCLUSION_CSV,
       "\nRun 4.4_wse_slope_domain_inclusion.R first (it writes this CSV next ",
       "to the shapefile).")
}

inclusion <- read_csv(INCLUSION_CSV, show_col_types = FALSE) %>%
  mutate(node_id = as_id_chr(node_id))

stopifnot("version_inclusion" %in% names(inclusion))

message("[4.6] inclusion table from 4.4:")
inclusion %>%
  count(version_inclusion) %>%
  mutate(meaning = case_when(
    version_inclusion == INCL_C_ONLY      ~ "C0 only",
    version_inclusion == INCL_BOTH        ~ "both",
    version_inclusion == INCL_D_ONLY      ~ "D0 only",
    version_inclusion == INCL_DOMAIN_ONLY ~ "in domain, no SWOT")) %>%
  print()

nodes_C_only <- inclusion %>% filter(version_inclusion == INCL_C_ONLY) %>% pull(node_id)
nodes_D_only <- inclusion %>% filter(version_inclusion == INCL_D_ONLY) %>% pull(node_id)
nodes_both   <- inclusion %>% filter(version_inclusion == INCL_BOTH)   %>% pull(node_id)

message(sprintf("[4.6] %d C0-only, %d D0-only, %d in both",
                length(nodes_C_only), length(nodes_D_only), length(nodes_both)))


# =============================================================================
# 2. THE ORIGINAL TIMESERIES, HARMONISED TO SWORD v17b
# =============================================================================
# Both files carry the same 31 data columns. The C0 file is SWORD v16, so its
# node ids are translated forward with the same strict 1:1 map used everywhere
# else in phase 4; the D0 file is already v17b and passes through.

node_translator <- read_csv(file.path(TRANSLATOR_DIR, "NA_NodeIDs_v17b_vs_v16.csv"),
                            show_col_types = FALSE)
node_lut <- chr_lut(build_id_lut(node_translator, "v16_node_id", "v17_node_id",
                                 "node translator"))

# --- 2a. read, tag, harmonise -------------------------------------------------
read_raw <- function(path, version) {
  raw <- read_csv(path, show_col_types = FALSE, guess_max = 100000)

  needed <- c("node_id", "reach_id", "time", "wse", "node_q", "node_q_b",
              "xtrk_dist", "dark_frac", "cycle_id", "pass_id")
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
    # 1.2 / 2.1 both start here
    distinct(node_id, time, wse, .keep_all = TRUE) %>%
    harmonise_ids("node_id", node_lut, version, label = paste("raw", version))
}

raw_C0 <- read_raw(RAW_C0, "C")
raw_D0 <- read_raw(RAW_D0, "D")

# --- 2b. apply the cascade, one flag per filter -------------------------------
# Each filter gets its own column so a node can be attributed to the filter that
# actually removed it, rather than only to whichever one happens to come first.
add_qc_flags <- function(df) {
  df %>%
    mutate(
      is_fill        = wse < FILL_WSE,
      pass_node_q    = node_q < NODE_Q_MAX,
      pass_xtrk_near = abs(xtrk_dist) >= XTRK_MIN,
      pass_xtrk_far  = abs(xtrk_dist) <= XTRK_MAX,
      pass_dark_080  = dark_frac <= DARK_FRAC_MATCH,
      pass_dark_050  = dark_frac <  DARK_FRAC_CMP,
      # the 1.2 / 2.1 matching stage
      pass_match     = pass_node_q & pass_xtrk_near & pass_xtrk_far & pass_dark_080,
      # the full cascade a row must survive to reach a 4.1 table
      pass_all       = pass_match & pass_dark_050,
      # first binding constraint, in cascade order (NA if the row passes)
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

raw_C0 <- add_qc_flags(raw_C0)
raw_D0 <- add_qc_flags(raw_D0)

message("[4.6] cascade survival, whole files:")
bind_rows(raw_C0, raw_D0) %>%
  summarise(rows = n(),
            fill_rows = sum(is_fill, na.rm = TRUE),
            pass_match = sum(pass_match, na.rm = TRUE),
            pass_all = sum(pass_all, na.rm = TRUE),
            .by = source) %>%
  print()

# Fill values must never survive the cascade (see the header note). If this
# stops, something changed in the products and the comparison files may contain
# fill WSEs.
fill_survivors <- bind_rows(raw_C0, raw_D0) %>% filter(is_fill, pass_all)
if (nrow(fill_survivors)) {
  print(count(fill_survivors, source, node_q))
  stop("fill-WSE rows survived the quality cascade. 1.2 / 2.1 assume node_q ",
       "removes them; that assumption no longer holds.")
}
message("[4.6] OK: no fill-WSE row survives the cascade in either version")


# =============================================================================
# 3. WHY IS EACH NODE VERSION-UNIQUE?
# =============================================================================
# For each version-unique node, look it up in the OTHER version's raw file.

# The set of v17b ids that have a v16 ancestor. A D0-only node outside this set
# is new in v17b and could not possibly appear in the C0 product.
v17b_with_v16_ancestor <- unique(node_lut$to_id)

classify_uniqueness <- function(unique_ids, other_raw, check_ancestor,
                                label = "") {

  present <- other_raw %>%
    filter(id_harmonised %in% unique_ids) %>%
    summarise(n_rows      = n(),
              n_pass_all  = sum(pass_all,   na.rm = TRUE),
              n_pass_match= sum(pass_match, na.rm = TRUE),
              n_fill      = sum(is_fill,    na.rm = TRUE),
              .by = id_harmonised)

  tibble(node_id = unique_ids) %>%
    left_join(present, by = c("node_id" = "id_harmonised")) %>%
    mutate(
      n_rows       = replace_na(n_rows, 0L),
      n_pass_all   = replace_na(n_pass_all, 0L),
      n_pass_match = replace_na(n_pass_match, 0L),
      n_fill       = replace_na(n_fill, 0L),
      has_ancestor = if (check_ancestor) node_id %in% v17b_with_v16_ancestor else TRUE,
      reason = case_when(
        !has_ancestor   ~ "no_sword_counterpart",
        n_rows == 0     ~ "absent_from_timeseries",
        n_pass_all == 0 ~ "failed_qc",
        TRUE            ~ "passed_qc_no_insitu_match"
      ),
      river = label_river(reach_id_from_node(node_id)),
      set   = label
    )
}

# C0-only nodes: why are they missing from D0? Their v17b id is what 4.1
# assigned them, so an ancestor check is not meaningful in this direction.
why_C_only <- classify_uniqueness(nodes_C_only, raw_D0,
                                  check_ancestor = FALSE, label = "C0_only")

# D0-only nodes: why are they missing from C0? Here a node may simply be new in
# SWORD v17b.
why_D_only <- classify_uniqueness(nodes_D_only, raw_C0,
                                  check_ancestor = TRUE,  label = "D0_only")

why_all <- bind_rows(why_C_only, why_D_only) %>%
  mutate(reason = factor(reason, levels = c("no_sword_counterpart",
                                            "absent_from_timeseries",
                                            "failed_qc",
                                            "passed_qc_no_insitu_match")))

message("[4.6] why each version-unique node is unique:")
table3_reasons <- why_all %>%
  count(set, reason) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = set)
print(table3_reasons, n = Inf)

message("[4.6] ... broken down by river:")
table3_reasons_river <- why_all %>% count(set, reason, river)
print(table3_reasons_river, n = Inf)

# --- 3b. within failed_qc, which filter did the removing? ---------------------
# Reported two ways, because they answer different questions:
#   binding   the first filter in cascade order that a row hit -- what would
#             have to change for the row to survive
#   any       how often each filter fails at all, ignoring order -- whether
#             failures are single-cause or piled up
failed_qc_nodes <- why_all %>% filter(reason == "failed_qc")

qc_rows <- bind_rows(
  raw_D0 %>% filter(id_harmonised %in% failed_qc_nodes$node_id[failed_qc_nodes$set == "C0_only"]) %>%
    mutate(set = "C0_only", looked_in = VERSION_D),
  raw_C0 %>% filter(id_harmonised %in% failed_qc_nodes$node_id[failed_qc_nodes$set == "D0_only"]) %>%
    mutate(set = "D0_only", looked_in = VERSION_C)
)

message("[4.6] failed_qc rows, first binding filter:")
table4_binding <- qc_rows %>%
  count(set, looked_in, first_fail) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(set, looked_in))
print(table4_binding, n = Inf)

message("[4.6] failed_qc rows, each filter failed at all (not mutually exclusive):")
table4_anyfail <- qc_rows %>%
  summarise(
    n_rows            = n(),
    `node_q >= 2`     = sum(!pass_node_q,    na.rm = TRUE),
    `xtrk < 10 km`    = sum(!pass_xtrk_near, na.rm = TRUE),
    `xtrk > 60 km`    = sum(!pass_xtrk_far,  na.rm = TRUE),
    `dark_frac > 0.80`= sum(!pass_dark_080,  na.rm = TRUE),
    `dark_frac >= 0.50` = sum(!pass_dark_050, na.rm = TRUE),
    fill_wse          = sum(is_fill,         na.rm = TRUE),
    .by = c(set, looked_in)) %>%
  pivot_longer(-c(set, looked_in, n_rows), names_to = "filter", values_to = "n_failing") %>%
  mutate(pct_of_rows = round(100 * n_failing / n_rows, 1))
print(table4_anyfail, n = Inf)


# =============================================================================
# 4. node_q_b BIT DECOMPOSITION
# =============================================================================
# The reference group is the nodes that produced a usable residual in BOTH
# versions, using their rows from the SAME file the failed_qc nodes were looked
# up in. That keeps the comparison within one product version, so a difference
# in bit prevalence is a difference between nodes, not between products.

ref_rows <- bind_rows(
  raw_D0 %>% filter(id_harmonised %in% nodes_both) %>%
    mutate(set = "C0_only", looked_in = VERSION_D),
  raw_C0 %>% filter(id_harmonised %in% nodes_both) %>%
    mutate(set = "D0_only", looked_in = VERSION_C)
) %>% mutate(group = "both_versions")

cmp_rows <- bind_rows(
  qc_rows %>% mutate(group = "failed_qc"),
  ref_rows
) %>%
  mutate(group = factor(group, levels = c("failed_qc", "both_versions")))

message(sprintf("[4.6] bit analysis on %d failed_qc rows vs %d reference rows",
                sum(cmp_rows$group == "failed_qc"),
                sum(cmp_rows$group == "both_versions")))

# --- 4a. prevalence of every bit, by set and group ---------------------------
bit_prevalence <- map_dfr(0:(N_QB_BITS - 1L), function(k) {
  cmp_rows %>%
    summarise(n_rows = n(),
              n_set  = sum(qb_bit(node_q_b, k), na.rm = TRUE),
              .by = c(set, looked_in, group)) %>%
    mutate(bit = k, .before = 1)
}) %>%
  mutate(prevalence = n_set / n_rows,
         flag = bit_label(bit))

# Drop bits that are never set anywhere -- there is nothing to say about them.
bits_used <- bit_prevalence %>%
  summarise(any_set = sum(n_set) > 0, .by = bit) %>%
  filter(any_set) %>% pull(bit)

bit_prevalence <- bit_prevalence %>% filter(bit %in% bits_used)
message(sprintf("[4.6] %d of %d bits are set at least once: %s",
                length(bits_used), N_QB_BITS, paste(bits_used, collapse = ", ")))

# --- 4b. the contrast: failed_qc minus both_versions --------------------------
# If one set has no failed_qc nodes at all, that column will not exist and the
# contrast is undefined -- say so rather than failing on a missing column.
table5_bits <- bit_prevalence %>%
  select(set, looked_in, bit, flag, group, prevalence) %>%
  pivot_wider(names_from = group, values_from = prevalence)

for (col in c("failed_qc", "both_versions")) {
  if (!col %in% names(table5_bits)) {
    warning("no '", col, "' rows to compare -- the bit contrast is undefined. ",
            "Check the reason counts in section 3.")
    table5_bits[[col]] <- NA_real_
  }
}

table5_bits <- table5_bits %>%
  mutate(
    diff        = round(failed_qc - both_versions, 3),
    ratio       = round(failed_qc / both_versions, 2),
    failed_qc   = round(failed_qc, 3),
    both_versions = round(both_versions, 3)
  ) %>%
  arrange(set, desc(abs(diff)))
print(table5_bits, n = Inf)

# --- 4c. whole node_q_b values and node_q classes -----------------------------
# Decoding happens after slice_max so only the top composites are expanded.
table5_qb_values <- cmp_rows %>%
  count(set, looked_in, group, node_q_b) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(set, looked_in, group)) %>%
  slice_max(order_by = n, n = 10, by = c(set, looked_in, group)) %>%
  arrange(set, group, desc(n)) %>%
  mutate(flags = decode_q_b(node_q_b))
print(table5_qb_values, n = Inf, width = Inf)

table5_node_q <- cmp_rows %>%
  count(set, looked_in, group, node_q) %>%
  mutate(pct = round(100 * n / sum(n), 1), .by = c(set, looked_in, group))
print(table5_node_q, n = Inf)


# =============================================================================
# 5. ATTRIBUTE COMPARISON
# =============================================================================
# Continuous attributes are summarised on NON-FILL rows only. A fill row carries
# -1e12 in wse and matching fill values elsewhere, and mixing those into a median
# would produce a number that describes the fill convention rather than the
# river. The fill fraction is reported separately, as its own diagnostic.

ATTR_QC <- ATTR_QC_COLS
ATTR_GEOM <- ATTR_GEOM_COLS

attr_cols <- intersect(c(ATTR_QC, ATTR_GEOM), names(cmp_rows))
missing_attrs <- setdiff(c(ATTR_QC, ATTR_GEOM), names(cmp_rows))
if (length(missing_attrs)) {
  message("[4.6] attribute(s) not in these files, skipped: ",
          paste(missing_attrs, collapse = ", "))
}

table6_attrs <- cmp_rows %>%
  filter(!is_fill) %>%
  pivot_longer(all_of(attr_cols), names_to = "attribute", values_to = "value") %>%
  filter(!is.na(value)) %>%
  summarise(
    n      = n(),
    p25    = round(quantile(value, 0.25), 3),
    median = round(median(value), 3),
    p75    = round(quantile(value, 0.75), 3),
    .by = c(set, looked_in, group, attribute)
  ) %>%
  pivot_wider(names_from = group, values_from = c(n, p25, median, p75)) %>%
  arrange(set, attribute)
print(table6_attrs, n = Inf)

table6_fill <- cmp_rows %>%
  summarise(n_rows = n(),
            n_fill = sum(is_fill, na.rm = TRUE),
            pct_fill = round(100 * mean(is_fill, na.rm = TRUE), 1),
            .by = c(set, looked_in, group))
print(table6_fill)

# --- 5b. spatial / temporal clustering ---------------------------------------
# Are version-unique nodes scattered, or do they sit on particular rivers,
# reaches, or overpasses?
table7_river <- why_all %>%
  count(set, reason, river) %>%
  pivot_wider(names_from = reason, values_from = n, values_fill = 0L) %>%
  arrange(set, river)
print(table7_river, n = Inf)

table7_reach <- why_all %>%
  mutate(reach_id_v17b = reach_id_from_node(node_id)) %>%
  filter(reason == "failed_qc") %>%
  count(set, river, reach_id_v17b, sort = TRUE)
print(table7_reach, n = Inf)

table7_overpass <- qc_rows %>%
  count(set, looked_in, cycle_id, pass_id, first_fail) %>%
  arrange(set, cycle_id, pass_id, desc(n))
print(table7_overpass, n = Inf)


# =============================================================================
# 6. ADJACENCY TEST
# =============================================================================
# WHAT THIS IS FOR
# -----------------------------------------------------------------------------
# Section 3 found that ~90% of version-unique nodes fall in
# `passed_qc_no_insitu_match`: the other version observed the node, the data
# passed every quality filter, and it still never produced a residual. The
# leading explanation is NOT a SWOT processing difference at all.
#
# The PT and GNSS node assignments were made against SWORD v16 for the C0 run
# and against SWORD v17b for the D0 run. Wherever SWORD renumbered or reshaped
# nodes between versions, the IN SITU observation moves to a neighbouring node
# even though SWOT observed both nodes perfectly well in both versions. The
# node then looks "unique" to whichever version still has the in situ pairing.
#
# If that is what is happening, version-unique nodes are not scattered at
# random: a C0-only node should sit next to a D0-only node, because the in situ
# observation moved from one to the other. If instead uniqueness reflects
# something real about SWOT's ability to see those nodes, the two classes have
# no reason to be adjacent.
#
# HOW IT IS TESTED
# -----------------------------------------------------------------------------
# Nodes are ordered along each river by p_dist_out (distance to outlet), which
# is monotonic downstream and, unlike node-number arithmetic, does not break at
# reach boundaries. For every C0-only node the distance to the nearest D0-only
# node is measured, in node positions and in metres, and vice versa.
#
# THE NULL MATTERS. With 173 C0-only and 241 D0-only nodes in a domain of a few
# thousand, plenty of them will be neighbours by chance, so a raw adjacency
# count proves nothing on its own. The observed statistics are therefore
# compared against a permutation null: the inclusion labels are shuffled among
# the nodes WITHIN each river (preserving how many of each label a river has,
# and the spatial layout of the nodes themselves), the statistic is recomputed,
# and the observed value is placed in that distribution.
#
# READING THE RESULT
#   observed adjacency much higher than null  -> uniqueness is a SWORD
#                                                reassignment artefact, and the
#                                                "more data in D0" claim is
#                                                largely bookkeeping
#   observed indistinguishable from null      -> uniqueness is spatially real
#                                                and something else explains it
#
# This tests a spatial signature, not the mechanism itself. A positive result is
# strong circumstantial evidence, not proof; the direct confirmation is to take
# a handful of the pairs this produces and check in the PT/GNSS node files
# whether the same physical instrument reading is assigned to the C0-only node
# under v16 and the D0-only node under v17b.
# -----------------------------------------------------------------------------

ADJ_N_PERM     <- 999L   # permutation replicates
ADJ_IMMEDIATE  <- 1L     # "adjacent" means within this many node positions
ADJ_SEED       <- 20260831L

# --- 6a. a v17b spatial frame -------------------------------------------------
# p_dist_out comes from the D0 timeseries because that file is already in v17b
# node space. A node observed on several overpasses has one prior distance, so
# any row will do; the median guards against a stray value.
node_positions <- raw_D0 %>%
  filter(!is.na(p_dist_out)) %>%
  summarise(p_dist_out = median(p_dist_out, na.rm = TRUE), .by = id_harmonised) %>%
  rename(node_id = id_harmonised)

adj_nodes <- inclusion %>%
  select(node_id, version_inclusion) %>%
  left_join(node_positions, by = "node_id") %>%
  mutate(river = label_river(reach_id_from_node(node_id)),
         reach_id_v17b = reach_id_from_node(node_id),
         status = case_when(
           version_inclusion == INCL_C_ONLY      ~ "C0_only",
           version_inclusion == INCL_D_ONLY      ~ "D0_only",
           version_inclusion == INCL_BOTH        ~ "both",
           version_inclusion == INCL_DOMAIN_ONLY ~ "domain_only"))

n_no_pos   <- sum(is.na(adj_nodes$p_dist_out))
n_no_river <- sum(is.na(adj_nodes$river))
message(sprintf("[adjacency] %d node(s) have no p_dist_out in the D0 timeseries and %d have no river label; both are excluded from the test",
                n_no_pos, n_no_river))
if (n_no_pos > 0) {
  adj_nodes %>% filter(is.na(p_dist_out)) %>% count(status) %>% print()
}

adj_nodes <- adj_nodes %>%
  filter(!is.na(p_dist_out), !is.na(river)) %>%
  arrange(river, p_dist_out) %>%
  mutate(rank = row_number(), .by = river)

message("[adjacency] nodes entering the test:")
print(count(adj_nodes, river, status))

# --- 6b. nearest-neighbour distance, vectorised -------------------------------
# For each value in `x`, the distance to the closest value in `targets`. The
# status classes are disjoint, so a node can never match itself.
nearest_dist <- function(x, targets) {
  if (!length(targets)) return(rep(NA_real_, length(x)))
  t   <- sort(targets)
  idx <- findInterval(x, t)                      # how many targets are <= x
  left  <- ifelse(idx >= 1L,        x - t[pmax(idx, 1L)],              Inf)
  right <- ifelse(idx < length(t),  t[pmin(idx + 1L, length(t))] - x,  Inf)
  pmin(left, right)
}

# Statistics computed on one labelling of the nodes. Returned as a named vector
# so the observed run and the permutation runs are guaranteed to be comparable.
adj_stats <- function(df, status_col = "status") {
  s <- df[[status_col]]
  out <- map_dfr(split(seq_along(s), df$river), function(ix) {
    r  <- df$rank[ix]
    tibble(status = s[ix],
           rank   = r,
           d_to_C = nearest_dist(r, r[s[ix] == "C0_only"]),
           d_to_D = nearest_dist(r, r[s[ix] == "D0_only"]))
  })

  c_only <- out %>% filter(status == "C0_only")
  d_only <- out %>% filter(status == "D0_only")
  both   <- out %>% filter(status == "both")

  c(
    # how often the opposite class is an immediate neighbour
    C_adj_D    = mean(c_only$d_to_D <= ADJ_IMMEDIATE, na.rm = TRUE),
    D_adj_C    = mean(d_only$d_to_C <= ADJ_IMMEDIATE, na.rm = TRUE),
    # how far away it is, typically
    C_med_dist = median(c_only$d_to_D, na.rm = TRUE),
    D_med_dist = median(d_only$d_to_C, na.rm = TRUE),
    # control: nodes usable in BOTH versions should show no such affinity
    both_adj_D = mean(both$d_to_D <= ADJ_IMMEDIATE, na.rm = TRUE)
  )
}

observed_stats <- adj_stats(adj_nodes)

# --- 6c. permutation null -----------------------------------------------------
set.seed(ADJ_SEED)
null_draws <- map_dfr(seq_len(ADJ_N_PERM), function(i) {
  shuffled <- adj_nodes %>% mutate(status = sample(status), .by = river)
  as_tibble_row(adj_stats(shuffled))
})

# NB: the vector is `observed_stats`, not `observed` -- the tibble has a column
# called `observed`, and inside mutate() a column always shadows an object of
# the same name in the calling environment.
table8_adjacency <- tibble(
  statistic = names(observed_stats),
  observed  = round(as.numeric(observed_stats), 4)
) %>%
  left_join(
    null_draws %>%
      pivot_longer(everything(), names_to = "statistic", values_to = "value") %>%
      summarise(null_mean = round(mean(value, na.rm = TRUE), 4),
                null_lo   = round(quantile(value, 0.025, na.rm = TRUE), 4),
                null_hi   = round(quantile(value, 0.975, na.rm = TRUE), 4),
                .by = statistic),
    by = "statistic") %>%
  mutate(
    # one-sided empirical p: adjacency statistics are extreme when HIGH,
    # distance statistics when LOW
    direction = if_else(str_detect(statistic, "adj"), "greater", "less"),
    p_value = map2_dbl(statistic, direction, function(st, dir) {
      nd <- null_draws[[st]]
      if (dir == "greater") (1 + sum(nd >= observed_stats[[st]], na.rm = TRUE)) / (ADJ_N_PERM + 1)
      else                  (1 + sum(nd <= observed_stats[[st]], na.rm = TRUE)) / (ADJ_N_PERM + 1)
    }),
    p_value = round(p_value, 4),
    verdict = case_when(
      p_value <= 0.01 ~ "far from chance",
      p_value <= 0.05 ~ "beyond chance",
      TRUE            ~ "consistent with chance")
  )

message("[adjacency] observed vs permutation null (", ADJ_N_PERM, " replicates):")
message("  C_adj_D    fraction of C0-only nodes with a D0-only node within ",
        ADJ_IMMEDIATE, " node position(s)")
message("  D_adj_C    the same, with the roles reversed")
message("  C_med_dist median node positions from a C0-only node to the nearest D0-only node")
message("  D_med_dist the same, reversed")
message("  both_adj_D CONTROL -- nodes usable in both versions should sit at the null")
print(table8_adjacency, n = Inf, width = Inf)

# --- 6d. the pairs themselves, for inspection in QGIS -------------------------
# Every version-unique node with its nearest opposite-class node, so the claim
# can be checked one pair at a time rather than taken on the p-value.
pos_lookup <- adj_nodes %>% select(node_id, river, rank, p_dist_out, reach_id_v17b, status)

nearest_partner <- function(from_status, to_status) {
  map_dfr(split(pos_lookup, pos_lookup$river), function(g) {
    a <- g %>% filter(status == from_status)
    b <- g %>% filter(status == to_status)
    if (!nrow(a) || !nrow(b)) return(NULL)
    j <- vapply(a$rank, function(r) which.min(abs(b$rank - r)), integer(1))
    tibble(
      river            = a$river,
      node_id          = a$node_id,
      node_status      = from_status,
      reach_id_v17b    = a$reach_id_v17b,
      partner_node_id  = b$node_id[j],
      partner_status   = to_status,
      partner_reach    = b$reach_id_v17b[j],
      same_reach       = a$reach_id_v17b == b$reach_id_v17b[j],
      sep_nodes        = abs(b$rank[j] - a$rank),
      sep_m            = round(abs(b$p_dist_out[j] - a$p_dist_out), 1)
    )
  })
}

table8_pairs <- bind_rows(
  nearest_partner("C0_only", "D0_only"),
  nearest_partner("D0_only", "C0_only")
) %>% arrange(node_status, sep_nodes, river)

message("[adjacency] separation between a version-unique node and its nearest opposite-class node:")
table8_pair_summary <- table8_pairs %>%
  summarise(n = n(),
            immediate_neighbour = sum(sep_nodes <= ADJ_IMMEDIATE),
            pct_immediate = round(100 * mean(sep_nodes <= ADJ_IMMEDIATE), 1),
            same_reach = sum(same_reach, na.rm = TRUE),
            median_sep_nodes = median(sep_nodes),
            median_sep_m = median(sep_m),
            .by = c(node_status, river))
print(table8_pair_summary, n = Inf)

message("[adjacency] closest 20 pairs -- these are the ones to check in the PT/GNSS node files:")
print(table8_pairs %>% slice_head(n = 20), n = Inf, width = Inf)


# =============================================================================
# 7. FIGURES
# =============================================================================

group_colours <- c(failed_qc = "#D55E00", both_versions = "#999999")

# --- Figure A: node_q_b bit prevalence, failed_qc vs both ---------------------
figA_data <- bit_prevalence %>%
  mutate(flag = fct_reorder(flag, bit))

ggplot(figA_data, aes(x = flag, y = prevalence, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ paste0(set, " (looked up in ", looked_in, ")"), ncol = 1) +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "node_q_b bit", y = "Rows with bit set",
       fill = NULL,
       title = "Quality bits: version-unique nodes that failed QC vs nodes usable in both") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
# export: 11 x 8 in

# --- Figure B: the contrast, sorted ------------------------------------------
figB_data <- table5_bits %>%
  filter(!is.na(diff)) %>%
  mutate(flag = fct_reorder(flag, diff))

ggplot(figB_data, aes(x = diff, y = flag, fill = diff > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ set, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(x = "Prevalence in failed_qc minus prevalence in both-version nodes",
       y = NULL,
       title = "Which quality bits separate the version-unique nodes") +
  theme_minimal(base_size = 18)
# export: 11 x 7 in

# --- Figure C: node_q class mix ----------------------------------------------
ggplot(table5_node_q, aes(x = factor(node_q), y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ set, ncol = 2) +
  scale_fill_manual(values = group_colours) +
  labs(x = "node_q  (0 good, 1 suspect, 2 degraded, 3 bad)",
       y = "Percent of rows", fill = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top")
# export: 9 x 5 in

# --- Figure D: attribute distributions ---------------------------------------
figD_data <- cmp_rows %>%
  filter(!is_fill) %>%
  select(set, group, all_of(intersect(ATTR_QC, attr_cols))) %>%
  pivot_longer(-c(set, group), names_to = "attribute", values_to = "value") %>%
  filter(!is.na(value))

ggplot(figD_data, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.8, colour = NA) +
  geom_boxplot(width = 0.18, fill = "white", outlier.size = 1) +
  facet_grid(attribute ~ set, scales = "free_y") +
  scale_fill_manual(values = group_colours, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Quality-driver attributes, non-fill rows") +
  theme_minimal(base_size = 16)
# export: 9 x 12 in

# --- Figure E: where the version-unique nodes are ----------------------------
figE_data <- why_all %>%
  filter(!is.na(river)) %>%
  count(set, reason, river) %>%
  mutate(river = factor(river, levels = river_levels))

ggplot(figE_data, aes(x = river, y = n, fill = reason)) +
  geom_col(width = 0.8) +
  facet_wrap(~ set, ncol = 2) +
  scale_x_discrete(breaks = river_levels, labels = river_labels) +
  labs(x = "River", y = "Version-unique nodes", fill = NULL) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 25, hjust = 0.9),
        legend.position = "top")
# export: 11 x 6 in


# --- Figure F: observed adjacency against the permutation null ----------------
figF_data <- null_draws %>%
  select(C_adj_D, D_adj_C, both_adj_D) %>%
  pivot_longer(everything(), names_to = "statistic", values_to = "value")

figF_obs <- table8_adjacency %>% filter(statistic %in% unique(figF_data$statistic))

ggplot(figF_data, aes(x = value)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(data = figF_obs, aes(xintercept = observed),
             colour = "#D55E00", linewidth = 1.2) +
  geom_text(data = figF_obs,
            aes(x = observed, y = Inf, label = paste0("observed = ", observed,
                                                      "\np = ", p_value)),
            hjust = -0.05, vjust = 1.4, size = 5, colour = "#D55E00") +
  facet_wrap(~ statistic, ncol = 1, scales = "free") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Fraction with an opposite-class node as immediate neighbour",
       y = "Permutation replicates",
       title = "Are version-unique nodes adjacent more often than chance?",
       subtitle = "Grey = label-shuffled null; both_adj_D is the control and should sit inside it") +
  theme_minimal(base_size = 17)
# export: 9 x 9 in

# --- Figure G: the nodes laid out along each river ----------------------------
# The persuasive picture. If uniqueness is a reassignment artefact, C0-only and
# D0-only points interleave in tight pairs rather than forming separate runs.
status_colours <- c(both = "#CCCCCC", domain_only = "#8C8C8C",
                    C0_only = version_colours[[VERSION_C]],
                    D0_only = version_colours[[VERSION_D]])

figG_data <- adj_nodes %>%
  mutate(river  = factor(river, levels = river_levels),
         status = factor(status, levels = c("both", "domain_only", "C0_only", "D0_only")))

ggplot(figG_data, aes(x = p_dist_out / 1000, y = status, colour = status)) +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_wrap(~ river, ncol = 1, scales = "free_x",
             labeller = labeller(river = setNames(river_labels, river_levels))) +
  scale_colour_manual(values = status_colours, guide = "none") +
  labs(x = "Distance to river outlet (km)", y = NULL,
       title = "Where the version-unique nodes sit along each river") +
  theme_minimal(base_size = 15) +
  theme(panel.grid.minor = element_blank())
# export: 12 x 14 in

# --- Figure H: separation to the nearest opposite-class node ------------------
ggplot(table8_pairs, aes(x = sep_nodes, fill = node_status)) +
  geom_histogram(binwidth = 1, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = c(C0_only = version_colours[[VERSION_C]],
                               D0_only = version_colours[[VERSION_D]])) +
  labs(x = "Node positions to the nearest opposite-class node",
       y = "Nodes", fill = NULL,
       title = "How close is each version-unique node to one unique to the other version?") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_minimal(base_size = 17) +
  theme(legend.position = "top")
# export: 9 x 5 in


# =============================================================================
# 8. EXPORT
# =============================================================================

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

write_csv(why_all,               file.path(OUT, "version_unique_node_reasons.csv"))
write_csv(table3_reasons,        file.path(OUT, "reason_counts.csv"))
write_csv(table3_reasons_river,  file.path(OUT, "reason_counts_by_river.csv"))
write_csv(table4_binding,        file.path(OUT, "failed_qc_binding_filter.csv"))
write_csv(table4_anyfail,        file.path(OUT, "failed_qc_filter_failures.csv"))
write_csv(table5_bits,           file.path(OUT, "node_q_b_bit_prevalence.csv"))
write_csv(table5_qb_values,      file.path(OUT, "node_q_b_top_values.csv"))
write_csv(table5_node_q,         file.path(OUT, "node_q_class_mix.csv"))
write_csv(table6_attrs,          file.path(OUT, "attribute_comparison.csv"))
write_csv(table6_fill,           file.path(OUT, "fill_fraction.csv"))
write_csv(table7_river,          file.path(OUT, "unique_nodes_by_river.csv"))
write_csv(table7_reach,          file.path(OUT, "failed_qc_by_reach.csv"))
write_csv(table7_overpass,       file.path(OUT, "failed_qc_by_overpass.csv"))
write_csv(table8_adjacency,      file.path(OUT, "adjacency_test.csv"))
write_csv(table8_pairs,          file.path(OUT, "adjacency_pairs.csv"))
write_csv(table8_pair_summary,   file.path(OUT, "adjacency_pair_summary.csv"))
write_csv(adj_nodes,             file.path(OUT, "adjacency_node_positions.csv"))

message("[4.6] wrote ", length(list.files(OUT, pattern = "\\.csv$")),
        " table(s) to ", OUT)
