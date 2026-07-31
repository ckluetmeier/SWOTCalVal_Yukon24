# =============================================================================
# Unique version flags
# -----------------------------------------------------------------------------

# =============================================================================

library(tidyverse)
library(lubridate)
library(dplyr)


# =============================================================================
# NODE LEVEL
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Read in node-level data
# -----------------------------------------------------------------------------

# --- PT ----------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
node_SWOT_PT_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
node_SWOT_PT_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# --- GNSS ---------------------------------------------------------------------

# Version C (SWORD v16 / RiverSP PIC0)
node_SWOT_GNSS_vC <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)

# Version D (SWORD v17b / RiverSP PGD0)
node_SWOT_GNSS_vD <- read_csv(
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)


# -----------------------------------------------------------------------------
# 2. Harmonize to SWORD v17b node IDs
# -----------------------------------------------------------------------------

# Translator: maps v16 node IDs to v17b
SWORD_translator <- read_csv(
  "/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv")

# Apply translation to Version C PT data
node_SWOT_PT_vC <- node_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# Apply translation to Version C GNSS data
node_SWOT_GNSS_vC <- node_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# Merge all four data frames; unify column names for WSE and time
node_SWOT_full_insitu <- bind_rows(
  node_SWOT_PT_vC,
  node_SWOT_PT_vD,
  node_SWOT_GNSS_vC,
  node_SWOT_GNSS_vD) %>%
  mutate(
    insitu_wse_m        = coalesce(pt_wse_m, mean_node_drift_wse_m),
    insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_node_drift_wse_no_bias_m),
    insitu_time_utc     = coalesce(pt_time_UTC, time_UTC))

# Compute version inclusion flag per node:
#   -1 = observed only in Version C
#    0 = observed in both versions
#    1 = observed only in Version D
all_nodes <- node_SWOT_full_insitu %>%
  distinct(node_id, source) %>%
  group_by(node_id) %>%
  summarise(
    has_PIC0   = any(source == "PIC0"),
    has_PGD0 = any(source == "PGD0"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_PIC0 & has_PGD0  ~  0L,
      has_PIC0 & !has_PGD0 ~ -1L,
      !has_PIC0 & has_PGD0 ~  1L,
      TRUE                         ~ NA_integer_
    )) %>%
  select(node_id, version_inclusion)

# Join version inclusion to full dataset
node_SWOT_full_insitu <- node_SWOT_full_insitu %>%
  left_join(all_nodes, by = "node_id")

# Create a separate PT and GNSS dataset
node_SWOT_PT   <- node_SWOT_full_insitu %>% filter(insitu_type == "PT")

node_SWOT_GNSS_vD_unique <- node_SWOT_full_insitu %>% 
  filter(insitu_type == "GNSS") %>% 
  filter(version_inclusion == 1)

node_SWOT_GNSS_vC_unique <- node_SWOT_full_insitu %>% 
  filter(insitu_type == "GNSS") %>% 
  filter(version_inclusion == -1)





