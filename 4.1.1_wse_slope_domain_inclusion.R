
# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP") %>%
  rename(old_reach_id = reach_id)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile")

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_SWOT_GNSS.csv")


# SWORD_translator
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")
YR_domain <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/YR_domain_v17b.csv") %>%
  distinct(v17b_reach_id) %>%
  rename(reach_id = v17b_reach_id)


reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(
    SWORD_translator %>% 
      select(v16_reach_id, v17_reach_id),
    by = c("old_reach_id" = "v16_reach_id")
  ) %>%
  rename(reach_id = v17_reach_id)



# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vC, reach_SWOT_PT_vD, reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_wse_m = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m))




all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source) %>%                 # one row per reach_id × source
  group_by(reach_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop"
  ) %>%
  mutate(
    version_inclusion = case_when(
      has_RiverSP & has_RiverTile ~ 0L,
      has_RiverSP & !has_RiverTile ~ -1L,
      !has_RiverSP & has_RiverTile ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(reach_id, version_inclusion)

all_YR_domain_reaches <- YR_domain %>%
  left_join(all_reaches, by = "reach_id") %>%
  mutate(version_inclusion = if_else(is.na(version_inclusion), 2L, version_inclusion))



library(sf)
library(dplyr)

# Read shapefile
sword_sf <- st_read("/Users/camryn/Desktop/SWORD_v17b/NA/na_sword_reaches_hb81_v17b.shp")

# Join to add version_inclusion
all_YR_domain_reaches_sf <- sword_sf %>%
  left_join(all_YR_domain_reaches, by = "reach_id")

# Subset to only YR-domain reaches
all_YR_domain_reaches_sf_subset <- all_YR_domain_reaches_sf %>%
  filter(!is.na(version_inclusion))

st_write(
  all_YR_domain_reaches_sf_subset,
  "/Users/camryn/Desktop/all_YR_domain_reaches_subset.shp",
  delete_layer = TRUE
)



