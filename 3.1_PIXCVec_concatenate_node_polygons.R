# =============================================================================
# PIXCVec Tile Concatenation and Node Polygon Generation
# -----------------------------------------------------------------------------
# Reads PIXCVec tiles (CSV format) from a directory, groups them by overpass
# timestamp, concatenates same-overpass tiles into a single sf point dataset,
# saves the merged outputs as CSV and shapefile, then builds Voronoi (Thiessen)
# node polygons from a selected subset of reaches and a single overpass date.
#
# Script sections:
#   1.  Read and concatenate PIXCVec tiles
#   2.  Save concatenated outputs (CSV and shapefile)
#   3.  Filter to AOI reaches and overpass timestamp
#   4.  Build Voronoi node polygons
#   5.  Save node and reach polygons
# =============================================================================

library(deldir)
library(dplyr)
library(tidyverse)
library(sf)


# =============================================================================
# 1. Read and concatenate PIXCVec tiles
# =============================================================================

# Directory containing PIXCVec tiles in CSV format
folder_path <- '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVec_v17b'
all_files   <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Extract unique overpass identifiers from filenames (characters 33–45 for RiverSP)
# Get the right substring range for different processing versions:
#   RiverSP (PIXCVec):        characters 33–47
#   RiverTile (PIXCVecRiver): characters 38–50
unique_ids <- unique(substring(basename(all_files), 33, 47))

# Pre-compile substitution pattern for creating short display names (15 chars starting at pos 30)
pattern <- "(.{29})(.{15}).*"

# Loop over unique overpass IDs: concatenate all tiles sharing the same ID,
# filter out zero-latitude rows, and convert to sf point objects
PIXCVec_list <- list()
for (id in unique_ids) {

  # Identify all tiles belonging to this overpass
  pair_files <- all_files[grepl(id, substring(basename(all_files), 33, 47))]  # RiverSP: 33,47 | RiverTile: 38,50

  # Short name derived from the first tile's filename (used as list key)
  short_name <- sub(pattern, "\\2", basename(pair_files[1]))

  # Read and row-bind all tiles; suppress scientific notation on node IDs
  concatenated_df <- bind_rows(lapply(pair_files, function(file) read.csv(file))) %>%
    mutate(node = format(node, scientific = FALSE, trim = TRUE))

  # Drop rows with no valid geographic position
  filtered_df <- concatenated_df %>% filter(latitude != 0)

  # Convert to sf; keep lon/lat columns after geometry creation
  # WGS84 (EPSG:4326) used as native SWOT CRS
  filtered_df <- st_as_sf(
    filtered_df,
    coords = c("longitude", "latitude"),
    crs    = 4326,
    remove = FALSE
  )

  PIXCVec_list[[id]] <- filtered_df
}


# =============================================================================
# 2. Save concatenated outputs (CSV and shapefile)
# =============================================================================

output_dir <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVec_v17b/merged"

# Save each overpass as a flat CSV (geometry dropped)
for (name in names(PIXCVec_list)) {
  df_to_save  <- st_drop_geometry(PIXCVec_list[[name]])
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(df_to_save, output_file)
}

# Save each overpass as a shapefile
for (name in names(PIXCVec_list)) {
  output_file <- file.path(output_dir, paste0(name, ".shp"))
  st_write(PIXCVec_list[[name]], output_file, delete_layer = TRUE)
}


# =============================================================================
# 3. Filter to AOI reaches and overpass timestamp
# =============================================================================

# Define the reaches of interest (small chunk; add buffer reach to avoid
# edge effects at polygon boundaries)
SWORD_reach <- data.frame(
  reach_id = c(81260300071, 81260300061, 81260300051, 81260300041,
               81260300031, 81260300211, 81260300201, 81260300191,
               81260300181)
)

# --- Reference: PT/ortho reach IDs by sub-basin ---
# SWORD v16
#   upper_YR:  81270501211, 81270501201, 81270501191, 81270501181, 81270501171
#   lower_YR:  81270100071, 81270100061, 81270100051, 81270100041, 81270100031
#   CD:        81250800091, 81250800081, 81250800071, 81250800061, 81250800051
#   SJ:        81260300251, 81260300241, 81260300231, 81260300061
#   lower_PR:  81260300081, 81260300071, 81260300051, 81260300041, 81260300031
#   CL:        81260401011, 81260401181, 81260401021
#   upper_PR:  81260500021, 81260500011, 81260300221, 81260300211, 81260300191
#
# SWORD v17b
#   upper_YR:  81270500181, 81270500171, 81270500161, 81270500151, 81270500141
#   lower_YR:  81270100071, 81270100061, 81270100051, 81270100041, 81270100031
#   CD:        81250800051, 81250800041, 81250800031, 81250800021, 81250800011
#   lower_PR:  81260300071, 81260300061, 81260300051, 81260300041, 81260300031
#   SJ:        81260300211, 81260300201, 81260300191, 81260300181, 81260300051
#   upper_PR:  81260500021, 81260500011, 81260300171, 81260300161, 81260300151
#   CL:        81260400031, 81260400021, 81260400011, 81260300171

# Combine all overpasses into one sf dataframe with a timestamp column
PIXCVec_points <- bind_rows(PIXCVec_list, .id = "timestamp")

# Spatial filter: keep only points inside the AOI reaches
PIXCVec_AOI_points <- semi_join(PIXCVec_points, SWORD_reach, by = c("reach" = "reach_id"))

# Temporal filter: keep only the overpass matching the orthomosaic collect date
#   upper_YR:       20240710T1809, 20240721T1632, 20240724T0746
#   upper_PR & CL:  20240710T1808, 20240716T0926
#   lower_YR:       20240716T0926, 20240722T1632
#   lower_PR & SJ:  20240710T1808, 20240726T0748
#   CD:             20240711T1809, 20240722T1632
PIXCVec_AOI_points <- PIXCVec_AOI_points %>%
  filter(timestamp == '20240726T0748')  # CHANGE TO MATCH

# Sanity check plot of filtered points
ggplot() +
  geom_sf(data = PIXCVec_AOI_points, color = "blue", size = 1) +
  theme_minimal() +
  ggtitle("PIXCVec points for Reach IDs in AOI")


# =============================================================================
# 4. Build Voronoi (Thiessen) node polygons from PIXCVec points
# =============================================================================

# Project to Web Mercator for Voronoi computation (must be a Cartesian CRS)
PIXCVec_AOI_points <- st_transform(PIXCVec_AOI_points, 3857)

# Extract coordinates and compute Voronoi tessellation
voronoi_points <- st_coordinates(PIXCVec_AOI_points)
voronoi_data   <- deldir(voronoi_points[, 1], voronoi_points[, 2])

# Convert Voronoi tiles to closed sf polygons
voronoi_sf <- tile.list(voronoi_data) %>%
  map(~{
    coords <- cbind(.x$x, .x$y)
    # Close the polygon ring if not already closed
    if (!identical(coords[1, ], coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
    st_polygon(list(coords))
  }) %>%
  st_sfc(crs = 3857) %>%
  st_sf(geometry = .)

# Assign reach and node IDs via nearest-feature spatial join
voronoi_sf <- st_join(voronoi_sf, PIXCVec_AOI_points["reach"], join = st_nearest_feature)
voronoi_sf <- st_join(voronoi_sf, PIXCVec_AOI_points["node"],  join = st_nearest_feature)

# Dissolve (union) polygons by reach and by node
voronoi_sf_reach <- voronoi_sf %>%
  group_by(reach) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

voronoi_sf_node <- voronoi_sf %>%
  group_by(node) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Reproject back to WGS84 for output
voronoi_sf_reach <- st_transform(voronoi_sf_reach, st_crs(4326))
voronoi_sf_node  <- st_transform(voronoi_sf_node,  st_crs(4326))

# Plot Thiessen polygons for visual QC
# # Reach polygons
# ggplot() +
#   geom_sf(data = voronoi_sf_reach, aes(fill = as.factor(reach)), alpha = 0.5) +
#   geom_sf(data = PIXCVec_AOI_points, color = "black", size = 1) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("Thiessen Polygons for PIXCVec Reach")
#
# # Node polygons
# ggplot() +
#   geom_sf(data = voronoi_sf_node, aes(fill = as.factor(node)), alpha = 0.5) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("Thiessen Polygons for PIXCVec Node")


# =============================================================================
# 5. Save node and reach polygons
# =============================================================================

st_write(
  voronoi_sf_node,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVec_v17b/polygons/PIXCVec_lowerPR_SJ_20240726T0748_node.shp")

st_write(
  voronoi_sf_reach,
  "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVec_v17b/polygons/PIXCVec_lowerPR_SJ_20240726T0748_reach.shp")
