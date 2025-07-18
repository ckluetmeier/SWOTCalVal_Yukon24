library(deldir)
library(dplyr)

# Concatenate PIXCVec tiles
# ---------------------------------------------------------------------------------------------------------------------------

# Import

# Directory with all pixcvec tiles in csv format
folder_path = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos'
all_files = list.files(folder_path, pattern="*.csv", full.names=TRUE)

# Extract unique identifiers from file names: start datetime only down to minutes of overpass
unique_ids <- unique(substring(basename(all_files), 33, 45)) #characters 33-47 in the filename string

# ---------------------------------------------------------------------------------------------------------------------------
# Concatenate
# ---------------------------------------------------------------------------------------------------------------------------

# Pre-compile the pattern for extracting short names (characters 33-47)
pattern <- "(.{29})(.{15}).*"

PIXCVec_list = list()
for (id in unique_ids) {
  
  # Identify files with the same unique ID
  pair_files = all_files[grepl(id, substring(basename(all_files), 33, 45))]
  
  # Create a short name for the sf object
  short_name <- sub(pattern, "\\2", basename(pair_files[1]))
  
  # Read and concatenate the grouped files
  concatenated_df = bind_rows(lapply(pair_files, function(file) {
    read_csv(file)})) %>% 
    mutate(node = format(node, scientific = FALSE, trim = TRUE))
  
  # Filter out rows where the latitude is 0
  filtered_df = concatenated_df %>% filter(latitude != 0)
  
  # Convert the concatenated dataframe to sf object
  filtered_df = st_as_sf(filtered_df,
                         coords = c("longitude", "latitude"),
                         crs = 4326, # doesn't work in a projected UTM (use global WGS84 or ITRF14)
                         remove = FALSE)
  
  # Store the sf object in the list with its short name
  PIXCVec_list[[id]] = filtered_df
}

# ---------------------------------------------------------------------------------------------------------------------------
# Save concatenate outputs as CSV or shapefiles
# ---------------------------------------------------------------------------------------------------------------------------

# Directory to save the files
output_dir <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/merged"

# Loop over PIXCVec_list and save each sf object as a CSV
for (name in names(PIXCVec_list)) {
  sf_object <- PIXCVec_list[[name]]
  
  # Drop geometry before saving
  df_to_save <- st_drop_geometry(sf_object)
  
  # Define the output filename
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write to CSV
  write_csv(df_to_save, output_file)
}

# Loop over PIXCVec_list and save each sf object as a Shapefile
for (name in names(PIXCVec_list)) {
  sf_object <- PIXCVec_list[[name]]
  
  # Define the output filename
  output_file <- file.path(output_dir, paste0(name, ".shp"))
  
  # Write to Shapefile
  st_write(sf_object, output_file, delete_layer = TRUE)
}

# ---------------------------------------------------------------------------------------------------------------------------
# Filter reaches and date for node polygon generation
# ---------------------------------------------------------------------------------------------------------------------------

# Define input data
SWORD_reach <- data.frame(reach_id = c()) # put small chunk of reaches here

# need to build node polygons in small chunks for code to run
# should come back and add a buffer reach so that reach/final node values are valid
# the final polygons close in a weird way

# All YR PT/ortho reaches:
# 81270100061, 81270100051, 81270100041,
# 81250800081, 81250800071, 81250800061,
# 81260200031, 81260200021, 81260200011, 81260300061,
# 81260300071, 81260300051, 81260300041,
# 81260401011, 81260401181,
# 81260500011, 81260300221, 81260300211

# Completed YR PT reaches
# one out of two days done:
# 81270501201, 81270501191, 81270501181
# both days done:
# 81260500011, 81260300221, 81260300211, 81260401011, 81260401181

# Combine all points from PIXCVec_list into one sf df
PIXCVec_points <- bind_rows(PIXCVec_list, .id = "timestamp")

# Filter PIXCVec_points to reach_id in SWORD_reach
PIXCVec_AOI_points <- semi_join(PIXCVec_points, SWORD_reach, by = c("reach" = "reach_id"))
# Filter to one day (since the orthos cover PT reaches twice)
PIXCVec_AOI_points <- PIXCVec_AOI_points %>%
  filter(timestamp == '20240716T0926')

# Plot to double check
ggplot() +
  geom_sf(data = PIXCVec_AOI_points, color = "blue", size = 1) +
  theme_minimal() +
  ggtitle("PIXCVec points for Reach IDs in AOI")

# ---------------------------------------------------------------------------------------------------------------------------
# Build node polygons from PIXCVec
# ---------------------------------------------------------------------------------------------------------------------------

# Ensure PIXCVec_AOI_points is projected (Voronoi works best in projected CRS)
PIXCVec_AOI_points <- st_transform(PIXCVec_AOI_points, 3857)  # Web Mercator projection

# Extract coordinates of points
voronoi_points <- st_coordinates(PIXCVec_AOI_points)
voronoi_data <- deldir(voronoi_points[,1], voronoi_points[,2]) # Compute Voronoi tessellation (Thiessen polygons)

# Convert Voronoi tiles to polygons and ensure they are closed
voronoi_sf <- tile.list(voronoi_data) %>%
  # Ensure polygon is closed by appending the first point at the end
  map(~{
    coords <- cbind(.x$x, .x$y)
    if (!identical(coords[1,], coords[nrow(coords),])) {coords <- rbind(coords, coords[1,])}
    st_polygon(list(coords))
  }) %>%
  st_sfc(crs = 3857) %>%
  st_sf(geometry = .)

# Add reach_id by spatially joining with original points
voronoi_sf <- st_join(voronoi_sf, PIXCVec_AOI_points["reach"], join = st_nearest_feature)
voronoi_sf <- st_join(voronoi_sf, PIXCVec_AOI_points["node"], join = st_nearest_feature)
voronoi_sf_reach <- voronoi_sf %>%
  group_by(reach) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")
voronoi_sf_node <- voronoi_sf %>%
  group_by(node) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Convert CRS
voronoi_sf_reach <- st_transform(voronoi_sf_reach, st_crs(4326)) 
voronoi_sf_node <- st_transform(voronoi_sf_node, st_crs(4326)) 

# # Check by plotting the Thiessen polygons
# # Reaches plot
# ggplot() +
#   geom_sf(data = voronoi_sf_reach, aes(fill = as.factor(reach)), alpha = 0.5) +
#   geom_sf(data = PIXCVec_AOI_points, color = "black", size = 1) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("Thiessen Polygons for PIXCVec Reach")
# 
# # Nodes plot
# ggplot() +
#   geom_sf(data = voronoi_sf_node, aes(fill = as.factor(node)), alpha = 0.5) +
#   #geom_sf(data = PIXCVec_AOI_points, color = "black", size = 1) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("Thiessen Polygons for PIXCVec Node")


# Save polygons
# ---------------------------------------------------------------------------------------------------------------------------

st_write(voronoi_sf_node, "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/polygons/PIXCVec_CL_upperPR_20240716T0926_node.shp")

st_write(voronoi_sf_reach, "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/polygons/voronoi_sf_reach2.shp")



