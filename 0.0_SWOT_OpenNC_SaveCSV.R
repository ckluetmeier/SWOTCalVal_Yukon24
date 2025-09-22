library(sf)
library(ncdf4)

# ---------------------------------------------------------------------------------------------------------------------------
# Read SWOT netCDFs of PIXC & PIXCVec and generate CSVs with selected variables
# ---------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# Set up
# -------------------------------------------------------------------------------------------

# Set working directory
wd = ("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVecRiver_v16")
setwd(wd)

# -------------------------------------------------------------------------------------------
# One file test / look at variables to include
# -------------------------------------------------------------------------------------------

# file path to one SWOT .nc
file_path = ("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVecRiver_v17b/SWOT_L2_HR_PIXCVecRiver_018_024_034L_20240710T180855_20240710T180906_DevPID0_01.nc")

# Open the NetCDF file
nc_data <- nc_open(file_path)
nc_layers <- names(nc_data$var)

# Print variable info

# more metadata then you ever wanted
print(nc_data)
# just the variable names
print(names(nc_data$var))

# Look at native CRS info
crs_info <- st_crs(nc_data)
print(crs_info) # prints at NA

# Close the NetCDF file
nc_close(nc_data)

# -------------------------------------------------------------------------------------
# Batch processing
# -------------------------------------------------------------------------------------

# List all files in the wd
file_list = list.files(wd, pattern='*.nc')

# Define variables (PIXCL) ----
# Loop to extract data of each
for (j in file_list){
  
  # Locate/open file
  file_path = j
  nc_data = nc_open(file_path)
  
  # Extract specific variables from the NetCDF file (PIX_CLOUD)
  if (grepl("_PIXC_", j)) {
    longitude = ncvar_get(nc_data, "pixel_cloud/longitude")
    latitude = ncvar_get(nc_data, "pixel_cloud/latitude")
    height = ncvar_get(nc_data, "pixel_cloud/height")
    phase_noise = ncvar_get(nc_data, "pixel_cloud/phase_noise_std")
    dheight = ncvar_get(nc_data, "pixel_cloud/dheight_dphase")
    #height_uncert = NA
    geoid = ncvar_get(nc_data, "pixel_cloud/geoid")
    solid_tide = ncvar_get(nc_data, "pixel_cloud/solid_earth_tide")
    load_tide = ncvar_get(nc_data, "pixel_cloud/load_tide_fes")
    pole_tide = ncvar_get(nc_data, "pixel_cloud/pole_tide")
    #geo_corr = NA
    #wse = NA
    class = ncvar_get(nc_data, "pixel_cloud/classification")
    classqual = ncvar_get(nc_data, "pixel_cloud/classification_qual")
    bright_land_flag = ncvar_get(nc_data, "pixel_cloud/bright_land_flag")
    ancillary_surf_class_flag = ncvar_get(nc_data, "pixel_cloud/ancillary_surface_classification_flag")
    waterfrac = ncvar_get(nc_data, "pixel_cloud/water_frac")
    waterfrac_uncert = ncvar_get(nc_data, "pixel_cloud/water_frac_uncert")
    prior_water_prob = ncvar_get(nc_data, "pixel_cloud/prior_water_prob")
    geolocqual = ncvar_get(nc_data, "pixel_cloud/geolocation_qual")
    sig0 = ncvar_get(nc_data, "pixel_cloud/sig0")
    sig0_uncert = ncvar_get(nc_data, "pixel_cloud/sig0_uncert")
    sig0_qual = ncvar_get(nc_data, "pixel_cloud/sig0_qual")
    crosstrack = ncvar_get(nc_data, "pixel_cloud/cross_track")
    pixel_area = ncvar_get(nc_data, "pixel_cloud/pixel_area")
    darea_dheight = ncvar_get(nc_data, "pixel_cloud/darea_dheight")
    #pixc_line_qual = ncvar_get(nc_data, "pixel_cloud/pixc_line_qual") 
    
    SWOT_Points = data.frame(longitude=longitude,
                           latitude=latitude,
                           height=height,
                           phase_noise=phase_noise,
                           dheight=dheight,
                           #height_uncert=phase_noise*dheight,
                           geoid=geoid,
                           solid_tide=solid_tide,
                           load_tide=load_tide,
                           pole_tide=pole_tide,
                           #geo_corr=geoid-solid_tide-load_tide-pole_tide,
                           #wse=height-geo_corr,
                           class=class,
                           classqual=classqual,
                           bright_land_flag=bright_land_flag,
                           ancillary_surf_class_flag=ancillary_surf_class_flag,
                           waterfrac=waterfrac,
                           waterfrac_uncert=waterfrac_uncert,
                           prior_water_prob=prior_water_prob,
                           geolocqual=geolocqual,
                           sig0=sig0,
                           sig0_uncert=sig0_uncert,
                           sig0_qual=sig0_qual,
                           crosstrack=crosstrack,
                           pixel_area=pixel_area,
                           darea_dheight=darea_dheight)
                           #pixc_line_qual=pixc_line_qual)
    SWOT_Points$height_uncert = SWOT_Points$phase_noise*SWOT_Points$dheight
    SWOT_Points$geo_corr = SWOT_Points$geoid-SWOT_Points$solid_tide-SWOT_Points$load_tide-SWOT_Points$pole_tide
    SWOT_Points$wse = SWOT_Points$height-SWOT_Points$geo_corr
    } # End PIXCloud processing
  
  # 2.2 Define variables (PIXCV) ----
  
  # Extract specific variables from the NetCDF file (PIXC_VECTOR)
  else if (grepl("_PIXCVec", j)) {
        latitude = ncvar_get(nc_data, "latitude_vectorproc")
        longitude = ncvar_get(nc_data, "longitude_vectorproc")
        height = ncvar_get(nc_data, "height_vectorproc")
        reach = ncvar_get(nc_data, "reach_id")
        node = ncvar_get(nc_data, "node_id")
        
        SWOT_Points = data.frame(longitude=longitude,
                             latitude=latitude,
                             height=height,
                             reach=reach,
                             node=node)
        } # End PIXC_VECTOR processing

  #Make data frames from the extracted data, which can then be output to .csv files
  SWOT_Points = na.omit(SWOT_Points)
  
  # Quality filtering
  # Optional filters to add!
  # SWOT_Points <- SWOT_Points %>%
  #   # phase unwrapping error flags
  #   filter(!geolocqual %in% c(
  #     4, 4101, 5, 6, 4100, 4102, 524292, 524293,
  #     524294, 524295, 528389, 528390, 7, 528388,
  #     16777220, 17301508, 17305604, 528391, 4103
  #   )) %>%
  #   # water near land, water, dark water
  #   filter(class %in% c(
  #     3, 4, 5
  #   )) %>%
  #   # Cross track between 10-60km
  #   filter(abs(crosstrack) >= 10000) %>%
  #   filter(abs(crosstrack) <= 60000)
  
  #Checks to make sure dimensions look good
  print(dim(SWOT_Points))
  
  # Save csv
  write.csv(SWOT_Points, paste0(sub("\\.nc$", "", j), ".csv"), row.names = TRUE)
  
  # Save shapefile
  # SWOT_sf <- st_as_sf(SWOT_Points, coords = c("longitude", "latitude"), crs = st_crs(4326))
  # 
  # st_write(SWOT_sf, gsub(".nc$", ".shp", j), driver = "ESRI Shapefile")
  # Writes shapefiles with data
  # st_write(SWOT_sf, dsn = paste0(j, ".shp"), driver = "ESRI Shapefile")
  
}
