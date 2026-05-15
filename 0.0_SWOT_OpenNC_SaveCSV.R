# =============================================================================
# Read SWOT NetCDFs (PIXC & PIXCVec) and Export CSVs
# -----------------------------------------------------------------------------
# Batch-reads SWOT Level 2 pixel-cloud (PIXC) and pixel-vector (PIXCVec)
# NetCDF files from a working directory and writes one CSV per file with a
# selected set of variables. Derived correction fields (geo_corr, wse,
# height_uncert) are computed from the raw NetCDF variables.
#
# Quality filters can be enabled before saving for the PIXC.
# =============================================================================

library(sf)
library(ncdf4)
library(tidyverse)


# =============================================================================
# 1. One-file inspection — check variable names and CRS
# =============================================================================

# Path to a single representative SWOT NetCDF (edit before running)
file_path <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVecRiver_v17b/SWOT_L2_HR_PIXCVecRiver_018_024_034L_20240710T180855_20240710T180906_DevPID0_01.nc"

# Open and inspect the file
nc_data  <- nc_open(file_path)
nc_layers <- names(nc_data$var)

# Full metadata dump
print(nc_data)

# Variable names only
print(names(nc_data$var))

# Native CRS — expected to return NA for most SWOT products
crs_info <- st_crs(nc_data)
print(crs_info)

nc_close(nc_data)


# =============================================================================
# 2. Batch processing — extract variables and write one CSV per file
# =============================================================================

# Working directory containing all .nc files to process
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/pixvec/for_orthos/PIXCVec_v17b"
setwd(wd)

# Enumerate all NetCDF files in the working directory
file_list <- list.files(wd, pattern = "*.nc")

for (j in file_list) {

  # Open the file
  file_path <- j
  nc_data   <- nc_open(file_path)

  # -------------------------------------------------------------------------
  # 2a. PIXC (pixel cloud) variables
  # -------------------------------------------------------------------------
  if (grepl("_PIXC_", j)) {

    # Extract all pixel cloud variables
    longitude                  <- ncvar_get(nc_data, "pixel_cloud/longitude")
    latitude                   <- ncvar_get(nc_data, "pixel_cloud/latitude")
    height                     <- ncvar_get(nc_data, "pixel_cloud/height")
    phase_noise                <- ncvar_get(nc_data, "pixel_cloud/phase_noise_std")
    dheight                    <- ncvar_get(nc_data, "pixel_cloud/dheight_dphase")
    geoid                      <- ncvar_get(nc_data, "pixel_cloud/geoid")
    solid_tide                 <- ncvar_get(nc_data, "pixel_cloud/solid_earth_tide")
    load_tide                  <- ncvar_get(nc_data, "pixel_cloud/load_tide_fes")
    pole_tide                  <- ncvar_get(nc_data, "pixel_cloud/pole_tide")
    class                      <- ncvar_get(nc_data, "pixel_cloud/classification")
    classqual                  <- ncvar_get(nc_data, "pixel_cloud/classification_qual")
    bright_land_flag           <- ncvar_get(nc_data, "pixel_cloud/bright_land_flag")
    ancillary_surf_class_flag  <- ncvar_get(nc_data, "pixel_cloud/ancillary_surface_classification_flag")
    waterfrac                  <- ncvar_get(nc_data, "pixel_cloud/water_frac")
    waterfrac_uncert           <- ncvar_get(nc_data, "pixel_cloud/water_frac_uncert")
    prior_water_prob           <- ncvar_get(nc_data, "pixel_cloud/prior_water_prob")
    geolocqual                 <- ncvar_get(nc_data, "pixel_cloud/geolocation_qual")
    sig0                       <- ncvar_get(nc_data, "pixel_cloud/sig0")
    sig0_uncert                <- ncvar_get(nc_data, "pixel_cloud/sig0_uncert")
    sig0_qual                  <- ncvar_get(nc_data, "pixel_cloud/sig0_qual")
    crosstrack                 <- ncvar_get(nc_data, "pixel_cloud/cross_track")
    pixel_area                 <- ncvar_get(nc_data, "pixel_cloud/pixel_area")
    darea_dheight              <- ncvar_get(nc_data, "pixel_cloud/darea_dheight")
    # pixc_line_qual           <- ncvar_get(nc_data, "pixel_cloud/pixc_line_qual")

    # Assemble data frame
    SWOT_Points <- data.frame(
      longitude                 = longitude,
      latitude                  = latitude,
      height                    = height,
      phase_noise               = phase_noise,
      dheight                   = dheight,
      geoid                     = geoid,
      solid_tide                = solid_tide,
      load_tide                 = load_tide,
      pole_tide                 = pole_tide,
      class                     = class,
      classqual                 = classqual,
      bright_land_flag          = bright_land_flag,
      ancillary_surf_class_flag = ancillary_surf_class_flag,
      waterfrac                 = waterfrac,
      waterfrac_uncert          = waterfrac_uncert,
      prior_water_prob          = prior_water_prob,
      geolocqual                = geolocqual,
      sig0                      = sig0,
      sig0_uncert               = sig0_uncert,
      sig0_qual                 = sig0_qual,
      crosstrack                = crosstrack,
      pixel_area                = pixel_area,
      darea_dheight             = darea_dheight
      # pixc_line_qual          = pixc_line_qual
    )

    # Compute uncertainty & geoid fields
    SWOT_Points$height_uncert <- SWOT_Points$phase_noise * SWOT_Points$dheight
    SWOT_Points$geo_corr      <- SWOT_Points$geoid - SWOT_Points$solid_tide -
                                  SWOT_Points$load_tide - SWOT_Points$pole_tide
    SWOT_Points$wse           <- SWOT_Points$height - SWOT_Points$geo_corr

  }

  # -------------------------------------------------------------------------
  # 2b. PIXCVec (pixel-vector) variables
  # -------------------------------------------------------------------------
  else if (grepl("_PIXCVec", j)) {

    latitude  <- ncvar_get(nc_data, "latitude_vectorproc")
    longitude <- ncvar_get(nc_data, "longitude_vectorproc")
    height    <- ncvar_get(nc_data, "height_vectorproc")
    reach     <- ncvar_get(nc_data, "reach_id")
    node      <- ncvar_get(nc_data, "node_id")

    SWOT_Points <- data.frame(
      longitude = longitude,
      latitude  = latitude,
      height    = height,
      reach     = reach,
      node      = node
    )

  }

  # -------------------------------------------------------------------------
  # 2c. Drop rows with any NA values
  # -------------------------------------------------------------------------
  SWOT_Points <- na.omit(SWOT_Points)

  # -------------------------------------------------------------------------
  # Optional quality filters
  # -------------------------------------------------------------------------
  # SWOT_Points <- SWOT_Points %>%
  #   # Phase unwrapping error flags
  #   filter(!geolocqual %in% c(
  #     4, 4101, 5, 6, 4100, 4102, 524292, 524293,
  #     524294, 524295, 528389, 528390, 7, 528388,
  #     16777220, 17301508, 17305604, 528391, 4103
  #   )) %>%
  #   # Water near land, water, dark water only
  #   filter(class %in% c(3, 4, 5)) %>%
  #   # Cross-track distance 10–60 km
  #   filter(abs(crosstrack) >= 10000) %>%
  #   filter(abs(crosstrack) <= 60000)

  # Sanity check: print dimensions before saving
  print(dim(SWOT_Points))

  # Write CSV (same name as input, .nc replaced with .csv)
  write.csv(SWOT_Points, paste0(sub("\\.nc$", "", j), ".csv"), row.names = TRUE)

  # Shapefile export
  # SWOT_sf <- st_as_sf(SWOT_Points, coords = c("longitude", "latitude"), crs = st_crs(4326))
  # st_write(SWOT_sf, gsub(".nc$", ".shp", j), driver = "ESRI Shapefile")

}
