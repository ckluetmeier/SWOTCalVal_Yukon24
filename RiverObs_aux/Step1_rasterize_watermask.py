
#!/usr/bin/env python
"""
=============================================================================
STEP 1 -- Rasterize the digitized water mask (and the orthomosaic footprint)
-----------------------------------------------------------------------------
RiverObs's CalVal path wants a pixel cloud, and the cleanest way to make one
from a digitized polygon is to burn the polygon onto a regular grid and hand
every water cell to the processor. That is what this script does. It writes
one raster per orthomosaic:
 
  <name>_water_<res>m.tif      1 = water, 0 = not water   -> the pixel cloud
 
That is the only file step 2 needs.
 
Optionally, pass --ortho-tif and it will ALSO write
 
  <name>_footprint_<res>m.tif  1 = imaged, 0 = not imaged
 
and cross-check the two. That raster is QC only -- nothing downstream reads it.
The cross-check is still worth having once per survey: it tells you whether any
of your digitized water lies outside the area the survey actually imaged, which
is a digitizing error worth knowing about before you review nodes by hand.
 
It does NOT substitute for the manual node review in step 3c. A nodata-derived
footprint cannot see cloud, because cloud is valid data.
 
WHY 3 m AND NOT 25 cm
  Node area is sum(cell area) over assigned cells, so the grid spacing only
  sets the quantization of that sum. At 3 m a 200 m x 200 m node holds ~4400
  cells, so area quantization is a fraction of a percent -- far below the
  digitization and SWORD-centerline error you already accept. Going to 25 cm
  multiplies the pixel cloud by 144x and makes SWOTRiverEstimator's
  connected-component segmentation allocate an image the size of the full
  orthomosaic grid, which will not fit in memory for a survey tens of km long.
 
Requires: geopandas, rasterio, numpy, shapely
=============================================================================
"""
 
PIPELINE_VERSION = '1.0.0'
 
 
import argparse
import os
 
import numpy as np
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.transform import from_origin
 
 
def snap_bounds(bounds, res, pad):
    """Expand bounds by `pad` metres, then snap outward to a multiple of res."""
    xmin, ymin, xmax, ymax = bounds
    xmin = np.floor((xmin - pad) / res) * res
    ymin = np.floor((ymin - pad) / res) * res
    xmax = np.ceil((xmax + pad) / res) * res
    ymax = np.ceil((ymax + pad) / res) * res
    return xmin, ymin, xmax, ymax
 
 
def burn(gdf, transform, shape, crs, out_path):
    """Burn polygon geometries to a uint8 1/0 raster."""
    arr = features.rasterize(
        ((geom, 1) for geom in gdf.geometry if geom is not None and not geom.is_empty),
        out_shape=shape,
        transform=transform,
        fill=0,
        all_touched=False,      # cell centre must fall inside -> unbiased area
        dtype='uint8',
    )
    profile = dict(
        driver='GTiff', height=shape[0], width=shape[1], count=1,
        dtype='uint8', crs=crs, transform=transform,
        compress='deflate', tiled=True, nodata=None,
    )
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(arr, 1)
    return int(arr.sum())
 
 
def footprint_from_ortho(ortho_path, target_crs):
    """
    Derive the imaged footprint from an orthomosaic's valid-data mask.
 
    Used only if you do not already have a footprint polygon. Reads the
    orthomosaic's internal mask, vectorizes it, and dissolves. For a large
    25 cm mosaic prefer passing --footprint-shp with a footprint you made once
    in QGIS (Raster > Extraction > Contour / Polygonize on the valid mask, or
    the mosaic outline Metashape exports).
    """
    from shapely.geometry import shape as shapely_shape
    from shapely.ops import unary_union
 
    with rasterio.open(ortho_path) as src:
        msk = src.dataset_mask()
        geoms = [
            shapely_shape(geom)
            for geom, val in features.shapes(msk, mask=msk > 0, transform=src.transform)
            if val > 0
        ]
        gdf = gpd.GeoDataFrame(geometry=[unary_union(geoms)], crs=src.crs)
    return gdf.to_crs(target_crs)
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('water_shp', help='digitized water mask shapefile')
    p.add_argument('out_dir', help='output directory')
    p.add_argument('--name', required=True,
                   help='output basename, e.g. upperYR_071024')
    p.add_argument('--res', type=float, default=3.0,
                   help='grid spacing in metres (default 3.0)')
    p.add_argument('--crs', default='EPSG:32606',
                   help='projected CRS for the grid (default EPSG:32606, UTM 6N)')
    p.add_argument('--pad', type=float, default=2000.0,
                   help='metres of padding around the water mask (default 2000). '
                        'Must exceed the widest SWORD max_width in the AOI so '
                        'RiverObs can see the full search corridor.')
    p.add_argument('--ortho-tif', default=None,
                   help='OPTIONAL. The orthomosaic for this survey. If given, '
                        'a footprint raster is written alongside the water '
                        'raster and the two are compared as a sanity check. '
                        'Nothing downstream reads the footprint raster; it is '
                        'a one-off check that your digitizing stayed inside '
                        'the imaged area.')
    p.add_argument('--footprint-shp', default=None,
                   help='OPTIONAL alternative to --ortho-tif, if you have a '
                        'footprint polygon (e.g. a Metashape mosaic outline).')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))
 
    os.makedirs(args.out_dir, exist_ok=True)
 
    water = gpd.read_file(args.water_shp).to_crs(args.crs)
    if water.empty:
        raise SystemExit('water mask is empty: ' + args.water_shp)
 
    foot = None
    if args.footprint_shp:
        foot = gpd.read_file(args.footprint_shp).to_crs(args.crs)
    elif args.ortho_tif:
        foot = footprint_from_ortho(args.ortho_tif, args.crs)
    else:
        print('no --ortho-tif or --footprint-shp given: writing the water '
              'raster only.')
        print('  That is fine -- nothing downstream reads the footprint '
              'raster. The orthomosaic is opened directly in QGIS for the '
              'manual node review.')
 
    # One grid for both rasters, snapped so cell area is exactly res^2.
    xmin, ymin, xmax, ymax = snap_bounds(water.total_bounds, args.res, args.pad)
    ncol = int(round((xmax - xmin) / args.res))
    nrow = int(round((ymax - ymin) / args.res))
    transform = from_origin(xmin, ymax, args.res, args.res)
 
    print('grid: {} rows x {} cols at {} m ({} cells)'.format(
        nrow, ncol, args.res, nrow * ncol))
    if nrow * ncol > 4e8:
        print('WARNING: grid is very large; consider a coarser --res or a '
              'smaller --pad')
 
    res_tag = ('{:g}'.format(args.res)).replace('.', 'p')
    water_tif = os.path.join(
        args.out_dir, '{}_water_{}m.tif'.format(args.name, res_tag))
    foot_tif = os.path.join(
        args.out_dir, '{}_footprint_{}m.tif'.format(args.name, res_tag))
 
    n_water = burn(water, transform, (nrow, ncol), args.crs, water_tif)
    cell_area = args.res ** 2
    print('water cells : {:>12,}  ({:,.0f} m2)'.format(n_water, n_water * cell_area))
    print('wrote', water_tif, '  <-- this is the file step 2 needs')
 
    if foot is not None:
        n_foot = burn(foot, transform, (nrow, ncol), args.crs, foot_tif)
        print('imaged cells: {:>12,}  ({:,.0f} m2)'.format(
            n_foot, n_foot * cell_area))
        if n_water > n_foot:
            print('WARNING: more water cells than imaged cells -- the water '
                  'mask extends beyond the footprint. Check both.')
        outside = int(((rasterio.open(water_tif).read(1) == 1) &
                       (rasterio.open(foot_tif).read(1) == 0)).sum())
        if outside:
            print('WARNING: {:,} water cells ({:.2f}% of the mask) fall '
                  'OUTSIDE the imaged footprint. Either the footprint is '
                  'wrong or the mask was digitized past the survey edge.'
                  .format(outside, 100.0 * outside / max(n_water, 1)))
        print('wrote', foot_tif, '  <-- QC only, not read downstream')
 
 
if __name__ == '__main__':
    main()