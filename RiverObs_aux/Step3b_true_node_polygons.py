#!/usr/bin/env python
"""
=============================================================================
STEP 3b -- Real node footprints, from RiverObs's own pixel assignment
-----------------------------------------------------------------------------
`--write-polygons` in step3_coverage_gate.py produces rectangles. Those are a
QC aid for the coverage gate and NOTHING ELSE.
 
This script produces the real thing. `calval2rivertile.py` writes a PIXCVec
file alongside the RiverTile, and it carries a `node_id` for every water cell
in your mask -- RiverObs's actual assignment. Dissolving the water cells by
that node_id gives the exact region each node's area and width were computed
from.
 
    python step3b_true_node_polygons.py \\
        /path/rasters/upperYR_071024_water_3m.tif \\
        /path/riverobs_out/v17b/ortho_upperYR_071024_v17b_pixcvec.nc \\
        /path/out/true_node_polygons_upperYR_071024.shp
 
Add --reaches to also write the reach-level dissolve.
 
WHY THIS IS THE HONEST FIGURE
RiverObs does not use polygons at all. Each water cell is assigned to the node
whose centerline point is nearest in along-track distance, in curvilinear
(along-reach, cross-reach) coordinates on a spline through the SWORD nodes,
subject to a cross-channel cutoff and a segmentation-based extension. The
equivalent region is a curvilinear cell, not a rectangle, and its lateral edges
follow the water mask rather than any geometric boundary. The only faithful way
to draw it is to draw what was actually assigned -- which is what this does.
 
These polygons are also the right basis for a redraw of Figures 9 and S1: same
idea as your current PIXCVec-points-over-node-polygons panels, but now the
polygons are the orthomosaic's own node footprints rather than Thiessen cells.
 
Requires: rasterio, geopandas, shapely, numpy, netCDF4, pandas
=============================================================================
"""
 
SCRIPT_VERSION = 'v1 2026-08-05 -- dissolves ortho water cells by RiverObs node_id'
 
import argparse
import os
 
import numpy as np
import netCDF4
import rasterio
from rasterio import features
import geopandas as gpd
from shapely.geometry import shape as shapely_shape
from shapely.ops import unary_union
 
 
def build(water_tif, pixcvec_nc, key='node_id'):
    with netCDF4.Dataset(pixcvec_nc, 'r') as ds:
        col = np.asarray(ds['range_index'][:]).astype(np.int64)
        row = np.asarray(ds['azimuth_index'][:]).astype(np.int64)
        ids = np.asarray(ds[key][:]).astype(np.int64)
 
    with rasterio.open(water_tif) as src:
        shape, transform, crs = (src.height, src.width), src.transform, src.crs
 
    if row.max() >= shape[0] or col.max() >= shape[1]:
        raise SystemExit(
            'PIXCVec indices exceed the raster shape -- the water .tif and the '
            'PIXCVec file are from different runs. They must be the same pair.')
 
    # rasterio.features.shapes cannot polygonize int64, so map the ids to a
    # dense int32 label and map back afterwards.
    uniq, inverse = np.unique(ids, return_inverse=True)
    label = np.zeros(shape, dtype=np.int32)          # 0 = unassigned
    label[row, col] = inverse.astype(np.int32) + 1
 
    print('{} water cells assigned to {} distinct {}s'.format(
        len(ids), len(uniq), key))
 
    geoms, values = [], []
    for geom, val in features.shapes(label, mask=label > 0, transform=transform):
        geoms.append(shapely_shape(geom))
        values.append(int(val))
    values = np.asarray(values)
 
    out_geom, out_id, out_area = [], [], []
    for v in np.unique(values):
        parts = [g for g, vv in zip(geoms, values) if vv == v]
        merged = unary_union(parts)
        out_geom.append(merged)
        out_id.append(int(uniq[v - 1]))
        out_area.append(merged.area)
 
    gdf = gpd.GeoDataFrame({key: out_id, 'area_m2': out_area},
                           geometry=out_geom, crs=crs)
    return gdf.sort_values(key).reset_index(drop=True)
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('water_tif', help='the water raster from step 1')
    p.add_argument('pixcvec_nc', help='the *_pixcvec.nc from the same run')
    p.add_argument('out_shp')
    p.add_argument('--reaches', action='store_true',
                   help='also write a reach-level dissolve alongside')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), SCRIPT_VERSION))
 
    gdf = build(args.water_tif, args.pixcvec_nc, 'node_id')
    gdf.to_file(args.out_shp)
    print('wrote {} node polygons to {}'.format(len(gdf), args.out_shp))
    print('  total assigned area: {:,.0f} m2'.format(gdf.area_m2.sum()))
    print('  per-node area: min {:,.0f}  median {:,.0f}  max {:,.0f} m2'.format(
        gdf.area_m2.min(), gdf.area_m2.median(), gdf.area_m2.max()))
 
    if args.reaches:
        rgdf = build(args.water_tif, args.pixcvec_nc, 'reach_id')
        rpath = args.out_shp.replace('.shp', '_reaches.shp')
        rgdf.to_file(rpath)
        print('wrote {} reach polygons to {}'.format(len(rgdf), rpath))
 
    print('\nThese are the regions RiverObs actually summed to get node area '
          'and width.\nSafe to put in a figure. The step3 rectangles are not.')
 
 
if __name__ == '__main__':
    main()