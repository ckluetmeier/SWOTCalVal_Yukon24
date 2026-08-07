#!/usr/bin/env python
"""
=============================================================================
STEP 3b -- Real node footprints, from RiverObs's own pixel assignment
-----------------------------------------------------------------------------
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
 
Two uses:
 
  1. Manual QC. Style the layer by `width` over the orthomosaic and nodes
     truncated by cloud or by the survey edge read as anomalously narrow
     against their neighbours. Mark them, then run step3c.
 
  2. Figures. These are the regions the reported areas and widths were summed
     over, so a panel drawn from them is a faithful depiction of the method.
 
Requires: rasterio, geopandas, shapely, numpy, netCDF4, pandas
=============================================================================
"""
 
PIPELINE_VERSION = '1.0.0'
 
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
 
    counts = {int(uniq[v - 1]): int((inverse == (v - 1)).sum())
              for v in np.unique(values)}
 
    gdf = gpd.GeoDataFrame({
        key: out_id,
        # String copy: shapefile DBF numeric fields are the usual place a
        # 14-digit id gets mangled. Round-trips fine in testing, but if the
        # numeric column ever looks wrong, this one is authoritative.
        key.replace('_id', '_str'): [str(i) for i in out_id],
        'n_cells': [counts.get(i, 0) for i in out_id],
        'area_m2': out_area,
        # Edited in QGIS. 1 = keep, 0 = drop. Read by step3c.
        'keep': [1] * len(out_id),
    }, geometry=out_geom, crs=crs)
    return gdf.sort_values(key).reset_index(drop=True)
 
 
def join_rivertile(gdf, rivertile_nc):
    """
    Attach width, p_length and p_dist_out from the RiverTile so the QGIS layer
    is directly judgeable: style by `width` and a node truncated by cloud or by
    the survey edge stands out as anomalously narrow against its neighbours.
    """
    import pandas as pd
    with netCDF4.Dataset(rivertile_nc, 'r') as ds:
        n = ds.groups['nodes']
        get = lambda k: np.ma.filled(n[k][:].astype('f8'), np.nan)
        rt = pd.DataFrame({
            'node_id': get('node_id').astype('int64'),
            'width': get('width'),
            'p_length': get('p_length'),
            'p_dist_out': get('p_dist_out'),
            'n_good_pix': get('n_good_pix'),
        })
    merged = gdf.merge(rt, on='node_id', how='left')
    return gpd.GeoDataFrame(merged, geometry=gdf.geometry.values, crs=gdf.crs)
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('water_tif', help='the water raster from step 1')
    p.add_argument('pixcvec_nc', help='the *_pixcvec.nc from the same run')
    p.add_argument('out_shp')
    p.add_argument('--rivertile', default=None,
                   help='the RiverTile .nc from the same run. Recommended: '
                        'adds width / p_length / p_dist_out so you can style '
                        'the layer by width during manual review.')
    p.add_argument('--reaches', action='store_true',
                   help='also write a reach-level dissolve alongside')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))
 
    gdf = build(args.water_tif, args.pixcvec_nc, 'node_id')
    if args.rivertile:
        gdf = join_rivertile(gdf, args.rivertile)
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
 
    if 'width' in gdf.columns:
        w = gdf['width'].dropna()
        if len(w):
            print('  width: min {:.1f}  median {:.1f}  max {:.1f} m'.format(
                w.min(), w.median(), w.max()))
            narrow = gdf[gdf['width'] < 0.6 * w.median()]
            if len(narrow):
                print('  {} node(s) narrower than 60% of the median -- worth a '
                      'look first in QGIS:'.format(len(narrow)))
                for _, r in narrow.sort_values('width').head(10).iterrows():
                    print('     {}  width {:.1f} m'.format(
                        int(r['node_id']), r['width']))
 
    print('\nNEXT: open in QGIS over the orthomosaic, set `keep` = 0 on nodes to')
    print('drop (cloud, survey edge), save, then run step3c_manual_node_qc.py.')
 
 
if __name__ == '__main__':
    main()