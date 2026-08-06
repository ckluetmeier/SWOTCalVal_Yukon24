#!/usr/bin/env python
"""
=============================================================================
STEP 3 -- Coverage gate: which nodes did the orthomosaic actually see in full?
-----------------------------------------------------------------------------
RiverObs has no way to know that your reference dataset stops at the edge of an
aerial survey. A SWORD node at the edge of an orthomosaic gets a real-looking
water area computed from only the imaged part of the node, which propagates
into an ortho width that is too small, which reads as SWOT overestimating
width. This is a false positive that would sit right on top of your headline
finding that SWOT is biased high.
 
Your current workflow suppresses this implicitly: the AOI clip in script 3.2
and the hand-curated `viable_reaches` list in script 3.4. This script makes it
explicit and reproducible, using the same node-polygon construction the
manuscript already describes for the GNSS comparison (section 3.3.3): a
rectangle centred on each SWORD node, node_length long in the flow direction
and scaled from SWORD max_width across it.
 
A node passes if its test rectangle is fully contained in the orthomosaic
footprint.
 
*** THE RECTANGLES ARE NOT NODE POLYGONS ***
The rectangle exists to answer one binary question: was this node's
neighbourhood fully inside the survey? It is a deliberate simplification and it
is intentionally a little more generous than RiverObs's real search corridor
(cross_frac 0.5 of max_width, versus RiverObs's max_width/3), so the gate errs
toward rejecting marginal nodes.
 
RiverObs does not use polygons at all. Each water cell goes to the node nearest
in along-track distance, in curvilinear (along-reach, cross-reach) coordinates
on a spline through the SWORD nodes. The region a node's area was actually
summed over is a curvilinear cell whose lateral edges follow the water mask.
To draw that -- for a figure, or to check attribution -- use
step3b_true_node_polygons.py, which dissolves the run's own PIXCVec node_id
assignment.
 
Output: a CSV of node_id plus `fully_imaged` (bool) and `imaged_frac` (the
fraction of the node polygon inside the footprint), to left-join onto the
step 2 node CSV.
 
Requires: geopandas, shapely, numpy, pandas
=============================================================================
"""
 
SCRIPT_VERSION = 'v4 2026-08-05 -- ortho .tif footprint; max-based guard; rectangles flagged as QC-only'
 
 
import argparse
 
import os
 
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
 
 
def load_footprint(path, crs):
    """
    Load the orthomosaic footprint from either a polygon shapefile or the
    orthomosaic raster itself.
 
    Passing the raster is usually easier -- no QGIS step -- but check the
    warning it prints. A gap-filled mosaic whose nodata value is not set will
    report every cell as valid, and the "footprint" degenerates to the full
    raster rectangle, which passes nodes the survey never actually imaged.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.shp', '.gpkg', '.geojson', '.json'):
        return gpd.read_file(path).to_crs(crs)
 
    import rasterio
    from rasterio import features
    from shapely.geometry import shape as shapely_shape, box as shapely_box
    from shapely.ops import unary_union
 
    with rasterio.open(path) as src:
        msk = src.dataset_mask()
        valid_frac = float((msk > 0).sum()) / msk.size
        geoms = [shapely_shape(g) for g, v in features.shapes(
            msk, mask=msk > 0, transform=src.transform) if v > 0]
        if not geoms:
            raise SystemExit('no valid data in ' + path)
        poly = unary_union(geoms)
        gdf = gpd.GeoDataFrame(geometry=[poly], crs=src.crs).to_crs(crs)
        full = shapely_box(*src.bounds)
        print('footprint derived from raster: {:.1f}% of the raster is valid '
              'data'.format(100 * valid_frac))
        if valid_frac > 0.995:
            print('  !! WARNING: essentially the whole raster reads as valid, '
                  'so the footprint is just its bounding rectangle.')
            print('     That happens when nodata is unset on a gap-filled '
                  'mosaic. If the survey outline is not a rectangle, this '
                  'footprint is wrong and will pass nodes you never imaged.')
            print('     Check the nodata value, or supply a footprint polygon.')
    return gdf
 
 
def node_polygons(nodes, cross_frac, min_half_width):
    """
    Build one rectangle per SWORD node.
 
    Flow direction at a node is taken from its neighbours within the same
    reach (central difference; one-sided at reach ends), after sorting by
    node_id -- SWORD node_ids increase monotonically along a reach.
 
    cross_frac scales SWORD max_width into a cross-channel half-width.
    RiverObs's own cross-channel search half-distance is max_width/3
    (RiverObs.get_ext_dist_threshold), so cross_frac=0.5 is deliberately a
    little more conservative than the processor: it asks for a slightly wider
    strip to be imaged than the processor will actually search.
    """
    polys, keep_idx = [], []
    for reach_id, grp in nodes.groupby('reach_id', sort=False):
        grp = grp.sort_values('node_id')
        x = grp.geometry.x.to_numpy()
        y = grp.geometry.y.to_numpy()
        n = len(grp)
        if n < 2:
            continue
 
        dx = np.gradient(x) if n > 2 else np.full(n, x[-1] - x[0])
        dy = np.gradient(y) if n > 2 else np.full(n, y[-1] - y[0])
        norm = np.hypot(dx, dy)
        norm[norm == 0] = 1.0
        ux, uy = dx / norm, dy / norm          # along flow
        px, py = -uy, ux                       # across flow
 
        node_len = grp['node_len'].to_numpy(dtype=float)
        max_w = grp['max_width'].to_numpy(dtype=float)
        half_w = np.maximum(max_w * cross_frac, min_half_width)
        half_l = node_len / 2.0
 
        for i in range(n):
            cx, cy = x[i], y[i]
            a = (ux[i] * half_l[i], uy[i] * half_l[i])
            b = (px[i] * half_w[i], py[i] * half_w[i])
            polys.append(Polygon([
                (cx - a[0] - b[0], cy - a[1] - b[1]),
                (cx + a[0] - b[0], cy + a[1] - b[1]),
                (cx + a[0] + b[0], cy + a[1] + b[1]),
                (cx - a[0] + b[0], cy - a[1] + b[1]),
            ]))
            keep_idx.append(grp.index[i])
 
    return gpd.GeoDataFrame(
        nodes.loc[keep_idx].drop(columns='geometry'),
        geometry=polys, crs=nodes.crs)
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('sword_nodes_shp',
                   help='e.g. na_sword_nodes_hb81_v17b.shp')
    p.add_argument('footprint',
                   help='orthomosaic footprint: either a polygon file (.shp/'
                        '.gpkg/.geojson) or the orthomosaic raster itself '
                        '(.tif), in which case the footprint is derived from '
                        'its valid-data mask')
    p.add_argument('out_csv')
    p.add_argument('--crs', default='EPSG:32606')
    p.add_argument('--cross-frac', type=float, default=0.5,
                   help='cross-channel half-width as a fraction of SWORD '
                        'max_width (default 0.5)')
    p.add_argument('--min-half-width', type=float, default=50.0,
                   help='floor on the cross-channel half-width in metres, for '
                        'nodes with a small or missing SWORD max_width '
                        '(default 50)')
    p.add_argument('--allow-watermask', action='store_true',
                   help='suppress the guard that refuses a footprint which '
                        'looks like a water mask rather than a survey outline')
    p.add_argument('--full-threshold', type=float, default=0.99,
                   help='imaged_frac at or above which a node counts as fully '
                        'imaged (default 0.99)')
    p.add_argument('--write-polygons', default=None,
                   help='optional shapefile of the coverage-test RECTANGLES, '
                        'for QGIS QC only. These are NOT the node footprints '
                        'RiverObs uses and must not be shown as such -- see '
                        'step3b_true_node_polygons.py for the real ones.')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), SCRIPT_VERSION))
 
    nodes = gpd.read_file(args.sword_nodes_shp).to_crs(args.crs)
    foot = load_footprint(args.footprint, args.crs)
    foot_union = foot.union_all() if hasattr(foot, 'union_all') else foot.unary_union
 
    # Only nodes anywhere near the footprint are worth testing.
    nodes = nodes[nodes.geometry.within(foot_union.buffer(5000))].copy()
    if nodes.empty:
        raise SystemExit('no SWORD nodes near this footprint -- check the CRS')
 
    missing = [c for c in ('node_id', 'reach_id', 'node_len', 'max_width')
               if c not in nodes.columns]
    if missing:
        raise SystemExit(
            'SWORD node shapefile is missing {}. Field names differ slightly '
            'between SWORD versions -- check with '
            'geopandas.read_file(...).columns and adjust.'.format(missing))
 
    poly = node_polygons(nodes, args.cross_frac, args.min_half_width)
    inter = poly.geometry.intersection(foot_union)
    poly['imaged_frac'] = (inter.area / poly.geometry.area).clip(0, 1)
    poly['fully_imaged'] = poly['imaged_frac'] >= args.full_threshold
 
    # Guard: a water mask passed as a footprint gives
    #   imaged_frac ~ channel_width / (2 * cross_frac * max_width)
    # for EVERY node -- the gate silently becomes a filter on the very width
    # being validated. Discriminate on the MAXIMUM, not the median: a real
    # survey footprint fully contains at least some interior nodes, so its max
    # is ~1.0, while a water mask can never reach 1.0 anywhere. Using the
    # median would wrongly flag a legitimate per-survey footprint tested
    # against a node set spanning the whole domain, where most nodes are
    # outside the survey and score 0.
    top = float(poly['imaged_frac'].max())
    overlapping = poly[poly['imaged_frac'] > 0]
    if top < 0.9 and not args.allow_watermask:
        raise SystemExit(
            'REFUSING TO CONTINUE.\n'
            '  No node is more than {:.1%} covered by this footprint.\n'
            '  A real survey footprint fully contains its interior nodes, so\n'
            '  the maximum should be ~1.0. A ceiling like this means the\n'
            '  footprint is a water mask or another thin feature, not the\n'
            '  imaged area. Using it would turn the coverage gate into a\n'
            '  filter on channel width -- the exact quantity you are\n'
            '  validating.\n'
            '  Pass the orthomosaic .tif instead, or a real footprint polygon.\n'
            '  Override with --allow-watermask only if you know why.'.format(top))
 
    print('{} nodes overlap the footprint at all; of those, {} are fully '
          'imaged'.format(len(overlapping),
                          int(overlapping['fully_imaged'].sum())))
    partial = overlapping[~overlapping['fully_imaged']]
    if len(partial):
        print('{} partially imaged (imaged_frac {:.2f}-{:.2f}) -- these are the '
              'ones that would bias widths low if kept'.format(
                  len(partial), partial['imaged_frac'].min(),
                  partial['imaged_frac'].max()))
 
    out = poly[['node_id', 'reach_id', 'imaged_frac', 'fully_imaged']].copy()
    out.to_csv(args.out_csv, index=False)
 
    n_full = int(out['fully_imaged'].sum())
    print('{} nodes tested, {} fully imaged ({:.1f}%)'.format(
        len(out), n_full, 100.0 * n_full / len(out)))
    print('wrote', args.out_csv)
 
    if args.write_polygons:
        poly.to_file(args.write_polygons)
        print('wrote', args.write_polygons,
              '-- overlay this on the orthomosaic in QGIS before trusting it')
 
 
if __name__ == '__main__':
    main()