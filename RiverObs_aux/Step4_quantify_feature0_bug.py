#!/usr/bin/env python
"""
=============================================================================
STEP 4a -- How much water did the ['features'][0] bug drop?
-----------------------------------------------------------------------------
Script 3.2 does:
 
    YR_water_geom = [json.loads(YR_manual_water_class.to_json())['features'][0]['geometry']]
 
`['features'][0]` takes only the FIRST row of the shapefile. Any survey whose
water mask is stored as more than one feature had everything after the first
silently set to nodata before the per-node pixel count, so its water area --
and therefore its orthomosaic width -- is too small.
 
You have two affected surveys: sheenjek_240725 (2 features) and
yukonUS_240710 (17 features).
 
This script quantifies it. For each mask it reports the area of feature 0, the
total area of all features, and what fraction was lost -- both overall and
restricted to the AOI, since water outside YR_AOI.shp was clipped away anyway
and never counted either way.
 
It changes nothing. It just tells you how big a correction you are looking at.
 
    python step4a_quantify_feature0_bug.py \\
        --aoi /Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/GIS/YR_AOI.shp \\
        /path/to/width_validation/*.shp
 
Requires: geopandas, shapely, pandas
=============================================================================
"""
 
import argparse
import os
 
import geopandas as gpd
import pandas as pd
 
 
def summarise(path, crs, aoi_geom):
    gdf = gpd.read_file(path).to_crs(crs)
    n = len(gdf)
 
    # Explode to count the actual polygon count, since one row may be a
    # MultiPolygon -- that distinction is exactly what makes n features
    # misleading on its own.
    exploded = gdf.explode(index_parts=False, ignore_index=True)
 
    area_all = float(gdf.geometry.area.sum())
    area_f0 = float(gdf.geometry.iloc[0].area)
 
    row = {
        'file': os.path.basename(path),
        'n_features': n,
        'n_polygons': len(exploded),
        'geom_type_f0': gdf.geometry.iloc[0].geom_type,
        'area_all_m2': area_all,
        'area_feature0_m2': area_f0,
        'area_dropped_m2': area_all - area_f0,
        'pct_dropped': 100.0 * (area_all - area_f0) / area_all if area_all else 0.0,
    }
 
    if aoi_geom is not None:
        clipped_all = gdf.geometry.intersection(aoi_geom)
        clipped_f0 = gdf.geometry.iloc[0].intersection(aoi_geom)
        a_all = float(clipped_all.area.sum())
        a_f0 = float(clipped_f0.area)
        row['area_all_in_aoi_m2'] = a_all
        row['area_f0_in_aoi_m2'] = a_f0
        row['pct_dropped_in_aoi'] = (
            100.0 * (a_all - a_f0) / a_all if a_all else 0.0)
 
    # Size of each dropped feature, largest first -- tells you whether the
    # missing water is one big side channel or a scatter of small ponds.
    if n > 1:
        others = gdf.geometry.iloc[1:].area.sort_values(ascending=False)
        row['dropped_feature_areas_m2'] = ', '.join(
            '{:.0f}'.format(a) for a in others.head(8))
        if n > 9:
            row['dropped_feature_areas_m2'] += ', ... ({} more)'.format(n - 9)
    return row
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('masks', nargs='+', help='water mask shapefiles')
    p.add_argument('--aoi', default=None,
                   help='YR_AOI.shp, to also report AOI-clipped areas')
    p.add_argument('--crs', default='EPSG:32606')
    p.add_argument('--out-csv', default=None)
    args = p.parse_args()
 
    aoi_geom = None
    if args.aoi:
        aoi = gpd.read_file(args.aoi).to_crs(args.crs)
        print('AOI: {} feature(s) in {}'.format(len(aoi), os.path.basename(args.aoi)))
        if len(aoi) > 1:
            print('  !! YR_AOI.shp has more than one feature, so script 3.2\'s')
            print('     YR_AOI_geom = [...][\'features\'][0] dropped some of the')
            print('     AOI as well. Every survey is affected, not just two.')
        aoi_geom = aoi.union_all() if hasattr(aoi, 'union_all') else aoi.unary_union
        # feature 0 only, matching what script 3.2 actually used
        aoi_geom_f0 = aoi.geometry.iloc[0]
        if len(aoi) > 1:
            print('  AOI area, all features: {:,.0f} m2'.format(aoi_geom.area))
            print('  AOI area, feature 0   : {:,.0f} m2'.format(aoi_geom_f0.area))
        aoi_geom = aoi_geom_f0  # what 3.2 used, so the comparison is like-for-like
        print()
 
    rows = [summarise(m, args.crs, aoi_geom) for m in args.masks]
    df = pd.DataFrame(rows)
 
    pd.set_option('display.width', 200)
    pd.set_option('display.max_colwidth', 60)
 
    cols = ['file', 'n_features', 'n_polygons', 'area_all_m2',
            'area_feature0_m2', 'pct_dropped']
    if 'pct_dropped_in_aoi' in df:
        cols.append('pct_dropped_in_aoi')
    print(df[cols].to_string(index=False,
                             float_format=lambda v: '{:,.1f}'.format(v)))
 
    affected = df[df['n_features'] > 1]
    print('\n' + '=' * 70)
    if affected.empty:
        print('No survey affected -- every mask is a single feature.')
    else:
        print('AFFECTED SURVEYS')
        for _, r in affected.iterrows():
            key = 'pct_dropped_in_aoi' if 'pct_dropped_in_aoi' in r else 'pct_dropped'
            print('\n  {}'.format(r['file']))
            print('    {} features, {:.1f}% of water area dropped{}'.format(
                int(r['n_features']), r[key],
                ' (within AOI)' if key.endswith('aoi') else ''))
            if 'dropped_feature_areas_m2' in r and pd.notna(
                    r.get('dropped_feature_areas_m2')):
                print('    dropped feature areas (m2, largest first): {}'.format(
                    r['dropped_feature_areas_m2']))
        print('\nInterpretation: dropped water means the orthomosaic width for')
        print('these surveys is too SMALL, which makes SWOT look like it')
        print('overestimates more than it does. Both affected rivers feed')
        print('results you report by river, so the correction is not cosmetic.')
    print('=' * 70)
 
    if args.out_csv:
        df.to_csv(args.out_csv, index=False)
        print('\nwrote', args.out_csv)
 
 
if __name__ == '__main__':
    main()