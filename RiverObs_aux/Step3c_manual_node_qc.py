
#!/usr/bin/env python
"""
=============================================================================
STEP 3c -- Turn your manual QGIS review into a node exclusion list
-----------------------------------------------------------------------------
There is no automated coverage gate in this workflow. A footprint derived from
a mosaic's nodata mask cannot see cloud -- cloud is valid data -- so anything
automated would pass exactly the nodes you most want to drop. You make the call
in QGIS instead; this script records it.
 
TWO WAYS TO WORK. The first is better.
 
--- A. Draw an exclusion polygon (recommended) --------------------------------
In QGIS, over the orthomosaic, digitize one polygon layer per survey covering
the areas you do not trust: cloud, haze, the ragged survey edge, anything else.
Then:
 
    python step3c_manual_node_qc.py \\
        --nodes  /path/true_node_polygons_upperYR_071024.shp \\
        --exclude-shp /path/badareas_upperYR_071024.shp \\
        --survey upperYR_071024 --sword-version v17b \\
        --out /path/node_qc_upperYR_071024_v17b.csv
 
Any node whose water polygon overlaps the exclusion area by more than
--max-overlap (default 5%) is dropped.
 
Why this is better: the exclusion is *spatial*, so one act of judgement per
survey serves EVERY prior-database version. Different SWORD versions place
nodes differently and give them different ids, but the cloud is in the same
place on the ground. Running this against each version's node polygons using
the SAME exclusion polygon filters them consistently. Reviewing each version's
node list separately by eye would not guarantee that, and an inconsistent
filter between versions contaminates any cross-version comparison.
 
It is also the reproducible artefact: a reader can see exactly which ground
area was excluded, rather than an unexplained list of node identifiers.
 
--- B. Edit the `keep` field ---------------------------------------------------
step3b writes a `keep` column, set to 1 everywhere. Open the layer in QGIS,
toggle editing, set `keep` = 0 on the nodes you want gone, save, then:
 
    python step3c_manual_node_qc.py \\
        --nodes /path/true_node_polygons_upperYR_071024.shp \\
        --survey upperYR_071024 --sword-version v17b \\
        --out /path/node_qc_upperYR_071024_v17b.csv
 
Deleting rows outright works too -- pass --original pointing at an untouched
copy so the script knows what the full set was.
 
--- Merging ---------------------------------------------------------------
Once every survey and version is done:
 
    python step3c_manual_node_qc.py --merge /path/node_qc_*.csv \\
        --out /path/node_qc_all.csv
 
The combined file is what downstream analysis joins against.
 
Requires: geopandas, pandas, shapely
=============================================================================
"""
 
PIPELINE_VERSION = '1.0.0'
 
import argparse
import glob
import os
 
import pandas as pd
import geopandas as gpd
 
 
def from_exclusion_polygon(nodes, exclude_shp, max_overlap):
    ex = gpd.read_file(exclude_shp).to_crs(nodes.crs)
    if ex.empty:
        raise SystemExit('exclusion layer is empty: ' + exclude_shp)
    ex_union = ex.union_all() if hasattr(ex, 'union_all') else ex.unary_union
 
    inter = nodes.geometry.intersection(ex_union)
    frac = (inter.area / nodes.geometry.area).fillna(0.0).clip(0, 1)
    nodes = nodes.copy()
    nodes['excl_frac'] = frac
    nodes['keep'] = (frac <= max_overlap).astype(int)
 
    touched = int((frac > 0).sum())
    dropped = int((nodes['keep'] == 0).sum())
    print('{} node(s) touch the exclusion area; {} exceed the {:.0%} threshold '
          'and are dropped'.format(touched, dropped, max_overlap))
    partial = nodes[(frac > 0) & (frac <= max_overlap)]
    if len(partial):
        print('{} node(s) overlap but stayed under the threshold -- check these '
              'if you want to be strict:'.format(len(partial)))
        for _, r in partial.sort_values('excl_frac', ascending=False).head(8).iterrows():
            print('   {}  {:.1%} inside the exclusion area'.format(
                int(r['node_id']), r['excl_frac']))
    return nodes
 
 
def from_keep_field(nodes, original_shp):
    if 'keep' not in nodes.columns:
        raise SystemExit(
            "no `keep` column in the node layer. Either regenerate it with "
            "step3b (v2 or later), or use --exclude-shp.")
    nodes = nodes.copy()
    nodes['keep'] = nodes['keep'].fillna(1).astype(int)
 
    if original_shp:
        orig = gpd.read_file(original_shp)
        missing = set(orig['node_id']) - set(nodes['node_id'])
        if missing:
            print('{} node(s) were deleted from the edited layer -- recording '
                  'them as dropped'.format(len(missing)))
            extra = pd.DataFrame({'node_id': sorted(missing), 'keep': 0})
            nodes = pd.concat(
                [pd.DataFrame(nodes.drop(columns='geometry')), extra],
                ignore_index=True)
    print('{} of {} nodes kept'.format(int((nodes['keep'] == 1).sum()),
                                       len(nodes)))
    return nodes
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--nodes', help='node polygon shapefile from step3b')
    p.add_argument('--exclude-shp', default=None,
                   help='polygon layer marking cloud / bad areas (route A)')
    p.add_argument('--original', default=None,
                   help='untouched copy of the node layer, if you deleted rows '
                        'rather than setting keep = 0 (route B)')
    p.add_argument('--exclude-ids', default=None,
                   help='comma-separated node_ids to drop, for small lists')
    p.add_argument('--max-overlap', type=float, default=0.05,
                   help='fraction of a node polygon that may fall inside the '
                        'exclusion area before it is dropped (default 0.05)')
    p.add_argument('--survey', default=None)
    p.add_argument('--sword-version', default=None,
                   help='prior-database version tag, e.g. v17b')
    p.add_argument('--out', required=True)
    p.add_argument('--merge', nargs='*', default=None,
                   help='combine per-survey CSVs into one and exit')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))
 
    if args.merge is not None:
        paths = []
        for m in args.merge:
            paths.extend(sorted(glob.glob(m)) or [m])
        frames = [pd.read_csv(f) for f in paths]
        out = pd.concat(frames, ignore_index=True)
        out.to_csv(args.out, index=False)
        print('merged {} file(s), {} rows -> {}'.format(
            len(paths), len(out), args.out))
        print('   kept {}, dropped {}'.format(
            int((out['keep'] == 1).sum()), int((out['keep'] == 0).sum())))
        return
 
    if not args.nodes:
        raise SystemExit('--nodes is required unless you are using --merge')
 
    nodes = gpd.read_file(args.nodes)
    if 'node_id' not in nodes.columns:
        raise SystemExit('no node_id column in ' + args.nodes)
 
    if args.exclude_shp:
        result = from_exclusion_polygon(nodes, args.exclude_shp, args.max_overlap)
    else:
        result = from_keep_field(nodes, args.original)
 
    result = pd.DataFrame(result).drop(columns='geometry', errors='ignore')
 
    if args.exclude_ids:
        ids = {int(x) for x in args.exclude_ids.replace(' ', '').split(',') if x}
        n_before = int((result['keep'] == 1).sum())
        result.loc[result['node_id'].isin(ids), 'keep'] = 0
        print('--exclude-ids dropped a further {} node(s)'.format(
            n_before - int((result['keep'] == 1).sum())))
 
    cols = ['node_id', 'keep']
    if 'excl_frac' in result.columns:
        cols.append('excl_frac')
    result = result[cols].copy()
    if args.survey:
        result.insert(0, 'survey', args.survey)
    if args.sword_version:
        result.insert(1, 'sword_version', args.sword_version)
 
    result.to_csv(args.out, index=False)
    print('wrote {} rows to {}  ({} kept, {} dropped)'.format(
        len(result), args.out,
        int((result['keep'] == 1).sum()), int((result['keep'] == 0).sum())))
 
 
if __name__ == '__main__':
    main()