#!/usr/bin/env python
"""Turn manual problem node review into a node exclusion list.

    python 3.1.4_manual_node_qc.py --batch \\
        --nodes-dir   /path/output/riverobs/node_qc \\
        --exclude-shp /path/digitizing/bad_ends_watermasks.shp \\
        --versions v17b v16 \\
        --out /path/output/riverobs/node_qc/node_qc_all.csv

Any node whose water polygon overlaps the exclusion area by more
than --max-overlap (default 5%) is dropped.

Requires: geopandas, pandas, shapely

-----------------------------------------------------------------------------
Script by:
Camryn Kluetmeier (camryn.kluetmeier@duke.edu)

Parts of this script were developed with assistance from Claude Code
(Anthropic) for debugging, documentation, and related editorial suggestions.

Last updated: 2026-09-13

"""

PIPELINE_VERSION = '1.1.0'

import argparse
import glob
import os
import re

import pandas as pd
import geopandas as gpd

NODE_GLOB = 'true_node_polygons_*_{version}.shp'
NODE_RE = r'^true_node_polygons_(.+)_{version}\.shp$'
ALL_TOKENS = {'', 'all', 'ALL', 'All', 'none', 'None'}


# Exclusion layer
def load_exclusion(path, field=None):
    """Read the exclusion layer once and report what is in it."""
    ex = gpd.read_file(path)
    if ex.empty:
        raise SystemExit('exclusion layer is empty: ' + path)
    attrs = [c for c in ex.columns if c != 'geometry']
    print('exclusion layer: {} feature(s) from {}'.format(len(ex), path))
    print('  fields: {}'.format(', '.join(attrs) if attrs else '(none)'))
    if field:
        if field not in ex.columns:
            raise SystemExit(
                "--exclude-field '{}' is not a field in {}. Available: {}"
                .format(field, path, ', '.join(attrs) or '(none)'))
        print("  per-survey mode: polygons are matched on '{}'".format(field))
    elif len(ex) > 1:
        print('  every polygon is applied to every survey (no --exclude-field). '
              'That is correct when the polygons mark ground that is bad in all '
              'surveys; see the header if two surveys overlap on the ground.')
    return ex


def exclusion_geom(ex, field=None, survey=None):
    """Union of the polygons that apply to `survey`. None if none apply."""
    sub = ex
    if field and survey is not None:
        vals = ex[field].astype('string').fillna('')
        sub = ex[(vals == survey) | (vals.isin(ALL_TOKENS))]
        if sub.empty:
            return None
    return sub.union_all() if hasattr(sub, 'union_all') else sub.unary_union


def from_exclusion_polygon(nodes, ex_union, max_overlap, quiet=False):
    nodes = nodes.copy()
    if ex_union is None:
        nodes['excl_frac'] = 0.0
        nodes['keep'] = 1
        if not quiet:
            print('  no exclusion polygon applies to this survey -- all nodes kept')
        return nodes

    inter = nodes.geometry.intersection(ex_union)
    frac = (inter.area / nodes.geometry.area).fillna(0.0).clip(0, 1)
    nodes['excl_frac'] = frac
    nodes['keep'] = (frac <= max_overlap).astype(int)

    touched = int((frac > 0).sum())
    dropped = int((nodes['keep'] == 0).sum())
    print('  {} node(s) touch the exclusion area; {} exceed the {:.0%} threshold '
          'and are dropped'.format(touched, dropped, max_overlap))
    partial = nodes[(frac > 0) & (frac <= max_overlap)]
    if len(partial) and not quiet:
        print('  {} node(s) overlap but stayed under the threshold -- check '
              'these if you want to be strict:'.format(len(partial)))
        for _, r in partial.sort_values(
                'excl_frac', ascending=False).head(8).iterrows():
            print('     {}  {:.1%} inside the exclusion area'.format(
                int(r['node_id']), r['excl_frac']))
    return nodes


def from_keep_field(nodes, original_shp):
    if 'keep' not in nodes.columns:
        raise SystemExit(
            "no `keep` column in the node layer. Either regenerate it with "
            "3.1.3_true_node_polygons.py, or use --exclude-shp.")
    nodes = nodes.copy()
    nodes['keep'] = nodes['keep'].fillna(1).astype(int)

    if original_shp:
        orig = gpd.read_file(original_shp)
        missing = set(orig['node_id']) - set(nodes['node_id'])
        if missing:
            print('  {} node(s) were deleted from the edited layer -- recording '
                  'them as dropped'.format(len(missing)))
            extra = pd.DataFrame({'node_id': sorted(missing), 'keep': 0})
            nodes = pd.concat(
                [pd.DataFrame(nodes.drop(columns='geometry')), extra],
                ignore_index=True)
    print('  {} of {} nodes kept'.format(int((nodes['keep'] == 1).sum()),
                                         len(nodes)))
    return nodes


# Shared
def tidy(result, survey, sword_version):
    result = pd.DataFrame(result).drop(columns='geometry', errors='ignore')
    cols = ['node_id', 'keep']
    if 'excl_frac' in result.columns:
        cols.append('excl_frac')
    result = result[cols].copy()
    if sword_version:
        result.insert(0, 'sword_version', sword_version)
    if survey:
        result.insert(0, 'survey', survey)
    return result


def discover(nodes_dir, versions, surveys=None):
    """Every (survey, version, shapefile) under nodes_dir, sorted."""
    found = []
    for v in versions:
        pattern = os.path.join(nodes_dir, NODE_GLOB.format(version=v))
        rx = re.compile(NODE_RE.format(version=re.escape(v)))
        for path in sorted(glob.glob(pattern)):
            m = rx.match(os.path.basename(path))
            if not m:
                continue
            name = m.group(1)
            if surveys and name not in surveys:
                continue
            found.append((name, v, path))
    return found


# Main
def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--nodes',
                   help='node polygon shapefile from '
                        '3.1.3_true_node_polygons.py')
    p.add_argument('--exclude-shp', default=None,
                   help='polygon layer marking cloud / bad ends / bad areas')
    p.add_argument('--exclude-field', default=None,
                   help='field in the exclusion layer naming the survey each '
                        'polygon belongs to; omit to apply every polygon to '
                        'every survey')
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
    p.add_argument('--out', required=True,
                   help='output CSV; in --batch mode this is the merged file')

    p.add_argument('--batch', action='store_true',
                   help='process every survey and version found in --nodes-dir')
    p.add_argument('--nodes-dir', default=None,
                   help='directory of true_node_polygons_*.shp written by '
                        '3.1.3_true_node_polygons.py (batch)')
    p.add_argument('--versions', nargs='+', default=['v17b'],
                   help='prior-database tags to process in batch mode')
    p.add_argument('--surveys', nargs='+', default=None,
                   help='restrict batch mode to these survey names')
    p.add_argument('--out-dir', default=None,
                   help='where the per-survey CSVs go (default: alongside --out)')

    p.add_argument('--merge', nargs='*', default=None,
                   help='combine existing per-survey CSVs into one and exit')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))

    # Merge only
    if args.merge is not None:
        paths = []
        for m in args.merge:
            paths.extend(sorted(glob.glob(m)) or [m])
        out = pd.concat([pd.read_csv(f) for f in paths], ignore_index=True)
        out.to_csv(args.out, index=False)
        print('merged {} file(s), {} rows -> {}'.format(
            len(paths), len(out), args.out))
        print('   kept {}, dropped {}'.format(
            int((out['keep'] == 1).sum()), int((out['keep'] == 0).sum())))
        return

    # Batch
    if args.batch:
        if not args.nodes_dir:
            raise SystemExit('--batch needs --nodes-dir')
        if not args.exclude_shp:
            raise SystemExit(
                '--batch needs --exclude-shp. The `keep`-field route is '
                'per-layer by nature and has no batch form.')

        jobs = discover(args.nodes_dir, args.versions, set(args.surveys or []))
        if not jobs:
            raise SystemExit(
                'no true_node_polygons_*_<version>.shp found in {} for '
                'version(s) {}. Run step 3b first.'.format(
                    args.nodes_dir, ' '.join(args.versions)))

        ex = load_exclusion(args.exclude_shp, args.exclude_field)
        out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.out))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        print('\n{} survey/version pair(s) to process'.format(len(jobs)))
        frames, summary, failed = [], [], []
        for survey, version, path in jobs:
            print('\n--- {} / {} ---'.format(survey, version))
            try:
                nodes = gpd.read_file(path)
                if 'node_id' not in nodes.columns:
                    raise SystemExit('no node_id column in ' + path)
                geom = exclusion_geom(ex.to_crs(nodes.crs),
                                      args.exclude_field, survey)
                res = tidy(from_exclusion_polygon(
                    nodes, geom, args.max_overlap), survey, version)
            except Exception as exc:          # keep going; report at the end
                print('  !!! FAILED: {}: {}'.format(type(exc).__name__, exc))
                failed.append('{}/{}'.format(survey, version))
                continue

            per = os.path.join(
                out_dir, 'node_qc_{}_{}.csv'.format(survey, version))
            res.to_csv(per, index=False)
            frames.append(res)
            summary.append((survey, version, len(res),
                            int((res['keep'] == 1).sum()),
                            int((res['keep'] == 0).sum())))
            print('  wrote {}'.format(per))

        if not frames:
            raise SystemExit('nothing succeeded -- see the FAILED lines above')

        merged = pd.concat(frames, ignore_index=True)
        merged.to_csv(args.out, index=False)

        print('\n' + '=' * 66)
        print('{:<20} {:<6} {:>7} {:>7} {:>8}'.format(
            'survey', 'sword', 'nodes', 'kept', 'dropped'))
        for row in summary:
            print('{:<20} {:<6} {:>7} {:>7} {:>8}'.format(*row))
        print('-' * 66)
        print('{:<20} {:<6} {:>7} {:>7} {:>8}'.format(
            'TOTAL', '', len(merged),
            int((merged['keep'] == 1).sum()),
            int((merged['keep'] == 0).sum())))
        print('=' * 66)
        print('merged -> {}'.format(args.out))
        if failed:
            print('FAILED: {}'.format(', '.join(failed)))
            raise SystemExit(1)
        return

    # Single survey/version
    if not args.nodes:
        raise SystemExit('--nodes is required unless you are using --batch or '
                         '--merge')

    nodes = gpd.read_file(args.nodes)
    if 'node_id' not in nodes.columns:
        raise SystemExit('no node_id column in ' + args.nodes)

    if args.exclude_shp:
        ex = load_exclusion(args.exclude_shp, args.exclude_field)
        geom = exclusion_geom(ex.to_crs(nodes.crs), args.exclude_field,
                              args.survey)
        result = from_exclusion_polygon(nodes, geom, args.max_overlap)
    else:
        result = from_keep_field(nodes, args.original)

    result = tidy(result, args.survey, args.sword_version)

    if args.exclude_ids:
        ids = {int(x) for x in args.exclude_ids.replace(' ', '').split(',') if x}
        n_before = int((result['keep'] == 1).sum())
        result.loc[result['node_id'].isin(ids), 'keep'] = 0
        print('--exclude-ids dropped a further {} node(s)'.format(
            n_before - int((result['keep'] == 1).sum())))

    result.to_csv(args.out, index=False)
    print('wrote {} rows to {}  ({} kept, {} dropped)'.format(
        len(result), args.out,
        int((result['keep'] == 1).sum()), int((result['keep'] == 0).sum())))


if __name__ == '__main__':
    main()