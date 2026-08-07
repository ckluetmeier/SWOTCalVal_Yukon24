
#!/usr/bin/env python
"""
=============================================================================
STEP 2 -- Run every orthomosaic water mask through RiverObs, both SWORD versions
-----------------------------------------------------------------------------
For each (water mask, SWORD version) pair this calls the patched
calval2rivertile.py, then flattens the resulting RiverTile netCDF into two
CSVs -- one node-level, one reach-level -- ready for downstream analysis.
 
Run it from inside the RiverObs conda environment. Keep it in the same folder
as run_calval2rivertile.py and riverobs_shim.py.
 
    python step2_run_riverobs.py --config surveys.yml     # or edit SURVEYS below
 
Requires: netCDF4, pandas, numpy (all already in the RiverObs environment)
=============================================================================
"""
 
PIPELINE_VERSION = '1.0.0'
 
 
import argparse
import os
import subprocess
import sys
 
import numpy as np
import netCDF4
import pandas as pd
 
 
# -----------------------------------------------------------------------------
# EDIT THIS BLOCK -- one entry per orthomosaic, mirroring the commented file
# `date` is the matching SWOT overpass date (UTC) and is carried through to the
# output CSV so the downstream temporal match is a join condition rather than a
# hand-edited constant.
# -----------------------------------------------------------------------------
SURVEYS = [
    {'name': 'CD_071024',        'date': '2024-07-11'},
    {'name': 'upperPR_CL_071024', 'date': '2024-07-10'},
    {'name': 'upperPR_CL_071624', 'date': '2024-07-16'},
    {'name': 'lowerPR_SJ_072624', 'date': '2024-07-26'},
    {'name': 'lowerYR_071624',   'date': '2024-07-16'},
    {'name': 'upperYR_071024',   'date': '2024-07-10'},
]
 
SWORD_VERSIONS = {
    'v16':  'riverobs_ortho_v16.rdf',   # pair with RiverSP version C / PIC0
    'v17b': 'riverobs_ortho_v17b.rdf',  # pair with RiverSP version D / PGD0
}
 
# Node fields worth keeping. `width` is the ortho width: RiverObs computed it
# as node area / SWORD p_length using the same RiverNode.width_area() call the
# operational processor uses for RiverSP.
NODE_FIELDS = [
    'reach_id', 'node_id', 'lat', 'lon',
    'width', 'width_u', 'area_total', 'area_detct',
    'n_good_pix', 'p_length', 'p_width', 'p_dist_out',
]
 
# Reach fields. `width` here is sum(node area) / sum(node p_length) over the
# observed nodes -- a length-weighted mean, not the unweighted mean of node
# widths. `obs_frac_n` is the fraction of the reach's SWORD nodes that received
# data, and is the defensible basis for a reach-completeness threshold.
REACH_FIELDS = [
    'reach_id', 'width', 'width_u', 'area_total', 'area_detct',
    'obs_frac_n', 'partial_f', 'n_good_nod', 'p_length', 'p_width',
]
 
 
def flatten_group(nc_path, group, fields):
    """Read one group of a RiverTile netCDF into a DataFrame."""
    out = {}
    with netCDF4.Dataset(nc_path, 'r') as ds:
        grp = ds.groups[group]
        for f in fields:
            if f not in grp.variables:
                print('    note: {} not in /{} -- skipping'.format(f, group))
                continue
            arr = grp.variables[f][:]
            if isinstance(arr, np.ma.MaskedArray):
                arr = arr.filled(np.nan) if arr.dtype.kind == 'f' else arr.filled(-999)
            out[f] = np.asarray(arr).ravel()
    return pd.DataFrame(out)
 
 
def run_one(survey, version, rdf, args):
    name = survey['name']
    res_tag = ('{:g}'.format(args.res)).replace('.', 'p')
    water_tif = os.path.join(
        args.raster_dir, '{}_water_{}m.tif'.format(name, res_tag))
    if not os.path.exists(water_tif):
        # FileNotFoundError, not SystemExit -- SystemExit does not inherit from
        # Exception, so --keep-going could not catch it and the whole batch died
        # on the first missing raster without writing any CSV.
        raise FileNotFoundError('missing water raster: ' + water_tif)
 
    out_dir = os.path.join(args.out_dir, version)
    os.makedirs(out_dir, exist_ok=True)
    rivertile_nc = os.path.join(out_dir, 'ortho_{}_{}.nc'.format(name, version))
    pixcvec_nc = os.path.join(out_dir, 'ortho_{}_{}_pixcvec.nc'.format(name, version))
    pixc_nc = os.path.join(out_dir, 'ortho_{}_{}_pixc.nc'.format(name, version))
 
    # run_calval2rivertile.py, not src/bin/calval2rivertile.py -- the wrapper
    # loads riverobs_shim first. See riverobs_shim.py for what it patches.
    cmd = [
        sys.executable,
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'run_calval2rivertile.py'),
        water_tif, 'airborne_watermask',
        rivertile_nc, pixcvec_nc, rdf, pixc_nc,
        '--riverobs-root', args.riverobs_root,
        '--log-level', args.log_level,
    ]
    print('\n=== {} / SWORD {} ==='.format(name, version))
    print(' '.join(cmd))
    subprocess.run(cmd, check=True)
 
    nodes = flatten_group(rivertile_nc, 'nodes', NODE_FIELDS)
    reaches = flatten_group(rivertile_nc, 'reaches', REACH_FIELDS)
 
    for df in (nodes, reaches):
        df.insert(0, 'survey', name)
        df.insert(1, 'SWOT_date', survey['date'])
        df.insert(2, 'sword_version', version)
 
    # RiverObs writes 'width'; rename so downstream code cannot confuse the
    # reference width with the SWOT width column.
    nodes = nodes.rename(columns={
        'width': 'ortho_width_m', 'width_u': 'ortho_width_u_m',
        'area_total': 'ortho_area_total_m2', 'area_detct': 'ortho_area_detct_m2'})
    reaches = reaches.rename(columns={
        'width': 'ortho_width_m', 'width_u': 'ortho_width_u_m',
        'area_total': 'ortho_area_total_m2', 'area_detct': 'ortho_area_detct_m2'})
 
    print('    {} nodes, {} reaches'.format(len(nodes), len(reaches)))
    return nodes, reaches
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--riverobs-root', required=True,
                   help='root of your RiverObs fork (the dir containing src/)')
    p.add_argument('--raster-dir', required=True,
                   help='output directory from step1_rasterize_watermask.py')
    p.add_argument('--rdf-dir', required=True,
                   help='directory holding riverobs_ortho_v16.rdf / _v17b.rdf')
    p.add_argument('--out-dir', required=True)
    p.add_argument('--res', type=float, default=3.0,
                   help='must match the --res used in step 1')
    p.add_argument('--versions', nargs='+', default=['v16', 'v17b'])
    p.add_argument('--survey', nargs='+', default=None,
                   help='run only these survey names (default: all in SURVEYS)')
    p.add_argument('--keep-going', action='store_true',
                   help='carry on if a survey fails, and still write the CSVs '
                        'for the ones that worked. Failures are listed at the '
                        'end and the exit status is non-zero.')
    p.add_argument('--log-level', default='info')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))
 
    surveys = SURVEYS
    if args.survey:
        wanted = set(args.survey)
        surveys = [s for s in SURVEYS if s['name'] in wanted]
        missing = wanted - {s['name'] for s in surveys}
        if missing:
            raise SystemExit('unknown survey name(s): {}. Known: {}'.format(
                ', '.join(sorted(missing)),
                ', '.join(s['name'] for s in SURVEYS)))
 
    all_nodes, all_reaches, failures = [], [], []
    for version in args.versions:
        rdf = os.path.join(args.rdf_dir, SWORD_VERSIONS[version])
        for survey in surveys:
            try:
                nodes, reaches = run_one(survey, version, rdf, args)
            except (Exception, SystemExit) as exc:
                if not args.keep_going:
                    raise
                failures.append((survey['name'], version, str(exc).split(chr(10))[0]))
                print('    FAILED: {} / {} -- {}'.format(
                    survey['name'], version, type(exc).__name__))
                continue
            all_nodes.append(nodes)
            all_reaches.append(reaches)
 
    if not all_nodes:
        raise SystemExit('\nevery run failed -- nothing to write')
 
    os.makedirs(args.out_dir, exist_ok=True)
    node_csv = os.path.join(args.out_dir, 'ortho_riverobs_nodes_all.csv')
    reach_csv = os.path.join(args.out_dir, 'ortho_riverobs_reaches_all.csv')
    pd.concat(all_nodes, ignore_index=True).to_csv(node_csv, index=False)
    pd.concat(all_reaches, ignore_index=True).to_csv(reach_csv, index=False)
    print('\nwrote', node_csv)
    print('wrote', reach_csv)
 
    if failures:
        print('\n{} run(s) FAILED and are absent from those CSVs:'.format(
            len(failures)))
        for name, version, msg in failures:
            print('   {} / {}'.format(name, version))
        print('Re-run just those with --survey <name> once fixed.')
        sys.exit(1)
 
 
if __name__ == '__main__':
    main()