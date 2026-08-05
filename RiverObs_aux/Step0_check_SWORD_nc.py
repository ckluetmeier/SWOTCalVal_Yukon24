
#!/usr/bin/env python
"""
=============================================================================
STEP 0 -- Check a SWORD netCDF is usable by RiverObs
-----------------------------------------------------------------------------
    python step0_check_sword_nc.py /path/na_sword_v17b.nc --verify
 
Run once per SWORD version. It reads only -- it no longer writes a padded copy,
because padding turned out to be unnecessary.
 
*** CORRECTION TO EARLIER VERSIONS OF THIS SCRIPT ***
v1 and v2 said a SWORD file missing variables that RiverObs declares would
raise AttributeError, and offered to pad it. That was wrong. I inferred it from
reading ReachDatabaseReaches.__call__()'s getattr loop without testing that
Product.__getattr__ intercepts first. It does:
 
    src/SWOTWater/products/product.py, Product.__getattr__
        if key in self.VARIABLES:
            ...
            return np.ma.masked_array(..., mask=np.ones(shape))
 
Any variable RiverObs declares but your file lacks comes back as a fully-masked
array of the right shape. The shape comes from RiverObs's own DIMENSIONS table,
so your file's dimension names are irrelevant too. Nothing raises.
 
Verified by running the complete pipeline against a SWORD file missing the same
35 node and 63 reach variables yours is, unpadded: node widths identical to the
padded run, reach_type correct.
 
/reaches/type self-heals the same way. Missing -> masked -> RiverObs's own
documented fallback fires, at ReachDatabase.py line 228:
 
    reach_type = this_reach['reaches']['type'][0]
    if reach_type is np.ma.masked:
        reach_type = reach_idx % 10       # last digit of the SWORD reach_id
 
So: use your ORIGINAL SWORD netCDFs. Delete any *_riverobs.nc padded copies.
No methods caveat needed.
 
*** WHY THIS SCRIPT STILL MATTERS ***
Because missing variables are silently masked rather than loudly fatal, a SWORD
file missing something RiverObs actually *uses* would not crash -- it would
quietly produce wrong widths. node_length is the divisor in
width = area / node_length; max_width and ext_dist_coef size the cross-channel
search corridor. Masked stand-ins for those give numbers that look plausible
and are not.
 
The ASSIGNMENT_CRITICAL check below is a guard against silent corruption, not
against a crash. --verify additionally proves the file loads and that
node_length carries real values.
 
Requires: netCDF4, numpy, plus RiverObs on PYTHONPATH
=============================================================================
"""
 
import argparse
import sys
 
import numpy as np
import netCDF4
 
from RiverObs.ReachDatabase import (
    ReachDatabaseNodes, ReachDatabaseReaches, ReachDatabaseCenterlines)
 
GROUPS = {
    'nodes': ReachDatabaseNodes,
    'reaches': ReachDatabaseReaches,
    'centerlines': ReachDatabaseCenterlines,
}
 
# Variables RiverObs reads to build the centerline, size the search corridor,
# and convert area to width. `type` is NOT here -- missing is fine, RiverObs
# derives it from reach_id.
ASSIGNMENT_CRITICAL = {
    'nodes': ['node_id', 'reach_id', 'x', 'y', 'node_length', 'width',
              'max_width', 'wth_coef', 'ext_dist_coef', 'dist_out'],
    'reaches': ['reach_id', 'x', 'y', 'x_min', 'x_max', 'y_min', 'y_max',
                'reach_length', 'n_nodes', 'width', 'max_width'],
    'centerlines': ['x', 'y', 'reach_id', 'node_id', 'cl_id'],
}
 
# Subgroups ReachDatabaseReaches.__call__ dereferences directly. These are NOT
# covered by the masked-array fallback -- .subset()/.__call__() get invoked on
# them, so they must genuinely exist.
REQUIRED_SUBGROUPS = ['area_fits', 'discharge_models']
 
# Missing but self-healing. Named so the report isn't alarming.
SELF_HEALING = {
    ('reaches', 'type'):
        'RiverObs derives it as reach_id % 10 (ReachDatabase.py:228)',
}
 
 
def report(nc_path):
    critical, subgroup_problem, healing = {}, [], []
 
    with netCDF4.Dataset(nc_path, 'r') as ds:
        print('file  :', nc_path)
        print('groups:', list(ds.groups))
 
        for gname, klass in GROUPS.items():
            if gname not in ds.groups:
                print('\n!! group /{} is absent -- not a SWORD netCDF RiverObs '
                      'can read'.format(gname))
                critical[gname] = ['<entire group>']
                continue
 
            have = set(ds.groups[gname].variables)
            want = set(klass.VARIABLES)
            miss = sorted(want - have)
            critical[gname] = [v for v in miss
                               if v in ASSIGNMENT_CRITICAL[gname]]
            healing += [(gname, v) for v in miss if (gname, v) in SELF_HEALING]
 
            n_crit = len(ASSIGNMENT_CRITICAL[gname])
            print('\n/{}: {} of {} declared variables present'.format(
                gname, len(want) - len(miss), len(want)))
            print('     {} of {} assignment-critical present'.format(
                n_crit - len(critical[gname]), n_crit))
            if miss:
                print('     {} missing -> masked arrays, harmless'.format(
                    len(miss)))
            for gn, v in healing:
                if gn == gname:
                    print('     /{}/{} missing, self-healing: {}'.format(
                        gn, v, SELF_HEALING[(gn, v)]))
            if critical[gname]:
                print('  *** MISSING AND CRITICAL: {} ***'.format(
                    ', '.join(critical[gname])))
 
        if 'reaches' in ds.groups:
            subgroup_problem = [
                sg for sg in REQUIRED_SUBGROUPS
                if sg not in ds.groups['reaches'].groups]
            if subgroup_problem:
                print('\n!! /reaches is missing subgroup(s): {}. Unlike '
                      'variables, these are NOT covered by the masked-array '
                      'fallback.'.format(', '.join(subgroup_problem)))
            else:
                print('\n/reaches subgroups present: {}'.format(
                    ', '.join(REQUIRED_SUBGROUPS)))
 
    return critical, subgroup_problem
 
 
def verify(path):
    """
    Prove the file loads: open it as a ReachDatabase and call it for one reach.
    That runs the getattr loop over all 105 reach variables plus the area_fits
    and discharge_models subsetting, so if it returns, RiverObs can read this
    file. Also confirms node_length carries real values, since that is the
    divisor in width = area / node_length.
 
    Reads the whole database into memory; on a continent-scale SWORD file this
    can take a minute or two.
    """
    import warnings
    from RiverObs.ReachDatabase import ReachDatabase
 
    print('\nverifying -- loading the whole database, may take a minute...')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        db = ReachDatabase.from_ncfile(path)
        reach_ids = np.asarray(db.reaches.reach_id[:]).ravel()
        if reach_ids.size == 0:
            print('  FAILED: no reaches in the file')
            return False
        rid = int(reach_ids[0])
        try:
            this = db(rid)
        except Exception as exc:
            print('  FAILED: {}: {}'.format(type(exc).__name__, exc))
            print('  Send me this output.')
            return False
 
        raw_type = this['reaches']['type'][0]
        node_length = this['nodes']['node_length']
        n_nodes = len(np.asarray(this['nodes']['x']).ravel())
 
    resolved = rid % 10 if raw_type is np.ma.masked else int(raw_type)
    print('  OK -- {} reaches in the database'.format(reach_ids.size))
    print('  reach {}: {} nodes, type {} ({})'.format(
        rid, n_nodes, resolved,
        'derived from reach_id' if raw_type is np.ma.masked
        else 'read from file'))
 
    if np.all(np.ma.getmaskarray(node_length)):
        print('  *** node_length is fully masked -- widths would be wrong ***')
        return False
    nl = np.asarray(node_length).ravel().astype(float)
    print('  node_length min {:.1f} max {:.1f} m -- real values, good '
          '(this is the width divisor)'.format(np.nanmin(nl), np.nanmax(nl)))
    return True
 
 
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('sword_nc')
    p.add_argument('--verify', action='store_true',
                   help='load the file through RiverObs and call one reach')
    p.add_argument('--write-padded', default=None, help=argparse.SUPPRESS)
    args = p.parse_args()
 
    if args.write_padded:
        print('NOTE: --write-padded has been removed. Padding is unnecessary '
              '-- RiverObs returns masked arrays for variables it declares but '
              'your file lacks. Use the original SWORD netCDF and delete any '
              '*_riverobs.nc copies. Continuing as a check only.\n')
 
    critical, subgroup_problem = report(args.sword_nc)
    blocked = bool(any(critical.values()) or subgroup_problem)
 
    print('\n' + '=' * 70)
    if subgroup_problem:
        print('NOT USABLE -- /reaches is missing a required subgroup.')
        print('Check you have the full SWORD netCDF, not a subset.')
    elif any(critical.values()):
        print('NOT USABLE -- a variable RiverObs actually uses is missing. It')
        print('would come back masked and your widths would be silently wrong.')
    else:
        print('READY -- point reach_db_path at THIS file.')
        print('No padding, no copies, no methods caveat.')
    print('=' * 70)
 
    if args.verify:
        if not verify(args.sword_nc):
            sys.exit(1)
    elif not blocked:
        print('\nRe-run with --verify to load it through RiverObs and confirm.')
 
    if blocked:
        sys.exit(1)
 
 
if __name__ == '__main__':
    main()