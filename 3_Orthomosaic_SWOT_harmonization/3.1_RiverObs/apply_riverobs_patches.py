#!/usr/bin/env python
"""Apply the RiverObs source patches.

Two files, four edits. Everything else the pipeline needs is applied
at runtime by riverobs_shim.py, which keys on class and function
names rather than on source text and is therefore robust across
RiverObs revisions.

Usage:
    python apply_riverobs_patches.py --riverobs-root /path/to/RiverObs

    --revert     restore the .orig_backup copies
    --status     report what is and isn't applied, and change nothing

The two files touched are:

    src/SWOTRiver/products/calval.py
        implement from_airborne_imagery, and declare the three
        pixel-cloud quality flags
    src/bin/calval2rivertile.py
        register the airborne_watermask format and its dispatch
        branch

Idempotent. The from_airborne_imagery edit replaces the whole method
whatever its current contents, and every other edit is skipped if it
is already present.
"""

PIPELINE_VERSION = '1.0.0'


import argparse
import os
import re
import shutil
import sys


NEW_METHOD = '''    @classmethod
    def from_airborne_imagery(cls, water_mask_file, water_value=1,
                              class_value=4):
        """
        Converter for a rasterized water mask derived from airborne imagery.
 
        Expects a single-band GeoTIFF on a regular grid in any projected CRS
        with metre units (e.g. EPSG:32606, UTM 6N), where cells equal to
        `water_value` are water and everything else is not.
 
        Every water cell becomes one pixel-cloud record carrying the grid's own
        cell area, so node area is sum(pixel_area) over the cells RiverObs
        assigns to that node, and node width is that area divided by the SWORD
        prior node_length -- the same definition the operational processor uses
        for RiverSP (RiverNode.width_area).
 
        `range_index` and `azimuth_index` are set to the raster column and row,
        which makes SWOTRiverEstimator's connected-component segmentation
        operate on the water mask's own 2D grid.
 
        Heights are zero and the quality flags are all-good: this converter is
        for area and width validation only. Do not read wse or slope from a
        RiverTile produced this way.
 
        Adapted from from_cnes_watermask() with two corrections:
          1. the CRS is read from the file rather than hard-coded to EPSG:32610,
          2. rasterio.transform.xy() returns flat sequences for 2D row/col
             input, so the coordinate arrays are reshaped before masking.
        """
        import rasterio
        from pyproj import Transformer
 
        with rasterio.open(water_mask_file) as src:
            if src.crs is None:
                raise ValueError(
                    "water mask {} has no CRS; RiverObs needs georeferenced "
                    "input".format(water_mask_file))
 
            arr = src.read(1)
            nrow, ncol = arr.shape
 
            # EPSG:4326 is lat-lon axis order, so transform() returns (lat, lon)
            transformer = Transformer.from_crs(src.crs, "EPSG:4326")
 
            cols, rows = np.meshgrid(np.arange(ncol), np.arange(nrow))
            xs, ys = rasterio.transform.xy(src.transform, rows, cols)
            xs = np.asarray(xs).reshape(nrow, ncol)
            ys = np.asarray(ys).reshape(nrow, ncol)
            lats, lons = transformer.transform(xs, ys)
 
            pixel_area = abs(src.transform[0] * src.transform[4])
            crs_str = str(src.crs)
 
        mask = (arr == water_value)
        npix = int(mask.sum())
        if npix == 0:
            raise ValueError(
                "no cells equal to water_value={} in {}".format(
                    water_value, water_mask_file))
 
        range_index, azimuth_index = np.meshgrid(
            np.arange(ncol, dtype=int), np.arange(nrow, dtype=int))
 
        klass = cls()
        klass.classification = np.full(npix, class_value, dtype='u1')
        klass.water_frac = np.ones(npix)
        klass.latitude = lats[mask]
        klass.longitude = lons[mask]
        klass.height = np.zeros(npix)
        klass.range_index = range_index[mask]
        klass.azimuth_index = azimuth_index[mask]
        klass.pixel_area = np.full(npix, pixel_area)
        klass.geolocation_qual = np.zeros(npix, dtype='u4')
        klass.classification_qual = np.zeros(npix, dtype='u4')
        klass.sig0_qual = np.zeros(npix, dtype='u4')
        klass.interferogram_size_azimuth = nrow
        klass.interferogram_size_range = ncol
 
        LOGGER.info(
            "from_airborne_imagery: %d water cells, %.4f m2 per cell, "
            "%.1f m2 total, CRS %s", npix, pixel_area,
            npix * pixel_area, crs_str)
 
        return klass
'''


# Quality flags on SimplePixelCloud
#
# Upstream has two of them commented out and misspells one as
# "classificaiton_qual"; SWOTRiverEstimator looks for
# "classification_qual", and sig0_qual is not listed at all. Without
# all three, process_node() raises
#   TypeError: 'NoneType' object is not subscriptable
QUAL_OLD = """        #['geolocation_qual', {'dtype': 'u4'}],
        #['classificaiton_qual', {'dtype': 'u4'}],
"""
QUAL_NEW = """        ['geolocation_qual', {'dtype': 'u4'}],
        ['classification_qual', {'dtype': 'u4'}],
        ['sig0_qual', {'dtype': 'u4'}],
"""
QUAL_MARKER = "['classification_qual', {'dtype': 'u4'}],"

DISPATCH_NEW = """
    elif args.format == 'airborne_watermask':
        pixc_simple = SimplePixelCloud.from_airborne_imagery(args.input_file)
"""
DISPATCH_ANCHOR = """    elif args.format == 'water_mask':
        pixc_simple = SimplePixelCloud.from_cnes_watermask(args.input_file)
"""


def backup(path):
    bak = path + '.orig_backup'
    if not os.path.exists(bak):
        shutil.copyfile(path, bak)


def replace_method(text, method_name, new_source):
    """
    Replace a whole class method, however it is currently written.

    Finds the `def <name>(` line, walks back over any decorator lines above it,
    then forward to the first line at the same indentation that starts a new
    decorator, def, or class. Robust to the method having been edited already.
    Returns (new_text, found).
    """
    lines = text.splitlines(keepends=True)
    def_re = re.compile(r'^(\s*)def\s+' + re.escape(method_name) + r'\s*\(')

    def_index = None
    indent = None
    for i, line in enumerate(lines):
        m = def_re.match(line)
        if m:
            indent = m.group(1)
            def_index = i
            break
    if def_index is None:
        return text, False

    # Walk back over decorators / comments attached to the def
    start = def_index
    while start > 0:
        prev = lines[start - 1].strip()
        if prev.startswith('@') or prev.startswith('#'):
            start -= 1
        else:
            break

    # Scan for the end starting AFTER the def line, not after `start` --
    # the def line itself matches the boundary pattern and would end the
    # block instantly.
    boundary = re.compile(r'^' + indent + r'(@|def\s|class\s)')
    end = len(lines)
    for j in range(def_index + 1, len(lines)):
        stripped = lines[j].strip()
        if not stripped:
            continue
        if boundary.match(lines[j]):
            end = j
            break
        # A line dedented past the method body ends it too
        leading = len(lines[j]) - len(lines[j].lstrip())
        if leading < len(indent) and stripped:
            end = j
            break

    return ''.join(lines[:start]) + new_source + ''.join(lines[end:]), True


def do_calval(path, status_only):
    results = []
    text = open(path).read()
    original = text

    # 1. from_airborne_imagery -- whole-method replacement
    if NEW_METHOD.strip() in text:
        results.append(('implement from_airborne_imagery', 'skip',
                        'already applied'))
    else:
        text, found = replace_method(text, 'from_airborne_imagery', NEW_METHOD)
        if found:
            results.append(('implement from_airborne_imagery', 'ok',
                            'replaced whole method'))
        else:
            results.append((
                'implement from_airborne_imagery', 'fail',
                'no "def from_airborne_imagery(" found in ' + path))

    # 2. quality flags
    if QUAL_MARKER in text:
        results.append(('declare pixel-cloud quality flags', 'skip',
                        'already applied'))
    elif QUAL_OLD in text:
        text = text.replace(QUAL_OLD, QUAL_NEW)
        results.append(('declare pixel-cloud quality flags', 'ok', ''))
    else:
        results.append((
            'declare pixel-cloud quality flags', 'fail',
            'could not find the commented-out quality flag lines in '
            'SimplePixelCloud.VARIABLES'))

    if text != original and not status_only:
        backup(path)
        open(path, 'w').write(text)
    return results


def do_calval2rivertile(path, status_only):
    results = []
    text = open(path).read()
    original = text

    # 3a. FORMATS
    if "'airborne_watermask'" in text.split('def main')[0]:
        results.append(('register airborne_watermask format', 'skip',
                        'already applied'))
    elif "'prior_water_mask']" in text:
        text = text.replace("'prior_water_mask']",
                            "'prior_water_mask', 'airborne_watermask']", 1)
        results.append(('register airborne_watermask format', 'ok', ''))
    else:
        results.append(('register airborne_watermask format', 'fail',
                        'could not find the FORMATS list'))

    # 3b. dispatch branch
    if 'from_airborne_imagery(args.input_file)' in text:
        results.append(('add airborne_watermask dispatch branch', 'skip',
                        'already applied'))
    elif DISPATCH_ANCHOR in text:
        text = text.replace(DISPATCH_ANCHOR, DISPATCH_ANCHOR + DISPATCH_NEW, 1)
        results.append(('add airborne_watermask dispatch branch', 'ok', ''))
    else:
        results.append((
            'add airborne_watermask dispatch branch', 'fail',
            "could not find the 'water_mask' elif branch to insert after"))

    if text != original and not status_only:
        backup(path)
        open(path, 'w').write(text)
    return results


FILES = ['src/SWOTRiver/products/calval.py', 'src/bin/calval2rivertile.py']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--riverobs-root', required=True)
    p.add_argument('--revert', action='store_true')
    p.add_argument('--status', action='store_true')
    args = p.parse_args()
    print('# {} {}'.format(os.path.basename(__file__), PIPELINE_VERSION))

    root = os.path.expanduser(args.riverobs_root)
    if not os.path.isdir(os.path.join(root, 'src')):
        sys.exit('no src/ under {} -- is that the right folder?'.format(root))

    if args.revert:
        for rel in FILES:
            bak = os.path.join(root, rel) + '.orig_backup'
            if os.path.exists(bak):
                shutil.copyfile(bak, os.path.join(root, rel))
                print('  reverted', rel)
        print('done. Backups left in place.')
        return

    calval = os.path.join(root, FILES[0])
    c2rt = os.path.join(root, FILES[1])
    for f in (calval, c2rt):
        if not os.path.exists(f):
            sys.exit('missing file: ' + f)

    results = do_calval(calval, args.status)
    results += do_calval2rivertile(c2rt, args.status)

    tag = {'ok': '[ok]  ', 'skip': '[skip]', 'fail': '[FAIL]'}
    for label, state, note in results:
        print('  {} {:<45s} {}'.format(tag[state], label, note))

    n_fail = sum(1 for _, s, _ in results if s == 'fail')
    print('\n{} applied, {} already in place, {} failed'.format(
        sum(1 for _, s, _ in results if s == 'ok'),
        sum(1 for _, s, _ in results if s == 'skip'), n_fail))

    if n_fail:
        print('\nIf an edit failed, your RiverObs revision differs from the')
        print('tested one. Record it with:')
        print('  cd "{}" && git log -1 --format=%H'.format(root))
        sys.exit(1)

    print('\nVerify with:')
    print('  python "{}/src/bin/calval2rivertile.py" --help'.format(root))
    print('  ...the format list should end with ",airborne_watermask}"')
    print('\nThe remaining fixes are applied at runtime by riverobs_shim.py.')
    print('Run jobs through run_calval2rivertile.py, which loads it.')


if __name__ == '__main__':
    main()