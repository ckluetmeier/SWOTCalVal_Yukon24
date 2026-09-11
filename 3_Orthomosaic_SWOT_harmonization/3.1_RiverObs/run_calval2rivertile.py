#!/usr/bin/env python
"""Run calval2rivertile.py with the runtime shim loaded.

Use this instead of calling src/bin/calval2rivertile.py directly.
Same arguments, same outputs; it imports riverobs_shim first.

Usage:
    python run_calval2rivertile.py \\
        upperYR_071024_water_3m.tif airborne_watermask \\
        out_rivertile.nc out_pixcvec.nc riverobs_ortho_v17b.rdf \\
        out_pixc.nc \\
        --riverobs-root "/path/to/RiverObs" --log-level info

`--riverobs-root` is optional if PYTHONPATH already points at
RiverObs/src. Keep this file in the same folder as riverobs_shim.py.

The search-corridor factors are set here, because this is where the
pixel-to-node assignment happens:

    --wth-coef-factor 3.0 --ext-dist-coef-factor 3.0

They are equivalent to exporting RIVEROBS_WTH_COEF_FACTOR /
RIVEROBS_EXT_DIST_COEF_FACTOR, and take precedence over those
variables if both are set. `1.0` is stock RiverObs. See
riverobs_shim.py for what the two gates are.

Whatever is used is written to the RiverTile and PIXCVec as the
global attributes `riverobs_wth_coef_factor` and
`riverobs_ext_dist_coef_factor`, so a product carries a record of the
corridor that produced it and 3.1.3_true_node_polygons.py can report
it back.
"""

PIPELINE_VERSION = '1.1.0'


import os
import sys

# --riverobs-root is ours, not calval2rivertile's -- pull it out of argv
# first.
_root = None
if '--riverobs-root' in sys.argv:
    i = sys.argv.index('--riverobs-root')
    _root = os.path.expanduser(sys.argv[i + 1])
    del sys.argv[i:i + 2]

# Likewise the corridor factors. These must be settled BEFORE
# riverobs_shim is imported, because the shim reads the environment at
# import time.
_factors = {}
for _flag, _env in (('--wth-coef-factor', 'RIVEROBS_WTH_COEF_FACTOR'),
                    ('--ext-dist-coef-factor', 'RIVEROBS_EXT_DIST_COEF_FACTOR')):
    if _flag in sys.argv:
        i = sys.argv.index(_flag)
        try:
            _val = float(sys.argv[i + 1])
        except (IndexError, ValueError):
            sys.exit('{} needs a number, e.g. {} 3.0'.format(_flag, _flag))
        if _val <= 0:
            sys.exit('{} must be positive (1.0 = stock RiverObs)'.format(_flag))
        del sys.argv[i:i + 2]
        os.environ[_env] = repr(_val)
        _factors[_env] = _val

if _root:
    src = os.path.join(_root, 'src')
    if not os.path.isdir(src):
        sys.exit('no src/ under {}'.format(_root))
    sys.path.insert(0, src)

# The shim must be importable; it lives next to this file.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import SWOTRiver  # noqa: F401
except ImportError:
    sys.exit(
        'Cannot import SWOTRiver. Either set PYTHONPATH to RiverObs/src or '
        'pass --riverobs-root /path/to/RiverObs')

import riverobs_shim  # noqa: F401  (installs the patches on import)

# calval2rivertile's main() is located and executed in this process,
# which riverobs_shim has already patched.
import importlib.util

_cv_path = None
for entry in sys.path:
    candidate = os.path.join(entry, '..', 'src', 'bin', 'calval2rivertile.py')
    candidate = os.path.normpath(candidate)
    if os.path.exists(candidate):
        _cv_path = candidate
        break
if _cv_path is None and _root:
    _cv_path = os.path.join(_root, 'src', 'bin', 'calval2rivertile.py')
if _cv_path is None or not os.path.exists(_cv_path):
    import SWOTRiver as _sr
    guess = os.path.normpath(os.path.join(
        os.path.dirname(_sr.__file__), '..', 'bin', 'calval2rivertile.py'))
    _cv_path = guess

if not os.path.exists(_cv_path):
    sys.exit('could not find calval2rivertile.py; pass --riverobs-root')

_spec = importlib.util.spec_from_file_location('calval2rivertile', _cv_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

def _stamp_corridor(paths):
    """Record the corridor factors on the outputs as global attributes.

    Purely additive, and never allowed to fail a run that has already produced
    its products -- a netCDF that cannot be reopened for append is worth a
    warning, not a lost result.
    """
    import netCDF4
    attrs = {
        'riverobs_wth_coef_factor': riverobs_shim.WTH_COEF_FACTOR,
        'riverobs_ext_dist_coef_factor': riverobs_shim.EXT_DIST_COEF_FACTOR,
        'riverobs_shim_version': riverobs_shim.PIPELINE_VERSION,
    }
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        try:
            with netCDF4.Dataset(path, 'a') as ds:
                for k, v in attrs.items():
                    setattr(ds, k, v)
        except Exception as exc:
            print('note: could not stamp corridor factors on {}: {}: {}'
                  .format(os.path.basename(path), type(exc).__name__, exc))


if __name__ == '__main__':
    # Positional order is fixed by calval2rivertile.py's own parser:
    #   input_file format out_riverobs_file out_pixcvec_file rdf_file [out_pixc_file]
    # Strip the one optional flag it takes so the positions hold no
    # matter where the caller put it.
    _pos, _skip = [], False
    for _a in sys.argv[1:]:
        if _skip:
            _skip = False
            continue
        if _a in ('-l', '--log-level'):
            _skip = True
            continue
        if _a.startswith('-') and len(_a) > 1:
            continue
        _pos.append(_a)
    _outputs = _pos[2:4]

    if _factors:
        print('corridor factors: wth_coef x{}, ext_dist_coef x{}'.format(
            riverobs_shim.WTH_COEF_FACTOR, riverobs_shim.EXT_DIST_COEF_FACTOR))

    _mod.main()
    _stamp_corridor(_outputs)
    absent = riverobs_shim.absent_variables()
    if absent:
        print('\nshim substituted missing_value for these pixel variables: {}'
              .format(', '.join(absent)))
        print('(time_from_prev_xover / time_to_next_xover are expected; '
              'anything else is worth a look)')
