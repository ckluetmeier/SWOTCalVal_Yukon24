#!/usr/bin/env python
"""
=============================================================================
run_calval2rivertile.py -- calval2rivertile.py with the runtime shim loaded
-----------------------------------------------------------------------------
Use this INSTEAD OF calling src/bin/calval2rivertile.py directly. Same
arguments, same outputs; it just imports riverobs_shim first (see that file for
what the three runtime patches do and why they don't affect widths).
 
    python run_calval2rivertile.py \\
        upperYR_071024_water_3m.tif airborne_watermask \\
        out_rivertile.nc out_pixcvec.nc riverobs_ortho_v17b.rdf out_pixc.nc \\
        --riverobs-root "/path/to/RiverObs" --log-level info
 
`--riverobs-root` is optional if PYTHONPATH already points at RiverObs/src.
Keep this file in the same folder as riverobs_shim.py.
=============================================================================
"""
 
import os
import sys
 
# --riverobs-root is ours, not calval2rivertile's -- pull it out of argv first.
_root = None
if '--riverobs-root' in sys.argv:
    i = sys.argv.index('--riverobs-root')
    _root = os.path.expanduser(sys.argv[i + 1])
    del sys.argv[i:i + 2]
 
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
 
# Locate and execute calval2rivertile's main() in this already-patched process.
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
 
if __name__ == '__main__':
    _mod.main()
    absent = riverobs_shim.absent_variables()
    if absent:
        print('\nshim substituted missing_value for these pixel variables: {}'
              .format(', '.join(absent)))
        print('(time_from_prev_xover / time_to_next_xover are expected; '
              'anything else is worth a look)')
 