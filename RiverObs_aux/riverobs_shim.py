"""
=============================================================================
riverobs_shim.py -- runtime fixes for the RiverObs CalVal path
-----------------------------------------------------------------------------
`import riverobs_shim` (or let run_calval2rivertile.py do it) BEFORE running a
cal/val job. It monkey-patches three things at import time and touches no
source files.
 
WHY A SHIM AND NOT SOURCE EDITS
These defects live in SWOTRiverEstimator.py, Estimate.py and ReachDatabase.py,
whose surrounding source text varies between RiverObs revisions -- a text patch
that matches one checkout can fail on another. Patching behaviour at runtime is
revision-independent: it keys on class and function names, which are stable.
 
WHAT IT CHANGES, AND WHAT IT DELIBERATELY DOES NOT
All four are guards and plumbing. None of them touch pixel-to-node assignment,
area aggregation, or the area-to-width conversion, so a RiverTile produced with
this shim loaded has the same node areas and widths as one produced without it
-- the difference is that without it, the run crashes.
 
  1. SWOTRiverEstimator.__init__
     CalValToRiverTile never passes reach_pct_good_sus_thresh (L2PixcToRiverTile
     does), so it arrives as None and get_reach_mask() raises
       TypeError: '<=' not supported between instances of 'float' and 'NoneType'
     The shim substitutes 0. Every water-mask cell is good, so the percentage
     being tested is always 100 and any value below that is neutral.
 
  2. RiverReach.__init__
     Under height_agg_method='orig', wse_s_u is a scalar constant rather than a
     per-node array, and the output packing loop raises
       ValueError: zero-dimensional arrays cannot be concatenated
     The shim broadcasts any 0-d numpy attribute to the reach's node count.
     Only genuinely 0-d arrays are touched; real per-node arrays are untouched.
 
  3. ReachExtractor.__init__
     SWORD reaches have different numbers of centerline points. The output
     packing loop does np.array([reach.metadata[k] for reach in ...]) over
     centerline_lon / centerline_lat, and since numpy 1.24 a ragged list raises
       ValueError: setting an array element with a sequence. The requested
       array has an inhomogeneous shape after 1 dimensions.
     The RiverTile product declares centerline_lat/lon as 2-D
     (reaches, centerlines), so a rectangular array is what it needs anyway.
     The shim pads every reach's centerline arrays to the longest one using the
     product's own fill value. rivertile.py already filters the padding out
     when it builds the reach LineString -- `is_valid = np.abs(lats) < 90`.
 
  4. RiverObs.get_node_stat
     time_from_prev_xover / time_to_next_xover are per-line PIXC variables a
     water mask does not have. SWOTRiverEstimator skips loading them, then asks
     for their node means anyway, raising
       AttributeError: 'RiverNode' object has no attribute 'time_from_prev_xover'
     The shim returns the missing-value sentinel for exactly these two and
     re-raises for anything else, so RiverObs's own AttributeError handlers
     (sig0, geoid, tides, layover_impact, bright_land_flag) still run as
     upstream wrote them.
=============================================================================
"""
 
PIPELINE_VERSION = '1.0.0'
 
 
import logging
 
import numpy as np
 
LOGGER = logging.getLogger('riverobs_shim')
 
_ABSENT_REPORTED = set()
 
# Only these are substituted. Every other missing variable is re-raised so that
# RiverObs's own try/except handlers (sig0, geoid, tides, layover_impact,
# bright_land_flag) run exactly as upstream designed them to. Intercepting those
# would replace upstream's NaN fill with the missing-value sentinel and quietly
# change what lands in the output product.
SUBSTITUTABLE = ('time_from_prev_xover', 'time_to_next_xover')
 
 
# The RiverTile product's declared fill for centerline_lat/lon. Chosen so that
# rivertile.py's own `np.abs(lats) < 90` test drops the padding.
MISSING_VALUE_FLT = -999999999999.0
 
 
def install():
    """Apply all four patches. Idempotent."""
    _patch_estimator_init()
    _patch_river_reach()
    _patch_reach_extractor()
    _patch_get_node_stat()
    LOGGER.info('riverobs_shim %s installed (4 runtime patches)',
                PIPELINE_VERSION)
 
 
# -- 3 -------------------------------------------------------------------------
def _patch_reach_extractor():
    from RiverObs import ReachDatabase as _RD
    klass = _RD.ReachExtractor
    if getattr(klass.__init__, '_shimmed', False):
        return
    original = klass.__init__
 
    KEYS = ('centerline_lon', 'centerline_lat')
 
    def __init__(self, *args, **kwargs):
        original(self, *args, **kwargs)
        reaches = getattr(self, 'reach', None) or []
        lengths = [len(np.atleast_1d(r.metadata[k]))
                   for r in reaches for k in KEYS if k in r.metadata]
        if not lengths or len(set(lengths)) == 1:
            return                      # already rectangular, nothing to do
        n_max = max(lengths)
        for r in reaches:
            for k in KEYS:
                if k not in r.metadata:
                    continue
                v = np.atleast_1d(r.metadata[k])
                v = np.ma.filled(np.ma.masked_invalid(
                    np.asarray(v, dtype='f8')), MISSING_VALUE_FLT)
                if len(v) < n_max:
                    v = np.concatenate(
                        [v, np.full(n_max - len(v), MISSING_VALUE_FLT)])
                r.metadata[k] = v
        LOGGER.info('shim: padded centerline_lon/lat across %d reaches to %d '
                    'points (was %d-%d) so the reach output packs to a '
                    'rectangular array', len(reaches), n_max,
                    min(lengths), max(lengths))
 
    __init__._shimmed = True
    klass.__init__ = __init__
 
 
# -- 1 -------------------------------------------------------------------------
def _patch_estimator_init():
    import SWOTRiver
    klass = SWOTRiver.SWOTRiverEstimator
    if getattr(klass.__init__, '_shimmed', False):
        return
    original = klass.__init__
 
    def __init__(self, *args, **kwargs):
        if kwargs.get('reach_pct_good_sus_thresh') is None:
            kwargs['reach_pct_good_sus_thresh'] = 0
            LOGGER.debug('shim: reach_pct_good_sus_thresh None -> 0')
        return original(self, *args, **kwargs)
 
    __init__._shimmed = True
    klass.__init__ = __init__
 
 
# -- 2 -------------------------------------------------------------------------
def _patch_river_reach():
    from RiverObs.RiverReach import RiverReach
    if getattr(RiverReach.__init__, '_shimmed', False):
        return
    original = RiverReach.__init__
 
    def __init__(self, **kwds):
        original(self, **kwds)
        lat = getattr(self, 'lat', None)
        if lat is None:
            return
        n_nodes = len(np.atleast_1d(lat))
        if n_nodes <= 1:
            return
        for key, value in list(self.__dict__.items()):
            if key in ('ds', 'metadata'):
                continue
            # wse_s_u is np.float64, which is np.generic and NOT np.ndarray --
            # an ndarray-only test misses it.
            is_zero_d = (
                (isinstance(value, np.ndarray) and value.ndim == 0) or
                isinstance(value, np.generic))
            if is_zero_d:
                setattr(self, key, np.repeat(np.asarray(value)[np.newaxis],
                                             n_nodes))
                LOGGER.debug('shim: broadcast scalar RiverReach.%s to %d nodes',
                             key, n_nodes)
 
    __init__._shimmed = True
    RiverReach.__init__ = __init__
 
 
# -- 4 -------------------------------------------------------------------------
def _patch_get_node_stat():
    from RiverObs.RiverObs import RiverObs
    if getattr(RiverObs.get_node_stat, '_shimmed', False):
        return
    original = RiverObs.get_node_stat
 
    def get_node_stat(self, stat, var, all_nodes=False, **kwargs):
        try:
            return original(self, stat, var, all_nodes=all_nodes, **kwargs)
        except AttributeError:
            if var not in SUBSTITUTABLE:
                raise
            if var not in _ABSENT_REPORTED:
                _ABSENT_REPORTED.add(var)
                LOGGER.warning(
                    'shim: pixel variable %r absent from this pixel cloud; '
                    'returning missing_value for its node statistics', var)
            n = (len(self.all_nodes) if all_nodes
                 else len(self.populated_nodes))
            return [self.missing_value] * n
 
    get_node_stat._shimmed = True
    RiverObs.get_node_stat = get_node_stat
 
 
def absent_variables():
    """Per-pixel variables the shim had to substitute for, this session."""
    return sorted(_ABSENT_REPORTED)
 
 
install()