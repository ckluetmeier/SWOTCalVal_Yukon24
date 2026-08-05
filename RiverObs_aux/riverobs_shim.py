"""
=============================================================================
riverobs_shim.py -- runtime fixes for the RiverObs CalVal path
-----------------------------------------------------------------------------
`import riverobs_shim` (or let run_calval2rivertile.py do it) BEFORE running a
cal/val job. It monkey-patches three things at import time and touches no
source files.
 
WHY A SHIM AND NOT SOURCE EDITS
These three bugs live in SWOTRiverEstimator.py and Estimate.py, and the exact
surrounding text varies between RiverObs revisions -- a text patch that matches
one clone fails on another. Patching the behaviour at runtime is revision-
independent: it keys on function and class names, which are stable, instead of
on whitespace.
 
WHAT IT CHANGES, AND WHAT IT DELIBERATELY DOES NOT
All three are guards and plumbing. None of them touch pixel-to-node assignment,
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
 
  3. RiverObs.get_node_stat
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
 
 
def install():
    """Apply all three patches. Idempotent."""
    _patch_estimator_init()
    _patch_river_reach()
    _patch_get_node_stat()
    LOGGER.info('riverobs_shim installed (3 runtime patches)')
 
 
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
 
 
# -- 3 -------------------------------------------------------------------------
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
 