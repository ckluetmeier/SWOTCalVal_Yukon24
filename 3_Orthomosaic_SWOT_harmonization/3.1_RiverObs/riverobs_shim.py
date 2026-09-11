"""Runtime fixes for the RiverObs CalVal path.

`import riverobs_shim` (or let run_calval2rivertile.py do it) before
running a cal/val job. It monkey-patches five things at import time
and touches no source files.

The defects it works around live in SWOTRiverEstimator.py,
Estimate.py and ReachDatabase.py, whose surrounding source text
varies between RiverObs revisions, so a text patch that matches one
checkout can fail on another. Patching behaviour at runtime is
revision-independent: it keys on class and function names, which are
stable.

Four of the five patches are guards and plumbing; the fifth
(corridor relaxation) is opt-in and off by default. None of them
touch pixel-to-node assignment, area aggregation, or the
area-to-width conversion, so a RiverTile produced with this shim
loaded has the same node areas and widths as one produced without
it. Without it, the run does not complete.

  1. SWOTRiverEstimator.__init__
     CalValToRiverTile never passes reach_pct_good_sus_thresh
     (L2PixcToRiverTile does), so it arrives as None and
     get_reach_mask() raises
       TypeError: '<=' not supported between instances of 'float'
       and 'NoneType'
     The shim substitutes 0. Every water-mask cell is good, so the
     percentage being tested is always 100 and any value below that
     is neutral.

  2. RiverReach.__init__
     Under height_agg_method='orig', wse_s_u is a scalar constant
     rather than a per-node array, and the output packing loop
     raises
       ValueError: zero-dimensional arrays cannot be concatenated
     A single-node reach hits the same wall from the other side,
     with a different message:
       ValueError: all the input arrays must have same number of
       dimensions, but the array at index 0 has 1 dimension(s) and
       the array at index N has 0 dimension(s)
     Estimate.py treats every RiverReach attribute except
     ds/metadata as a NODE variable and feeds it to np.concatenate,
     which rejects anything 0-d no matter how many nodes the reach
     has. The shim broadcasts any scalar attribute -- 0-d ndarray,
     np.float64 and friends, or a plain Python number -- to the
     reach's node count, including when that count is 1. Real
     per-node arrays and string attributes are untouched.

  3. ReachExtractor.__init__
     SWORD reaches have different numbers of centerline points. The
     output packing loop does
     np.array([reach.metadata[k] for reach in ...]) over
     centerline_lon / centerline_lat, and since numpy 1.24 a ragged
     list raises
       ValueError: setting an array element with a sequence. The
       requested array has an inhomogeneous shape after 1
       dimensions.
     The RiverTile product declares centerline_lat/lon as 2-D
     (reaches, centerlines), so a rectangular array is what it needs
     anyway. The shim pads every reach's centerline arrays to the
     longest one using the product's own fill value. rivertile.py
     already filters the padding out when it builds the reach
     LineString -- `is_valid = np.abs(lats) < 90`.

  4. ReachExtractor.__init__ -- corridor relaxation (opt-in,
     default off)
     RiverObs sizes its cross-channel search corridor from the prior
     database's channel width. Two gates apply, both scaled by the
     per-node wth_coef and ext_dist_coef:

       gate 1  SWOTRiverEstimator.assign_reaches sets
                   search_width = 2 * prior_width * wth_coef
               and RiverObs.get_ext_dist_threshold takes
                   half-corridor = search_width / 3
                                 = 0.667 * prior_width * wth_coef

       gate 2  SWOTRiverEstimator.assign_reaches_ext_dist_coef
               re-masks with
                   extreme_dist = ext_dist_coef
                                  * max(node_spacing,
                                        max(prior_max_width,
                                            prior_width)
                                        * wth_coef)

     The effective corridor is the tighter of the two. Because both
     are keyed to the PRIOR channel width, water that lies well
     outside the prior channel -- anabranches, secondary threads of
     a braidplain, wide side channels -- is excluded even when it is
     unambiguously part of the river. For a reference dataset where
     every digitized polygon is known-good river water, that
     exclusion is a loss, not a filter.

     Setting RIVEROBS_WTH_COEF_FACTOR /
     RIVEROBS_EXT_DIST_COEF_FACTOR (or calling
     set_corridor_factors()) multiplies those two coefficients,
     which widens both gates together. Nothing else is touched: the
     coefficients are corridor knobs only and are not written to any
     output product.

  5. RiverObs.get_node_stat
     time_from_prev_xover / time_to_next_xover are per-line PIXC
     variables a water mask does not have. SWOTRiverEstimator skips
     loading them, then asks for their node means anyway, raising
       AttributeError: 'RiverNode' object has no attribute
       'time_from_prev_xover'
     The shim returns the missing-value sentinel for exactly these
     two and re-raises for anything else, so RiverObs's own
     AttributeError handlers (sig0, geoid, tides, layover_impact,
     bright_land_flag) still run as upstream wrote them.
"""

PIPELINE_VERSION = '1.0.1'


import logging
import os

import numpy as np

LOGGER = logging.getLogger('riverobs_shim')

_ABSENT_REPORTED = set()

# Only these are substituted. Every other missing variable is re-raised
# so that RiverObs's own try/except handlers (sig0, geoid, tides,
# layover_impact, bright_land_flag) run exactly as upstream designed
# them to. Intercepting those would replace upstream's NaN fill with the
# missing-value sentinel and quietly change what lands in the output
# product.
SUBSTITUTABLE = ('time_from_prev_xover', 'time_to_next_xover')

# Corridor relaxation. 1.0 = stock RiverObs behaviour. Set with the
# environment variables below, or by calling set_corridor_factors()
# before processing.
#
#   RIVEROBS_WTH_COEF_FACTOR       multiplies wth_coef      (both gates)
#   RIVEROBS_EXT_DIST_COEF_FACTOR  multiplies ext_dist_coef (gate 2, and the
#                                  dominant-label extension in gate 1)
#
# Raising these admits more water per node. A corridor wider than needed
# lets adjacent reaches claim the same water (cross-reach double
# counting), which is the failure mode to watch for.
WTH_COEF_FACTOR = float(os.environ.get('RIVEROBS_WTH_COEF_FACTOR', 1.0))
EXT_DIST_COEF_FACTOR = float(os.environ.get('RIVEROBS_EXT_DIST_COEF_FACTOR', 1.0))


def set_corridor_factors(wth=None, ext_dist=None):
    """Override the corridor factors programmatically. Call before processing."""
    global WTH_COEF_FACTOR, EXT_DIST_COEF_FACTOR
    if wth is not None:
        WTH_COEF_FACTOR = float(wth)
    if ext_dist is not None:
        EXT_DIST_COEF_FACTOR = float(ext_dist)


# The RiverTile product's declared fill for centerline_lat/lon. Chosen
# so that rivertile.py's own `np.abs(lats) < 90` test drops the padding.
MISSING_VALUE_FLT = -999999999999.0


def install():
    """Apply all patches. Idempotent."""
    _patch_estimator_init()
    _patch_river_reach()
    _patch_reach_extractor()
    _patch_get_node_stat()
    LOGGER.info('riverobs_shim %s installed (5 runtime patches)',
                PIPELINE_VERSION)


# Patches 3 and 4: ReachExtractor.__init__
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

        _apply_corridor_factors(reaches)

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


# Patch 1: SWOTRiverEstimator.__init__
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


# Patch 2: RiverReach.__init__
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
        if n_nodes < 1:
            return
        for key, value in list(self.__dict__.items()):
            if key in ('ds', 'metadata'):
                continue
            # Every attribute except ds/metadata is treated by
            # Estimate.py as a NODE variable and fed to np.concatenate,
            # which rejects anything 0-dimensional -- including a
            # single-node reach. Three shapes have to be caught:
            #   * a 0-d ndarray
            #   * np.float64 and friends, which are np.generic, NOT np.ndarray
            #   * a plain Python float/int
            # Strings and bytes are excluded: broadcasting one would be
            # wrong, and a string attribute would fail further
            # downstream anyway.
            is_zero_d = (
                (isinstance(value, np.ndarray) and value.ndim == 0) or
                (isinstance(value, np.generic)
                 and not isinstance(value, (np.str_, np.bytes_))) or
                isinstance(value, (float, int, bool)))
            if is_zero_d:
                setattr(self, key, np.repeat(np.asarray(value)[np.newaxis],
                                             n_nodes))
                LOGGER.debug('shim: broadcast scalar RiverReach.%s to %d nodes',
                             key, n_nodes)

    __init__._shimmed = True
    RiverReach.__init__ = __init__


# Patch 5: RiverObs.get_node_stat
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


def _apply_corridor_factors(reaches):
    """
    Multiply the per-node corridor coefficients, and report the resulting
    cross-channel half-corridor in metres so the effect is visible rather than
    implicit.
    """
    if WTH_COEF_FACTOR == 1.0 and EXT_DIST_COEF_FACTOR == 1.0:
        return
    for r in reaches:
        for attr, factor in (('wth_coef', WTH_COEF_FACTOR),
                             ('ext_dist_coef', EXT_DIST_COEF_FACTOR)):
            v = getattr(r, attr, None)
            if v is None:
                continue
            setattr(r, attr, np.asarray(v, dtype='f8') * factor)

    # Report what this means on the ground, for the first few reaches.
    for r in reaches[:5]:
        try:
            width = np.ma.filled(np.asarray(r.width, dtype='f8'), np.nan)
            maxw = np.ma.filled(np.asarray(r.max_width, dtype='f8'), np.nan)
            wth = np.asarray(r.wth_coef, dtype='f8')
            ext = np.asarray(r.ext_dist_coef, dtype='f8')
            node_len = np.asarray(r.node_length, dtype='f8')
            gate1 = np.nanmedian(2.0 * width * wth / 3.0)
            gate2 = np.nanmedian(ext * np.maximum(
                node_len, np.fmax(maxw, width) * wth))
            LOGGER.info('shim: corridor reach %s -- prior width %.0f m, '
                        'half-corridor gate1 %.0f m, gate2 %.0f m '
                        '(binding: %.0f m)',
                        getattr(r, 'reach_index', '?'), np.nanmedian(width),
                        gate1, gate2, min(gate1, gate2))
        except Exception:
            pass
    LOGGER.warning('shim: corridor relaxed -- wth_coef x%.2f, ext_dist_coef '
                   'x%.2f. Node assignment may pick up water from adjacent '
                   'reaches (cross-reach double counting).',
                   WTH_COEF_FACTOR, EXT_DIST_COEF_FACTOR)


def absent_variables():
    """Per-pixel variables the shim had to substitute for, this session."""
    return sorted(_ABSENT_REPORTED)


install()