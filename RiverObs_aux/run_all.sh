
#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- steps 1, 2 and 3b for all six orthomosaics
# -----------------------------------------------------------------------------
#   conda activate RiverObs
#   cd .../YR2024_scripts/SWOTCalVal_Yukon24/RiverObs_aux
#   bash run_all.sh                  # everything
#   bash run_all.sh 1                # step 1 only
#   bash run_all.sh 2 3b             # steps 2 and 3b only
#
# Written for bash 3.2, which is what /bin/bash is on macOS -- no associative
# arrays, no mapfile.
#
# After this finishes, go to step 5.4 in the protocol: open the node polygons in
# QGIS over each orthomosaic and mark what to drop. Nothing here excludes any
# node.
# =============================================================================
 
set -uo pipefail
 
# ---------------------------------------------------------------- paths -----
BASE=/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats
SCRIPTS=/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/RiverObs_aux
RIVEROBS=/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/RiverObs
 
MASK_DIR=$BASE/Orthomosaics/width_validation
ORTHO_DIR=$BASE/Orthomosaics/orthos_resampled_3m/gap_filled
RASTER_DIR=$BASE/Orthomosaics/width_validation/rasterized_water_masks
OUT_DIR=$BASE/Orthomosaics/width_validation/RiverObs_output
QC_DIR=$BASE/Orthomosaics/width_validation/RiverObs_output/node_qc
 
RES=3
VERSIONS="v17b v16"          # e.g. "v17b v16" once you are ready for both
LOG_LEVEL=info
 
# --------------------------------------------------- survey definitions -----
# name | water mask basename | orthomosaic basename
# Names must match the SURVEYS list at the top of step2_run_riverobs.py.
SURVEYS="
CD_071024|chandalar_240710_ortho_25cm_watermask|chandalar_240710_ortho_3m_epsg32606
upperPR_CL_071024|coleen_240710_ortho_25cm_watermask|coleen_240710_ortho_3m_epsg32606
upperPR_CL_071624|coleen_240716_ortho_25cm_watermask|coleen_240716_ortho_3m_epsg32606
lowerPR_SJ_072624|sheenjek_240725_ortho_25cm_watermask|sheenjek_240725_ortho_3m_epsg32606
lowerYR_071624|yukonDS_240716_ortho_25cm_watermask|yukonDS_240716_ortho_3m_epsg32606
upperYR_071024|yukonUS_240710_ortho_25cm_watermask|yukonUS_240710_ortho_3m_epsg32606
"
 
STEPS="${*:-1 2 3b}"
want () { case " $STEPS " in *" $1 "*) return 0;; *) return 1;; esac; }
 
mkdir -p "$RASTER_DIR" "$OUT_DIR" "$QC_DIR"
cd "$SCRIPTS" || exit 1
 
FAILED=""
 
# =============================================================================
# STEP 1 -- rasterize each water mask to a 3 m grid
# =============================================================================
if want 1; then
echo
echo "################ STEP 1: rasterize water masks ################"
echo "$SURVEYS" | while IFS='|' read -r NAME MASK ORTHO; do
  [ -z "$NAME" ] && continue
  echo
  echo "--- $NAME ---"
  python step1_rasterize_watermask.py \
      "$MASK_DIR/$MASK.shp" \
      "$RASTER_DIR" \
      --name "$NAME" \
      --res "$RES" \
      --ortho-tif "$ORTHO_DIR/$ORTHO.tif" \
    || echo "!!! step1 FAILED for $NAME"
done
fi
 
# =============================================================================
# STEP 2 -- run RiverObs. This script loops over SURVEYS and versions itself.
# --keep-going means one bad survey does not cost you the CSVs for the rest.
# =============================================================================
if want 2; then
echo
echo "################ STEP 2: RiverObs ################"
python step2_run_riverobs.py \
    --riverobs-root "$RIVEROBS" \
    --raster-dir "$RASTER_DIR" \
    --rdf-dir . \
    --out-dir "$OUT_DIR" \
    --res "$RES" \
    --versions $VERSIONS \
    --keep-going \
    --log-level "$LOG_LEVEL"
STEP2_RC=$?
[ $STEP2_RC -ne 0 ] && FAILED="$FAILED step2(rc=$STEP2_RC)"
fi
 
# =============================================================================
# STEP 3b -- node polygons for the manual QGIS review, per survey per version
# =============================================================================
if want 3b; then
echo
echo "################ STEP 3b: node polygons ################"
for V in $VERSIONS; do
  echo "$SURVEYS" | while IFS='|' read -r NAME MASK ORTHO; do
    [ -z "$NAME" ] && continue
    WATER="$RASTER_DIR/${NAME}_water_${RES}m.tif"
    PIXCVEC="$OUT_DIR/$V/ortho_${NAME}_${V}_pixcvec.nc"
    RIVERTILE="$OUT_DIR/$V/ortho_${NAME}_${V}.nc"
    OUT_SHP="$QC_DIR/true_node_polygons_${NAME}_${V}.shp"
 
    if [ ! -f "$RIVERTILE" ]; then
      echo "--- $NAME / $V : SKIPPED, no RiverTile (step 2 failed for this one)"
      continue
    fi
    echo
    echo "--- $NAME / $V ---"
    python step3b_true_node_polygons.py \
        "$WATER" "$PIXCVEC" "$OUT_SHP" \
        --rivertile "$RIVERTILE" \
      || echo "!!! step3b FAILED for $NAME / $V"
  done
done
fi
 
# =============================================================================
echo
echo "############################################################"
if [ -n "$FAILED" ]; then
  echo "Finished with problems:$FAILED"
  echo "Scroll up for '!!!' and 'FAILED' lines."
else
  echo "Finished."
fi
echo
echo "Outputs:"
echo "  water rasters   $RASTER_DIR"
echo "  RiverObs        $OUT_DIR"
echo "     combined CSVs: ortho_riverobs_nodes_all.csv, ortho_riverobs_reaches_all.csv"
echo "  node polygons   $QC_DIR"
echo
echo "NEXT (step 5.4, manual): open each true_node_polygons_*.shp in QGIS over"
echo "its orthomosaic, style graduated on 'width', and mark the nodes to drop."
echo "Then run step3c_manual_node_qc.py. Nothing above excludes any node."
echo "############################################################"