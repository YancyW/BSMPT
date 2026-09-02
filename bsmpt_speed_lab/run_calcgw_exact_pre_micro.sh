#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
binary="${BSMPT_CALCGW_BINARY:-${script_dir}/build-conda-control/bin/CalcGW}"

# Exact control immediately before the three most recent micro-optimizations:
# low-order fermion-series unrolling, pair-result reserve, and cached profiler
# gate are deliberately absent.  All earlier accepted exact switches remain.
exec env \
  BSMPT_USE_PATH_GEOMETRY_JET=1 \
  BSMPT_USE_SPLINE_INTERVAL_HINT=1 \
  BSMPT_USE_COMBINED_NUMERICAL_DERIVATIVES=1 \
  BSMPT_USE_R2HDM_HIGGS_PAIR=1 \
  BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF=1 \
  BSMPT_USE_R2HDM_HIGGS_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_GAUGE_PAIR=1 \
  BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_QUARK_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_QUARK_FIXED12_DIFF0=1 \
  BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0=1 \
  BSMPT_USE_V1LOOP_THERMAL_CONTEXT=1 \
  BSMPT_USE_R2HDM_VTREE_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_COUNTERTERM_FLAT_CACHE=1 \
  "${binary}" "$@"
