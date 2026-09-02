#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
binary="${BSMPT_CALCGW_BINARY:-${script_dir}/build-conda-control/bin/CalcGW}"

# Same isolated lab binary, with every switch used by the exact-fast wrapper
# explicitly removed.  This is the corrected dense numerical control, not the
# original install tree.
exec env \
  -u BSMPT_USE_PATH_GEOMETRY_JET \
  -u BSMPT_USE_SPLINE_INTERVAL_HINT \
  -u BSMPT_USE_COMBINED_NUMERICAL_DERIVATIVES \
  -u BSMPT_USE_R2HDM_HIGGS_PAIR \
  -u BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF \
  -u BSMPT_USE_R2HDM_HIGGS_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_GAUGE_PAIR \
  -u BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_QUARK_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_QUARK_FIXED12_DIFF0 \
  -u BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0 \
  -u BSMPT_USE_V1LOOP_THERMAL_CONTEXT \
  -u BSMPT_USE_R2HDM_VTREE_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE \
  -u BSMPT_USE_R2HDM_COUNTERTERM_FLAT_CACHE \
  -u BSMPT_USE_FERMION_LOW4_EXACT \
  -u BSMPT_USE_V1LOOP_MASS_VECTOR_RESERVE \
  -u BSMPT_USE_CACHED_PROFILE_GATE \
  "${binary}" "$@"
