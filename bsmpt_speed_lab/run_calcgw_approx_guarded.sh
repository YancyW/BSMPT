#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
binary="${BSMPT_CALCGW_BINARY:-${script_dir}/build-conda-control/bin/CalcGW}"

# Research-only approximate profile.  It starts from the validated exact-fast
# configuration and replaces only the diff==0 R2HDM lepton eigenspectrum with
# its analytic block spectrum.  The C++ implementation checks the complete
# sparse/symmetric structure and falls back to Eigen if that structure is not
# present.  Scan-level boundary fallback is handled by the validation policy;
# this wrapper is not a replacement for final exact confirmation.
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
  BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0=1 \
  BSMPT_USE_V1LOOP_THERMAL_CONTEXT=1 \
  BSMPT_USE_R2HDM_VTREE_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE=1 \
  BSMPT_USE_R2HDM_COUNTERTERM_FLAT_CACHE=1 \
  BSMPT_USE_FERMION_LOW4_EXACT=1 \
  BSMPT_USE_V1LOOP_MASS_VECTOR_RESERVE=1 \
  BSMPT_USE_CACHED_PROFILE_GATE=1 \
  "${binary}" "$@"
