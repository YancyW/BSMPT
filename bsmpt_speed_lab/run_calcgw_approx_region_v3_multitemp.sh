#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Experimental only.  Keep v1/v2 defaults frozen and inject the v3 components
# solely for this invocation.
exec env \
  BSMPT_APPROX_EXACT_ANCHORS="${script_dir}/approx_exact_direct_anchors_v3.tsv" \
  BSMPT_APPROX_TWO_NEIGHBOR_ANCHORS="${script_dir}/approx_two_neighbor_safe_anchors_v3.tsv" \
  BSMPT_APPROX_RUNNER="${script_dir}/run_calcgw_approx_multitemp_shadow_diagnostic.sh" \
  BSMPT_APPROX_DECISION="${script_dir}/region_v3_multitemp_decision.py" \
  "${script_dir}/run_calcgw_approx_region_v2_nlo_cert.sh" "$@"
