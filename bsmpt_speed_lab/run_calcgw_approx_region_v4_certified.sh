#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Certified-domain router: every approximate terminal decision is either
# backed by the exact NLO prefix or restricted to an explicit validated table.
exec env \
  BSMPT_APPROX_ALL_THREE_ANCHORS="${script_dir}/approx_all_three_safe_anchors_v4.tsv" \
  BSMPT_APPROX_COEX_ANCHORS="${script_dir}/approx_coex_safe_anchors_v4.tsv" \
  BSMPT_REQUIRE_EXACT_NLO_PREFIX=1 \
  "${script_dir}/run_calcgw_approx_region_v3_multitemp.sh" "$@"
