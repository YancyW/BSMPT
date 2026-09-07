#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Region-v5 keeps every v4 guard and only expands the positive certificate to
# independently validated Yukawa type-2 and type-4 neighborhoods.
exec env \
  BSMPT_APPROX_ALL_THREE_ANCHORS="${script_dir}/approx_all_three_safe_anchors_v5.tsv" \
  BSMPT_APPROX_COEX_ANCHORS="${script_dir}/approx_coex_safe_anchors_v4.tsv" \
  BSMPT_REQUIRE_EXACT_NLO_PREFIX=1 \
  "${script_dir}/run_calcgw_approx_region_v3_multitemp.sh" "$@"
