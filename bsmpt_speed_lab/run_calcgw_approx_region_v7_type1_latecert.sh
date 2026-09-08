#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Type-I late-failure certificate. BOUNCE_CERT only records quantities already
# computed by the approximate solver; no extra action or bounce evaluation.
exec env \
  BSMPT_EMIT_BOUNCE_CERTIFICATE=1 \
  BSMPT_APPROX_LATE_FAILURE_ANCHORS="${script_dir}/approx_late_failure_safe_anchors_v7.tsv" \
  "${script_dir}/run_calcgw_approx_region_v6_type1_safe.sh" "$@"
