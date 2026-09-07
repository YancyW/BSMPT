#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Diagnostic only: production region-v1 does not set this variable and
# therefore retains its frozen single-temperature shadow behavior.
exec env BSMPT_LEGACY_SHADOW_TEMPERATURES="${BSMPT_LEGACY_SHADOW_TEMPERATURES:-3}" \
  "${script_dir}/run_calcgw_approx_legacy_shadow.sh" "$@"
