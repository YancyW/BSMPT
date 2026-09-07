#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Lab-only C2 prototype. The shadow hook is opt-in and exists only in the
# isolated approximate build; the strict executable and original tree remain
# untouched.
exec env BSMPT_RUN_LEGACY_GRADIENT_SHADOW=1 \
  "${script_dir}/run_calcgw_approx_c2_adaptive_r500_thermal_fast.sh" "$@"
