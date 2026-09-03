#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
binary="${BSMPT_APPROX_THERMAL_BINARY:-${script_dir}/build-approx-thermal/bin/CalcGW}"

if [[ ! -x "$binary" ]]; then
  echo "thermal-fast approximate binary not found: ${binary}" >&2
  echo "build the isolated bsmpt_speed_lab/build-approx-thermal tree first" >&2
  exit 2
fi

# Non-strict guarded first pass only.  The exact fallback is launched by the
# parent safe runner after this process exits, so this binary override and the
# approximate thermal arithmetic cannot leak into the strict calculation.
exec env BSMPT_CALCGW_BINARY="$binary" BSMPT_USE_THERMAL_FAST_POWERS=1 \
  "${script_dir}/run_calcgw_approx_c2_adaptive_raster500.sh" "$@"
