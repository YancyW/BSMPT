#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Experimental guarded-first-pass candidate.  Reduce the bounce dV/dl raster
# from 1000 to 500 intervals on top of central2 + adaptive threshold.
exec env BSMPT_BOUNCE_RASTER_INTERVALS=500 \
  "${script_dir}/run_calcgw_approx_central2_adaptive.sh" "$@"
