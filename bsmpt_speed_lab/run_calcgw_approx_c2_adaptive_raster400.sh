#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Experimental midpoint between the validated raster500 default and the
# rejected lower-density raster250 ablation.  Guarded validation is required.
exec env BSMPT_BOUNCE_RASTER_INTERVALS=400 \
  "${script_dir}/run_calcgw_approx_central2_adaptive.sh" "$@"
