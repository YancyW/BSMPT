#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Minimum density accepted by the experimental raster implementation.  This
# is an exploratory rejection/acceptance candidate, never the strict path.
exec env BSMPT_BOUNCE_RASTER_INTERVALS=100 \
  "${script_dir}/run_calcgw_approx_central2_adaptive.sh" "$@"
