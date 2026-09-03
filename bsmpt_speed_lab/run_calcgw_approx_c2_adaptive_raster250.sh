#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Next experimental density step.  This is not the guarded default unless it
# separately passes the full acceptance matrix.
exec env BSMPT_BOUNCE_RASTER_INTERVALS=250 \
  "${script_dir}/run_calcgw_approx_central2_adaptive.sh" "$@"
