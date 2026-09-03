#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Exploratory guarded-first-pass alternative to central2.  The analytic
# gradient has separate historical error behavior and must pass independently.
exec env BSMPT_USE_ANALYTIC_GRADIENT=1 BSMPT_BOUNCE_RASTER_INTERVALS=500 \
  BSMPT_ADAPTIVE_THRESHOLD=1 BSMPT_ADAPTIVE_THRESHOLD_GRID=64 \
  "${script_dir}/run_calcgw_approx_guarded.sh" "$@"
