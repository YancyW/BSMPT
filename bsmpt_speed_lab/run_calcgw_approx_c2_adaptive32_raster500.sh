#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Exploratory guarded-first-pass profile.  Keep the current central2 and
# raster-500 choices, but halve the adaptive exact-solution threshold grid.
# This wrapper is not an accepted-output entry point; use only for validation.
exec env BSMPT_ADAPTIVE_THRESHOLD=1 BSMPT_ADAPTIVE_THRESHOLD_GRID=32 \
  BSMPT_BOUNCE_RASTER_INTERVALS=500 \
  "${script_dir}/run_calcgw_approx_central2.sh" "$@"
