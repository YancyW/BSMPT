#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Experimental guarded-first-pass candidate.  Add the 64-point adaptive exact
# solution threshold search to central2 + analytic lepton.  Select only through
# run_calcgw_approx_safe.sh after validation; direct output is not accepted.
exec env BSMPT_ADAPTIVE_THRESHOLD=1 BSMPT_ADAPTIVE_THRESHOLD_GRID=64 \
  "${script_dir}/run_calcgw_approx_central2.sh" "$@"
