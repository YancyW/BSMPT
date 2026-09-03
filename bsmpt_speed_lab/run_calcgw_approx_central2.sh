#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# More aggressive research profile: retain the guarded analytic lepton first
# pass and replace the numerical gradient stencil with central two-point
# differences.  Never use this wrapper directly for accepted scan output;
# select it through run_calcgw_approx_safe.sh so risky rows fall back exactly.
exec env BSMPT_USE_CENTRAL2_GRADIENT=1 \
  "${script_dir}/run_calcgw_approx_guarded.sh" "$@"
