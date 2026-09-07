#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Type-I safety revision: the broad-C row3 neighborhood remains qualitatively
# positive but changes transition schema, so route it directly to strict.
exec env \
  BSMPT_APPROX_EXACT_ANCHORS="${script_dir}/approx_exact_direct_anchors_v6.tsv" \
  "${script_dir}/run_calcgw_approx_region_v5_certified.sh" "$@"
