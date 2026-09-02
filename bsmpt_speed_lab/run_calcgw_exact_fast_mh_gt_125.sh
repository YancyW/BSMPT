#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Kept as a compatibility entry point.  The formerly branch-specific flat
# CounterTerm and fixed fermion paths are now validated and beneficial in the
# generic exact wrapper as well.
exec "${script_dir}/run_calcgw_exact_fast.sh" "$@"
