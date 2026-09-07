#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
args=("$@")
output=""
output_index=-1
for i in "${!args[@]}"; do
  case "${args[$i]}" in
    --output=*) output="${args[$i]#--output=}"; output_index="$i" ;;
  esac
done
if [[ -z "$output" && ${#args[@]} -ge 3 && "${args[0]}" != --* ]]; then
  output="${args[2]}"; output_index=2
fi
if [[ -z "$output" || $output_index -lt 0 ]]; then
  echo "region-v2 requires a positional output or --output=PATH" >&2
  exit 2
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/bsmpt-region-v2.XXXXXX")"
trap 'rm -rf -- "$tmpdir"' EXIT
nlo_output="${tmpdir}/exact-nlo.tsv"
nlo_args=("${args[@]}")
if [[ "${nlo_args[$output_index]}" == --output=* ]]; then
  nlo_args[$output_index]="--output=${nlo_output}"
else
  nlo_args[$output_index]="$nlo_output"
fi

BSMPT_NLO_ONLY=1 "${script_dir}/run_calcgw_exact_fast.sh" "${nlo_args[@]}"
if reason="$(python3 "${script_dir}/nlo_only_decision.py" "$nlo_output")"; then
  echo "region-v2: ${reason}" >&2
  mv -- "$nlo_output" "$output"
  exit 0
fi

echo "region-v2: ${reason}:continue_region_v1" >&2
exec "${script_dir}/run_calcgw_approx_region_v1.sh" "$@"
