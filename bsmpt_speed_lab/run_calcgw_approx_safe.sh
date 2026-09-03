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
  output="${args[2]}"
  output_index=2
fi
if [[ -z "$output" || $output_index -lt 0 ]]; then
  echo "run_calcgw_approx_safe.sh requires a positional output or --output=PATH" >&2
  exit 2
fi

anchors="${BSMPT_APPROX_EXACT_ANCHORS:-${script_dir}/approx_exact_direct_anchors.tsv}"
if [[ -f "$anchors" ]]; then
  if reason="$(python3 "${script_dir}/approx_input_prefilter.py" \
      --anchors "$anchors" -- "${args[@]}")"; then
    echo "approx-safe: direct exact: ${reason}" >&2
    exec "${script_dir}/run_calcgw_exact_fast.sh" "$@"
  fi
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/bsmpt-approx-safe.XXXXXX")"
trap 'rm -rf -- "$tmpdir"' EXIT
approx_output="${tmpdir}/approx.tsv"
if [[ "${args[$output_index]}" == --output=* ]]; then
  args[$output_index]="--output=${approx_output}"
else
  args[$output_index]="$approx_output"
fi

first_pass="${BSMPT_APPROX_FIRST_PASS:-${script_dir}/run_calcgw_approx_c2_adaptive_raster500.sh}"
case "$first_pass" in
  "${script_dir}"/*) ;;
  *)
    echo "BSMPT_APPROX_FIRST_PASS must name a wrapper inside ${script_dir}" >&2
    exit 2
    ;;
esac
"${first_pass}" "${args[@]}"

check_args=("$approx_output" --snr-floor "${BSMPT_APPROX_SNR_FLOOR:-1e-20}")
if [[ -n "${BSMPT_APPROX_SNR_CUTS:-}" ]]; then
  IFS=',' read -r -a cuts <<< "${BSMPT_APPROX_SNR_CUTS}"
  for cut in "${cuts[@]}"; do check_args+=(--snr-cut "$cut"); done
fi
check_args+=(--cut-margin "${BSMPT_APPROX_CUT_MARGIN:-0.15}")

if reason="$(python3 "${script_dir}/approx_needs_fallback.py" "${check_args[@]}")"; then
  echo "approx-safe: exact fallback: ${reason}" >&2
  "${script_dir}/run_calcgw_exact_fast.sh" "$@"
else
  echo "approx-safe: accepted approximate result" >&2
  mv -- "$approx_output" "$output"
fi
