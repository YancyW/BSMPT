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
  echo "region-v1 requires a positional output or --output=PATH" >&2
  exit 2
fi

anchors="${BSMPT_APPROX_EXACT_ANCHORS:-${script_dir}/approx_exact_direct_anchors.tsv}"
approx_runner="${BSMPT_APPROX_RUNNER:-${script_dir}/run_calcgw_approx_legacy_shadow.sh}"
decision_script="${BSMPT_APPROX_DECISION:-${script_dir}/region_v1_decision.py}"
two_neighbor_anchors="${BSMPT_APPROX_TWO_NEIGHBOR_ANCHORS:-}"
all_three_anchors="${BSMPT_APPROX_ALL_THREE_ANCHORS:-}"
coex_anchors="${BSMPT_APPROX_COEX_ANCHORS:-}"
late_failure_anchors="${BSMPT_APPROX_LATE_FAILURE_ANCHORS:-}"
if reason="$(python3 "${script_dir}/approx_input_prefilter.py" --anchors "$anchors" -- "${args[@]}")"; then
  echo "region-v1: direct exact: ${reason}" >&2
  exec "${script_dir}/run_calcgw_exact_fast.sh" "$@"
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/bsmpt-region-v1.XXXXXX")"
trap 'rm -rf -- "$tmpdir"' EXIT
approx_output="${tmpdir}/approx.tsv"
diagnostics="${tmpdir}/shadow.log"
if [[ "${args[$output_index]}" == --output=* ]]; then
  args[$output_index]="--output=${approx_output}"
else
  args[$output_index]="$approx_output"
fi

"$approx_runner" "${args[@]}" 2>"$diagnostics"
allow_two_neighbor=0
if [[ -n "$two_neighbor_anchors" ]] &&
   python3 "${script_dir}/approx_input_prefilter.py" \
     --anchors "$two_neighbor_anchors" -- "${args[@]}" >/dev/null; then
  allow_two_neighbor=1
fi
allow_all_three=1
if [[ -n "$all_three_anchors" ]]; then
  allow_all_three=0
  if python3 "${script_dir}/approx_input_prefilter.py" \
       --anchors "$all_three_anchors" -- "${args[@]}" >/dev/null; then
    allow_all_three=1
  fi
fi
allow_coex=1
if [[ -n "$coex_anchors" ]]; then
  allow_coex=0
  if python3 "${script_dir}/approx_input_prefilter.py" \
       --anchors "$coex_anchors" -- "${args[@]}" >/dev/null; then
    allow_coex=1
  fi
fi
allow_late_failure=0
if [[ -n "$late_failure_anchors" ]] &&
   python3 "${script_dir}/approx_input_prefilter.py" \
     --anchors "$late_failure_anchors" -- "${args[@]}" >/dev/null; then
  allow_late_failure=1
fi
if reason="$(BSMPT_ALLOW_TWO_NEIGHBOR="$allow_two_neighbor" \
  BSMPT_ALLOW_ALL_THREE="$allow_all_three" BSMPT_ALLOW_COEX="$allow_coex" \
  BSMPT_ALLOW_LATE_FAILURE="$allow_late_failure" \
  python3 "$decision_script" "$approx_output" "$diagnostics")"; then
  echo "region-v1: ${reason}" >&2
  mv -- "$approx_output" "$output"
else
  echo "region-v1: exact fallback: ${reason}" >&2
  "${script_dir}/run_calcgw_exact_fast.sh" "$@"
fi
