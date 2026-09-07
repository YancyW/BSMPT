#!/usr/bin/env python3
"""Compare strict, raw approximate and optional guarded CalcGW outcomes.

The primary contract is qualitative and bidirectional: finite SNR > 0 versus
failure/non-positive.  Amplitude differences are secondary diagnostics only.
"""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


INPUT_FIELDS = {"yuktype", "L1", "L2", "L3", "L4", "L5", "m12sq", "tbeta"}
EARLY_STAGES = (
    ("status_nlo_stability", "nlo"),
    ("status_tracing", "tracing"),
    ("status_coex_pairs", "coexistence"),
)
TRANSITION_STAGES = (
    ("status_crit_{}", "critical"),
    ("status_bounce_sol_{}", "bounce"),
    ("status_nucl_approx_{}", "nucleation_approx"),
    ("status_nucl_{}", "nucleation"),
    ("status_perc_{}", "percolation"),
    ("status_compl_{}", "completion"),
    ("status_gw_{}", "gw"),
)
INACTIVE = {"", "nan", "not_set"}


def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def transition_ids(row):
    ids = set()
    for key in row:
        if key.startswith("status_gw_"):
            suffix = key.removeprefix("status_gw_")
            if suffix.isdigit():
                ids.add(int(suffix))
    return sorted(ids)


def total_snr(row, index):
    for key in (f"SNR(LISA-3yrs)_{index}", f"SNR_total_{index}"):
        if key in row:
            return finite(row.get(key))
    if index == 0:
        for key in ("SNR_total", "SNR"):
            if key in row:
                return finite(row.get(key))
    return None


def outcome(row):
    positives = []
    for index in transition_ids(row):
        snr = total_snr(row, index)
        if row.get(f"status_gw_{index}", "") == "success" and snr is not None and snr > 0:
            positives.append((index, snr))
    if positives:
        return "positive", "positive", positives

    for key, stage in EARLY_STAGES:
        value = (row.get(key) or "").strip()
        if key == "status_nlo_stability" and value not in {"success", *INACTIVE}:
            return "fail_or_nonpositive", f"{stage}:{value}", []
        if key != "status_nlo_stability" and value not in {"success", *INACTIVE}:
            return "fail_or_nonpositive", f"{stage}:{value}", []

    for index in transition_ids(row):
        for pattern, stage in TRANSITION_STAGES:
            key = pattern.format(index)
            value = (row.get(key) or "").strip()
            if value not in {"success", *INACTIVE}:
                return "fail_or_nonpositive", f"{stage}_{index}:{value}", []
    return "fail_or_nonpositive", "no_finite_positive_snr", []


def secondary_snr_error(strict_row, approx_row):
    if transition_ids(strict_row) != transition_ids(approx_row):
        return 0.0, "incompatible_transition_ids", 0
    if strict_row.get("transition_history") != approx_row.get("transition_history"):
        return 0.0, "incompatible_transition_history", 0
    maximum = 0.0
    field = ""
    compared = 0
    for key, ref_text in strict_row.items():
        if not key.startswith("SNR(") or key not in approx_row:
            continue
        ref, cand = finite(ref_text), finite(approx_row.get(key))
        if ref is None or cand is None or ref <= 0:
            continue
        compared += 1
        error = abs(cand - ref) / abs(ref)
        if error > maximum:
            maximum, field = error, key
    return maximum, field, compared


def same_nonruntime(left, right):
    common = set(left) & set(right) - INPUT_FIELDS - {"runtime"}
    return bool(common) and all(left.get(key, "") == right.get(key, "") for key in common)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("strict")
    parser.add_argument("approx")
    parser.add_argument("--guarded")
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    strict_rows, approx_rows = load(args.strict), load(args.approx)
    guarded_rows = load(args.guarded) if args.guarded else None
    lengths = [len(strict_rows), len(approx_rows)] + ([len(guarded_rows)] if guarded_rows is not None else [])
    if len(set(lengths)) != 1:
        raise SystemExit(f"row-count mismatch: {lengths}")

    details = []
    raw_matrix = Counter()
    guarded_matrix = Counter()
    for index, (strict_row, approx_row) in enumerate(zip(strict_rows, approx_rows), 1):
        mismatched_inputs = [
            key for key in INPUT_FIELDS
            if strict_row.get(key) != approx_row.get(key)
        ]
        if mismatched_inputs:
            raise SystemExit(
                f"row {index} input mismatch: {','.join(sorted(mismatched_inputs))}"
            )
        strict_class, strict_stage, _ = outcome(strict_row)
        approx_class, approx_stage, _ = outcome(approx_row)
        raw_matrix[(strict_class, approx_class)] += 1
        if strict_class == "positive" and approx_class != "positive":
            error_kind = "false_negative"
        elif strict_class != "positive" and approx_class == "positive":
            error_kind = "false_positive"
        else:
            error_kind = "none"
        max_error, max_field, compared = secondary_snr_error(strict_row, approx_row)
        record = {
            "row": index,
            "strict_class": strict_class,
            "raw_approx_class": approx_class,
            "raw_error_kind": error_kind,
            "strict_failure_stage": strict_stage,
            "raw_approx_failure_stage": approx_stage,
            "raw_qualitative_match": str(error_kind == "none").lower(),
            "secondary_max_snr_relative_error": max_error,
            "secondary_max_snr_field": max_field,
            "secondary_snr_fields_compared": compared,
        }
        if guarded_rows is not None:
            guarded_row = guarded_rows[index - 1]
            guarded_class, guarded_stage, _ = outcome(guarded_row)
            guarded_matrix[(strict_class, guarded_class)] += 1
            record.update({
                "guarded_class": guarded_class,
                "guarded_failure_stage": guarded_stage,
                "guarded_qualitative_match": str(guarded_class == strict_class).lower(),
                "guarded_equals_strict_nonruntime": str(same_nonruntime(strict_row, guarded_row)).lower(),
                "guarded_equals_raw_nonruntime": str(same_nonruntime(approx_row, guarded_row)).lower(),
                "fallback_inferred": str(
                    same_nonruntime(strict_row, guarded_row)
                    and not same_nonruntime(approx_row, guarded_row)
                ).lower(),
            })
        details.append(record)

    def matrix_dict(matrix):
        return {
            "true_positive": matrix[("positive", "positive")],
            "false_negative": matrix[("positive", "fail_or_nonpositive")],
            "false_positive": matrix[("fail_or_nonpositive", "positive")],
            "true_negative": matrix[("fail_or_nonpositive", "fail_or_nonpositive")],
        }

    summary = {
        "rows": len(details),
        "primary_contract": "strict positive iff approximate positive",
        "raw_approx_confusion": matrix_dict(raw_matrix),
        "raw_qualitative_errors": sum(d["raw_error_kind"] != "none" for d in details),
        "secondary_snr_amplitude_over_10pct": sum(
            d["strict_class"] == d["raw_approx_class"] == "positive"
            and d["secondary_max_snr_relative_error"] > 0.10
            for d in details
        ),
        "strict_sha256": hashlib.sha256(Path(args.strict).read_bytes()).hexdigest(),
        "approx_sha256": hashlib.sha256(Path(args.approx).read_bytes()).hexdigest(),
        "secondary_incompatible_transition_schema": sum(
            d["secondary_max_snr_field"].startswith("incompatible_transition_")
            for d in details
        ),
    }
    if guarded_rows is not None:
        summary["guarded_confusion"] = matrix_dict(guarded_matrix)
        summary["guarded_qualitative_errors"] = sum(d["guarded_qualitative_match"] != "true" for d in details)
        summary["fallback_inferred"] = sum(d["fallback_inferred"] == "true" for d in details)

    prefix = Path(args.output_prefix)
    detail_path = prefix.with_name(prefix.name + "_details.tsv")
    summary_path = prefix.with_name(prefix.name + "_summary.json")
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(details[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(details)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
