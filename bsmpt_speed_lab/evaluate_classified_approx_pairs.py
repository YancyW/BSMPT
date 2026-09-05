#!/usr/bin/env python3
"""Evaluate stored strict classified rows against newly computed approximations."""

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


INPUT_NAMES = {"yuktype", "L1", "L2", "L3", "L4", "L5", "m12sq", "m12sq_calgw", "tbeta"}


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata")
    parser.add_argument("strict")
    parser.add_argument("approx")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--snr-limit", type=float, default=0.10)
    parser.add_argument("--snr-floor", type=float, default=1e-20)
    args = parser.parse_args()

    def load(path):
        with open(path, newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    metadata, strict, approx = load(args.metadata), load(args.strict), load(args.approx)
    if not (len(metadata) == len(strict) == len(approx)):
        raise SystemExit(f"row mismatch metadata={len(metadata)} strict={len(strict)} approx={len(approx)}")

    common = [k for k in approx[0] if k in strict[0] and k not in INPUT_NAMES]
    status_fields = [k for k in common if k.startswith("status_")] + (["transition_history"] if "transition_history" in common else [])
    numeric_fields = [k for k in common if k not in status_fields and k != "runtime"]
    snr_fields = [k for k in numeric_fields if k.startswith("SNR(")]
    results = []
    for i, (meta, ref, cand) in enumerate(zip(metadata, strict, approx), 1):
        reasons = []
        changed_status = [k for k in status_fields if ref.get(k, "") != cand.get(k, "")]
        if changed_status:
            reasons.append("status_or_history_mismatch")
        mask_changes = []
        max_snr_error = 0.0
        max_snr_field = ""
        informative_snr = 0
        for key in numeric_fields:
            rv, cv = number(ref.get(key)), number(cand.get(key))
            if (rv is None) != (cv is None):
                mask_changes.append(key)
            if key not in snr_fields or rv is None or cv is None:
                continue
            if abs(rv) < args.snr_floor:
                if abs(cv) >= args.snr_floor:
                    max_snr_error, max_snr_field = math.inf, key
                continue
            informative_snr += 1
            error = abs(cv - rv) / abs(rv)
            if error > max_snr_error:
                max_snr_error, max_snr_field = error, key
        if mask_changes:
            reasons.append("finite_mask_mismatch")
        if max_snr_error > args.snr_limit:
            reasons.append("snr_error_gt_limit")
        strict_total = number(meta.get("strict_snr"))
        weak = strict_total is None or abs(strict_total) < args.snr_floor or informative_snr == 0
        if weak:
            reasons.append("weak_or_no_strict_snr")
        pair_exact_enough = not any(r in reasons for r in ("status_or_history_mismatch", "finite_mask_mismatch", "snr_error_gt_limit"))
        safe_for_routing = pair_exact_enough and not weak
        results.append({
            **meta,
            "pair_exact_enough": str(pair_exact_enough).lower(),
            "safe_for_routing": str(safe_for_routing).lower(),
            "max_snr_relative_error": max_snr_error,
            "max_snr_field": max_snr_field,
            "status_mismatch_count": len(changed_status),
            "finite_mask_mismatch_count": len(mask_changes),
            "reasons": ";".join(reasons) if reasons else "ok",
        })

    prefix = Path(args.output_prefix)
    detail_path = prefix.with_name(prefix.name + "_details.tsv")
    summary_path = prefix.with_name(prefix.name + "_summary.json")
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    overall = Counter(r["reasons"] for r in results)
    grouped = defaultdict(Counter)
    for r in results:
        key = "|".join((r["validation_stratum"], r["branch"], r["machine"]))
        grouped[key]["n"] += 1
        grouped[key]["pair_exact_enough"] += r["pair_exact_enough"] == "true"
        grouped[key]["safe_for_routing"] += r["safe_for_routing"] == "true"
    summary = {
        "rows": len(results),
        "snr_limit": args.snr_limit,
        "snr_floor": args.snr_floor,
        "pair_exact_enough": sum(r["pair_exact_enough"] == "true" for r in results),
        "safe_for_routing": sum(r["safe_for_routing"] == "true" for r in results),
        "reason_combinations": dict(overall),
        "groups": {k: dict(v) for k, v in sorted(grouped.items())},
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
