#!/usr/bin/env python3
"""Aggregate N1 qualitative detail TSVs without discarding group identity."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("details", nargs="+")
parser.add_argument("--output-prefix", required=True)
args = parser.parse_args()

rows = []
for filename in args.details:
    group = Path(filename).name.removesuffix("_details.tsv")
    with open(filename, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({"group": group, **row})

matrix = Counter((r["strict_class"], r["raw_approx_class"]) for r in rows)
errors = Counter(r["raw_error_kind"] for r in rows)
strict_stages = Counter(r["strict_failure_stage"] for r in rows if r["strict_class"] != "positive")
approx_stages = Counter(r["raw_approx_failure_stage"] for r in rows if r["raw_approx_class"] != "positive")
summary = {
    "rows": len(rows),
    "raw_approx_confusion": {
        "true_positive": matrix[("positive", "positive")],
        "false_negative": matrix[("positive", "fail_or_nonpositive")],
        "false_positive": matrix[("fail_or_nonpositive", "positive")],
        "true_negative": matrix[("fail_or_nonpositive", "fail_or_nonpositive")],
    },
    "qualitative_error_kinds": dict(errors),
    "raw_qualitative_accuracy": sum(r["raw_error_kind"] == "none" for r in rows) / len(rows),
    "strict_failure_stages": dict(strict_stages),
    "raw_approx_failure_stages": dict(approx_stages),
    "warning": "sample-specific; do not generalize to full BSMPT parameter space",
}

prefix = Path(args.output_prefix)
with prefix.with_name(prefix.name + "_details.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
prefix.with_name(prefix.name + "_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, sort_keys=True))
