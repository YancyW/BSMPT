#!/usr/bin/env python3
"""Audit the low-temperature-count late-failure certificate against labels."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--case", nargs=3, action="append", metavar=("LABEL", "DETAILS", "DIAGNOSTICS"), required=True)
parser.add_argument("--output-prefix", required=True)
parser.add_argument("--max-temperatures", type=int, default=6)
args = parser.parse_args()

rows = []
for label, details_path, diagnostics_path in args.case:
    temperatures = defaultdict(set)
    with open(diagnostics_path, newline="") as stream:
        for record in csv.DictReader(stream, delimiter="\t"):
            parts = [record["diagnostic"]] + (record.get(None) or [])
            if not parts or parts[0] != "BSMPT_BOUNCE_CERT":
                continue
            fields = dict(item.split("=", 1) for item in parts[2:] if "=" in item)
            if "T" in fields:
                temperatures[int(record["input_row"])].add(fields["T"])
    with open(details_path, newline="") as stream:
        for detail in csv.DictReader(stream, delimiter="\t"):
            row_number = int(detail["row"])
            approximate_failure = detail["raw_approx_class"] != "positive"
            count = len(temperatures[row_number])
            accepted = approximate_failure and 0 < count <= args.max_temperatures
            unsafe = accepted and detail["strict_class"] == "positive"
            rows.append({
                "case": label,
                "row": row_number,
                "strict_class": detail["strict_class"],
                "raw_approx_class": detail["raw_approx_class"],
                "raw_error_kind": detail["raw_error_kind"],
                "temperature_count": count,
                "would_accept": str(accepted).lower(),
                "unsafe_accept": str(unsafe).lower(),
            })

prefix = Path(args.output_prefix)
with prefix.with_name(prefix.name + "_details.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(
        stream, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
summary = {
    "rows": len(rows),
    "would_accept": sum(row["would_accept"] == "true" for row in rows),
    "unsafe_accept": sum(row["unsafe_accept"] == "true" for row in rows),
    "max_temperatures": args.max_temperatures,
    "warning": "certificate is valid only together with explicit type-I anchor gating",
}
prefix.with_name(prefix.name + "_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, sort_keys=True))
