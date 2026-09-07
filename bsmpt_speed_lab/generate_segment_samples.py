#!/usr/bin/env python3
"""Generate deterministic interior samples between two BSMPT TSV rows."""

import argparse
import csv


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--row-a", type=int, required=True, help="1-based data row")
parser.add_argument("--row-b", type=int, required=True, help="1-based data row")
parser.add_argument("--count", type=int, default=10)
args = parser.parse_args()

with open(args.input, newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames

if args.count < 1:
    parser.error("--count must be positive")
left, right = rows[args.row_a - 1], rows[args.row_b - 1]
generated = []
for index in range(1, args.count + 1):
    fraction = index / (args.count + 1)
    row = {}
    for name in fields:
        if name == "yuktype":
            if left[name] != right[name]:
                parser.error("segment endpoints have different yuktype")
            row[name] = left[name]
        else:
            value = (1.0 - fraction) * float(left[name]) + fraction * float(right[name])
            row[name] = format(value, ".17g")
    generated.append(row)

with open(args.output, "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(generated)
