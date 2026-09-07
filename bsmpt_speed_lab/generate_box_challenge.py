#!/usr/bin/env python3
"""Generate deterministic multi-axis corners and target-directed challenges."""

import argparse
import csv
import math


PARAMETERS = ("L1", "L2", "L3", "L4", "L5", "m12sq", "tbeta")

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--anchor-row", type=int, required=True)
parser.add_argument("--target-row", type=int, required=True)
parser.add_argument("--radius", type=float, required=True)
args = parser.parse_args()

with open(args.input, newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames
anchor = rows[args.anchor_row - 1]
target = rows[args.target_row - 1]

generated = []
# Eight balanced Walsh-like sign patterns exercise simultaneous parameter
# extremes without pretending to cover all 2^7 corners.
for index in range(8):
    row = anchor.copy()
    for dimension, name in enumerate(PARAMETERS):
        sign = 1.0 if ((index >> (dimension % 3)) & 1) else -1.0
        value = float(anchor[name])
        row[name] = format(value + sign * args.radius * abs(value), ".17g")
    generated.append(row)

# Four points move toward a known strict-positive target, normalized so the
# largest component reaches the requested fraction of the certified box.
deltas = {
    name: (float(target[name]) - float(anchor[name])) / max(abs(float(anchor[name])), 1e-12)
    for name in PARAMETERS
}
maximum = max(abs(value) for value in deltas.values())
for fraction in (0.25, 0.5, 0.75, 1.0):
    row = anchor.copy()
    for name in PARAMETERS:
        value = float(anchor[name])
        relative = fraction * args.radius * deltas[name] / maximum
        row[name] = format(value + relative * abs(value), ".17g")
    generated.append(row)

with open(args.output, "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(generated)
