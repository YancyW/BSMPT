#!/usr/bin/env python3
"""Generate the frozen type-I v7 A9/B10 edge challenge panel."""

import argparse
import csv
import math
import random


parser = argparse.ArgumentParser()
parser.add_argument("a_input")
parser.add_argument("b_input")
parser.add_argument("output")
args = parser.parse_args()


def load_row(path, number):
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames[:8], rows[number - 1]


fields, a9 = load_row(args.a_input, 9)
b_fields, b10 = load_row(args.b_input, 10)
if fields != b_fields or fields != ["yuktype", "L1", "L2", "L3", "L4", "L5", "m12sq", "tbeta"]:
    raise SystemExit("unexpected input schema")
if a9["yuktype"] != "1" or b10["yuktype"] != "1":
    raise SystemExit("challenge anchors must be Yukawa type-I")


def displaced(anchor, direction, radius):
    row = {name: anchor[name] for name in fields}
    norm = max(abs(value) for value in direction)
    for name, component in zip(fields[1:], direction):
        value = float(anchor[name])
        scale = max(abs(value), 1e-12)
        row[name] = format(value + radius * scale * component / norm, ".17g")
    return row


generated = []
for radius in (4e-5, 6e-5, 8e-5, 9.5e-5):
    direction = [0.0] * 7
    direction[3] = -1.0  # L4
    generated.append(displaced(a9, direction, radius))

rng_a = random.Random(9080801)
for _ in range(2):
    direction = [rng_a.gauss(0.0, 1.0) for _ in range(7)]
    direction[3] = -abs(direction[3])
    direction[6] = -abs(direction[6])
    generated.append(displaced(a9, direction, 1.05e-4))
    generated.append(displaced(a9, [-value for value in direction], 1.05e-4))

rng_b = random.Random(9080802)
for _ in range(2):
    direction = [rng_b.gauss(0.0, 1.0) for _ in range(7)]
    generated.append(displaced(b10, direction, 3.0e-4))
    generated.append(displaced(b10, [-value for value in direction], 3.0e-4))

with open(args.output, "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(generated)
