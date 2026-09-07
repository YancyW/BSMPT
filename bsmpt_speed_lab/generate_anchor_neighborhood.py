#!/usr/bin/env python3
"""Generate deterministic local BSMPT input perturbations around one TSV row."""

import argparse
import csv
import math
import random


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--row", type=int, required=True, help="1-based data row")
parser.add_argument("--delta", type=float, default=1e-4)
parser.add_argument("--random-directions", type=int, default=8)
parser.add_argument(
    "--max-points",
    type=int,
    help="Deterministically down-sample the generated axis/random directions",
)
parser.add_argument("--seed", type=int, default=20260905)
args = parser.parse_args()

with open(args.input, newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames
anchor = rows[args.row - 1]
vary = [name for name in fields if name != "yuktype"]

generated = []
for name in vary:
    for sign in (-1.0, 1.0):
        row = anchor.copy()
        value = float(row[name])
        scale = max(abs(value), 1.0)
        row[name] = format(value + sign * args.delta * scale, ".17g")
        generated.append(row)

rng = random.Random(args.seed)
for _ in range(args.random_directions):
    direction = [rng.gauss(0.0, 1.0) for _ in vary]
    norm = math.sqrt(sum(value * value for value in direction))
    row = anchor.copy()
    for name, component in zip(vary, direction):
        value = float(row[name])
        scale = max(abs(value), 1.0)
        row[name] = format(value + args.delta * scale * component / norm, ".17g")
    generated.append(row)

if args.max_points is not None:
    if args.max_points < 1:
        parser.error("--max-points must be positive")
    if args.max_points < len(generated):
        # Cover the complete ordered direction list instead of retaining only
        # the first parameters' axis directions.
        indices = [
            round(i * (len(generated) - 1) / (args.max_points - 1))
            for i in range(args.max_points)
        ] if args.max_points > 1 else [0]
        generated = [generated[index] for index in indices]

with open(args.output, "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(generated)
