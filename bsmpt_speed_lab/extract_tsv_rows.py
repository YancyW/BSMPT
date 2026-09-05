#!/usr/bin/env python3
"""Extract 1-based data rows from a TSV while preserving LF line endings."""
import argparse
import csv


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("rows", nargs="+", type=int)
args = parser.parse_args()
wanted = set(args.rows)
with open(args.input, newline="") as source, open(args.output, "w", newline="") as target:
    reader = csv.reader(source, delimiter="\t")
    writer = csv.writer(target, delimiter="\t", lineterminator="\n")
    writer.writerow(next(reader))
    for index, row in enumerate(reader, 1):
        if index in wanted:
            writer.writerow(row)
