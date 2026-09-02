#!/usr/bin/env python3
"""Numerically compare two CalcGW TSV files, ignoring runtime by default."""

import argparse
import csv
import math


def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()
    reference, candidate = load(args.reference), load(args.candidate)
    if len(reference) != len(candidate):
        raise SystemExit(f"row count differs: {len(reference)} != {len(candidate)}")

    failures = []
    for row_index, (left, right) in enumerate(zip(reference, candidate), start=2):
        if left.keys() != right.keys():
            failures.append((row_index, "header", "different", "different"))
            continue
        for field in left:
            if field == "runtime" or left[field] == right[field]:
                continue
            try:
                x, y = float(left[field]), float(right[field])
                equal = math.isclose(x, y, rel_tol=args.rtol, abs_tol=args.atol)
            except ValueError:
                equal = False
            if not equal:
                failures.append((row_index, field, left[field], right[field]))

    if failures:
        for failure in failures[:50]:
            print("row=%s field=%s reference=%s candidate=%s" % failure)
        raise SystemExit(f"FAIL: {len(failures)} values outside tolerance")
    print(f"PASS: {len(reference)} rows agree (runtime ignored)")


if __name__ == "__main__":
    main()
