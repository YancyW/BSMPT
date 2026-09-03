#!/usr/bin/env python3
"""Return 0 when a single CalcGW input matches a known exact-direct region."""

import argparse
import csv
import math
from pathlib import Path


PARAMETERS = ("L1", "L2", "L3", "L4", "L5", "m12sq", "tbeta")


def calcgw_selection(arguments):
    named = {}
    positional = []
    for value in arguments:
        if value.startswith("--") and "=" in value:
            key, item = value[2:].split("=", 1)
            named[key.lower()] = item
        elif not value.startswith("--"):
            positional.append(value)
    input_path = named.get("input")
    first = named.get("firstline")
    last = named.get("lastline")
    if input_path is None and len(positional) >= 2:
        input_path = positional[1]
    if first is None and len(positional) >= 4:
        first = positional[3]
    if last is None and len(positional) >= 5:
        last = positional[4]
    if input_path is None or first is None or last is None:
        return None
    try:
        first_number, last_number = int(first), int(last)
    except ValueError:
        return None
    if first_number != last_number or first_number < 2:
        return None
    return Path(input_path), first_number


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchors", required=True)
    parser.add_argument("calcgw_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    raw = args.calcgw_args
    if raw and raw[0] == "--":
        raw = raw[1:]
    selection = calcgw_selection(raw)
    if selection is None:
        raise SystemExit(1)
    input_path, line_number = selection
    with input_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row_index = line_number - 2
    if row_index < 0 or row_index >= len(rows):
        raise SystemExit(1)
    point = rows[row_index]
    with open(args.anchors, newline="") as handle:
        anchors = list(csv.DictReader(handle, delimiter="\t"))
    for anchor in anchors:
        if int(point["yuktype"]) != int(anchor["yuktype"]):
            continue
        radius = float(anchor["rel_radius"])
        matched = True
        for key in PARAMETERS:
            actual, reference = float(point[key]), float(anchor[key])
            tolerance = radius * max(abs(reference), 1e-12)
            if not math.isclose(actual, reference, rel_tol=0.0,
                                abs_tol=tolerance):
                matched = False
                break
        if matched:
            print(f"known_exact_region:{anchor['reason']}")
            return
    raise SystemExit(1)


if __name__ == "__main__":
    main()
