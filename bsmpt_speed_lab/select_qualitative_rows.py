#!/usr/bin/env python3
"""Select input rows whose paired CalcGW output is qualitatively positive."""

import argparse
import csv
import math

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output_tsv")
parser.add_argument("selected_input")
args = parser.parse_args()

with open(args.input, newline="") as stream:
    inputs = list(csv.reader(stream, delimiter="\t"))
with open(args.output_tsv, newline="") as stream:
    outputs = list(csv.DictReader(stream, delimiter="\t"))

selected = []
for source, result in zip(inputs[1:], outputs):
    positive = False
    for key, value in result.items():
        if not key.startswith("status_gw_") or value != "success":
            continue
        suffix = key.removeprefix("status_gw_")
        try:
            snr = float(result.get(f"SNR(LISA-3yrs)_{suffix}", "nan"))
        except ValueError:
            continue
        positive |= math.isfinite(snr) and snr > 0
    if positive:
        selected.append(source)

with open(args.selected_input, "w", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(inputs[0])
    writer.writerows(selected)
