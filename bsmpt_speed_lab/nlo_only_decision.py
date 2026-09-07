#!/usr/bin/env python3
"""Exit zero when an exact NLO-only CalcGW result rejects the point."""

import csv
import sys


with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
if len(rows) != 1:
    raise SystemExit(2)
status = (rows[0].get("status_nlo_stability") or "").strip()
runtime = (rows[0].get("runtime") or "nan").strip()
if status == "success":
    print(f"exact_nlo_success:runtime={runtime}")
    raise SystemExit(1)
if status in {"", "nan", "not_set"}:
    raise SystemExit(2)
print(f"exact_nlo_reject:{status}:runtime={runtime}")
