#!/usr/bin/env python3
"""Exit 0 only for outcomes accepted by the frozen region-v1 router."""

import argparse
import re

from compare_qualitative_outcomes import load, outcome

parser = argparse.ArgumentParser()
parser.add_argument("output")
parser.add_argument("diagnostics")
args = parser.parse_args()

rows = load(args.output)
if len(rows) != 1:
    print("unsupported_row_count")
    raise SystemExit(1)

classification, stage, _ = outcome(rows[0])
if classification == "fail_or_nonpositive":
    if stage.startswith("nlo:") or stage.startswith("coexistence:"):
        print(f"accepted_early_failure:{stage}")
        raise SystemExit(0)
    print(f"fallback_late_failure:{stage}")
    raise SystemExit(1)

shadow_lines = [
    line for line in open(args.diagnostics, errors="replace")
    if line.startswith("BSMPT_BOUNCE_SHADOW")
]
if not shadow_lines:
    print("fallback_missing_shadow")
    raise SystemExit(1)
if all(re.search(r"(?:^|\t)success=1(?:\t|$)", line) for line in shadow_lines):
    print("accepted_positive_shadow")
    raise SystemExit(0)
print("fallback_shadow_disagreement")
raise SystemExit(1)
