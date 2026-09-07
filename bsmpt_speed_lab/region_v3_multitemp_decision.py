#!/usr/bin/env python3
"""Experimental positive certificate using three legacy-shadow temperatures."""

import argparse
import os
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
    if stage.startswith("nlo:"):
        if os.environ.get("BSMPT_REQUIRE_EXACT_NLO_PREFIX") == "1":
            print(f"fallback_approx_nlo_disagrees_with_exact_prefix:{stage}")
            raise SystemExit(1)
        print(f"accepted_early_failure:{stage}")
        raise SystemExit(0)
    if stage.startswith("coexistence:") and os.environ.get("BSMPT_ALLOW_COEX", "1") == "1":
        print(f"accepted_early_failure:{stage}")
        raise SystemExit(0)
    if stage.startswith("coexistence:"):
        print(f"fallback_coexistence_outside_certified_region:{stage}")
        raise SystemExit(1)
    print(f"fallback_late_failure:{stage}")
    raise SystemExit(1)

rank_success = {}
for line in open(args.diagnostics, errors="replace"):
    if not line.startswith("BSMPT_BOUNCE_SHADOW"):
        continue
    rank = re.search(r"(?:^|\t)rank=(\d+)(?:\t|$)", line)
    success = re.search(r"(?:^|\t)success=([01])(?:\t|$)", line)
    if rank and success:
        rank_success[int(rank.group(1))] = success.group(1) == "1"

if set(rank_success) != {0, 1, 2}:
    print("fallback_incomplete_three_temperature_shadow")
    raise SystemExit(1)
if all(rank_success.values()) and os.environ.get("BSMPT_ALLOW_ALL_THREE", "1") == "1":
    print("accepted_positive_all_three_shadow")
    raise SystemExit(0)
if all(rank_success.values()):
    print("fallback_all_three_outside_certified_region")
    raise SystemExit(1)
if (os.environ.get("BSMPT_ALLOW_TWO_NEIGHBOR") == "1"
        and not rank_success[0] and rank_success[1] and rank_success[2]):
    print("accepted_positive_two_neighbor_shadow")
    raise SystemExit(0)
if not rank_success[0] and rank_success[1] and rank_success[2]:
    print("fallback_two_neighbor_pattern_outside_certified_region")
    raise SystemExit(1)
print("fallback_three_temperature_shadow_pattern")
raise SystemExit(1)
