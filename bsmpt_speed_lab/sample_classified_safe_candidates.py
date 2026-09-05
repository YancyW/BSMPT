#!/usr/bin/env python3
"""Stream classified CSVs and build a deterministic approximation-validation set.

The source trees are read only.  Crucially, thdmTools-rejected rows are not
treated as CalcGW failures: the candidate population is conditioned on
``th_passed == true`` and on a complete CalcGW input record.
"""

import argparse
import csv
import hashlib
import heapq
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


INPUT_FIELDS = ("yuktype", "L1", "L2", "L3", "L4", "L5", "m12sq_calgw", "tbeta")
KEY_FIELDS = ("mH", "mA", "mHc", "m12sq", "sba", "tb")
DEFAULT_ROOTS = (
    "/home/yancy/Software/Bin/MathRelated/MultiNest/data_office",
    "/home/yancy/Software/Bin/MathRelated/MultiNest/data_qiushi",
    "/home/yancy/Software/Bin/MathRelated/MultiNest/data_yuanyuan",
)


def finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def truth(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def traceable_source(value):
    value = (value or "").strip()
    return bool(value) and value.endswith(".md")


def branch(row):
    mh, sba = finite(row.get("mH")), finite(row.get("sba"))
    if mh is None or sba is None:
        return "outside"
    if 62 <= mh < 125 and 0.06 <= sba <= 0.21:
        return "light"
    if 125 <= mh <= 500 and 0.99 <= sba <= 0.999:
        return "heavy"
    return "outside"


def complete_input(row):
    return all(finite(row.get(k)) is not None for k in INPUT_FIELDS)


def total_snr(row):
    for key in ("SNR_total", "SNR", "SNR(LISA-3yrs)_0"):
        value = finite(row.get(key))
        if value is not None:
            return value
    return None


def active_transition_is_risky(row):
    """Check routing-relevant states without misclassifying ew_sym_res.

    Inactive transition blocks are stored as nan/not_set.  EWSR also uses a
    physical label rather than the word success, so a blanket status scan would
    incorrectly mark essentially every valid row as a numerical failure.
    """
    if (row.get("status_tracing") or "") != "success":
        return True
    active = [i for i in range(3) if (row.get(f"status_crit_{i}") or "") == "success"]
    if not active:
        return True
    return any((row.get(f"status_gw_{i}") or "") not in {"success", "nan", "not_set"} for i in active)


def edge_score(row, br):
    mh, ma, mhc, sba = (finite(row.get(k)) for k in ("mH", "mA", "mHc", "sba"))
    if None in (mh, ma, mhc, sba):
        return 1.0
    if br == "light":
        distances = ((mh - 62) / 63, (125 - mh) / 63, (sba - 0.06) / 0.15, (0.21 - sba) / 0.15)
    elif br == "heavy":
        distances = ((mh - 125) / 375, (500 - mh) / 375, (sba - 0.99) / 0.009, (0.999 - sba) / 0.009)
    else:
        return 1.0
    branch_edge = max(0.0, 1.0 - 8.0 * min(distances))
    splitting_edge = min(1.0, abs(ma - mhc) / 60.0)
    return max(branch_edge, splitting_edge)


def classify(row, category, br):
    multi = (row.get("multistep_status") or "").strip()
    snr = total_snr(row)
    if category == "thdmtool_only":
        return "strict_negative"
    if multi not in {"single_step", ""} or active_transition_is_risky(row):
        return "numeric_boundary"
    if snr is None or snr <= 1.0 or edge_score(row, br) >= 0.65:
        return "physical_boundary"
    return "interior"


def stable_rank(machine, source, row):
    values = [machine, source] + [row.get(k, "") for k in KEY_FIELDS]
    return int.from_bytes(hashlib.sha256("|".join(values).encode()).digest()[:8], "big")


def add_reservoir(reservoirs, group, cap, rank, payload):
    heap = reservoirs[group]
    item = (-rank, payload)
    if len(heap) < cap:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", action="append", default=[])
    parser.add_argument("--output-prefix", default="bsmpt_speed_lab/classified_safe_candidates_240")
    parser.add_argument("--per-stratum", type=int, default=12)
    parser.add_argument("--limit", type=int, default=240)
    args = parser.parse_args()

    roots = [Path(p) for p in (args.root or DEFAULT_ROOTS)]
    reservoirs = defaultdict(list)
    counts, exclusions = Counter(), Counter()
    seen_keys = set()

    for root in roots:
        machine = root.name.removeprefix("data_")
        for category in ("both", "thdmtool_only"):
            paths = sorted((root / "classified").glob(f"{category}*.csv"))
            for path in paths:
                with path.open(newline="") as handle:
                    for line_number, row in enumerate(csv.DictReader(handle), 2):
                        counts[f"read:{machine}:{category}"] += 1
                        if not truth(row.get("th_passed")):
                            exclusions["thdm_not_passed"] += 1
                            continue
                        if not complete_input(row):
                            exclusions["incomplete_calcgw_input"] += 1
                            continue
                        source = (row.get("source") or "").strip()
                        if not traceable_source(source):
                            exclusions["untraceable_source"] += 1
                            continue
                        br = branch(row)
                        if br == "outside":
                            exclusions["outside_known_thdm_branch"] += 1
                            continue
                        key = tuple(round(float(row[k]), 6) for k in KEY_FIELDS)
                        if key in seen_keys:
                            exclusions["duplicate_round6"] += 1
                            continue
                        seen_keys.add(key)
                        mode = classify(row, category, br)
                        yuk = str(int(float(row["yuktype"])))
                        group = (mode, br, yuk, machine)
                        rank = stable_rank(machine, source, row)
                        payload = {
                            "machine": machine,
                            "category": category,
                            "file": str(path),
                            "line": line_number,
                            "source": source,
                            "branch": br,
                            "validation_stratum": mode,
                            "edge_score": edge_score(row, br),
                            "strict_snr": total_snr(row),
                            "row": row,
                        }
                        add_reservoir(reservoirs, group, args.per_stratum, rank, payload)
                        counts[f"eligible:{mode}:{br}:{yuk}:{machine}"] += 1

    candidates = []
    for group in sorted(reservoirs):
        candidates.extend(payload for _, payload in sorted(reservoirs[group], reverse=True))
    candidates.sort(key=lambda p: (p["validation_stratum"], p["branch"], p["machine"], p["source"], p["line"]))
    if len(candidates) > args.limit:
        # A second stable ranking keeps the final cap unbiased across populated strata.
        candidates = sorted(candidates, key=lambda p: stable_rank(p["machine"], p["source"], p["row"]))[: args.limit]
        candidates.sort(key=lambda p: (p["validation_stratum"], p["branch"], p["machine"], p["line"]))

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    input_path = prefix.with_name(prefix.name + "_input.tsv")
    metadata_path = prefix.with_name(prefix.name + "_metadata.tsv")
    reference_path = prefix.with_name(prefix.name + "_strict_reference.tsv")
    summary_path = prefix.with_name(prefix.name + "_summary.json")

    with input_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("yuktype", "L1", "L2", "L3", "L4", "L5", "m12sq", "tbeta"))
        for p in candidates:
            r = p["row"]
            writer.writerow((r["yuktype"], r["L1"], r["L2"], r["L3"], r["L4"], r["L5"], r["m12sq_calgw"], r["tbeta"]))

    meta_fields = ("sample_id", "machine", "category", "branch", "validation_stratum", "edge_score", "strict_snr", "source", "file", "line") + KEY_FIELDS
    with metadata_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=meta_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, p in enumerate(candidates, 1):
            record = {k: p[k] for k in meta_fields if k in p}
            record.update({k: p["row"].get(k, "") for k in KEY_FIELDS})
            record["sample_id"] = i
            writer.writerow(record)

    reference_fields = list(candidates[0]["row"]) if candidates else []
    with reference_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reference_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(p["row"] for p in candidates)

    selected = Counter((p["validation_stratum"], p["branch"], p["machine"], str(int(float(p["row"]["yuktype"])))) for p in candidates)
    summary = {
        "selection_semantics": "conditioned_on_thdmTools_pass; neither/gw_only excluded",
        "selected_rows": len(candidates),
        "selected_by_stratum_branch_machine_yuk": {"|".join(k): v for k, v in sorted(selected.items())},
        "stream_counts": dict(sorted(counts.items())),
        "exclusions": dict(sorted(exclusions.items())),
        "roots": [str(p) for p in roots],
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"selected": len(candidates), "input": str(input_path), "metadata": str(metadata_path), "reference": str(reference_path), "summary": str(summary_path)}))


if __name__ == "__main__":
    main()
