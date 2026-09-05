#!/usr/bin/env python3
"""Stream classified CSVs and build a deterministic approximation-validation set.

The source trees are read only.  Selection uses evidence that CalcGW actually
ran, never thdmTools outcomes.  The resulting archive sample is biased and is
only E1 auxiliary evidence, not a representation of the full BSMPT space.
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


def traceable_source(value):
    value = (value or "").strip()
    return bool(value) and value.endswith(".md")


def parameter_cell(row):
    mh, sba = finite(row.get("mH")), finite(row.get("sba"))
    if mh is None or sba is None:
        return "unknown"
    mass_edges = (-math.inf, 62, 125, 220, 500, math.inf)
    sba_edges = (-math.inf, 0.06, 0.21, 0.9, 0.99, 0.999, math.inf)
    mi = next(i for i in range(len(mass_edges) - 1) if mass_edges[i] <= mh < mass_edges[i + 1])
    si = next(i for i in range(len(sba_edges) - 1) if sba_edges[i] <= sba < sba_edges[i + 1])
    return f"m{mi}_s{si}"


def complete_input(row):
    return all(finite(row.get(k)) is not None for k in INPUT_FIELDS)


def calcgw_actually_ran(row):
    """Use only BSMPT execution evidence; never inspect th_* fields."""
    runtime = finite(row.get("runtime"))
    nlo = (row.get("status_nlo_stability") or "").strip()
    tracing = (row.get("status_tracing") or "").strip()
    return runtime is not None and runtime >= 0 and bool(nlo or tracing)


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


def classify(row):
    multi = (row.get("multistep_status") or "").strip()
    snr = total_snr(row)
    if multi not in {"single_step", ""} or active_transition_is_risky(row):
        return "numeric_boundary"
    if snr is None or abs(snr) < 1e-20:
        return "strict_no_signal"
    if snr <= 1.0 or any(abs(snr - cut) <= 0.2 * cut for cut in (10.0, 100.0)):
        return "snr_boundary"
    return "resolved_signal"


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
    parser.add_argument("--category", action="append", choices=("both", "gw_only", "thdmtool_only", "neither"), default=[])
    parser.add_argument("--output-prefix", default="bsmpt_speed_lab/classified_safe_candidates_240")
    parser.add_argument("--per-stratum", type=int, default=12)
    parser.add_argument("--limit", type=int, default=240)
    args = parser.parse_args()

    roots = [Path(p) for p in (args.root or DEFAULT_ROOTS)]
    categories = args.category or ["both", "gw_only", "thdmtool_only"]
    reservoirs = defaultdict(list)
    counts, exclusions = Counter(), Counter()
    seen_keys = set()

    for root in roots:
        machine = root.name.removeprefix("data_")
        for category in categories:
            paths = sorted((root / "classified").glob(f"{category}*.csv"))
            for path in paths:
                with path.open(newline="") as handle:
                    for line_number, row in enumerate(csv.DictReader(handle), 2):
                        counts[f"read:{machine}:{category}"] += 1
                        if not complete_input(row):
                            exclusions["incomplete_calcgw_input"] += 1
                            continue
                        if not calcgw_actually_ran(row):
                            exclusions["calcgw_not_run_or_unverifiable"] += 1
                            continue
                        source = (row.get("source") or "").strip()
                        if not traceable_source(source):
                            exclusions["untraceable_source"] += 1
                            continue
                        cell = parameter_cell(row)
                        key = tuple(round(float(row[k]), 6) for k in KEY_FIELDS)
                        if key in seen_keys:
                            exclusions["duplicate_round6"] += 1
                            continue
                        seen_keys.add(key)
                        mode = classify(row)
                        yuk = str(int(float(row["yuktype"])))
                        group = (mode, cell, yuk, machine)
                        rank = stable_rank(machine, source, row)
                        payload = {
                            "machine": machine,
                            "category": category,
                            "file": str(path),
                            "line": line_number,
                            "source": source,
                            "parameter_cell": cell,
                            "validation_stratum": mode,
                            "strict_snr": total_snr(row),
                            "row": row,
                        }
                        add_reservoir(reservoirs, group, args.per_stratum, rank, payload)
                        counts[f"eligible:{mode}:{cell}:{yuk}:{machine}"] += 1

    queues = {group: [payload for _, payload in sorted(heap, reverse=True)] for group, heap in reservoirs.items()}
    candidates = []
    while len(candidates) < args.limit and any(queues.values()):
        for group in sorted(queues):
            if queues[group] and len(candidates) < args.limit:
                candidates.append(queues[group].pop(0))
    candidates.sort(key=lambda p: (p["validation_stratum"], p["parameter_cell"], p["machine"], p["source"], p["line"]))

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

    meta_fields = ("sample_id", "machine", "archive_category", "parameter_cell", "validation_stratum", "strict_snr", "source", "file", "line") + KEY_FIELDS
    with metadata_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=meta_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, p in enumerate(candidates, 1):
            record = {k: p[k] for k in meta_fields if k in p}
            record["archive_category"] = p["category"]
            record.update({k: p["row"].get(k, "") for k in KEY_FIELDS})
            record["sample_id"] = i
            writer.writerow(record)

    reference_fields = list(candidates[0]["row"]) if candidates else []
    with reference_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reference_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(p["row"] for p in candidates)

    selected = Counter((p["validation_stratum"], p["parameter_cell"], p["machine"], str(int(float(p["row"]["yuktype"])))) for p in candidates)
    summary = {
        "evidence_level": "E1_biased_archive_auxiliary_only",
        "selection_semantics": "CalcGW execution evidence and complete BSMPT input; th_* fields ignored",
        "categories_scanned_for_file_location_only": categories,
        "coverage_warning": "biased archive; must not define full-space safety or production error rates",
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
