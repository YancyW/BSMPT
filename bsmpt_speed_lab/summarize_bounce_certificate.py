#!/usr/bin/env python3
"""Summarize opt-in BSMPT_BOUNCE_CERT sidecar lines by temperature."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_fields(line: str) -> tuple[str, dict[str, str]] | None:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 2 or parts[0] != "BSMPT_BOUNCE_CERT":
        return None
    return parts[1], dict(part.split("=", 1) for part in parts[2:] if "=" in part)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    by_temperature: dict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    with args.log.open(errors="replace") as stream:
        for line in stream:
            parsed = parse_fields(line)
            if parsed is not None:
                kind, fields = parsed
                by_temperature[fields.get("T", "nan")].append((kind, fields))

    fieldnames = [
        "T",
        "path_checks",
        "deformations",
        "converged_deformations",
        "stalled_deformations",
        "min_best_relative_error",
        "max_best_relative_error",
        "certificate_stable",
    ]
    rows = []
    for temperature, records in sorted(
        by_temperature.items(), key=lambda item: float(item[0])
    ):
        deformations = [fields for kind, fields in records if kind == "path_deformation"]
        errors = [float(fields["best_relative_error"]) for fields in deformations]
        converged = sum(fields.get("converged_without_1d") == "1" for fields in deformations)
        stalled = sum(
            fields.get("converged_without_1d") != "1"
            and int(fields.get("iterations", "0")) >= 21
            for fields in deformations
        )
        # This is a deliberately conservative candidate, not a validated router:
        # every observed deformation at this temperature must reach the solver's
        # own 0.05 convergence threshold.
        stable = bool(deformations) and converged == len(deformations)
        rows.append(
            {
                "T": temperature,
                "path_checks": sum(kind == "path_check" for kind, _ in records),
                "deformations": len(deformations),
                "converged_deformations": converged,
                "stalled_deformations": stalled,
                "min_best_relative_error": min(errors) if errors else "",
                "max_best_relative_error": max(errors) if errors else "",
                "certificate_stable": str(stable).lower(),
            }
        )

    output = args.output.open("w", newline="") if args.output else None
    stream = output or __import__("sys").stdout
    try:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if output:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
