#!/usr/bin/env python3
"""Run independent CalcGW points concurrently without changing BSMPT numerics."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
from pathlib import Path
import subprocess
import tempfile


def run_point(task):
    index, header, row, binary, extra, temp_root = task
    point_dir = Path(temp_root) / f"point_{index:08d}"
    point_dir.mkdir()
    input_path = point_dir / "input.tsv"
    output_path = point_dir / "output.tsv"
    input_path.write_text(header + "\n" + row + "\n")

    command = [
        binary,
        "--model=r2hdm",
        f"--input={input_path}",
        f"--output={output_path}",
        "--firstline=2",
        "--lastline=2",
        *extra,
    ]
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        env=environment,
        check=False,
    )
    if completed.returncode != 0 or not output_path.exists():
        message = completed.stderr.strip() or "CalcGW did not create output"
        raise RuntimeError(f"point {index}: {message}")

    lines = output_path.read_text().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"point {index}: incomplete CalcGW output")
    return index, lines[0], lines[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to bin/CalcGW")
    parser.add_argument("--input", required=True, help="BSMPT input TSV")
    parser.add_argument("--output", required=True, help="Merged output TSV")
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument(
        "--extra",
        nargs="*",
        default=["--thigh=400", "--checknlo=on", "--checkewsr=on", "--multistepmode=default"],
        help="Additional CalcGW arguments",
    )
    args = parser.parse_args()

    source_lines = Path(args.input).read_text().splitlines()
    if len(source_lines) < 2:
        parser.error("input must contain a header and at least one point")
    header, rows = source_lines[0], [line for line in source_lines[1:] if line.strip()]
    jobs = min(max(1, args.jobs), len(rows))
    binary = str(Path(args.binary).resolve())

    with tempfile.TemporaryDirectory(prefix="bsmpt-parallel-") as temp_root:
        tasks = [
            (i, header, row, binary, args.extra, temp_root)
            for i, row in enumerate(rows)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(run_point, tasks))

    results.sort(key=lambda item: item[0])
    output_headers = {item[1] for item in results}
    if len(output_headers) != 1:
        raise RuntimeError("CalcGW workers produced inconsistent output headers")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(results[0][1] + "\n" + "\n".join(item[2] for item in results) + "\n")


if __name__ == "__main__":
    main()
