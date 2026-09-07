#!/usr/bin/env python3
"""Run independent CalcGW points concurrently without changing BSMPT numerics."""

from __future__ import annotations

import argparse
import concurrent.futures
import fcntl
import os
from pathlib import Path
import subprocess
import tempfile


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


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
    diagnostics = [
        line
        for line in completed.stderr.splitlines()
        if line.startswith("BSMPT_BOUNCE_SHADOW")
        or line.startswith("BSMPT_BOUNCE_CERT")
        or line.startswith("region-v1:")
        or line.startswith("region-v2:")
    ]
    return index, lines[0], lines[1], diagnostics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, help="Path to bin/CalcGW")
    parser.add_argument("--input", required=True, help="BSMPT input TSV")
    parser.add_argument("--output", required=True, help="Merged output TSV")
    parser.add_argument(
        "--diagnostic-output",
        help="Optional line-numbered bounce certificate and region decision sidecar",
    )
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument(
        "--extra",
        nargs="*",
        default=["--thigh=400", "--checknlo=on", "--checkewsr=on", "--multistepmode=default"],
        help="Additional CalcGW arguments",
    )
    args = parser.parse_args()

    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output.with_name(output.name + ".lock")
    lock_stream = lock_path.open("w")
    try:
        fcntl.flock(lock_stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error(f"another runner is already writing {output}")
    lock_stream.write(f"pid={os.getpid()}\n")
    lock_stream.flush()

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
    header_fields = [item[1].split("\t") for item in results]
    widest_header = max(header_fields, key=len)
    if any(not set(fields).issubset(widest_header) for fields in header_fields):
        raise RuntimeError("CalcGW workers produced incompatible output headers")
    normalized_rows = []
    for _, header, row, _ in results:
        values = dict(zip(header.split("\t"), row.split("\t")))
        normalized_rows.append("\t".join(values.get(name, "nan") for name in widest_header))
    atomic_write(output, "\t".join(widest_header) + "\n" + "\n".join(normalized_rows) + "\n")
    if args.diagnostic_output:
        diagnostic = Path(args.diagnostic_output)
        diagnostic.parent.mkdir(parents=True, exist_ok=True)
        atomic_write(
            diagnostic,
            "input_row\tdiagnostic\n"
            + "".join(
                f"{index + 1}\t{line}\n"
                for index, _, _, lines in results
                for line in lines
            ),
        )


if __name__ == "__main__":
    main()
