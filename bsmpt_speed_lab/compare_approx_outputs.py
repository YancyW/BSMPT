#!/usr/bin/env python3
import argparse
import csv
import math


def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite_number(value):
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def relative_error(reference, candidate):
    scale = max(abs(reference), 1e-300)
    return abs(candidate - reference) / scale


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--snr-limit", type=float, default=0.10)
    args = parser.parse_args()
    ref_rows, cand_rows = load(args.reference), load(args.candidate)
    if len(ref_rows) != len(cand_rows):
        raise SystemExit(f"FAIL: row count {len(ref_rows)} != {len(cand_rows)}")

    status_changes = []
    max_snr = (0.0, "", 0)
    max_numeric = (0.0, "", 0)
    for row_index, (ref, cand) in enumerate(zip(ref_rows, cand_rows), start=1):
        if ref.keys() != cand.keys():
            raise SystemExit(f"FAIL: columns differ at row {row_index}")
        for key in ref:
            if key == "runtime":
                continue
            rnum, cnum = finite_number(ref[key]), finite_number(cand[key])
            categorical = key.startswith("status_") or key == "transition_history"
            if categorical or rnum is None or cnum is None:
                if ref[key] != cand[key]:
                    status_changes.append((row_index, key, ref[key], cand[key]))
                continue
            error = relative_error(rnum, cnum)
            if error > max_numeric[0]:
                max_numeric = (error, key, row_index)
            if key.startswith("SNR(") and error > max_snr[0]:
                max_snr = (error, key, row_index)

    if status_changes:
        print(f"FAIL: {len(status_changes)} categorical/non-finite changes")
        for item in status_changes[:20]:
            print("  row=%d field=%s reference=%r candidate=%r" % item)
        raise SystemExit(1)
    print(f"PASS: {len(ref_rows)} rows preserve every status/history field")
    print(f"max_snr_relative_error={max_snr[0]:.9g} field={max_snr[1]} row={max_snr[2]}")
    print(f"max_numeric_relative_error={max_numeric[0]:.9g} field={max_numeric[1]} row={max_numeric[2]}")
    if max_snr[0] > args.snr_limit:
        raise SystemExit(f"FAIL: SNR error exceeds {args.snr_limit:.1%}")


if __name__ == "__main__":
    main()
