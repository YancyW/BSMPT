#!/usr/bin/env python3
"""Return 0 when an approximate CalcGW row must be recomputed exactly."""
import argparse
import csv
import math


RISK_STATES = {"failure", "nan", "not_met", "not_set"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--snr-floor", type=float, default=1e-20)
    parser.add_argument("--snr-cut", action="append", type=float, default=[])
    parser.add_argument("--cut-margin", type=float, default=0.15)
    args = parser.parse_args()

    with open(args.output, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    reasons = []
    if len(rows) != 1:
        reasons.append(f"expected_one_row:{len(rows)}")
    for row in rows:
        for key, value in row.items():
            if key.startswith("status_") and value in RISK_STATES:
                reasons.append(f"risk_state:{key}={value}")
        totals = [
            (key, value)
            for key, value in row.items()
            if key.startswith("SNR(") and key.rsplit("_", 1)[-1].isdigit()
            and not any(part in key for part in ("_col_", "_sw_", "_turb_"))
        ]
        if not totals:
            reasons.append("missing_total_snr")
        for key, value in totals:
            try:
                snr = float(value)
            except ValueError:
                snr = math.nan
            if not math.isfinite(snr):
                reasons.append(f"nonfinite_snr:{key}")
                continue
            if abs(snr) < args.snr_floor:
                reasons.append(f"low_snr:{key}={snr:.9g}")
            for cut in args.snr_cut:
                scale = max(abs(cut), args.snr_floor)
                if abs(snr - cut) <= args.cut_margin * scale:
                    reasons.append(f"near_cut:{key}={snr:.9g},cut={cut:.9g}")

    if reasons:
        print(";".join(dict.fromkeys(reasons)))
        return
    raise SystemExit(1)


if __name__ == "__main__":
    main()
