#!/usr/bin/env python3
"""Generate the expanded comparison table by column name."""
import csv, math
from pathlib import Path

root = Path(__file__).resolve().parent
pairs = (
    ("nlo_boundary", "expanded_v314_nlo_boundary_7.tsv", "expanded_current_nlo_boundary_7.tsv"),
    ("yukawa_families_2_3", "expanded_v314_yukawa_families_2_3.tsv", "expanded_current_yukawa_families_2_3.tsv"),
)
statuses = ("status_nlo_stability", "status_ewsr", "status_tracing", "status_coex_pairs",
            "status_crit_0", "status_bounce_sol_0", "status_nucl_approx_0", "status_nucl_0",
            "status_perc_0", "status_compl_0", "status_gw_0")

def rel(old, new):
    try:
        old, new = float(old), float(new)
        return f"{100*(new-old)/abs(old):.12g}" if old and math.isfinite(old) and math.isfinite(new) else "nan"
    except ValueError:
        return "nan"

cols = ("group row yuktype L1 m12sq status_flip history_flip status_gw_v314 "
        "status_gw_current alpha_rel_pct betaH_rel_pct snr_rel_pct snr_v314 "
        "snr_current runtime_v314_s runtime_current_s").split()
with (root / "expanded_result_summary.tsv").open("w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=cols)
    writer.writeheader()
    for group, old_file, new_file in pairs:
        old = list(csv.DictReader((root / old_file).open(), delimiter="\t"))
        new = list(csv.DictReader((root / new_file).open(), delimiter="\t"))
        if len(old) != len(new): raise RuntimeError(f"row mismatch: {group}")
        for row, (a, b) in enumerate(zip(old, new), 1):
            writer.writerow(dict(group=group, row=row, yuktype=a["yuktype"], L1=a["L1"],
                m12sq=a["m12sq"], status_flip=int(any(a.get(k) != b.get(k) for k in statuses)),
                history_flip=int(a["transition_history"] != b["transition_history"]),
                status_gw_v314=a["status_gw_0"], status_gw_current=b["status_gw_0"],
                alpha_rel_pct=rel(a["alpha_PT_0"], b["alpha_PT_0"]),
                betaH_rel_pct=rel(a["beta/H_0"], b["beta/H_0"]),
                snr_rel_pct=rel(a["SNR(LISA-3yrs)_0"], b["SNR(LISA-3yrs)_0"]),
                snr_v314=a["SNR(LISA-3yrs)_0"], snr_current=b["SNR(LISA-3yrs)_0"],
                runtime_v314_s=a["runtime"], runtime_current_s=b["runtime"]))
