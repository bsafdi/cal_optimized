"""Production (alpha, beta) sweep.  Writes production_sweep.csv."""

from __future__ import annotations
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from bem_solver import run_one, run_one_har


ALPHAS = [0.5, 1.0, 2.0]
BETAS = [0.20, 0.12, 0.08, 0.05, 0.035, 0.025]


def main(out_csv: str, n_t: int = 15, w: float = 0.03,
         wire_type: str = "production"):
    rows = []
    t0_all = time.time()
    print(f"Sweep: wire_type={wire_type}, n_t={n_t}, w={w}", flush=True)
    for alpha in ALPHAS:
        for beta in BETAS:
            t0 = time.time()
            if wire_type == "harmonic":
                res, _ = run_one_har(alpha=alpha, beta=beta, w=w,
                                     n_t_per_seg=n_t)
            else:
                res, _ = run_one(alpha=alpha, beta=beta, w=w, n_t_per_seg=n_t,
                                 wire_type=wire_type)
            rows.append({
                "alpha": alpha,
                "beta": beta,
                "F": res.F,
                "ell": res.ell,
                "w_nominal": res.w_nominal,
                "w_eff": res.w_eff,
                "n_elements": res.n_elements,
                "cond": res.matrix_cond,
                "wall_s": res.wall_time_s,
            })
            print(f"  alpha={alpha:.2f} beta={beta:.4f}  F={res.F:+.4e}  ell={res.ell:.4e}"
                  f"  cond={res.matrix_cond:.2e}  t={time.time()-t0:.1f}s",
                  flush=True)
    total = time.time() - t0_all
    print(f"\nSweep wall time: {total:.1f}s", flush=True)

    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"Wrote {out_csv}.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/production_sweep.csv")
    ap.add_argument("--nt", type=int, default=15)
    ap.add_argument("--w", type=float, default=0.03)
    ap.add_argument("--wire", default="production",
                     choices=["production", "linked", "none", "harmonic"])
    args = ap.parse_args()
    main(args.out, n_t=args.nt, w=args.w, wire_type=args.wire)
