#!/usr/bin/env python
"""Run geomux v5 (adaptive LOR, the published method) on every sim dataset.
Reads each dataset dir's ext_counts.mtx (guides x cells) + ext_guides.txt +
ext_cells.txt, writes geomux_calls.csv (guide,cell) of assigned pairs.
Run with .venv_geomux_v5/bin/python from the folder root."""
import os, glob, sys, numpy as np, scipy.io, geomux

DS = "results/sim_framework/datasets"

def run_one(d):
    m = scipy.io.mmread(os.path.join(d, "ext_counts.mtx")).tocsr()      # guides x cells
    guides = np.array([l.strip() for l in open(os.path.join(d, "ext_guides.txt"))])
    cells  = np.array([l.strip() for l in open(os.path.join(d, "ext_cells.txt"))])
    X = m.T.tocsr()                                                     # cells x guides
    res = geomux.geomux(X, cell_names=cells, guide_names=guides, fdr_threshold=0.05)  # adaptive LOR (default)
    df = res.to_pandas()
    rows = []
    for _, r in df.iterrows():
        gio = r.get("guide_ids_original", "")
        if isinstance(gio, str) and gio:
            for gi in gio.split("|"):
                rows.append((guides[int(gi)], r["cell"]))
    import pandas as pd
    pd.DataFrame(rows, columns=["guide", "cell"]).to_csv(os.path.join(d, "geomux_calls.csv"), index=False)
    return len(rows)

if __name__ == "__main__":
    dirs = sorted(d for d in glob.glob(os.path.join(DS, "*")) if os.path.exists(os.path.join(d, "ext_counts.mtx")))
    for i, d in enumerate(dirs):
        try:
            n = run_one(d)
            print(f"[{i+1}/{len(dirs)}] {os.path.basename(d)}: {n} geomux calls", flush=True)
        except Exception as e:
            print(f"[{i+1}/{len(dirs)}] {os.path.basename(d)}: ERROR {e}", flush=True)
