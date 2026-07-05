#!/usr/bin/env python
"""Run crispat (Poisson-Gaussian, the lab's representative crispat method) on sim
datasets.  Reads each dataset's ext_counts.mtx; writes crispat_calls.csv
(guide,cell).  Diagnostic plotting is patched out (Agg + no-op savefig) -- the
assignment logic is untouched, only the per-guide PNGs are skipped for speed.
Shardable for parallelism:  run_crispat_sims.py <shard_idx> <n_shards>
Run with .venv_crispat/bin/python from the folder root."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt; plt.savefig = lambda *a, **k: None
import os, glob, sys, shutil, numpy as np, scipy.io, anndata as ad, pandas as pd
import crispat

DS = "results/sim_framework/datasets"

def run_one(d):
    m = scipy.io.mmread(os.path.join(d, "ext_counts.mtx")).tocsr()      # guides x cells
    guides = [l.strip() for l in open(os.path.join(d, "ext_guides.txt"))]
    cells  = [l.strip() for l in open(os.path.join(d, "ext_cells.txt"))]
    A = ad.AnnData(X=m.T.tocsr().astype("float32"))                    # cells x guides
    A.var_names = guides; A.obs_names = cells
    tmp_h5 = os.path.join(d, "_crispat_in.h5ad")
    outdir = os.path.join(d, "_crispat_out")
    try:
        A.write(tmp_h5)
        shutil.rmtree(outdir, ignore_errors=True)                      # crispat creates plot subdirs; start clean
        crispat.ga_poisson_gauss(tmp_h5, outdir + "/")
        asg_path = os.path.join(outdir, "assignments.csv")
        try:
            asg = pd.read_csv(asg_path)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            asg = pd.DataFrame(columns=["cell", "gRNA", "UMI_counts"]) # legitimate empty (weak-signal regimes)
        out_calls = os.path.join(d, "crispat_calls.csv")
        if len(asg):
            asg[["gRNA", "cell"]].rename(columns={"gRNA": "guide"}).to_csv(out_calls, index=False)
        else:
            pd.DataFrame(columns=["guide", "cell"]).to_csv(out_calls, index=False)
        return len(asg)
    finally:
        if os.path.exists(tmp_h5): os.remove(tmp_h5)
        shutil.rmtree(outdir, ignore_errors=True)

if __name__ == "__main__":
    shard = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    n_shards = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    dirs = sorted(d for d in glob.glob(os.path.join(DS, "*")) if os.path.exists(os.path.join(d, "ext_counts.mtx")))
    dirs = [d for d in dirs if not os.path.exists(os.path.join(d, "crispat_calls.csv"))]  # idempotent: skip done
    dirs = [d for i, d in enumerate(dirs) if i % n_shards == shard]
    for i, d in enumerate(dirs):
        try:
            n = run_one(d)
            print(f"[shard{shard} {i+1}/{len(dirs)}] {os.path.basename(d)}: {n} crispat calls", flush=True)
        except Exception as e:
            print(f"[shard{shard} {i+1}/{len(dirs)}] {os.path.basename(d)}: ERROR {e}", flush=True)
