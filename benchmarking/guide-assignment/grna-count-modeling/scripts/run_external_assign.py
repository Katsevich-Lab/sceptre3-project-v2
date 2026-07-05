#!/usr/bin/env python
# External assignment runner: crispat (Poisson-Gaussian) or cleanser (MCMC).
# Reads a guides x cells matrix (mtx + guides.txt + cells.txt), runs the method,
# writes calls.csv with columns guide,cell (the assigned (g,c) pairs).
# Usage: python run_external_assign.py <method> <mtx> <guides.txt> <cells.txt> <out_calls.csv> [n_iter]
import sys, os, shutil, tempfile
import numpy as np, scipy.io, pandas as pd

method   = sys.argv[1]
mtx      = sys.argv[2]
guides   = [l.strip() for l in open(sys.argv[3])]
cells    = [l.strip() for l in open(sys.argv[4])]
out_csv  = sys.argv[5]
n_iter   = int(sys.argv[6]) if len(sys.argv) > 6 else 500

m = scipy.io.mmread(mtx).tocsr()          # guides x cells
assert m.shape[0] == len(guides) and m.shape[1] == len(cells), (m.shape, len(guides), len(cells))

def write_calls(pairs):                    # pairs: list of (guide, cell)
    pd.DataFrame(pairs, columns=["guide","cell"]).to_csv(out_csv, index=False)

if method == "crispat":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt; plt.savefig = lambda *a, **k: None
    import anndata as ad, crispat
    A = ad.AnnData(X=m.T.tocsr().astype("float32"))     # cells x guides
    A.var_names = guides; A.obs_names = cells
    td = tempfile.mkdtemp(); tmp_h5 = os.path.join(td,"in.h5ad"); outdir = os.path.join(td,"out")+"/"
    try:
        A.write(tmp_h5)
        crispat.ga_poisson_gauss(tmp_h5, outdir, n_iter=n_iter)
        try:
            asg = pd.read_csv(os.path.join(outdir,"assignments.csv"))
        except (pd.errors.EmptyDataError, FileNotFoundError):
            asg = pd.DataFrame(columns=["cell","gRNA"])
        pairs = list(zip(asg["gRNA"], asg["cell"])) if len(asg) else []
        write_calls(pairs)
        print(f"crispat: {len(pairs)} calls", flush=True)
    finally:
        shutil.rmtree(td, ignore_errors=True)

elif method == "cleanser":
    # CLEANSER (Gersbach lab): per-guide Stan mixture (ambient vs native) with a
    # signal-free cell exposure. CLI: cleanser -i mtx -o posteriors --dc/--cs.
    # Posteriors output: dims header line, then  guide_idx \t cell_idx \t P(native).
    import subprocess
    td = tempfile.mkdtemp()
    int_mtx = os.path.join(td, "in.mtx"); post = os.path.join(td, "post.txt")
    scipy.io.mmwrite(int_mtx, m.astype(int), field="integer")   # guides x cells, integer
    cli = os.path.join(os.path.dirname(sys.executable), "cleanser")
    model = os.environ.get("CLNS_MODEL", "--dc")                # direct-capture (NB) default
    args = [cli, "-i", int_mtx, "-o", post, model,
            "-n", os.environ.get("CLNS_SAMPLES","400"),
            "-w", os.environ.get("CLNS_WARMUP","200"),
            "-c", os.environ.get("CLNS_CHAINS","2"),
            "-p", os.environ.get("CLNS_PARALLEL","8")]
    print("running:", " ".join(args), flush=True)
    subprocess.run(args, check=True)
    pairs = []
    with open(post) as fh:
        first = True
        for line in fh:
            if line.startswith("%"): continue
            if first: first = False; continue          # dims header line
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            gi, ci, p = int(parts[0]), int(parts[1]), float(parts[2])
            if p >= 0.5: pairs.append((guides[gi-1], cells[ci-1]))
    write_calls(pairs)
    shutil.rmtree(td, ignore_errors=True)
    print(f"cleanser: {len(pairs)} calls", flush=True)
else:
    raise SystemExit(f"unknown method {method}")
