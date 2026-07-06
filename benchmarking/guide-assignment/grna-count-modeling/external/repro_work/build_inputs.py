#!/usr/bin/env python
"""
Faithfully replicate fishash_analysis 21_process_barnyard_data.R for the 4
species-mix samples, producing:
  - <sample>_grna.h5ad : AnnData(cells x 200 guides), counts as int  (input to geomux)
  - <sample>_grna_counts.mtx : guides x cells grna counts (input to fishash / convert)
  - <sample>_meta.csv : per-cell GEX QC fields needed for Table-2 scoring

guide var names are nt_1..nt_200 in matrix row order; guide_type homo for nt_1..100.
No QC / no species filtering at this stage (matches the authors).
"""
import gzip
import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

def _config_path(var):
    # mirror ~/.Rprofile's .get_config_path: resolve a ~/.research_config variable
    import subprocess
    return subprocess.run(["sh", "-c", f"source ~/.research_config; echo ${var}"],
                          capture_output=True, text=True).stdout.strip()

BARN = os.path.join(_config_path("LOCAL_EXTERNAL_DATA_DIR"), "liu-2025-cleanser/GSE272457")
OUT = os.path.dirname(os.path.abspath(__file__))

SAMPLES = {
    "mix0hr_Cropseq":        "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
    "mix72hr_Cropseq":       "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix",
    "mix0hr_DirectCapture":  "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
    "mix72hr_DirectCapture": "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
}


def read_features(stub):
    rows = []
    with gzip.open(os.path.join(BARN, stub + "_features.tsv.gz"), "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append(p[:3])
    return pd.DataFrame(rows, columns=["ID", "Symbol", "type"])


def read_barcodes(stub):
    with gzip.open(os.path.join(BARN, stub + "_barcodes.tsv.gz"), "rt") as f:
        return [l.strip() for l in f]


for name, stub in SAMPLES.items():
    feat = read_features(stub)
    bc = read_barcodes(stub)
    m = scipy.io.mmread(os.path.join(BARN, stub + "_matrix.mtx.gz")).tocsr()
    assert m.shape[0] == len(feat)
    assert m.shape[1] == len(bc)

    is_guide = (feat["type"] == "CRISPR Guide Capture").values
    is_gex = (feat["type"] == "Gene Expression").values
    assert (is_guide | is_gex).all()

    ref = np.where(feat["Symbol"].str.startswith("GRCh38"), "homo",
                   np.where(feat["Symbol"].str.startswith("mm10"), "mus", "??"))
    # authors only require this for the gene-expression rows
    assert (ref[is_gex] != "??").all()

    gex = m[is_gex]
    grna = m[is_guide]  # guides x cells
    gex_ref = ref[is_gex]

    homo_sum_gex = np.asarray(gex[gex_ref == "homo"].sum(axis=0)).ravel()
    mus_sum_gex = np.asarray(gex[gex_ref == "mus"].sum(axis=0)).ravel()

    # mito + features for QC
    sym = feat["Symbol"].values[is_gex]
    is_mito = np.array([s.startswith("GRCh38_MT-") or s.startswith("mm10___mt-") for s in sym])
    mito_sum = np.asarray(gex[is_mito].sum(axis=0)).ravel()
    features_gex = np.asarray((gex > 0).sum(axis=0)).ravel()

    # guide var names: authors use rownames make.unique(Symbol); for guides the
    # Symbol IS nt_1..nt_200. guide_type homo for nt_1..100.
    guide_sym = feat["Symbol"].values[is_guide]
    # sanity: should be nt_1 .. nt_200
    gidx = np.array([int(s.replace("nt_", "")) for s in guide_sym])
    guide_type = np.where(np.isin(guide_sym, ["nt_%d" % i for i in range(1, 101)]),
                          "homo_guide", "mus_guide")

    # build AnnData: cells x guides, int counts (matches convert_to_anndata.R)
    X = grna.T.tocsr().astype(np.int64)
    adata = ad.AnnData(X)
    adata.obs_names = [b for b in bc]
    adata.var_names = [s for s in guide_sym]
    adata.obs["batch"] = 0
    adata.write(os.path.join(OUT, name + "_grna.h5ad"))

    # save grna counts mtx (guides x cells) for fishash
    scipy.io.mmwrite(os.path.join(OUT, name + "_grna_counts.mtx"), grna.tocoo())

    meta = pd.DataFrame({
        "barcode": bc,
        "homo_sum_gex": homo_sum_gex,
        "mus_sum_gex": mus_sum_gex,
        "mito_sum": mito_sum,
        "features_gex": features_gex,
    })
    meta.to_csv(os.path.join(OUT, name + "_meta.csv"), index=False)

    # guide order/type file
    pd.DataFrame({"guide": guide_sym, "gidx": gidx, "guide_type": guide_type}).to_csv(
        os.path.join(OUT, name + "_guides.csv"), index=False)

    print(f"[{name}] cells={len(bc)} guides={int(is_guide.sum())} gex={int(is_gex.sum())} "
          f"homo_guides={int((guide_type=='homo_guide').sum())}")
