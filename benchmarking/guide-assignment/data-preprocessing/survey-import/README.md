# survey-import — parse code for the perturbseq-survey datasets

Version-controlled copies of the gRNA-extraction code that produced the survey `grna_matrix.rds`
files, moved **out of the data directory** (`~/data/external/perturbseq-survey/`) per the lab
convention that data dirs hold no code. Full provenance:
[`../../grna-count-modeling/DATA.md`](../../grna-count-modeling/DATA.md).

- `parse_10x_h5_guides.R <input.h5> <out_dir> [regex]` — extract the CRISPR Guide Capture modality from a
  CellRanger `.h5` → `grna_matrix.rds` + `guide_features.csv`. Used for: **a549** (GSE236304),
  **erythroid_multiome** (GSE274113), **cd4_tcell** (GSE314342).
- `parse_10x_mtx_guides.R <mtx_dir> <out_dir> [regex]` — same, from a 10x mtx triplet. Used for:
  **cd8_tcell** (GSE279498), **dctap** (GSE303901), **endoc_t2d** (GSE273677, per-rep untar + rbind).
- `build_grna_matrix_{gastric_organoid,invivo_cortex,ipsc}.R` — per-dataset custom builders (non-10x
  formats: cell-identities CSV + guide library, legacy Seurat `.qs`, HipSci figshare UMI-count CSV).

Verified faithful: re-running the relevant parser on the local raw reproduces the committed
`grna_matrix.rds` **byte-identically** (checked for cd8_tcell; the others share the identical code path).

These scripts are the parse step of the per-study `import-<firstauthor>-<year>` repos (see DATA.md);
each import repo pairs a `download_raw` step (GEO/figshare) with the matching parser here.
