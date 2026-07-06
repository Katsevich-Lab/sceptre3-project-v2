# survey-import — shared 10x gRNA-extraction utilities

Shared, argument-driven parsers for the standard-10x survey datasets (no hardcoded paths). They are the
parse step that each dataset's `import-<firstauthor>-<year>` repo self-contains a port of, and they double
as the tool for the full-registry reproducibility audit in
[`../../grna-count-modeling/DATA.md`](../../grna-count-modeling/DATA.md).

- `parse_10x_h5_guides.R <input.h5> <out_dir> [regex]` — extract the CRISPR Guide Capture modality from a
  CellRanger `.h5` → `grna_matrix.rds` + `guide_features.csv`. Verifies: **a549** (GSE236304),
  **erythroid_multiome** (GSE274113), **cd4_tcell** (GSE314342).
- `parse_10x_mtx_guides.R <mtx_dir> <out_dir> [regex]` — same, from a 10x mtx triplet. Verifies:
  **cd8_tcell** (GSE279498, HKC045), **dctap** (GSE303901), **endoc_t2d** (GSE273677, GSM8434996),
  **mccutcheon** (GSE218988, GSM6761464).

Re-running the relevant parser on the local raw reproduces the committed `grna_matrix.rds` **byte-identically**.

The **non-10x** datasets (gastric = cell-identities CSV, invivo = legacy Seurat `.qs`, ipsc = HipSci figshare
CSV) have their bespoke builders in their own import repos — `import-chen-2025`, `import-zheng-2023`,
`import-feng-2025` — not here (those were dataset-specific + path-hardcoded; the config-based repo versions
supersede them).
