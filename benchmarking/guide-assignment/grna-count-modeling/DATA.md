# Dataset provenance & acquisition

Provenance for every dataset in `scripts/datasets.R` (`dataset_paths()`), and how to obtain each the
lab-standard way. Data locations resolve through `~/.research_config` (via `.get_config_path`), never
hardcoded paths. Processed artifacts travel between machines **only via HPCC** (`hpcc pull`/`hpcc push`),
per the [lab wiki](https://github.com/Katsevich-Lab/lab-resources/wiki/Synchronizing-Data-and-Code).

## Two provenance classes

1. **Lab-native** (`gasperini`, `replogle-rd7`): built by `../data-preprocessing/make_{gasperini,replogle-rd7}.R`
   from the lab's canonical import repos ([import-gasperini-2019-v3](https://github.com/Katsevich-Lab/import-gasperini-2019-v3),
   [import-replogle-2022](https://github.com/Katsevich-Lab/import-replogle-2022)) → external ODMs
   (`LOCAL_GASPERINI_2019_V3_DATA_DIR`, `LOCAL_REPLOGLE_2022_DATA_DIR`). Nothing to do — already reproducible.

2. **Survey collection** (`~/data/external/perturbseq-survey/`): a curation assembled for this project by
   downloading each study's **GEO/figshare processed output** and extracting the gRNA-count matrix. These are
   **guide-count matrices** (`grna_matrix.rds`), distinct from any from-raw lab reprocessing (see the
   a549/dctap note below). Each now gets a lab-standard `import-<firstauthor>-<year>` repo that reproduces
   this exact acquisition.

## Per-dataset provenance

| registry key | study | accession | my source file | parse | reproducibility | import repo |
|---|---|---|---|---|---|---|
| gasperini | Gasperini 2019 | GSE120861 | lab ODM | `make_gasperini.R` | lab-native | import-gasperini-2019-v3 |
| replogle-rd7 | Replogle 2022 | figshare (Weissman) | lab ODM | `make_replogle-rd7.R` | lab-native | import-replogle-2022 |
| barnyard ×4 | Liu 2025 (CLEANSER) | GSE272457 | mtx triplets (GEO suppl) | import raw guides (all cells) → `make_barnyard.R` species-purity QC | full | **import-liu-2025** ✅ |
| a549 | Sakellaropoulos 2024 | GSE236304 | `filtered_feature_bc_matrix.h5` | `read_10x_h5` | full (raw h5 local) | see note ↓ (import-Sakellaropoulos-2024 exists, SRA path) |
| dctap ×2 | Ray 2025 (ENCODE DC-TAP) | GSE303901 | mtx triplet / tar.gz | `read_10x_mtx` | full (raw local) | see note ↓ (`ray-2025` on HPCC, different processing) |
| cd8_tcell | Chung 2024 | GSE279498 | mtx triplet (GEO suppl) | `read_10x_mtx` | **full** (raw triplet local) | **import-chung-2024** ✅ |
| endoc_t2d | Cao 2024 | GSE273677 | `GSM8434996` GWAS-rep1 triplet (GEO) | `read_10x_mtx` (single rep) | **full** (raw local) | **import-cao-2024** ✅ |
| erythroid_multiome | Caulier 2025 | GSE274113 | `..._filtered_feature_bc_matrix_1.h5` | `read_10x_h5` | **full** (raw h5 local) | **import-caulier-2025** ✅ |
| gastric_organoid | Chen 2025 | GSE280506 | `cell_identities.csv.gz` | `build_grna_matrix.R` | partial (gRNA only; GEX = GEO RAW.tar) | **import-chen-2025** ✅ |
| invivo_cortex | Zheng 2023 | GSE249416 | `Perturb_sg.qs.gz` (Seurat) | `build_grna_matrix.R` (R 4.4 + qs) | full (legacy `.qs`) | **import-zheng-2023** ✅ |
| ipsc | Feng 2025 (HipSci) | figshare 27989294 | `Guide-UMI-Counts.csv.gz` + metadata | `build_grna_matrix.R` | gRNA only (GEX on ENA/SRA) | **import-feng-2025** ✅ |
| cd4_tcell | Zhu (unpublished) | GSE314342 | `GSM9393828..._filtered_feature_bc_matrix.h5` | `read_10x_h5` | partial (1 of ~10 GEX reps; unpublished) | **import-zhu-2025** (private) |
| mccutcheon | McCutcheon 2023 | GSE218988 | `GSM6761464_CRISPRi_D1.tar.gz` | not versioned | partial (raw tar local) | **none yet** — one-off local object, no import repo |

GEO supplementary URLs follow `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE<XXX>nnn/GSE<XXXXX>/suppl/<file>`
(e.g. `GSE273nnn/GSE273677/suppl/GSE273677_RAW.tar`).

## Note: a549 & dctap vs the lab-native reprocessings (investigated)

The lab has these two studies on HPCC, but as **different objects** from my guide-count matrices — do **not**
swap (it would change the data under committed writeups):
- **a549** — my matrix (253 guides × 25,220 cells, `Non_Targeting_Human_CRi_*`) is from GEO's authors-processed
  `filtered_feature_bc_matrix.h5`. `sakellaropoulos-2024-sra/` (on HPCC; `hpcc pull`) — from
  [import-Sakellaropoulos-2024](https://github.com/Katsevich-Lab/import-Sakellaropoulos-2024) — is a
  **from-SRA Cell Ranger** reprocessing (per-SRR outputs, no assembled guide-count matrix).
- **dctap** — my matrix (110 guides × 10,932 cells, `safe_*`) is from GEO mtx. `ray-2025/` (on HPCC)
  is a **Cell Ranger → sceptre** reprocessing whose guides are annotated differently
  (`K562_Random_Screen_Crop_*`) — a distinct curation, not a reformat of mine.

Recommendation: keep my GEO-derived versions (preserves results) and document the correspondence, as done here.

## "No code under `data/`" (lab convention)

The survey parse scripts that previously lived in `~/data/external/perturbseq-survey/` (`parse_10x_h5_guides.R`,
`parse_10x_mtx_guides.R`, per-dataset `build_grna_matrix.R`) violated the convention that data dirs hold no
code. They have been **moved** into `../data-preprocessing/survey-import/` (versioned) and ported into the
per-study `import-<author>-<year>` repos (below); the data-dir originals were deleted.

## Import-repo status

Each survey dataset's `import-<author>-<year>` repo downloads the exact GEO/figshare source I used and
runs the versioned parse to reproduce the `grna_matrix.rds` the benchmarking consumes.
Config var per repo: `LOCAL_<AUTHOR>_<YEAR>_DATA_DIR=$LOCAL_EXTERNAL_DATA_DIR"<author>-<year>/"`.

**Import = raw, no QC.** Import repos only *reformat* the published data (select the CRISPR-Guide-Capture
modality / extract counts) and keep **all cells** — cell QC belongs downstream in this project. barnyard is
the worked example: `import-liu-2025` outputs raw all-cell matrices, and `../data-preprocessing/make_barnyard.R`
applies the species-purity gate. (One borderline case: `import-chen-2025` builds from the authors'
cell-identities table and so drops unassigned cells — structural to that data form, not a quality gate.)

| repo | status | verified |
|---|---|---|
| [import-caulier-2025](https://github.com/Katsevich-Lab/import-caulier-2025) (erythroid) | ✅ live, public | byte-identical (63 × 14,525) |
| [import-zhu-2025](https://github.com/Katsevich-Lab/import-zhu-2025) (cd4_tcell) | ✅ live, private (unpublished) | byte-identical (27,272 × 137,716) |
| [import-chung-2024](https://github.com/Katsevich-Lab/import-chung-2024) (cd8_tcell, HKC045) | ✅ live, public | byte-identical (44 × 2,426) |
| [import-chen-2025](https://github.com/Katsevich-Lab/import-chen-2025) (gastric) | ✅ live, public | byte-identical (96 × 42,324) |
| [import-feng-2025](https://github.com/Katsevich-Lab/import-feng-2025) (ipsc) | ✅ live, public | byte-identical (6,824 × 322,746) |
| [import-cao-2024](https://github.com/Katsevich-Lab/import-cao-2024) (endoc_t2d) | ✅ live, public | byte-identical (225 × 8,329) |
| [import-zheng-2023](https://github.com/Katsevich-Lab/import-zheng-2023) (invivo_cortex) | ✅ live, public | byte-identical (17 × 11,688) — requires R 4.4 + `qs` |
| [import-liu-2025](https://github.com/Katsevich-Lab/import-liu-2025) (barnyard ×4) | ✅ live, public | raw matrices (all cells; no QC). The species-purity QC (`frac_homo > 0.9 \| < 0.1`, recovered from the 2026-06-25 ingest session) lives in `../data-preprocessing/make_barnyard.R` and reproduces `input_data/barnyard_*/sceptre/grna_matrix.rds` byte-identically (8,265→7,453 / 8,959→8,793 / 11,512→8,902 / 8,207→8,057). |

Per-dataset parse code lives in `../data-preprocessing/survey-import/`; each repo ports the relevant one.
Partial datasets (gastric/invivo/ipsc) fetch GEX **count** matrices from GEO/figshare where available (not SRA);
GEX is skipped only where no count matrix is published.

## Obtaining the data (collaborator quick-start)

```sh
# lab-native + already-imported survey datasets are on HPCC; pull the benchmarking inputs:
hpcc pull SCEPTRE3 benchmarking/guide_assignment/input_data      # barnyard, gasperini, replogle-rd7
# survey datasets: reproduce via the per-study import repo, e.g.
git clone https://github.com/Katsevich-Lab/import-caulier-2025 && cd import-caulier-2025 && ./run_all.sh
```
