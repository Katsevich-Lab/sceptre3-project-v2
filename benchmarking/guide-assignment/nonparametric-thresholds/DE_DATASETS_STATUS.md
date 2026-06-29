# DE pushthrough dataset status

For each dataset, status of (GEX matrix, guide matrix) availability and the
location of the sceptre-ready triple (`response_matrix.rds`,
`grna_matrix_aligned.rds`, derived from the local raw data).

| dataset | n_cells | n_genes | n_guides | sceptre-ready | location |
|---|---|---|---|---|---|
| **gasperini** | 207,324 | 13,135 | 13,077 | ✓ ondisc | `~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini/sceptre-pipeline/` |
| **replogle-rd7** | 616,184 | ~16,000 | 2,666 | ✓ ondisc | `~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline/` |
| **a549** | 25,220 | 36,601 | 253 | ✓ in-memory | `~/data/external/perturbseq-survey/a549_crispri_perturbseq_GSE236304/sceptre/` |
| **cd8_tcell** | 2,426 | 32,285 | 44 | ✓ in-memory | `~/data/external/perturbseq-survey/cd8_tcell_perturbcite_GSE279498/sceptre/` |
| **dctap_highmoi** | 10,932 | **93** | 110 | ✓ in-memory | `~/data/external/perturbseq-survey/dctap_k562_highmoi_GSE303901/sceptre/` |
| **dctap_lowmoi** | 5,397 | **93** | 110 | ✓ in-memory | `~/data/external/perturbseq-survey/dctap_k562_lowmoi_GSE303901/sceptre/` |
| **endoc** | 41,642 | 36,601 | 225 | ✓ in-memory | `~/data/external/perturbseq-survey/endoc_t2d_perturbseq_GSE273677/sceptre/` (7 GWAS reps combined) |
| **multiome_erythroid** | 14,525 | 36,601 | 63 | ✓ in-memory | `~/data/external/perturbseq-survey/perturb_multiome_erythroid_GSE274113/sceptre/` |
| **tcell_cd4** | 137,716 | 18,130 | **27,272** | ✓ in-memory | `~/data/external/perturbseq-survey/tcell_cd4_perturbseq_GSE314342/sceptre/` (1 rep only on disk) |
| ~~gastric_organoid~~ | — | — | 19,168 | ✗ gRNA only | `cell_identities.csv.gz` (preprocessed); GEX needs GEO download (GSM8602xxx_RAW.tar) |
| ~~invivo_cortex~~ | — | — | 17 | ✗ gRNA only | qs Seurat object has GEX assay but requires **R 4.4 + qs v1** to load (legacy format) |
| ~~ipsc~~ | — | — | 6,824 | ✗ gRNA only | figshare has only guide-UMI counts + cell metadata; GEX is on ENA/SRA (~30 GB raw) |
| barnyard ×4 | — | — | 200 NT | n/a | no real perturbations (species-mixing test); not useful for DE |

**Total useable for DE pushthrough: 9 datasets** spanning 2,426 → 616,184 cells
and 44 → 27,272 guides, including both the lab's reference datasets and 7
diverse survey datasets.

## Per-dataset caveats

- **dctap (high/lowMOI)**: targeted GEX panel — only **93 genes** profiled
  (not whole transcriptome). DE pushthrough is informative for the panel genes
  only; absolute power numbers will look different from whole-transcriptome
  screens.
- **tcell_cd4**: only 1 GEX rep on disk (`GSM9393828_CD4i_R1L01_D1_Rest`); the
  full study has ~10 reps in GEO. The 27,272 guides came from a combined parse
  the lab did. To match cells across the full grna_matrix.rds with the loaded
  rep we may need to either (a) re-parse the matching rep of grna only,
  (b) download the remaining reps, or (c) just intersect cell barcodes.
- **endoc**: 7 GWAS reps × ~6k cells each = 41,642 cells, prefixed barcodes
  per rep to avoid collisions. The 2 RQC reps use a different library and were
  excluded.
- **multiome_erythroid**: GEX is the RNA modality; ATAC is dropped.

## To complete the remaining 3 datasets

- **gastric_organoid**: download the GSE280506 mtx tar (per-cell expression).
  Likely ~few GB. CRISPRa = activation; expect on-target = upregulation.
- **invivo_cortex**: install R 4.4 + qs (v1.5.x) in a side venv:
  `R 4.4 -e 'install.packages("qs", repos="https://cloud.r-project.org")'`,
  then re-run `scripts/sim_de_extract_survey.R` with an extractor for the qs
  Seurat object. The gRNA was already extracted this way.
- **ipsc**: download GEX from ENA project PRJEB75103 (the HipSci scRNA-seq).
  Pricey (~30 GB raw + ~10 GB processed). Likely defer.
