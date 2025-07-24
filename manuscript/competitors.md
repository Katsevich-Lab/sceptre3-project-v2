# SCEPTRE Competitor Methods

The following competitor methods to SCEPTRE are available in the `../../sceptre3-competitors` directory:

## Identified Competitors

Here are the competitors that are recent and/or have some parallel computation:

1. **[FR-Perturb](https://github.com/douglasyao/FR-Perturb)** (*Nature Biotechnology* 2023) - A Python (recommended) and R method for analyzing functional readouts in perturbation experiments  
2. **[crispat](https://github.com/velten-group/crispat)** (*Bioinformatics* 2024) - A Python package for CRISPR perturbation analysis with multiple statistical approaches
3. **[CLEANSER](https://github.com/Gersbachlab-Bioinformatics/CLEANSER)** (*Cell Genomics* 2025) - A Python method for analyzing single-cell CRISPR data with guide mixture modeling
4. **[Mixscale](https://satijalab.github.io/Mixscale/)** (*Nature Cell Biology* 2025) - An R package for perturbation scoring and differential expression analysis via Seurat
5. **[pySpade](https://pypi.org/project/pySpade/)** (*Genome Biology* 2025) - A Python package for single-cell perturbation analysis and differential expression
6. **[pertpy](https://pertpy.readthedocs.io/)** (*Nature Methods* minor revision 2025+) - A comprehensive Python toolkit for perturbation analysis in single-cell data

Not prioritized for benchmarking: MUSIC (*Nature Communications* 2019), scMAGeCK (*Genome Biology* 2020), GSFA (*Nature Methods* 2023).

## Comparative Summary Tables

### Data Import

| Software | Summary of Functionality | Data Structure |
|----------|-------------------------|----------------|
| pertpy | Import from AnnData, sparse matrices | AnnData objects, scipy.sparse matrices |
| CLEANSER | Import from Cell Ranger, Matrix Market files | Dictionaries, tuples, sparse matrices (CSR/DOK) |
| FR-Perturb | Import from AnnData (.h5ad) files | AnnData objects, numpy arrays |
| Mixscale | Import via Seurat objects | Seurat objects, sparse matrices |
| crispat | Import from Cell Ranger, CSV files | AnnData objects, sparse CSR matrices |
| pySpade | Import from Cell Ranger (.h5 files) | HDF5 files, numpy arrays |

### Guide Assignment

| Software | Summary of Functionality | Mode of Parallelization |
|----------|-------------------------|------------------------|
| pertpy | Most frequent guide, mixture models, threshold-based | GPU (JAX/CUDA/ROCm) |
| CLEANSER | Mixture model for ambient vs native gRNA presence | CPU (configurable parallel guide models) |
| crispat | 11 different assignment methods (GLM, mixture models, thresholds) | CPU (Dask parallelization) |

### Differential Expression Testing

| Software | Summary of Functionality | Mode of Parallelization |
|----------|-------------------------|------------------------|
| pertpy | PyDESeq2, EdgeR, simple tests (t-test, Wilcoxon) | CPU (PyDESeq2 multithreading only) |
| FR-Perturb | Factorize-Recover algorithm, elastic net regression, permutation-based significance | Manual (batch processing for HPC/clusters) |
| Mixscale | Weighted negative binomial GLM, perturbation scoring | None (single-threaded R) |
| pySpade | Hypergeometric test, t-test, gamma distribution significance | CPU (configurable threads) |

## Detailed Analysis: pertpy

### Data Import and Data Structures
- Uses AnnData objects as primary data structure
- Provides dataset loading functions (e.g., `pt.dt.papalexi_2021()`)
- Supports sparse matrix formats (scipy.sparse.csr_matrix)
- Requires guide RNA data as separate AnnData object with counts per cell
- No direct Cell Ranger import utilities found in documentation

### Guide Assignment
- **Methods**: Most frequent guide, mixture models, threshold-based assignment
- **CPU Parallelization**: Not explicitly documented
- **GPU Acceleration**: ✅ Uses JAX for GPU acceleration in mixture models
- **Distributed Computing**: ❌ Not implemented
- **Technology**: JAX arrays, NumPyro for MCMC sampling, automatic differentiation

### Differential Expression Testing
- **Methods**: PyDESeq2, EdgeR (via rpy2), simple tests (t-test, Wilcoxon)
- **CPU Parallelization**: ✅ Limited - only PyDESeq2 uses multithreading (n_cpus parameter, defaults to all CPUs)
- **GPU Acceleration**: ❌ Not implemented for differential expression methods
- **Distributed Computing**: ❌ Not implemented
- **Large Dataset Support**: Claims sparse and memory-efficient implementations for millions of cells

### Special Features for Large Datasets
- Sparse matrix support throughout pipeline
- Memory-efficient implementations
- JAX-based GPU acceleration (primarily for guide assignment)
- Numba acceleration for some operations
- Designed to handle datasets from thousands to millions of cells

## Detailed Analysis: CLEANSER

### Data Import and Data Structures
- **Cell Ranger Support**: ✅ Includes `cr2cleanser` utility to convert Cell Ranger outputs
- **Data Format**: Matrix Market file format with columns: "Guide ID", "Cell ID", "Guide Count"
- **Input Requirements**: Filters out gene expression data, keeps only CRISPR Guide Capture entries
- **Data Structure**: Uses CmdStan for statistical modeling, probabilistic approach

### Guide Assignment
- **Method**: Mixture of two distinct distributions to model ambient vs native gRNA presence
- **Algorithm**: Accounts for gRNA-specific and cell-specific biases
- **Output**: Generates probability values for whether gRNA is expressed natively in each cell
- **Experiment Types**: Supports both direct capture and CROP-seq experiments
- **CPU Parallelization**: ✅ Configurable parallel processing of guide models (`-p` parameter)
- **GPU Acceleration**: ❌ Uses CmdStan (CPU-based MCMC sampling)
- **Distributed Computing**: ❌ Not implemented

### Differential Expression Testing
- **Capability**: ❌ Not implemented - CLEANSER is specifically for guide assignment only

### Special Features for Large Datasets
- **CPU Parallelization**: Multiple guide models can run in parallel
- **MCMC Configuration**: Configurable chains, samples, and warmup iterations
- **Quality Control**: Comprehensive QC metrics (MOI, coverage, sample statistics)
- **Normalization**: Configurable guide count normalization with low-pass filtering

## Detailed Analysis: FR-Perturb

### Data Import and Data Structures
- **Data Format**: AnnData (.h5ad) files for expression matrices
- **Perturbation Tracking**: Binary indicator matrix or metadata column in AnnData
- **Data Structure**: AnnData objects, numpy arrays for matrix operations
- **Input Requirements**: Expression data with perturbation status information

### Guide Assignment
- **Capability**: ❌ Not implemented - focuses on downstream differential expression analysis
- **Approach**: Assumes guide assignment already completed upstream

### Differential Expression Testing
- **Primary Method**: Factorize-Recover algorithm (Sharan et al. 2019) for gene expression effect sizes
- **Alternative Method**: Elastic net regression (`--elastic-net` flag)
- **Output**: Log-fold changes with permutation-based p-values
- **CPU Parallelization**: ❌ Manual implementation required
- **GPU Acceleration**: ❌ Not implemented
- **Distributed Computing**: ✅ Designed for HPC/cluster environments but requires manual orchestration

### Special Features for Large Datasets
- **Memory Efficiency**: Randomized partial SVD for datasets with millions of cells
- **Batch Processing**: `--temp-out` flag enables p-value computation in separate batches
- **Manual Parallelization**: Three-step process requiring user implementation:
  1. Main script with temporary output
  2. Parallel batch jobs (`compute_pvalues_batched.py`)
  3. Results combination
- **HPC Integration**: Designed for compute clusters via separate job submission
- **Flexibility**: Users control parallelization strategy for their specific computing environment

## Detailed Analysis: Mixscale

### Data Import and Data Structures
- **Data Format**: Seurat objects (primary data structure)
- **Import Methods**: Standard Seurat import functions for count matrices
- **Data Structure**: Seurat objects with sparse matrices, requires perturbation metadata
- **Requirements**: Guide RNA assignments and target gene labels in metadata

### Guide Assignment
- **Capability**: ✅ Supports per-gRNA analysis via "fine mode"
- **Methods**: Target gene grouping, non-targeting control handling
- **Multi-guide Support**: Can handle multiple gRNAs per target gene
- **Perturbation Scoring**: Uses RunMixscale() to calculate perturbation strength scores

### Differential Expression Testing
- **Primary Method**: Weighted negative binomial GLM using glmGamPoi
- **Weighted Regression**: Incorporates perturbation scores as weights (Run_wmvRegDE)
- **Standard DE**: Falls back to unweighted analysis when scores unavailable
- **CPU Parallelization**: ❌ Single-threaded processing
- **GPU Acceleration**: ❌ Not implemented
- **Distributed Computing**: ❌ Limited by R memory constraints

### Special Features for Large Datasets
- **Memory Limitations**: Constrained by standard R memory capacity
- **Single-core Processing**: No built-in parallelization
- **Advanced Analysis**: PCA permutation testing, hierarchical clustering, multi-CCA
- **Enrichment Testing**: Fisher's exact test, rank biased overlap (RBO)
- **Seurat Integration**: Full compatibility with Seurat ecosystem

## Detailed Analysis: crispat

### Data Import and Data Structures
- **Cell Ranger Support**: ✅ Direct import via create_anndata_from_cellranger()
- **Data Formats**: Cell Ranger output, CSV files, AnnData objects
- **Data Structure**: AnnData objects with sparse CSR matrices
- **Multi-batch Support**: Handles batch information in adata.obs['batch']

### Guide Assignment
- **Methods**: 11 different assignment methods across 4 categories:
  - Independent: UMI threshold
  - Across gRNAs: Maximum, ratio threshold
  - Across cells: Gaussian mixture, Poisson-Gaussian mixture
  - Across both: Beta mixtures (2&3-component), GLMs (Poisson/NB/Binomial), quantiles
- **CPU Parallelization**: ✅ Dask-based parallelization for computationally intensive methods
- **GPU Acceleration**: ❌ CPU-only processing
- **Distributed Computing**: ✅ HPC-friendly with gRNA subset processing

### Differential Expression Testing
- **Capability**: ❌ Not implemented - focuses on guide assignment only
- **Integration**: Provides tutorials for downstream tools (SCEPTRE, Seurat, scanpy)

### Special Features for Large Datasets
- **Dask Parallelization**: Automatic parallelization across CPU cores
- **Memory Efficiency**: Sparse matrix storage, optional downsampling
- **HPC Support**: Configurable for cluster computing environments
- **Method Comparison**: Tools for comparing assignment strategies

## Detailed Analysis: pySpade

### Data Import and Data Structures
- **Cell Ranger Support**: ✅ Direct import from Cell Ranger outs folder (.h5 files)
- **Data Formats**: Cell Ranger filtered_feature_bc_matrix.h5, CSV, pickle files
- **Data Structure**: HDF5 files for storage, numpy arrays for processing
- **Guide Integration**: Requires separate sgRNA matrix with cell barcode consistency

### Guide Assignment
- **Capability**: ❌ Not implemented - assumes guide assignment completed upstream
- **Data Processing**: Removes experimental doublets and sgRNA outlier cells
- **Quality Control**: Filters based on sgRNA noise and cell multiplex information

### Differential Expression Testing
- **Primary Method**: Hypergeometric test for genome-wide differential expression
- **Statistical Approach**: Compares observed vs random cell selection backgrounds
- **Significance Testing**: Gamma distribution approximation for adjusted p-values
- **Analysis Types**: Local hits (±2 Mb) and global hits (genome-wide)
- **CPU Parallelization**: ✅ Configurable thread support for DE analysis
- **GPU Acceleration**: ❌ Not implemented
- **Distributed Computing**: ❌ Single-machine processing

### Special Features for Large Datasets
- **Threading Support**: Configurable number of threads for barcode comparisons
- **Batch Processing**: Supports iteration-based random background generation
- **Memory Efficiency**: HDF5 storage with optional compression
- **Comprehensive Pipeline**: End-to-end from Cell Ranger to Manhattan plots
- **Human Genome Focus**: Currently supports human genome only