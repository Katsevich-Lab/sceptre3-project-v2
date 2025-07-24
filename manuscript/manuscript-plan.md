# Manuscript outline

## Description and demonstration of new functionality (Tim)

- Package
- Pipeline

## Benchmarking (Louis)

- Data import
- gRNA assignment
- DE testing

# Submission venue

I suggest a [Genome Biology Short Report](https://genomebiology.biomedcentral.com/submission-guidelines/preparing-your-manuscript/short-report): 

> Genome Biology publishes Short Reports that are concise studies of high quality and broad interest. Short Reports can present new research findings, or can present a new method or software. For a paper to be suitable as a Short Report, the main results (or description, for a method or software paper) should be able to be clearly shown in a **maximum of 2 figures or tables.** The main text should be around 1000-1500 words, excluding abstract, references, and figure legends, and should contain no headings.
>
> **Short Report computational method articles may be new versions of previously published methods, and they should be benchmarked against the original version as well as any other similar methods published in the interim.**  A Short Report may also be a visualization tool, software, or web-based application that would be widely used in the field and examples of how the tool will be used should be provided. Short Report methods may also be wet lab methods that must be shown to be of broad utility. Short methods and software do not need to include novel biological findings, but it should be clear that the method represents a significant advance.

## Exemplary manuscripts

In recent years, two manuscripts have reported updates to previously published softwares as Genome Biology Short Reports: [Kraken 2 (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) and [JACUSA2 (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0).

### Kraken 2

A k‑mer–based classifier that tags metagenomic sequencing reads with taxonomic labels on the fly. Kraken 2 replaces the original hash scheme with a compact probabilistic table, slashing RAM needs by ≈85 %, boosting speed about five‑fold, and adding a translated‑search mode that improves sensitivity to highly divergent (e.g., viral) sequences.

#### Benchmarking conducted

| Functionality        | Competitors                                   | Datasets         | Evaluation metrics         |
|----------------------|-----------------------------------------------|------------------|----------------------------|
| Nucleotide search    | Kraken 1, KrakenUniq, CLARK, Centrifuge       | 1 simulated, 1 real | accuracy, runtime, memory |
| Translated search    | Kraken 2X, Kaiju                              | 1 simulated, 1 real | accuracy, runtime, memory |

#### Manuscript structure

Figure 1: Illustration of algorithmic differences between Kraken 1 and Kraken 2, and a preview of benchmarking of Kraken 1 versus Kraken 2. 

Figure 2: Benchmarking with more competitor methods, evaluation metrics, and datasets.

Supplementary figures: More benchmarking and other details

### JACUSA2

A read‑signature framework for pinpointing RNA modifications from sequencing data. JACUSA 2 broadens support to both Illumina and Nanopore platforms, captures richer signatures (substitutions, indels, read‑arrest events), integrates replicate/condition information, and runs far faster than earlier JACUSA 1.x versions.

#### Benchmarking conducted

| Functionality                           | Competitors                                | Datasets | Evaluation metrics  |
|-----------------------------------------|--------------------------------------------|----------|---------------------|
| Illumina m⁶A assays                     | MAZTER‑mine, JACUSA 1.0, JACUSA 1.3         | 2 real   | accuracy, runtime   |
| Nanopore direct‑RNA m⁶A detection       | ELIGOS2, xpore, Epinano 1.2                | 2 real   | accuracy, runtime   |
| Nanopore rRNA pseudouridine detection   | nanopseudo_U                               | 1 real   | accuracy, runtime   |


#### Manuscript structure

Figures 1-3: Demonstration of three kinds of software functionality.

Supplementary figures: Benchmarking