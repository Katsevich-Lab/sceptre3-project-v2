# Bioconductor submission cover letter — `sceptre`

To the Bioconductor reviewers,

We submit `sceptre` for inclusion in Bioconductor. `sceptre` is an R package for statistically rigorous, scalable, and user-friendly analysis of single-cell CRISPR screen data. The package implements an end-to-end pipeline: data import from 10X CellRanger, Parse CRISPR Detect, or R objects; gRNA-to-cell assignment via one of three principled methods; quality control on cells and gRNAs; a calibration check that verifies false-positive control on the dataset under analysis; a power check; and a discovery analysis identifying high-confidence perturbation-gene links. The methodology relies on several novel statistical and computational algorithms designed for high statistical accuracy and computational throughput on large single-cell CRISPR perturbation datasets.

## Why Bioconductor

> **Placeholder — Tim to draft.** Justify hosting `sceptre` on Bioconductor despite the package currently having no direct Bioconductor software dependencies (its `Imports` are all on CRAN). The core of the argument: **explain why the `ondisc` data structure is better suited to `sceptre`'s needs than `SingleCellExperiment`** — single-cell CRISPR screens are large enough that an on-disk-backed representation is essential, and `ondisc` is purpose-built for the access patterns `sceptre` actually performs (large sparse expression matrices, paired gRNA-by-cell matrices, per-pair iteration). Brief supporting points worth including:
>
> - The single-cell CRISPR screen domain sits squarely within Bioconductor's single-cell ecosystem.
> - `sceptre`'s user community overlaps substantially with users of `SingleCellExperiment`, `scran`, and related Bioc infrastructure, even though `sceptre` does not currently depend on those classes.
> - **Kasper Hansen will back up this argument.**
>
> A brief paragraph (5–8 sentences) is typically sufficient.

## Notes on BiocCheck output

A small number of BiocCheck findings warrant brief acknowledgment:

### `cat()` / `print()` outside `show()` methods — 2 intentional sites

Two flagged sites remain after the previous cleanup pass; both are intentional:

- **`R/s4_helpers.R:146`** — this `cat()` is the body of a real S4 `show()` method for the `sceptre_object` class. BiocCheck's static parser does not recognize the surrounding `setMethod("show", ...)` and reports the `cat()` anyway. Replacing it with `message()` would break the `show()` contract (which writes to `stdout`, not to a condition).
- **`R/s4_analysis_functs_2.R:501`** — this `print()` is wrapped in a `sink()/print()/sink()` block that writes the analysis summary to disk. We considered rewriting via `writeLines(capture.output(print(x)))`, but BiocCheck flags the bare `print` token identically in either form, and the `sink()` idiom is clearer and standard for this kind of disk capture.

We believe both are BiocCheck false positives in spirit and have left them in place.

### `\dontrun{}` / `\donttest{}` rate (3%)

`man/sceptre-package.Rd` carries a `\donttest{}` block for the package-level orientation example. BiocCheck reports a 3% rate across all man pages, which is at or below typical community norms.

## Subscription

The maintainer (Timothy Barry, <tbarry@hsph.harvard.edu>) is subscribed to the Bioc-Devel mailing list at <https://stat.ethz.ch/mailman/listinfo/bioc-devel>.

## Acknowledgments

This work was supported by Wharton Analytics (Eugene Katsevich) and by U.S. National Science Foundation grants DMS-2113072 and DMS-2310654 (Eugene Katsevich).

Thank you for considering `sceptre` for inclusion in Bioconductor.

Sincerely,

The sceptre authors
