[03:06, 02/11/2026] Nisrine Hassini: # EpiRegion

[![R-CMD-check](https://github.com/hassininisrine05/EpiRegion/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hassininisrine05/EpiRegion/actions)

*Region-level statistical methods for episignature detection from long-read methylation data — a simulation-based benchmark.*


## Question

Les épisignatures reposent sur des différences de méthylation …
[03:09, 02/11/2026] Nisrine Hassini: # EpiRegion

*An exploratory R toolkit for regional DNA methylation analysis using long-read sequencing data.*

## Overview

EpiRegion is an open-source bioinformatics project exploring statistical methods for detecting subtle DNA methylation differences across genomic regions.

The project investigates whether aggregating methylated-read counts across neighboring CpG sites can improve the detection of small methylation changes compared with single-CpG analysis, particularly under moderate sequencing coverage.

The broader objective is to explore computational strategies relevant to DNA methylation profiling and episignature research using long-read sequencing technologies

*Project status:* Research prototype based on simulated data. Performance has not yet been validated on real sequencing datasets or clinical samples.

## Objectives

- Simulate DNA methylation count data for case and control groups.
- Aggregate CpG-level counts into genomic regions.
- Explore statistical methods for differential methylation analysis.
- Investigate sample classification using regional methylation features.
- Benchmark site-level and region-level approaches using simulated datasets.
- Support exploratory analysis of methylation data exported in bedMethyl format.

## Features

| Component | Description |
|---|---|
| Data simulation | Generate methylation counts under configurable coverage and effect sizes. |
| Regional aggregation | Pool methylated counts and total coverage across CpGs within each region. |
| Differential analysis | Explore statistical testing and multiple-testing correction. |
| Classification | Investigate penalized logistic regression for methylation-based classification. |
| Genomic tiling | Assign CpG positions to fixed-width genomic windows. |
| Data import | Read supported modkit bedMethyl output formats. |
| Benchmarking | Compare analytical strategies under controlled simulation scenarios. |

These features are experimental and may evolve as the project develops.

## Methodological approach

### 1. Methylation count aggregation

For each genomic region, methylation is calculated by pooling counts across its CpG sites:

\[
\beta_{\mathrm{region}} =
\frac{\sum_i \mathrm{methylated}_i}
{\sum_i \mathrm{coverage}_i}
\]

This approach accounts for differences in sequencing coverage across CpG sites.

### 2. Differential methylation analysis

The project explores statistical tests for identifying regions with differences between case and control groups, followed by multiple-testing correction using the Benjamini–Hochberg procedure.

### 3. Classification

Penalized logistic regression is explored as a way to evaluate whether regional methylation features can distinguish simulated sample groups.

### 4. Simulation-based benchmarking

Site-level and region-level strategies can be compared under different coverage levels and methylation effect sizes. Benchmark results should be interpreted in the context of the assumptions used to generate the simulated data.

## Installation

Once the repository is available on GitHub and the package metadata has been updated:

r
install.packages("remotes")
remotes::install_github("hassininisrine05/EpiRegion")


Load the package:

r
library(epiregion)


Note: the repository name and the R package name are separate. The package must be named epiregion in its DESCRIPTION file for the installation and loading commands above to work as written.

## Example

Generate a simulated dataset:

r
sim <- simulate_methylation(
  n_cases = 10,
  n_controls = 10,
  n_regions = 20,
  n_signal = 4,
  cpgs_per_region = 5,
  coverage = 30,
  delta = 0.10,
  seed = 1
)


Inspect the simulated data:

r
dim(sim$meth)
head(sim$cpg)
head(sim$pheno)


The example illustrates the simulation workflow. It does not demonstrate performance on real biological samples.

## Limitations

- The current evaluation relies on simulated data.
- The simulation may not capture the full technical and biological complexity of PacBio HiFi methylation measurements.
- Spatial dependence between neighboring CpGs, batch effects and cell-type composition may require additional modelling.
- The performance of regional aggregation depends on the simulation assumptions and requires independent evaluation.
- No clinical episignature or diagnostic application has been validated.
- Support for importing modification calls does not, by itself, establish their biological validity.

EpiRegion is intended for research and methodological exploration, not for clinical diagnosis.

## Future directions

- Explore spatial correlation between neighboring CpGs.
- Evaluate partial and heterogeneous methylation changes within regions.
- Investigate models incorporating relevant covariates.
- Compare fixed-width windows with alternative region definitions.
- Evaluate robustness to coverage variability, noise and batch effects.
- Benchmark methods using suitable public long-read methylation datasets.
- Compare results with established differential methylation methods.

## Technologies

- R
- Statistical modelling
- DNA methylation analysis
- Long-read sequencing data
- Reproducible computational research

## License

MIT License. See the LICENSE file for details.
