# Multi-Tissue Single-Cell Analysis of Cellular Senescence

This repository contains the source code and computational workflows developed for bioinformatic studies of cellular senescence using single-cell RNA sequencing (scRNA-seq) data.

The repository includes two related research projects addressing complementary aspects of senescence analysis: the characterization of senescent cell heterogeneity across human tissues and the comparative benchmarking of computational methods for senescent cell identification.

To keep the analyses and their respective documentation independent, each project is maintained in a separate Git branch.

## Repository branches

### [`sen_heterogeneity`](https://github.com/estakathersun/sen_heterogeneity_paper/tree/sen_heterogeneity)

**Senescent Cell Heterogeneity: Tissue-Dependent Signatures and Age Dynamics Revealed by Human scRNA-Seq Data**

This branch contains the source code and computational workflows for the study of senescent cell heterogeneity across human tissues using publicly available scRNA-seq datasets.

The analysis focuses on the identification of cell type-specific senescence signatures, characterization of tissue-dependent transcriptional profiles, age-associated changes in senescent cell abundance, and interactions between senescent cells and their tissue microenvironment.

>**Associated publication:**
>
>Matveeva K, Shevyrev D. *Senescent Cell Heterogeneity: Tissue-Dependent Signatures and Age Dynamics Revealed by Human scRNA-Seq Data.* International Journal of Molecular Sciences. 2026;27(18):8416.
>
>[https://doi.org/10.3390/ijms27188416](https://doi.org/10.3390/ijms27188416)

### [`benchmark`](https://github.com/estakathersun/sen_heterogeneity_paper/tree/benchmark)

**Benchmarking Computational Methods for Senescent Cell Identification**

This branch contains the source code and computational workflows for a comparative evaluation of computational methods for identifying senescent cells in scRNA-seq data.

The study benchmarks three algorithms — SenePy, SenCID, and DeepScence — across heterogeneous datasets derived from independent research groups. The analysis evaluates their performance in terms of age-associated senescent cell accumulation, concordance with established senescence gene signatures, inter-method agreement, and robustness to data sparsity, dropout, and subsampling. Different levels of data aggregation are also considered to assess their impact on senescence identification.

### `SCaging` Database

The results of both analyses are integrated into an interactive database, providing access to tissue- and cell type-specific gene expression profiles, age-associated expression patterns, and cellular senescence scores. The database is available at [scaging.online](https://scaging.online).


## Repository structure and usage

Each branch is maintained independently and contains its own `README.md` with a detailed description of the corresponding study, data sources, computational methods, software dependencies, and instructions for reproducing the analyses.

To access the code and documentation for a particular study, switch to the corresponding branch or follow the links provided above.

```bash
# Clone the repository
git clone https://github.com/estakathersun/sen_heterogeneity_paper.git

# Switch to the desired branch
git checkout sen_heterogeneity
# or
git checkout benchmark
```
