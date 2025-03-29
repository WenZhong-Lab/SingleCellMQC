**SingleCellMQC**: A comprehensive quality control workflow for single-cell multi-omics in human health and disease monitoring

## Installation

You can install the development version of `SingleCellMQC` from
[GitHub](https://github.com/DaihanJi/SingleCellMQC) with:

```R
install.packages('devtools')
devtools::install_github('DaihanJi/SingleCellMQC')
library(SingleCellMQC)
```

The SingleCellMQC pipeline is developed in R and accepts input files from scRNA-seq, surface protein seqencing, scTCR-seq, and scBCR-seq to perform QC analysis. The pipeline utilizes the Seurat and BPCells objects to store and analyze large data. The pipeline consists of four major QC modules: (i) sample QC, (ii) cell QC, (iii) feature QC, and (iv) batch QC. Each module encompasses a range of functions, including QC metrics assessment, outlier sample detection, abnormal cell identification, background noise detection, batch effect evaluation, etc. Data visualization for each module encompasses a variety of elements, including static graphs, interactive graphs, and interactive tables. In particular, each module supports the generation of HTML reports for multi-omics quality control (QC) based on R Markdown, which contains QC results and potential warnings.

