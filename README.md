
## SingleCellMQC

**SingleCellMQC**: A comprehensive quality control workflow for single-cell multi-omics in human health and disease monitoring. The SingleCellMQC pipeline is developed in R and accepts input files from scRNA-seq, surface protein seqencing, scTCR-seq, and scBCR-seq to perform QC analysis. The pipeline utilizes the **Seurat** and **BPCells** objects to store and analyze large data. The pipeline consists of four major QC modules:

## R Installation

Before installation, we recommend installing these dependencies first:

-   [Seurat](https://github.com/satijalab/seurat) (single-cell analysis toolkit)
-   [BPCells](https://github.com/bnprks/BPCells) (memory-efficient single-cell data processing)

You can install the development version of `SingleCellMQC` :
```         
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")  
}
devtools::install_github('WenZhong-Lab/SingleCellMQC') 
library(SingleCellMQC)
```
## Note
Version 0.7.0 is scheduled for release in September. This will be an entirely new version, featuring numerous optimizations, bug fixes, and added functions. The logical error present in the RunLQ_MAD function of Version 0.6 will also be fixed in 0.7.0. We highly recommend using our future new version, as Version 0.6 will subsequently be deprecated.


### Integrated QC Tools in SingleCellMQC

[scater](https://bioconductor.org/packages/release/bioc/html/scater.html), [Seurat](https://satijalab.org/seurat/), [miQC](https://github.com/greenelab/miQC), [ddqc](https://github.com/ayshwaryas/ddqc), [decontX](https://github.com/campbio/decontX), [scCDC](https://github.com/ZJU-UoE-CCW-LAB/scCDC), [scDblFinder](https://github.com/plger/scDblFinder), [scds](https://github.com/kostkalab/scds), [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder), [DropletUtils](https://github.com/MarioniLab/DropletUtils), [scRepertoire](https://github.com/BorchLab/scRepertoire), [GTEs](https://github.com/yzhou1999/GTEs/).

**Other Key Packages For QC**: [BPCells](https://github.com/bnprks/BPCells), [COSG](https://github.com/genecell/COSGR), [ScType](https://github.com/IanevskiAleksandr/sc-type#readme), [variancePartition](https://github.com/GabrielHoffman/variancePartition), [dreamlet](https://github.com/GabrielHoffman/dreamlet), [reactable](https://glin.github.io/reactable/index.html)



