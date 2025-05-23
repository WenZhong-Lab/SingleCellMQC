% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Read_data.R
\name{Read10XData}
\alias{Read10XData}
\title{Load multi-omics data from 10x}
\usage{
Read10XData(
  dir_GEX = NULL,
  dir_TCR = NULL,
  dir_BCR = NULL,
  sample = NULL,
  saveBPCells = F,
  dir_BPCells = "./BPCellData",
  gene.column = 2,
  return.type = "Seurat"
)
}
\arguments{
\item{dir_GEX}{Directory for Gene/Protein Expression data provided by 10X. Directory containing the matrix.mtx, features.tsv, and barcodes.tsv files provided by 10X. e.g. "/data/SingleCellMQC/CellRanger/TP1/sample_filtered_feature_bc_matrix/".
A vector or named vector can be given in order to load several data directories.}

\item{dir_TCR}{Directory for TCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.}

\item{dir_BCR}{Directory for BCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.}

\item{sample}{Optional; sample identifiers. If not provided, default names will be generated.}

\item{saveBPCells}{Logical; whether to save the GEX matrix using BPCells for more efficient data handling. Defaults to \code{FALSE}.}

\item{dir_BPCells}{Directory to save BPCells data. Used if \code{saveBPCells} is set to \code{TRUE}. Defaults to \code{./BPCellData}.}

\item{gene.column}{Integer; the column of the features.tsv file that contains gene names. Defaults to \code{2}.}

\item{return.type}{Character; specifies the format of the returned GEX data. Options are \code{"Seurat"} (returns a Seurat object) or \code{"matrix"} (returns a dgCMatrix or IterableMatrix ).
Defaults to \code{"Seurat"}.}
}
\value{
A list containing the following components:
\itemize{
\item \code{GEX}: Gene/Protein Expression data. If \code{return.type} is \code{"Seurat"}, this is a merged Seurat object; if \code{"matrix"}, it is a merged matrix.
\item \code{TCR}: T-cell receptor data, returned as a list of data frames.
\item \code{BCR}: B-cell receptor data, returned as a list of data frames.
}

If only one component is loaded (e.g., only GEX data), the function returns that component directly instead of a list.
}
\description{
Read gene/protein expression (GEX), T-cell receptor (TCR), and B-cell receptor (BCR) data provided by 10X genomics.
}
