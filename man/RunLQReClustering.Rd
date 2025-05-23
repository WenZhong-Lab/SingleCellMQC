% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cell_QC.R
\name{RunLQReClustering}
\alias{RunLQReClustering}
\title{LQ cells re-clustering}
\usage{
RunLQReClustering(
  object,
  method_columns = NULL,
  cluster_resolution = 0.5,
  preprocess = "rna.umap"
)
}
\arguments{
\item{object}{A Seurat object. The input must be a valid Seurat object.}

\item{method_columns}{Character. The name of the metadata column used for specifying the method.}

\item{cluster_resolution}{Numeric. The resolution used for clustering the low-quality cells. Default is 0.5.}

\item{preprocess}{Character. The name of the preprocessing method to use. Default is "rna.umap".}
}
\value{
A Seurat object containing only the low-quality cells that were re-clustered, with updated cluster assignments.
}
\description{
This function performs re-clustering on low quality cells.
Cells marked as 'Fail' in any of the LQ columns are extracted and preprocessed.
}
