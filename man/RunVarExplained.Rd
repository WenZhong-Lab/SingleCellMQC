% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Batch_effect.R, R/Generic.R
\name{RunVarExplained.Seurat}
\alias{RunVarExplained.Seurat}
\alias{RunVarExplained.IterableMatrix}
\alias{RunVarExplained.default}
\alias{RunVarExplained}
\title{RunVarExplained Per-Feature Variance Explained by Variables}
\usage{
\method{RunVarExplained}{Seurat}(object, assay = "RNA", variables = NULL, ...)

\method{RunVarExplained}{IterableMatrix}(object, metadata, variables = NULL, ...)

\method{RunVarExplained}{default}(object, metadata, variables = NULL, ...)

RunVarExplained(object, ...)
}
\arguments{
\item{object}{A Seurat object, a sparse matrix (\code{dgCMatrix}), \code{IterableMatrix}, or a dense matrix (\code{matrix}). If a Seurat object is provided, the function will extract the specified assay data (RNA or ADT) and normalize it if necessary.}

\item{assay}{Character string specifying the assay to use. Options are \code{"RNA"} (default) or \code{"ADT"}. If \code{"RNA"}, the function will normalize the data using \code{Seurat::NormalizeData}. If \code{"ADT"}, the function will normalize the data using CLR normalization.}

\item{variables}{A character vector specifying the variables of interest from the metadata. These variables will be used to explain the variance for each feature.}

\item{...}{Additional arguments passed to the underlying variance calculation function.}

\item{metadata}{A data frame or matrix containing the metadata. If \code{object} is a Seurat object, this parameter is optional, and the metadata will be extracted from \code{object@meta.data}. If \code{object} is a matrix-like, this parameter is required.}
}
\value{
A matrix containing the percentage of variance explained for each feature by the specified variables. The rows correspond to features, and the columns correspond to the variables of interest.
}
\description{
This function calculates the percentage of variance explained by one or more variables of interest for each feature. It supports both RNA and ADT assays and can handle sparse matrices efficiently. The function is designed to work with Seurat objects and can also accept raw matrices or metadata directly.
}
\details{
The function first checks the type of input object (\code{Seurat}, \code{dgCMatrix}, \code{IterableMatrix} or \code{matrix}) and processes it accordingly. For Seurat objects, it extracts the specified assay data and normalizes it if necessary. For sparse or dense matrices, it converts the data into an \code{IterableMatrix} format for efficient computation.
}
