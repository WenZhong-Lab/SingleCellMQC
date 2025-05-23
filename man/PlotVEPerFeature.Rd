% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature_QC.R
\name{PlotVEPerFeature}
\alias{PlotVEPerFeature}
\title{Plot variance explained per feature for specified variables}
\usage{
PlotVEPerFeature(
  object,
  assay = "RNA",
  variables = NULL,
  plot.type = "density",
  color.density = NULL,
  color.bar = "#56B4E9",
  ntop = 10,
  return.type = "plot",
  csv.name = paste0(assay, "_variance_explained")
)
}
\arguments{
\item{object}{A Seurat object or a data frame. If a Seurat object is provided, the variance explained for each feature will be computed using the \code{RunVarExplained} function. If a data frame is provided, it should contain precomputed variance explained values.}

\item{assay}{Character string specifying which assay to use. Options include \code{"RNA"} (default) or other assays present in the Seurat object.}

\item{variables}{Character vector specifying the variables of interest from the metadata (e.g., \code{"condition"}, \code{"batch"}). These variables will be used to explain the variance for each feature.}

\item{plot.type}{Character string specifying the type of plot to generate. Options are \code{"density"} to visualize the distribution of variance explained or \code{"bar"} to visualize the top features. Default is \code{"density"}.}

\item{color.density}{Optional vector of colors for density plots. If \code{NULL}, default colors are generated.}

\item{color.bar}{Character string specifying the color of bars in the bar plot. Default is \code{"#56B4E9"}.}

\item{ntop}{Numeric value indicating the number of top features to highlight in the bar plot based on their variance explained. Default is \code{10}.}

\item{return.type}{Character string or vector indicating the type of output to return. Options are \code{"plot"}, \code{"interactive_table"}, or both. Default is \code{"plot"}.}

\item{csv.name}{Character. Name for interactive_table csv file.}
}
\value{
The function returns either a \code{ggplot} object, an interactive table, or both, depending on the \code{return.type} specified.
\itemize{
\item If \code{"plot"} is specified, a density plot or bar plot (or both) showing the variance explained per feature is returned.
\item If \code{"interactive_table"} is specified, an interactive table summarizing the variance explained per feature is returned.
\item If both are specified, a list containing both the plot and the interactive table is returned.
}
}
\description{
This function visualizes the percentage of variance explained per feature (e.g., gene) for specified variables in a Seurat object.
It supports creating density plots and bar plots to represent the distribution or highlight key features with high variance explained by the selected variables.
}
