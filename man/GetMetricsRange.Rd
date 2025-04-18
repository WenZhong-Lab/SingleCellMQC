% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotCell.R
\name{GetMetricsRange}
\alias{GetMetricsRange}
\title{Get range of quality control metrics for different tissues}
\usage{
GetMetricsRange(
  tissue,
  metrics = c("nCount_min", "nCount_max", "nFeature_min", "nFeature_max", "nCell_min",
    "MT", "RB", "HB")
)
}
\arguments{
\item{tissue}{Character vector specifying the tissue type(s) to analyze.}

\item{metrics}{Character vector specifying the quality control metrics of interest (\code{"nCount_min"}, \code{"nCount_max"}, \code{"nFeature_min"}, \code{"nFeature_max"}, \code{"nCell_min"}, \code{"MT"}, \code{"RB"}, \code{"HB"}).}
}
\value{
A list of \code{ggplot} objects, each representing a bar plot for the specified metric.
}
\description{
This function retrieves and visualizes the range of specified quality control metrics for different tissue types. It generates bar plots for each metric, showing the distribution of values.
}
