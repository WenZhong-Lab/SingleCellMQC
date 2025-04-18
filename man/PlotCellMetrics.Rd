% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotCell.R
\name{PlotCellMetrics}
\alias{PlotCellMetrics}
\title{Plot Cell Metrics}
\usage{
PlotCellMetrics(
  object,
  group.by = "orig.ident",
  metrics,
  color.by = group.by,
  split.by = NULL,
  lwd.box = 0.2,
  color = NULL,
  alpha = 0.6,
  log.y = T,
  combine = F,
  ncol = 1,
  log.mad = T,
  nmads = 3,
  type = "lower",
  guides = NULL,
  linewidth = 0.5,
  add.mad = F,
  linetable = NULL,
  ...
)
}
\arguments{
\item{object}{A Seurat object or a data frame containing metrics.}

\item{group.by}{A string indicating the column name in \code{object} to use for sample identification.}

\item{metrics}{A character vector specifying the metrics to plot. Each metric should correspond to
a column in \code{object}.}

\item{color.by}{(Optional) A string indicating the column name to group data points in the plot.
Defaults to \code{color.by} if not specified.}

\item{split.by}{(Optional) A string indicating the column name to split the data into separate plots.}

\item{lwd.box}{Line width for the box in violin plots. Default is 0.2.}

\item{color}{Color to use for plotting. Can be a single color or a vector of colors corresponding to groups.}

\item{alpha}{Transparency level for the plots. Default is 0.6.}

\item{log.y}{Logical indicating whether to log-transform the Y-axis. Default is TRUE.}

\item{combine}{Logical indicating whether to combine plots into a single plot. Not implemented if TRUE. Default is FALSE.}

\item{ncol}{Number of columns to arrange the combined plots into, if \code{combine} is TRUE. Default is 1.}

\item{log.mad}{Logical indicating whether to log-transform the MAD calculation, applicable only when \code{add.mad} is TRUE. Default is TRUE.}

\item{nmads}{Number of MADs to use for determining outliers, applicable only when \code{add.mad} is TRUE. Default is 3.}

\item{type}{Specifies the type of threshold for outlier detection, applicable only when \code{add.mad} is TRUE. Can be "lower", "upper", or both.}

\item{guides}{A string specifying how guides should be treated in the layout. Details see \code{\link[patchwork]{wrap_plots}} function.}

\item{linewidth}{Line width for the MAD-based outlier highlighting, applicable only when \code{add.mad} is TRUE. Default is 0.5.}

\item{add.mad}{Logical indicating whether to add MAD-based outlier highlighting. Default is TRUE.}

\item{linetable}{An optional data frame containing predefined lines (thresholds) to add to the plots.
Must contain a column matching \code{split.by} and a column with threshold values.}

\item{...}{Additional arguments.}
}
\value{
A list of ggplot objects, each representing a plot for a metric specified in \code{metrics}.
If \code{combine} is TRUE, returns a single ggplot object combining all plots.
}
\description{
This function visualizes various metrics for samples.
It supports plotting violin plots of specified metrics, with options to log-transform the Y-axis,
combine plots, adjust aesthetics, and highlight outliers based on Median Absolute Deviation (MAD).
}
\examples{
\dontrun{
PlotCellMetrics(object = seurat_obj,
                   group.by = "orig.ident",
                   metrics = c("nFeature_RNA", "nCount_RNA"),
                   color = dittoSeq::dittoColors() )
}

}
