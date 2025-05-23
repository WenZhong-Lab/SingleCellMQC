% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generic.R, R/Read_data.R
\name{AddSampleMeta}
\alias{AddSampleMeta}
\alias{AddSampleMeta.Seurat}
\alias{AddSampleMeta.default}
\title{Add Sample metadata to Object}
\usage{
AddSampleMeta(
  object,
  merge_by_object = "orig.ident",
  SampleMeta,
  merge_by_meta,
  ...
)

\method{AddSampleMeta}{Seurat}(
  object,
  merge_by_object = "orig.ident",
  SampleMeta,
  merge_by_meta,
  ...
)

\method{AddSampleMeta}{default}(
  object,
  merge_by_object = "orig.ident",
  SampleMeta,
  merge_by_meta,
  ...
)
}
\arguments{
\item{object}{An object}

\item{merge_by_object}{Object column name to merge by}

\item{SampleMeta}{Data Frame containing Sample metadata}

\item{merge_by_meta}{Sample metadata column name to merge by}

\item{...}{Additional arguments passed to methods}
}
\value{
Object with Sample metadata added
}
\description{
Add Sample metadata to an object
}
