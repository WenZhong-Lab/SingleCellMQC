#' @title Show Database Tissue
#' @description This function returns the unique tissue types in the ScType or Cell_Taxonomy database.
#'
#' @param database A character string specifying the database to use. Supported values include:
#' \itemize{
#' \item "ScType": Use the ScType database.
#' \item "Cell_Taxonomy": Use the Cell_Taxonomy database.
#' }
#' @return A character vector of unique tissue types.
#' @export
#' @examples
#' ShowDatabaseTissue("ScType")
#' ShowDatabaseTissue("Cell_Taxonomy")
#'

ShowDatabaseTissue <- function(database="ScType"){
  Tissue_unique <- switch (database,
                           "ScType" = unique(marker_all$ScType$pos$Tissue_standard),
                           "Cell_Taxonomy" = unique(marker_all$Cell_Taxonomy$pos$Tissue_standard),
                           stop("Invalid database")
  )
  return(Tissue_unique)
}


#' Show Common PCT Tissue
#'
#' @param object NULL
#'
#' @return A character vector containing the available tissue types.
#' @export
ShowCommonPCTTissue <- function(object=NULL){
  return(unique(pct_stat_out$tissue))
}



#' @title Show available cell quality control metrics
#'
#' @description This function retrieves and displays the names of available cell-level quality control metrics from a Seurat object.
#'
#' @param object A Seurat object containing the per-cell quality control metrics. The input must be a valid Seurat object to extract the quality control data.
#'
#' @return A character vector containing the names of the available cell quality control metrics in the Seurat object.
#'
#' @export
ShowCellMetricsName <- function(object){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }
  QC_misc <- GetSingleCellMQCData(object)
  out <- colnames(QC_misc$perQCMetrics$perCell)[-1]
  return(out)
}




#' @title Show Sample Metrics
#' @description This function show the available metrics for the samples in the Seurat object.
#'
#' @param object A Seurat object.
#' @param type The type of metrics to show. Supported values include "count", "summary", and "Metrics_10x". Default is "count".
#'
#' @return A character vector containing the available metrics for the samples in the Seurat object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object:
#' # Show the available sample count metrics for the samples in the Seurat object
#' ShowSampleMetrics(seurat_obj, type = "count")
#' }
#'
ShowSampleMetricsName <- function(object, type="count"){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }
  QC_misc <- GetSingleCellMQCData(object)

  out <-  switch (type,
                  "count" = setdiff(colnames(QC_misc$perQCMetrics$perSample$count), "sample"),
                  "summary" = names(QC_misc$perQCMetrics$perSample$summary),
                  "Metrics_10x" = setdiff(colnames(QC_misc$perQCMetrics$perSample$Metrics_10x), "sample"),
                  stop("Invalid `type`, only `count`, `summary`, `Metrics_10x`. ")
  )
  return(out)
}

