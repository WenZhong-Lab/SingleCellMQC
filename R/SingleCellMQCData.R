

#' Retrieve Single Cell Multi-omics QC Data from a Seurat Object
#'
#' Extracts the Single Cell Multi-omics Quality Control (MQC) data stored in the miscellaneous slot of a Seurat object. This function is useful for accessing pre-computed metrics related to quality control of single-cell data.
#'
#' @param object A Seurat object.
#' @return A list containing the SingleCellMQC data stored within the Seurat object.
#'
#' @seealso \code{\link{AddSingleCellMQCData}} for adding Single Cell MQC data into a Seurat object.
#' @examples
#' \dontrun{
#' # Assuming `seuratObj` is a Seurat object with MQC data stored
#' mqcData <- GetSingleCellMQCData(seuratObj)
#'}
#' @export

GetSingleCellMQCData <- function(object){
  out <- Seurat::Misc(object)$SingleCellMQC
  return(out)
}


#' Add Single Cell Multi-omics QC Data to a Seurat Object
#'
#' This function adds or updates metadata in the Seurat object's miscellaneous (`misc`) slot, specifically targeting the Single Cell MQC (Quality Control) data.
#'
#' @param object A Seurat object.
#' @param listname The name of the list within the `SingleCellMQC` slot of the `misc` object where the metadata should be stored.
#' @param metadata A data frame containing the new metadata to add. Each row corresponds to a cell, and columns represent different metadata attributes.The column name in `metadata` used to match rows with cells in the Seurat object.
#' @param sample.by Defaults to "orig.ident".
#'
#' @return The Seurat object with the updated `SingleCellMQC` data in its `misc`.
#'
#' @seealso \code{\link{GetSingleCellMQCData}} for retrieving Single Cell MQC data from a Seurat object.
#' @export
AddSingleCellMQCData <- function(object, listname, metadata, sample.by="orig.ident"){
  # name <- colnames(metadata)
  # object@meta.data[, match(name, colnames(object@meta.data)) ] <- NULL
  object <- SeuratObject::AddMetaData(object = object, metadata = metadata)
  if( is.null(object@misc$SingleCellMQC[[listname]])  ){
    col_name <- colnames(metadata)
  }else{
    col_name <- unique(c(colnames(object@misc$SingleCellMQC[[listname]]),colnames(metadata)) )
  }
  metadata <- object@meta.data[, unique(c(sample.by, col_name))]
  object@misc$SingleCellMQC[[listname]] <-metadata
  return(object)
}


