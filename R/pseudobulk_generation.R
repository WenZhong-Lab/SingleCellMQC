#' Create Pseudobulk Data
#'
#' @description
#' This function generates a pseudobulk object consisting of expression matrices,
#' corresponding metadata, and cell count information. It aggregates single-cell
#' expression data into pseudobulk samples based on specified grouping variables.
#'
#' @param object A single-cell object (e.g., Seurat object) from which
#'   expression data and metadata will be extracted.
#' @param assay A character string specifying the assay to use (e.g., "RNA").
#'   Defaults to "RNA".
#' @param slot A character string specifying the slot to use for expression
#'   data (e.g., "counts", "data"). Defaults to "counts".
#' @param sample.by A character string specifying the metadata column representing
#'   the samples (e.g., patient ID, donor ID). This defines the primary unit
#'   for pseudobulk.
#' @param cluster.by (Optional) A character string specifying the metadata column
#'   representing cell clusters or cell types. If provided, pseudobulk matrices
#'   will be generated separately for each unique cluster, resulting in a list
#'   of matrices. Defaults to `NULL`.
#' @param other_cols_contain (Optional) A character vector of additional metadata
#'   column names to include in the pseudobulk metadata. Defaults to `NULL`.
#' @param method A character string indicating the aggregation method for pseudobulking.
#'   Accepted values depend on the underlying `BPCells::pseudobulk_matrix` function,
#'   common options include "sum" (default), "mean".
#'
#' @return A list containing three elements:
#'   \itemize{
#'     \item `pseudobulk_mat`: A matrix or a named list of matrices of pseudobulk
#'       expression values (genes x pseudobulk samples). If `cluster.by` is
#'       specified, it's a list where names are clusters.
#'     \item `pseudobulk_metadata`: A data frame of metadata for the pseudobulk
#'       samples. Rows correspond to pseudobulk samples, and columns include
#'       `sample.by` and any `other_cols_contain`.
#'     \item `ncells`: A data frame summarizing the number of cells contributing
#'       to each pseudobulk sample. Includes columns for `sample.by`, `cluster.by`
#'       (if specified), and "ncells".
#'   }
#' @details
#' This function leverages `BPCells::pseudobulk_matrix` for efficient aggregation.
#' It constructs unique keys for aggregation from `sample.by` and `cluster.by`
#' (if provided). The pseudobulk metadata is constructed by aggregating unique or
#' mean values of specified metadata columns for each pseudobulk sample.
#'
#' @importFrom BPCells pseudobulk_matrix
#' @importFrom stats aggregate
#' @export
RunPseudobulkData <- function(object, assay = "RNA",slot = "counts",
                              sample.by = NULL, cluster.by = NULL,
                              other_cols_contain=NULL, method="sum"){
  mat <- getMatrix(object, assay = assay, slot = slot)
  cell_meta <- getMetaData(object)
  pse_bulk <- .getPseudobulkMatrix(mat, cell_meta, sample.by = sample.by, group.by = cluster.by, method=method)
  pse_meta <- .getPseudobulkMetadata(cell_meta, sample.by = sample.by, variable.by = unique(c(cluster.by, other_cols_contain) ))
  rownames(pse_meta) <- pse_meta[[sample.by]]
  if(is.null(cluster.by)){
    ncells <- .getCounts(cell_meta[[sample.by]])
    colnames(ncells) <- c(sample.by, "ncells")
    ncells[[sample.by]] <- as.character(ncells[[sample.by]])
  }else{
    ncells <- .getCounts(cell_meta[[sample.by]], cell_meta[[cluster.by]])
    colnames(ncells) <- c(sample.by, cluster.by, "ncells")
    ncells[[sample.by]] <- as.character(ncells[[sample.by]])
    ncells[[cluster.by]] <- as.character(ncells[[cluster.by]])
  }
  return( list(pseudobulk_mat = pse_bulk, pseudobulk_metadata = pse_meta, ncells=ncells) )
}


.getCounts <- function(x, y = NULL){
  if (is.null(y) || length(y) == 0) {
    tab <- table(x)
    result <- as.data.frame(tab)
    colnames(result) <- c("x", "Freq")  #
  } else {
    tab <- table(x, y)
    result <- as.data.frame(tab)
    # colnames(result)
  }
  return(result)
}


.getPseudobulkMatrix <- function(object, metadata, sample.by = NULL, group.by = NULL, method="sum"){

  # check IterableMatrix
  if(!is(object, "IterableMatrix")) {
    object <- as(object, "IterableMatrix")
  }

  ## check sample.by & group.by
  if (is.null(sample.by) && is.null(group.by)) {
    stop("Both 'sample.by' and 'group.by' cannot be NULL. Please provide at least one grouping column.")
  }
  cols <- c(group.by, sample.by)
  cols <- cols[!sapply(cols, is.null)]
  if (length(cols) == 1) {
    cell_groups <- metadata[[cols]]
  } else {
    cell_groups <- paste0(metadata[[cols[1]]], "_SCMQC_", metadata[[cols[2]]])
  }

  ## whole pseudobulk_matrix
  mat <- BPCells::pseudobulk_matrix(object, cell_groups = cell_groups, method = method)

  if( !(is.null(sample.by)) && !(is.null(group.by)) ){
    split_vec <- strsplit(colnames(mat), "_SCMQC_")  #
    split_rbind <- do.call(rbind, split_vec)
    split_index_list <- split(1:dim(mat)[2], split_rbind[,1] )

    count_list <- lapply(split_index_list, function(x){
      out <- mat[, x, drop=F]
      colnames(out) <- split_rbind[x,2]
      return(out)
    })
    names(count_list) <- names(split_index_list)
    return(count_list)
  }else{
    return(mat)
  }
}

.getPseudobulkMetadata <- function(object, sample.by, variable.by){
  if(is.null(sample.by)) stop("sample.by cannot be NULL. Please specify a grouping variable.")
  object[[sample.by]] <- as.character(object[[sample.by]])
  variable.by <- setdiff(variable.by, sample.by)
  out_list <- lapply(variable.by, function(x){
    if (!is.numeric(object[[x]])) {
      tab <- stats::aggregate(object[[x]], list(object[[sample.by]]),
                              function(vals) paste(unique(vals), collapse = "|"))
      colnames(tab) <- c(sample.by, x)
    } else {
      # 检查所有组是否都只有一个唯一值
      nuniq <- stats::aggregate(object[[x]], list(object[[sample.by]]), function(vals) length(unique(vals[!is.na(vals)])))
      all_unique <- all(nuniq$x == 1)

      if (all_unique) {
        tab <- stats::aggregate(object[[x]], list(object[[sample.by]]), function(vals) unique(vals[!is.na(vals)]))
        colnames(tab) <- c(sample.by, x)
      } else {
        tab <- stats::aggregate(object[[x]], list(object[[sample.by]]), mean, na.rm = TRUE)
        colnames(tab) <- c(sample.by, paste0("mean_", x))
      }
    }
    tab
  })

  result <- Reduce(function(x, y) merge(x, y, by = sample.by, all = TRUE, sort = FALSE), out_list)
  rownames(result) <- NULL
  return(result)
}
