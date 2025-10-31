
#' Filter Cells in a Seurat Object Based on Multiple Criteria
#'
#' This function filters cells in a Seurat object based on multiple criteria columns.
#' Each column in the Seurat object's metadata should contain "Pass" or "Fail" values,
#' and optionally NA values. The function can filter cells based on two logical conditions:
#' 1. "and": Filter cells where all specified columns are "Fail".
#' 2. "or": Filter cells where at least one of the specified columns is "Fail".
#'
#' @param object A Seurat object to be filtered.
#' @param filter_columns A character vector of column names in the Seurat object's metadata
#'                      that contain "Pass" or "Fail" values.
#' @param filter_logic A string specifying the filtering logic. Must be either "and" or "or".  Default is "or".
#'                     - `and`: Filter cells where all specified columns are "Fail".
#'                     - `or`: Filter cells where at least one of the specified columns is "Fail".
#' @param split.by (Optional) A column name in the Seurat object's metadata to split the data by.
#'                The function will report the number of cells before and after filtering for each group.
#'                Default is "orig.ident".
#' @param return.table Logical. If TRUE, returns a data.frame summarizing the number of cells
#'                     before and after filtering, including cases with and without split.by.
#'                     If FALSE, returns the filtered Seurat object. Default is FALSE.
#' @param verbose A logical value. If `TRUE` (default), the message will be printed.
#'              If `FALSE` , the message will be suppressed (not printed).
#'
#' @return If `return.table` is TRUE, returns a data.frame summarizing the number of cells
#'         before and after filtering. If `return.table` is FALSE, returns the filtered Seurat object.
#'
#' @export
FilterCells <- function(object,
                        filter_columns = NULL,
                        filter_logic = "or",
                        split.by = "orig.ident",
                        return.table = FALSE,
                        verbose=TRUE){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!all(filter_columns %in% colnames(object@meta.data))){
    stop("Error: The filter_columns must be in the Seurat object's metadata!!")
  }

  if(!filter_logic %in% c("and", "or")){
    stop("Error: The filter_logic must be 'and' or 'or'!!")
  }

  if(!split.by %in% colnames(object@meta.data)){
    stop("Error: The split.by column must be in the Seurat object's metadata!!")
  }


  filter_data <- object@meta.data[, filter_columns, drop = FALSE]
  filter_data[is.na(filter_data)] <- "Pass"

  if(filter_logic == "and"){
    # Filter cells where all specified columns are "Fail"
    out <- !(rowSums(filter_data == "Fail") == length(filter_columns))
  } else {
    # Filter cells where at least one of the specified columns is "Fail"
    out <- !(rowSums(filter_data == "Fail") > 0)
  }

  summary_table <- data.frame(
    Sample = character(),
    Before_Filtering = integer(),
    After_Filtering = integer(),
    stringsAsFactors = FALSE
  )
  # Print filtering stats for the entire object
  .log_cat(paste0("Total cells before filtering:", ncol(object), "\n"), verbose = verbose )
  .log_cat(paste0("Total cells after filtering:", sum(out), "\n"), verbose = verbose )

  summary_table <- rbind(summary_table, data.frame(
    Sample = "Total",
    Before_Filtering = ncol(object),
    After_Filtering = sum(out),
    stringsAsFactors = FALSE
  ))

  if(!is.null(split.by)){
    split_groups <- unique(object@meta.data[[split.by]])
    for(group in split_groups){
      group_cells <- object@meta.data[[split.by]] == group
      .log_cat(paste0(">>>> Sample:", group, "\n"), verbose = verbose )
      .log_cat(paste0("  Cells before filtering:", sum(group_cells), "\n"), verbose = verbose )
      .log_cat(paste0("  Cells after filtering:", sum(group_cells & out), "\n"), verbose = verbose )
      summary_table <- rbind(summary_table, data.frame(
        Sample = group,
        Before_Filtering = sum(group_cells),
        After_Filtering = sum(group_cells & out),
        stringsAsFactors = FALSE
      ))
    }
  }

  if(return.table){
    return(summary_table)
  } else {
    # Return the filtered Seurat object
    object <- object[, out]
    return(object)
  }
}

