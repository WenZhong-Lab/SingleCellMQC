
#' @title Generate Statistics for Cell Filtering
#'
#' @description This function generates statistics for cell filtering based on low-quality (lq) or doublet (db) detection.
#' It calculates the number of cells before and after filtering, as well as the percentage of cells filtered out.
#' The results can be returned as a data table or an interactive table.
#'
#' @param object A Seurat object or a data frame containing cell metadata.
#' @param type.detection Type of detection to analyze. Options are "lq" (low-quality cells) or "db" (doublets). Default is "lq".
#' @param return.type Type of output to return. Options are "table" (data frame) or "interactive_table" (interactive HTML table). Default is "table".
#'
#' @return A data frame or an interactive table containing the following columns:
#'   - `Method`: The filtering method used.
#'   - `Pre`: Number of cells before filtering.
#'   - `Post`: Number of cells after filtering.
#'   - `Post%`: Percentage of cells remaining after filtering.
#'   - `Filtered`: Number of cells filtered out.
#'   - `Filtered%`: Percentage of cells filtered out.
#'
#' @export
#'
StatsCellFilter <- function(object,
                            type.detection = c("lq", "db"),
                            return.type = c("table", "interactive_table")) {
  type.detection <- match.arg(type.detection)
  return.type <- match.arg(return.type)
  object <- getMetaData(object)

  pattern <- switch(type.detection,
                    "lq" = "^lq_(?!.*score$)",
                    "db" = "^db_(?!.*score$)")

  grep_col <- grep(pattern, colnames(object), perl = TRUE, value = TRUE)

  if (length(grep_col) == 0) {
    stop(sprintf("No columns matching pattern '%s' found", pattern))
  }

  result <- data.table::rbindlist(lapply(grep_col, function(col) {
    data.table::data.table(
      Method = col,
      Pre = nrow(object),
      Post = sum( !(object[[col]] %in% "Fail"), na.rm = TRUE),
      Filtered = sum(object[[col]] %in% "Fail", na.rm = TRUE)
    )
  }))

  result[, `:=`(
    'Post%' = round(Post / Pre * 100, 3),
    'Filtered%' = round(Filtered / Pre * 100, 3)
  )]

  data.table::setcolorder(result, c("Method", "Pre", "Post", "Post%",
                                    "Filtered", "Filtered%"))

  if (return.type == "interactive_table") {
    return(.re_table(result,
                     maxWidth = NULL,
                     csv.name = paste0("CellPrePostTable_", type.detection),
                     elementId = paste0("CellPrePostTable_", type.detection), subtitle = "Pre and Post Filtering Statistics") )
  }

  return(as.data.frame(result))
}


#' @title Calculate Feature Percentage Across Groups
#'
#' @description This function calculates the percentage of cells expressing a given feature (e.g., a gene) across specified groups.
#' The results can be returned as a data table or an interactive table.
#'
#' @param object A Seurat object containing single-cell data.
#' @param feature A character vector specifying the features (e.g., genes) to analyze. Must be present in the Seurat object.
#' @param group.by A character string specifying the column in the metadata to group cells by. Default is NULL.
#' @param split.by A character string specifying the column in the metadata to split cells by. Default is "orig.ident".
#' @param assay The assay to use for feature extraction. Default is "RNA".
#' @param slot The slot to use for feature extraction. Default is "counts".
#' @param return.type Type of output to return. Options are "table" (data frame) or "interactive_table" (interactive HTML table). Default is "table".
#'
#' @return A data frame or an interactive table containing the following columns:
#'   - `Feature`: The feature (e.g., gene) analyzed.
#'   - `Group`: The group or split category (if specified).
#'   - `PCT`: The percentage of cells expressing the feature in the group.
#'
#' @export
#'
StatsFeaturePCT <- function(object, feature=NULL, group.by = NULL, split.by = "orig.ident", assay = "RNA", slot="counts", return.type="table") {
  if(!is(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  if(is.null(feature)) {
    stop("feature must be specified")
  }
  metadata <- object@meta.data

  object <- Seurat::GetAssayData(object, assay=assay, slot=slot)

  if (!is.null(feature) && length(intersect(feature ,rownames(object)))<length(feature) ) {
    stop("feature not found in object")
  }
  object <- object[rownames(object) %in% feature, , drop = FALSE]


  if (is.null(group.by) && is.null(split.by)) {
    stop("At least one of group.by or split.by must be specified")
  }

  if (!is.null(group.by) && !group.by %in% colnames(metadata)) {
    stop("group.by column not found in metadata")
  }
  if (!is.null(split.by) && !split.by %in% colnames(metadata)) {
    stop("split.by column not found in metadata")
  }

  if (is.null(split.by) && !is.null(group.by)) {
    metadata$group_combo <- as.character(metadata[[group.by]])
  } else if (!is.null(split.by) && is.null(group.by)) {
    metadata$group_combo <- as.character(metadata[[split.by]])
  } else {
    metadata <- dplyr::mutate(metadata,
                              group_combo = paste(!!rlang::sym(group.by), !!rlang::sym(split.by), sep = "++++")
    )
  }
  group_indices <- split(seq_along(colnames(object)), metadata$group_combo)
  result_list <- lapply(names(group_indices), function(group_name) {
    cols <- group_indices[[group_name]]
    proportions <- BPCells::rowSums(object[, cols,drop=F] > 0) / length(cols)*100

    if (is.null(split.by)) {
      data.table::data.table(row_id = rownames(object), group = group_name, proportion = proportions)
    } else if (is.null(group.by)) {
      data.table::data.table(split = group_name, row_id = rownames(object), PCT = proportions)
    } else {
      data.table::data.table(row_id = rownames(object), group_combo = group_name, proportion = proportions)
    }
  })

  result_dt <- data.table::rbindlist(result_list)

  if (!is.null(split.by) && !is.null(group.by)) {
    result_dt[, c(group.by, split.by) := data.table::tstrsplit(group_combo, "\\+\\+\\+\\+")]
    result_dt[, group_combo := NULL]
    result_wide <- data.table::dcast(result_dt, get(split.by) + row_id ~ get(group.by),
                                     value.var = "proportion")
  } else if (is.null(split.by)) {
    result_wide <- data.table::dcast(result_dt, row_id ~ group, value.var = "proportion")
  } else {
    setnames(result_dt, "split", split.by)
    result_wide <- result_dt[, .(get(split.by), row_id, PCT)]
    setnames(result_wide, c("V1", "row_id", "PCT"), c(split.by, "row_name", "PCT"))
    return(result_wide[])
  }

  data.table::setnames(result_wide, "row_id", "Feature")

  out <- list()
  if (return.type == "table") {
    out$table <- result_wide[]
  }

  if (return.type == "interactive_table") {
    out$interactive_table <- .re_table(result_wide,
                                       maxWidth = NULL,
                                       csv.name = paste0("FeaturePCTTable_",group.by,"_",split.by),
                                       elementId = paste0("FeaturePCTTable_",group.by,"_",split.by), subtitle = "Feature PCT Table")
  }

  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}






