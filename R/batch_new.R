




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




RunBatchGTEs <- function(object, assay = "RNA", slot = "data", group.by=NULL, batch.by, block.size = 500) {
  mat <- getMatrix(object,assay=assay,  slot=slot)
  mat_meta <- getMetaData(object)
  colnames(mat) <- NULL
  row_index_list <- split(1:nrow(mat), ceiling(seq_along(1:nrow(mat)) / block.size))
  out_list <- lapply(row_index_list, function(idx) {
    mat1 <- as(mat[idx, , drop=FALSE], "dgCMatrix")
    out <- suppressMessages(GTEs::Run.GroupTechEffects(mat1, mat_meta, g_factor = group.by, b_factor = batch.by))
    return(out)
  })
  names(out_list) <- NULL
  out <- do.call(c, lapply(out_list, function(x) x$OverallTechEffects))
  out <- data.frame(Feature=names(out), GTE=out)
  colnames(out)[2] <- batch.by
  return(out)
}







#' @title Visualize Reduced Dimensional Data
#' @description This function generates scatter plots of reduced dimensional data, either at the cell level or pseudobulk level.
#'
#' @param object A Seurat object or a list of pseudobulk data.
#' @param group.by A character string specifying the metadata column used for grouping cells (default: `"seurat_clusters"`).
#' @param reduction A character string specifying the type of reduction to plot (default: `"rna.umap"`).
#' @param split.by A character string specifying the metadata column used to split the data into subsets (default: `NULL`).
#' @param ncol An integer specifying the number of columns for the facet plot (default: `NULL`).
#' @param ggside Logical, whether to use ggside for facet plots (default: `FALSE`).
#' @param color A vector of colors to use for visualizations (default: `NULL`).
#' @param guide.nrow An integer specifying the number of rows for the legend (default: `10`).
#' @param raster.cutoff An integer specifying the cutoff for rasterization (default: `100000`).
#'
#' @return A plot based on the specified `plot.type`.
#' @export
#'
PlotReducedDim <- function(object,
                           group.by = "seurat_clusters",
                           process.pseudobulk="pca",
                           reduction="rna.umap",
                           split.by= NULL,
                           ncol=NULL,
                           ggside=F,
                           color = NULL,
                           guide.nrow=10,
                           raster.cutoff=100000,
                           size=0.3){
  if( !("Seurat" %in% is(object)) & process.pseudobulk=="pca" ){
    mat <- object$pseudobulk_mat
    meta_data <- object$pseudobulk_metadata
    if(is.list(mat)){
      group_name <- names(mat)
      plot_list <- lapply(group_name, function(x){
        message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- runPseudobulkPCA: ", x))
        out <- runPseudobulkPCA(mat[[x]], metadata = meta_data)
        p <- plotReducedDim(object=out, group.by = group.by, reduction = "pca" , split.by = split.by, ncol = ncol, ggside = ggside,
                             color = color, guide.nrow = guide.nrow, raster.cutoff = raster.cutoff, size = size)
        return(p)
      })
      names(plot_list) <- group_name
      return(plot_list)
    }else{
      out <- runPseudobulkPCA(mat, metadata = meta_data)
      p <- plotReducedDim(object=out, group.by = group.by, reduction = "pca" , split.by = split.by, ncol = ncol, ggside = ggside,
                           color = color, guide.nrow = guide.nrow, raster.cutoff = raster.cutoff, size = size)
      return(p)
    }

  }else{
    p <- plotReducedDim(object=object, group.by = group.by, reduction = reduction , split.by = split.by, ncol = ncol, ggside = ggside,
                         color = color, guide.nrow = guide.nrow, raster.cutoff = raster.cutoff, size = size)
    return(p)
  }
}


plotReducedDim <- function(object,
                            group.by = "seurat_clusters",
                            reduction="rna.umap",
                            split.by= NULL,
                            ncol=NULL,
                            ggside=F,
                            color = NULL,
                            guide.nrow=10,
                            raster.cutoff=100000,
                            size=0.3
){
  cluster_data <- SeuratObject::Embeddings(object = object, reduction= reduction)
  index <- union(group.by, split.by)
  cluster_data <- data.frame(cluster_data, object@meta.data[, index, drop=F])
  if( is.null(color)){
    color <- get_colors( length(unique(cluster_data[[group.by]])))
  }
  sample_num <- length(unique(object@meta.data[,split.by]))
  if(is.null(ncol)){
    ncol = ceiling(sqrt(sample_num))
  }
  cluster_data[[group.by]] <- factor(cluster_data[[group.by]], levels = sort(unique(cluster_data[[group.by]])) )

  p <- plotScatter(cluster_data, x= colnames(cluster_data)[1] , y= colnames(cluster_data)[2], group.by = group.by, color = color,size=size,
                   log.x = F,log.y = F,ggside = ggside,split.by = split.by, ncol = ncol,guide.nrow=guide.nrow, raster.cutoff=raster.cutoff)

  p[["facet"]][["params"]][["ncol"]] <- ncol
  p[["facet"]][["params"]][["nrow"]] <- ceiling(sample_num / ncol)
  return(p)

}


runPseudobulkPCA <- function(object, metadata){
  suppressWarnings(seu <- Seurat::CreateSeuratObject(counts = object, meta.data = metadata))
  seu <- Seurat::NormalizeData(object = seu, verbose =F)
  seu <- Seurat::FindVariableFeatures(object = seu, verbose =F)
  seu <- Seurat::ScaleData(object = seu, verbose =F)
  seu <- Seurat::RunPCA(object = seu, verbose =F)
  #seu <- Seurat::RunUMAP(object = seu, verbose =F)
  return(seu)

}







PlotGTEBar <- function(object, color=NULL){
  object <- data.table::data.table(object)
  rsquared_long <- data.table::melt(object,id.vars="Feature")
  rsquared_long <- rsquared_long[order(rsquared_long$value, decreasing = TRUE)]
  rsquared_long$Feature <- factor(rsquared_long$Feature, levels =  rev(unique(rsquared_long$Feature))  )
  ggplot2::ggplot(data = rsquared_long) +
    ggplot2::geom_col(ggplot2::aes(x=Feature, y = .data[["value"]]),alpha=0.7,fill=color)+
    ggplot2::coord_flip()+
    ggplot2::theme_classic(base_size = 13)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text = ggplot2::element_text(color = "black"),
                   axis.title.y = ggplot2::element_blank()
    )+
    # ggplot2::theme(panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank())+
    ggplot2::labs( y = "GTE score")
}



