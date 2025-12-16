
#' Quantify Feature-level Technical Effects (GTE) for Single-cell Data
#'
#' @description
#' This function calculates feature-level Technical Effects (GTE) scores for
#' individual genes (or other features) in single-cell expression data,
#' measuring the influence of user-defined technical factors (batch variables).
#' It is an adapted and extended implementation of the methodology originally
#' presented in the `GTEs::Run.GroupTechEffects` function, designed to offer
#' broader compatibility with various single-cell data formats, such as
#' `BPCells` matrices, `Seurat` object.
#' The core GTE principle is based on the work by Zhou et al. (2025).
#'
#' @param object A single-cell object.
#' @param assay A character string specifying the assay to use (e.g., "RNA")
#'   if `object` is a single-cell object. Defaults to "RNA".
#' @param slot A character string specifying the slot to use for expression
#'   data (e.g., "data", "counts") if `object` is a single-cell object.
#'   Defaults to "data".
#' @param group.by (Optional) A character string specifying the metadata column
#'   representing biological groups (e.g., "cell_type", "condition") to be
#'   controlled for in the GTE calculation. This helps isolate batch effects
#'   from true biological variation.
#' @param batch.by A character string specifying the metadata column representing
#'   the technical factor (e.g., "batch_id", "sequencing_run") whose influence
#'   is to be quantified.
#'
#' @return A `data.frame`
#'   \itemize{
#'     \item \code{Feature}: Gene or protein names.
#'     \item \code{GTE_Overall}: The overall GTE score for each gene, quantifying
#'       the total influence of the `batch.by` variable.
#'     \item Additional columns: For each subgroup defined by `group.by`, if
#'       provided, a column (e.g., "CellTypeA", "CellTypeB") will represent
#'       the GTE score specific to that subgroup.
#'   }
#'
#' @seealso \code{\link[GTEs]{Run.GroupTechEffects}}
#' @references
#' Zhou Y, Sheng Q, Wang G, Xu L, Jin S. Quantifying batch effects
#' for individual genes in single-cell data. Nat Comput Sci. 2025 Aug;5(8):612-660.
#' doi: 10.1038/s43588-025-00824-7. Epub 2025 Jun 27. PMID: 40579473.
#'
#' @export
RunBatchGTE <- function(object, assay = "RNA", slot = "data", group.by, batch.by) {
  mat <- getMatrix(object,assay=assay,  slot=slot)
  mat_meta <- getMetaData(object)

  #name
  g_name = unique(mat_meta[[group.by]])
  b_name = unique(mat_meta[[batch.by]])

  #cellcount
  cellcount <- as.data.frame(table(mat_meta[[group.by]], mat_meta[[batch.by]]))
  colnames(cellcount) <- c(group.by, batch.by, "count" )

  # var
  gb <- .getPseudobulkMatrix(mat, metadata = mat_meta, sample.by=batch.by, group.by = group.by, method = "variance")
  g <- .getPseudobulkMatrix(mat, metadata = mat_meta, sample.by=group.by, method = "variance")

  # GroupTechEffects
  GroupTechEffects <- lapply(g_name, function(x){
    cellcount_temp = cellcount[ cellcount[, group.by] %in% x, , drop=F]
    g_count = sum(cellcount_temp[,"count" ])
    ssg = g[, match(x, colnames(g))] * (g_count-1)
    index = match( colnames(gb[[x]]) , cellcount_temp[, batch.by] )
    ssgb = t(gb[[x]]) * (cellcount_temp[index, "count"]-1)
    ssgb = colSums(ssgb)
    return((ssg - ssgb)/g_count)
  })
  names(GroupTechEffects) <- g_name
  GroupTechEffects <- do.call(cbind, GroupTechEffects)
  scale_factor = length(which(colSums(GroupTechEffects) > 0))
  GroupTechEffects = GroupTechEffects / scale_factor
  GTE_Overall <- rowSums(GroupTechEffects)
  out <- data.frame(Feature=rownames(GroupTechEffects), GTE_Overall=GTE_Overall,GroupTechEffects, check.names = F)
  return(out)
}



#' Visualize Reduced Dimensional Data
#'
#' @description
#' This function generates scatter plots of reduced dimensional data, either
#' at the cell level (from a single-cell object) or at the pseudobulk PCA level.
#'
#' @param object A Seurat object directly, or a list/object (typically output
#'   from `RunPseudobulkData`)
#'   for pseudobulk-level PCA plotting.
#' @param group.by A character string specifying the metadata column used
#'   for coloring points (e.g., "seurat_clusters").
#' @param process.pseudobulk A character string. If "pca" and `object` is
#'   a pseudobulk data list, it will first run PCA on the pseudobulk matrix
#'   and then plot its PCA embeddings. If `object` is a Seurat object, this
#'   parameter is ignored. Defaults to "pca".
#' @param reduction A character string specifying the name of the dimensionality
#'   reduction to plot (e.g., "umap", "pca"). For single-cell Seurat object,
#'   this typically refers to existing reductions. For pseudobulk data,
#'   if `process.pseudobulk="pca"`, this will automatically be "pca".
#'   Defaults to "rna.umap".
#' @param split.by (Optional) A character string specifying a metadata column
#'   by which to facet the plots. Defaults to `NULL`.
#' @param ncol (Optional) An integer specifying the number of columns for the
#'   facet plot. If `NULL`, it will be automatically determined.
#' @param ggside A logical value. If `TRUE`, it attempts to add marginal plots
#'   using `ggside` (if installed and loaded). Defaults to `FALSE`.
#' @param color (Optional) A character vector of colors to use for `group.by`.
#'   If `NULL`, default colors are assigned.
#' @param guide.nrow An integer specifying the number of rows for the legend guide.
#'   Defaults to 10.
#' @param raster.cutoff An integer specifying the number of points above which
#'   data will be rasterized to improve plot rendering performance. Points below
#'   this threshold will be plotted as individual points. Defaults to 100000.
#' @param size A numeric value specifying the size of the plotted points.
#'   Defaults to 0.3.
#' @param pca.label A character string specifying the metadata column used for label.
#'
#' @return A `ggplot` object, or a named list of `ggplot` objects if `object$pseudobulk_mat`
#'   is a list (e.g., pseudobulk PCA plots per cluster).
#'
#'
#' @export
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
                           pca.label=NULL,
                           size=0.3){
  if( !("Seurat" %in% is(object)) & process.pseudobulk=="pca" ){
    mat <- object$pseudobulk_mat
    meta_data <- object$pseudobulk_metadata
    assay=object$assay
    if(is.list(mat)){
      group_name <- names(mat)
      plot_list <- lapply(group_name, function(x){
        message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- runPseudobulkPCA: ", x))
        out <- runPseudobulkPCA(mat[[x]], metadata = meta_data)
        p <- plotReducedDim(object=out, group.by = group.by, reduction = "pca" , split.by = split.by, ncol = ncol, ggside = ggside,
                             color = color, guide.nrow = guide.nrow, raster.cutoff = raster.cutoff, size = size,
                            label =pca.label)
        return(p)
      })
      names(plot_list) <- group_name
      return(plot_list)
    }else{
      out <- runPseudobulkPCA(mat, metadata = meta_data, assay=assay)
      p <- plotReducedDim(object=out, group.by = group.by, reduction = "pca" ,
                          split.by = split.by,
                          ncol = ncol,
                          ggside = ggside,
                           color = color,
                          guide.nrow = guide.nrow,
                          raster.cutoff = raster.cutoff,
                          size = size,
                          label =pca.label)
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
                           label=NULL,
                            size=0.3
){
  cluster_data <- SeuratObject::Embeddings(object = object, reduction= reduction)
  index <- union( union(group.by, split.by), label)
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
                   log.x = F,log.y = F,ggside = ggside,split.by = split.by, ncol = ncol,guide.nrow=guide.nrow,
                   raster.cutoff=raster.cutoff, label=label)

  p[["facet"]][["params"]][["ncol"]] <- ncol
  p[["facet"]][["params"]][["nrow"]] <- ceiling(sample_num / ncol)
  return(p)

}


runPseudobulkPCA <- function(object, metadata, assay){
  suppressWarnings(seu <- Seurat::CreateSeuratObject(counts = object, meta.data = metadata))
  if(assay=="RNA"){
    seu <- Seurat::NormalizeData(object = seu, verbose =F)
  }else if(assay=="ADT"){
    seu <- Seurat::NormalizeData(object = seu, verbose =F, normalization.method = 'CLR', margin = 2)
  } else {
    stop("Unsupported assay type provided. Only 'RNA' and 'ADT' are supported.")
  }
  seu <- Seurat::FindVariableFeatures(object = seu, verbose =F)
  seu <- Seurat::ScaleData(object = seu, verbose =F)
  npcs = min(20, nrow(metadata)-1)
  seu <- Seurat::RunPCA(object = seu, verbose =F, npcs= npcs)
  #seu <- Seurat::RunUMAP(object = seu, verbose =F)
  return(seu)

}




#' Plot GTE Scores as a Bar Plot
#'
#' @description
#' This function generates a bar plot to visualize the top genes ranked by
#' their GTE scores, indicating their association with
#' a batch effect.
#'
#' @param object A data frame, typically the output from `RunBatchGTEs`. It
#'   should contain at least two columns: 'Feature' (gene names) and a column
#'   representing the GTE scores (usually named after the `batch.by` variable
#'   used in `RunBatchGTE`).
#' @param color (Optional) A character string specifying the color for the bars.
#'   If `NULL`, a default color is used.
#' @param ntop An integer specifying the number of top features (genes) to display
#'   in the plot. If `NULL` or less than 1, all features will be plotted.
#'   Defaults to 10.
#' @param column_name A character string specifying the name of the column in
#'   `object` that contains the GTE scores (or desired metric) to be plotted.
#'   This column must be numeric. Default is "GTE_Overall".
#'
#' @return A `ggplot` object representing the bar plot of GTE scores.
#'
#' @details
#' The function sorts the features by their GTE scores in descending order and
#' selects the top `ntop` features. It then creates a horizontal bar plot where
#' each bar's length corresponds to the GTE score of a feature.
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot geom_col aes labs theme element_text element_blank coord_flip
#' @export
PlotGTEBar <- function(object, column_name = "GTE_Overall", color = NULL, ntop = 10){
  if(is.null(color)){
    color = get_colors(2)[2]
  }
  object <- object[,c("Feature", column_name) ]
  object <- data.table::data.table(object)
  rsquared_long <- data.table::melt(object, id.vars = "Feature")
  rsquared_long <- rsquared_long[order(rsquared_long$value, decreasing = TRUE)]

  if (!is.null(ntop) && is.numeric(ntop) && ntop > 0) {
    if (ntop > nrow(rsquared_long)) {
      warning("ntop (", ntop, ") is greater than the number of available features (", nrow(rsquared_long), "). Plotting all available features.")
    }
    rsquared_long <- utils::head(rsquared_long, ntop)
  }

  rsquared_long$Feature <- factor(rsquared_long$Feature, levels =  rev(unique(rsquared_long$Feature)))

  p <- ggplot2::ggplot(data = rsquared_long) +
    ggplot2::geom_col(ggplot2::aes(x = Feature, y = .data[["value"]]), alpha = 0.7, fill = color) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
                   axis.text = ggplot2::element_text(color = "black"),
                   axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = "GTE score", subtitle = column_name)

  return(p)
}

#' @title Run Batch Test
#' @description This function performs a batch test on a Seurat object to evaluate the presence of batch effects. It calculates the Kendall W statistic, pairwise Kendall correlations, and the number of negative correlations between samples.
#' The function also checks for severe batch effects based on the results.
#'
#' @param object A Seurat object containing the single-cell data to be analyzed.
#' @param sample.by A character string specifying the metadata column used for defining samples (default: `"orig.ident"`).
#' @param cluster.by A character string specifying the metadata column used for grouping cells (default: `"seurat_clusters"`).
#' @param k An integer specifying the number of clusters to use for batch testing (default: `3`).
#' @param seed An integer specifying the seed for random sampling (default: `1`).
#' @param n An integer specifying the number of samples to use for batch testing (default: `100`).
#' @return A list containing the correlation matrix, Kendall W results, and the qbinom results.
#' @export
RunBatchTest <- function(object, sample.by="orig.ident", cluster.by="seurat_clusters", k=3, seed=1, n=100 ){
  test_data <- data.table(object@meta.data)[, .(count = .N), by = .(sample = get(sample.by), cluster = get(cluster.by))]
  test_data <- dcast(test_data, cluster ~ sample, value.var = "count", fill = 0)
  test_data <- as.data.frame(test_data)
  rownames(test_data) <- test_data[,1]
  test_data <- test_data[,-1]
  cor_data <- stats::cor(test_data,method ="kendall")
  dist_matrix <- stats::as.dist(1-cor_data)
  hc <- stats::hclust(dist_matrix)
  clusters <- stats::cutree(hc, k = k)

  k_list <- lapply(1:k, function(x){
    test_data[,colnames(test_data) %in% names(which(clusters==x)), drop=F]
  })

  sampled_list <- lapply(1:k, function(x){
    set.seed(seed)
    sampled_chars <- sample(colnames(k_list[[x]]), size = n, replace = TRUE)
    test_data[, match(sampled_chars, colnames(test_data)), drop=F]
  })

  kendallW_result <- lapply(1:n, function(x){
    merge_data <- lapply(sampled_list, function(y){
      y[, x,drop=F]
    })
    merge_data <- do.call(cbind, merge_data)
    merge_data <- merge_data[rowSums(merge_data) > 0, ]
    kendallW<-irr::kendall(merge_data)

    cor_data <- apply(utils::combn(dim(merge_data)[2], 2), 2, function(x){
      cor_temp <- merge_data[, x]
      cor_temp <- cor_temp[rowSums(cor_temp) > 0, ]
      stats::cor( as.numeric(cor_temp[,1]),  as.numeric(cor_temp[,2]), method ="kendall")
    })
    c(kendallW$value, kendallW$p.value, sum(cor_data<0) )
  })
  kendallW_result <- do.call(rbind, kendallW_result)
  colnames(kendallW_result) <- c("W", "pvalue", "InterNegNum")
  kendallW_result <- as.data.frame(kendallW_result)
  kendallW_result$Warning <- unlist(apply(kendallW_result, 1, function(x){
    (x[2] > 0.1) & (x[3] >= 1)
  }))

  if(sum(kendallW_result$Warning) >= stats::qbinom(0.05,n, 0.5)){
    message("The dataset may have severe clustering differences, which may be batch effect !!!")

  }
  return(list(cor_matrix=cor_data, kendallW_result=kendallW_result, qbinom_result=data.frame(sum=sum(kendallW_result$Warning), ntimes=n, qbinom_min= stats::qbinom(0.05,n, 0.5) ) ))
}










