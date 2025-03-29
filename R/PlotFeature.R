# PlotFeature -------------------------------------------------------------

#' Plot Feature Metrics
#'
#' This function generates plots for visualizing the distribution of specific feature metrics.
#'
#' @param object A Seurat object processed through the \code{\link{CalculateMetrics}} or \code{\link{CalculateMetricsPerFeature}} function.
#' @param assay A character string specifying the assay to analyze. The default is `"RNA"`.
#' @param type A character string specifying the type of feature metric to visualize.
#' Valid options include:"pct", "mean", "nCell", "variance", "variance.standardized". Default is "pct".
#'
#' @param sample Optional. Character string specifying a sample name to subset the data. If `NULL`, the entire dataset is used.
#' @param log.x Logical. If `TRUE`, the x-axis will be scaled logarithmically. Default is `TRUE`.
#' @param color Character string specifying the color of the plot (in hex or other color formats). Default is `"#0072B2"`.
#' @param ntop Numeric. Specifies the number of top features to highlight in the bar plot. Default is `10`.
#' @param ylab Character string specifying the y-axis label for the top bar plot. Default is the same as the `type` parameter.
#' @param return.type Character string or vector indicating the type of output to return. Options are `"plot"`, `"interactive_table"`, or both. Default is `"plot"`.
#'
#' @return The function returns either a `ggplot` object, an interactive table, or both, depending on the `return.type` specified.
#'   \itemize{
#'     \item If `"plot"` is specified, a `ggplot` object showing the distribution of feature metrics is returned.
#'     \item If `"interactive_table"` is specified, an interactive table showing the feature metrics is returned.
#'     \item If both are specified, a list containing both the plot and the interactive table is returned.
#'   }
#'
#' @export
#'
PlotFeatureMetrics <- function(object, assay = "RNA", type="pct", sample=NULL, log.x=T, color="#0072B2", ntop=10, ylab=type, return.type="plot" ){
  if( "Seurat" %in% is(object)){
    misc_data <- GetSingleCellMQCData(object)$perQCMetrics$perFeature
    if(is.null(misc_data)){
      stop("Error: Please run `CalculateMetricsPerFeature` first. ")
    }else{
      exp_bind <- GetSingleCellMQCData(object)$perQCMetrics$perFeature[[assay]]
    }
  }else{
    exp_bind <- object[[assay]]
  }
  if( length(setdiff(return.type, c("plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` or/and `interactive_table` ")
  }
  if (!(type %in% c("Feature", "nCell", "pct", "mean", "variance", "variance.standardized") )) {
    stop(sprintf("Error: 'type' '%s' is not allowed! Allowed types are: %s",
                 type, paste(c("Feature", "nCell", "pct", "mean", "variance", "variance.standardized"), collapse = ", ")))
  }

  if( !is.null(sample) ){
    exp_bind <- exp_bind[sample]
  }

  out <- list()

  if("interactive_table" %in% return.type){

    out$interactive_table <- lapply(exp_bind, function(x){
      rownames(x) <- NULL
      reactable::reactable(x)
    })
    names(out$interactive_table) <- names(exp_bind)
    if(length(out$interactive_table)==1){
      out$interactive_table <- out$interactive_table[[1]]
    }
  }

  if("plot" %in% return.type){
    plot_list <- lapply(exp_bind, function(x){
      p1 <- ggplot2::ggplot(x, ggplot2::aes_string(x = type)) +
        ggplot2::geom_line(stat = "density",alpha = 0.8, linewidth = 1) +
        # scale_color_manual(values =color )+
        ggplot2::theme_bw(base_size = 14)+
        ggplot2::theme(#text = element_text(face = "bold"),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),

          axis.text = ggplot2::element_text(color = "black"),
          axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,vjust = 1)
        )+
        ggplot2::labs(x=paste0(assay, " ", type)  )
      if(log.x){
        p1 <- p1+ggplot2::scale_x_log10()
      }
      p2 <- .gene_top_bar(x, type, color = color, subtitle ="", ntop=ntop, ylab=ylab)

      return(p1+p2)
    })
    names(plot_list) <- names(exp_bind)
    if(length(plot_list)==1){
      plot_list <- plot_list[[1]]
    }
    out$plot <- plot_list
  }


  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}

.gene_top_bar <- function(object, metrics, color, subtitle=NULL, ntop=10, ylab="% variance explained"){
  object <- object[, c(colnames(object)[1],metrics)]
  top10 <- object[order(-object[[metrics]]),1:2][1:ntop,]
  top10[,1] <- factor(top10[,1], levels = top10[,1][length(top10[,1]):1])
  colnames(top10)[1] <- c("gene")
  p1 <- ggplot2::ggplot(data = top10) +
    ggplot2::geom_col(ggplot2::aes(x=gene, y = .data[[metrics]]),alpha=0.8,fill=color)+
    ggplot2::labs(y = ylab)+
    # geom_text(aes(x = gene, y = .data[[metrics]], label = round(.data[[metrics]],1)),hjust = 1.2,color="black",fontface="bold") +
    ggplot2::coord_flip()+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),

          axis.text = ggplot2::element_text(color = "black"),
          axis.title.y = ggplot2::element_blank()
    )+
    ggplot2::labs(subtitle = subtitle)
  return(p1)
}

#' Plot scatter plot of feature metrics in single cell data
#'
#' This function generates scatter plots to visualize relationships between different feature metrics in single cell data.
#' It can create either 2D scatter plots or 3D scatter plots based on the number of metrics provided.
#'
#' @param object A Seurat object processed through the \code{\link{CalculateMetrics}} or \code{\link{CalculateMetricsPerFeature}} function.
#' @param group.by Character string specifying a column in the metadata to use for grouping. If provided, points in the plot will be colored by the groups specified in this column.
#' @param metrics.by Character vector specifying the metrics to plot. Must contain either 2 or 3 metric names for generating 2D or 3D scatter plots, respectively.
#' @param color Optional. A vector of colors used to distinguish different groups.
#' @param assay Character string specifying which assay to use. Default is `"RNA"`.
#' @param sample Optional. Character string specifying a sample to subset the data. If `NULL`, all samples are used.
#' @param ggside Logical. If `TRUE`, includes marginal density plots on the side of the scatter plot. Default is `TRUE`.
#' @param log.x Logical. If `TRUE`, applies a log10 transformation to the x-axis. Default is `TRUE`.
#' @param log.y Logical. If `TRUE`, applies a log10 transformation to the y-axis. Default is `TRUE`.
#' @param log.z Logical. If `TRUE`, applies a log10 transformation to the z-axis in 3D plots. Default is `TRUE`.
#' @param raster.cutoff Numeric. The number of points at which rasterization is applied to the plot to improve performance. Default is `100000`.
#' @param size Numeric. Point size for non-rasterized points. Default is `0.3`.
#' @param size.raster Numeric. Point size for rasterized points. Default is `1`.
#' @param ticks Numeric vector specifying tick positions for the 3D scatter plot. Optional.
#'
#' @return A ggplot object or list of ggplot objects representing the scatter plot(s).
#' If multiple assays or samples are provided, a list of plots is returned. For a single assay/sample, a single plot is returned.
#'
#' @export
#'
PlotFeatureMetricsScatter <- function(object, group.by=NULL, metrics.by = NULL, color=NULL,assay="RNA",sample=NULL,
                                      ggside=T, log.x = T, log.y = T,log.z=T, raster.cutoff=100000, size =0.3, size.raster= 1, ticks=NULL){

  if( "Seurat" %in% is(object)){
    misc_data <- GetSingleCellMQCData(object)$perQCMetrics$perFeature
    if(is.null(misc_data)){
      stop("Error: Please run `CalculateMetricsPerFeature` first. ")
    }else{
      exp_bind <- GetSingleCellMQCData(object)$perQCMetrics$perFeature[[assay]]
    }
  }else{
    exp_bind <- object[[assay]]
  }
  if( !is.null(sample) ){
    exp_bind <- exp_bind[sample]
  }

  object <- exp_bind
  if(length(metrics.by) <2 & length(metrics.by) >3){
    stop("Error: metrics.by must contain 2 or 3 metrics!")
  }

  plot_list <- lapply(object, function(x){
    if(!is.null(group.by)){
      if(is.null(color)){
        color <- get_colors( length(unique(x[[group.by]])) )
      }
    }else{
      if(is.null(color)){
        color <- get_colors(2)[2:1]
      }
    }

    if(length(metrics.by)==2){
      plot_out <- plotScatter(object=x, x = metrics.by[1] ,y = metrics.by[2], group.by = group.by, color=color,
                              ggside=ggside, log.x = log.x, log.y = log.y, raster.cutoff=raster.cutoff, size =size, size.raster= size.raster)
    }else{
      plot_out <- plotScatter3D(object=x, x = metrics.by[1] ,y = metrics.by[2], z=metrics.by[3], group.by = group.by, color=color,
                                log.x = log.x, log.y = log.y, log.z=log.z, ticks=ticks, text="Feature")
    }

    return(plot_out)
  })
  names(plot_list) <- names(object)
  if(length(plot_list)==1){
    plot_list <-  plot_list[[1]]
  }
  return(plot_list)
}



#' @title Plot cosine similarity between RNA and ADT features in CITE-Seq data
#'
#' @description This function calculates and plots the cosine similarity between specified RNA and ADT features in a Seurat object from CITE-Seq data.
#'
#' @param object A Seurat object containing RNA and ADT data. This object is required for extracting the specified features and calculating cosine similarity.
#' @param RNA_name A character vector of RNA feature names. These features will be compared with ADT features for cosine similarity calculation.
#' @param ADT_name A character vector of ADT feature names corresponding to the RNA features specified in `RNA_name`.
#' @param slot Character string specifying which slot to use for feature data (e.g., `"data"`, `"counts"`, or `"scale.data"`). Default is `"data"`.
#' @param return.type Character string or vector indicating the type of output to return. Options are `"plot"`, `"interactive_table"`, or both. Default is `"plot"`.
#' @param color Character string specifying the color used in the plot. Default is `"#8B658B"`.
#' @param ntop Numeric value indicating the number of top feature pairs to highlight in the bar plot based on their cosine similarity. Default is `15`.
#'
#' @return The function returns either a `ggplot` object, an interactive table, or both, depending on the `return.type` specified.
#'   \itemize{
#'     \item If `"plot"` is specified, a combined `ggplot` object showing a bar plot of the top `ntop` cosine similarities and a histogram of all cosine similarities is returned.
#'     \item If `"interactive_table"` is specified, an interactive table summarizing the cosine similarity between RNA and ADT features is returned.
#'     \item If both are specified, a list containing both the plot and the interactive table is returned.
#'   }
#'
#' @export
#'
PlotFeatureCosineCiteSeq <- function(object, RNA_name=NULL, ADT_name=NULL, slot="data", return.type="plot", color= "#8B658B", ntop=15){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }
  if( length(setdiff(return.type, c("plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` or/and `interactive_table` ")
  }
  ##

  map_id <- data.frame(RNA=RNA_name, ADT= ADT_name)
  if(slot %in% "data"){
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(object, assay = "RNA", slot = "counts")), "package") ) {
      expRNA <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
      expRNA <- BPCells::multiply_cols(expRNA, 1/Matrix::colSums(expRNA))
    }else{
      Seurat::DefaultAssay(object) <- "RNA"
      object <- Seurat::NormalizeData(object)
      expRNA <- Seurat::GetAssayData(object, assay = "RNA", slot = "data")
    }

    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(object, assay = "ADT", slot = "counts")), "package") ) {
      expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = "counts")
      log_sums <- BPCells::colSums(log1p(expADT ), na.rm = TRUE)
      scale_factors <- exp(log_sums / nrow(expADT))
      expADT <- log1p(  BPCells::multiply_cols(expADT ,1/ scale_factors) )
    }else{
      Seurat::DefaultAssay(object) <- "ADT"
      Seurat::VariableFeatures(object) <- rownames(object[["ADT"]])
      object <- Seurat::NormalizeData(object, normalization.method = 'CLR', margin = 2)
      expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = "data")
    }

  }else{
    expRNA <- Seurat::GetAssayData(object, assay = "RNA", slot = slot)
    expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = slot)
  }

  map_id1 <- map_id[!is.na(match(map_id$RNA, rownames(expRNA))),]
  expRNA <- expRNA[map_id1$RNA, ,drop=F]
  expADT <- expADT[map_id1$ADT, ,drop=F]
  cosine_data<-proxyC::simil(expRNA, expADT, method = "cosine")
  cos_data <- data.frame(map_id1, pair=paste0(map_id1$ADT, " - ", map_id1$RNA), cosine_similarity=diag( as.matrix(cosine_data)))
  out <- list()
  if("plot" %in% return.type){
    if(ntop > nrow(cos_data)){
      ntop <- nrow(cos_data)
    }
    p1 <- .gene_top_bar(cos_data, metrics = "cosine_similarity",color = color, ylab = "Cosine similarity", ntop = ntop)
    p2 <- ggplot2::ggplot(data = cos_data) +
      ggplot2::geom_histogram(ggplot2::aes(x = cosine_similarity),
                     fill = color,
                     color = "black") +
      ggplot2::theme_bw(base_size = 14)+ggplot2::labs(x="Cosine similarity")
    out$plot <- p1+p2
  }

  if("interactive_table" %in% return.type){
    interactive_table <- .re_table(cos_data[,c(2,1,3,4)], csv.name = "CosineCiteSeq",elementId = "CosineCiteSeq-table", right.sparkline = F, down.sparkline = F, subtitle = "Cosine similarity", maxWidth = 100)
    out$interactive_table <- interactive_table
  }

  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}


