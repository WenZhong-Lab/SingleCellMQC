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




