
#' Plot Feature Metrics
#'
#' This function generates plots for visualizing the distribution of specific feature metrics.
#'
#' @param object A Seurat object processed through the \code{\link{CalculateMetrics}} or \code{\link{CalculateMetricsPerFeature}} function.
#' @param assay A character string specifying the assay to analyze. The default is `"RNA"`.
#' @param metric A character string specifying the feature metric to visualize. Default is "pct".
#'
#' @param sample Optional. Character string specifying a sample name to subset the data. If `NULL`, the entire dataset is used.
#' @param log.x Logical. If `TRUE`, the x-axis will be scaled logarithmically. Default is `TRUE`.
#' @param color Character string specifying the color of the plot (in hex or other color formats). Default is `"#0072B2"`.
#' @param ntop Numeric. Specifies the number of top features to highlight in the bar plot. Default is `10`.
#' @param ylab Character string specifying the y-axis label for the top bar plot. Default is the same as the `metric` parameter.
#'
#' @return The function returns a `ggplot` object.
#'
#' @export
#'
PlotFeatureMetrics <- function(object, assay = "RNA", metric="pct", sample=NULL,
                               log.x=F, color="#0072B2", ntop=10, ylab=metric){
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

  plot_list <- lapply(exp_bind, function(x) {
    if(metric=="pct"){
      data_df <- data.frame(Proportion = x$pct)
      breaks <- c(0,0.01, seq(0.05, 1, by=0.05))
      num_zero_proportion = sum(data_df$Proportion<=0)
      num_zero_proportion_0.1 = sum(data_df$Proportion<=0.001)
      num_all <- dim(data_df)[1]
      p1 <- ggplot2::ggplot(data_df, ggplot2::aes(x = Proportion)) +
        ggplot2::geom_histogram(
          breaks = breaks,
          fill = "#A6CEE3"
        ) +
        ggplot2::geom_text(
          stat = "bin",
          breaks = breaks,
          ggplot2::aes(label = ggplot2::after_stat(count), y = ggplot2::after_stat(count)),
          vjust = -0.5, color = "black", size=3
        ) +
        ggplot2::scale_x_continuous(
          breaks = breaks,
          limits = c(0, 1),
          labels = paste0( c(0,1,seq(5, 100, by=5)) )
        ) +
        ggplot2::labs(subtitle = "Expression proportion distribution (5% Bins)",
             x = "Feature expression proportion (%)",
             y = "Frequency")+
        ggplot2::annotate("text",
                 x = 0.1,
                 y = Inf,
                 label = paste0(num_zero_proportion, "/",num_all,  " features with <= 0% proportion"),
                 hjust = -0.05, vjust = 2.5, #
                 color = "red", size = 4) +
        ggplot2::annotate("text",
                 x = 0.1,
                 y = Inf,
                 label = paste0(num_zero_proportion_0.1,"/",num_all,  " features with <= 0.1% proportion"),
                 hjust = -0.05, vjust = 5, #
                 color = "red", size = 4) +
        ggplot2::theme_classic(base_size = 13)+
        ggplot2::theme(#text = element_text(face = "bold"),
          axis.text = ggplot2::element_text(color = "black")
        )
    }else{
      p1 <- ggplot2::ggplot(x, ggplot2::aes_string(x = metric)) +
        ggplot2::geom_line(stat = "density",
                           alpha = 0.8,
                           linewidth = 1) +
        # scale_color_manual(values =color )+
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          axis.text = ggplot2::element_text(color = "black"),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
        ) +
        ggplot2::labs(x = paste0(assay, " ", metric))
      if (log.x) {
        p1 <- p1 + ggplot2::scale_x_log10()
      }
    }

    p2 <- .gene_top_bar(
      x,
      metric,
      color = color,
      subtitle = "",
      ntop = ntop,
      ylab = ylab
    )
    return(patchwork::wrap_plots(list(p1,p2), widths =c(2,1) ))
  })
  names(plot_list) <- names(exp_bind)
  if (length(plot_list) == 1) {
    plot_list <- plot_list[[1]]
  }

  return(plot_list)
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
    ggplot2::theme_classic(base_size = 14)+
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
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
                                      ggside=T, log.x = F, log.y = F,log.z=F,
                                      raster.cutoff=100000, size =0.3, size.raster= 1, ticks=NULL){

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
                              ggside=ggside, log.x = log.x, log.y = log.y, raster.cutoff=raster.cutoff,
                              size =size, size.raster= size.raster)
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


#' Add feature plot labels
#'
#' This function takes a ggplot object, sorts its data by a specified metric
#' column in descending order, and then adds labels for the top 'ntop'
#' features using `ggrepel::geom_text_repel`. This helps highlight
#' the most prominent features on the plot.
#'
#' @param p A ggplot object. It's expected that the data within the ggplot
#'   object contains a 'Feature' column.
#' @param metric_order_col A string, the name of the column in the plot data
#'   to be used for sorting features (e.g., 'variance'). Features with
#'   higher values in this column will be prioritized.
#' @param ntop An integer, the number of top features to label. Defaults to 10.
#' @return A modified ggplot object with added text labels for the top
#'   'ntop' features.
#' @export
AddFeaturePlotLabel <- function(p, metric_order_col, ntop = 10) {
  if (!inherits(p, "ggplot")) {
    stop("Input 'p' must be a ggplot object.")
  }
  feature_col_name = "Feature"
  plot_data <- p$data
  if (!feature_col_name %in% colnames(plot_data)) {
    stop(paste0("Column '", feature_col_name, "' not found in the plot data."))
  }
  if (!metric_order_col %in% colnames(plot_data)) {
    stop(paste0("Column '", metric_order_col, "' not found in the plot data for sorting."))
  }
  sorted_features <- plot_data %>%
    dplyr::arrange(dplyr::desc(.data[[metric_order_col]])) %>%
    head(ntop)

  label_data <- sorted_features
  p_labeled <- p +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = .data[[feature_col_name]]),color = "red"
    )
  return(p_labeled)
}



#' Plot mean-variance plot by sample
#'
#' For each sample within an assay, this function extracts feature-level
#' mean and variance statistics and plots the mean-variance relationship
#' using a GAM smoother. Optionally, the ggplot object can be converted
#' to an interactive plotly plot.
#'
#' @param object A list object gets from `CalculateMetricsPerFeature` function.
#' @param assay Character string specifying which assay to use. Defaults to \code{"RNA"}.
#' @param mean_col Column name that stores the per-feature means. Defaults to \code{"mean_lognorm"}.
#' @param var_col Column name that stores the per-feature variances. Defaults to \code{"variance_lognorm"}.
#' @param return.type Character string specifying the type of object to return.
#'   Can be \code{"plot"} (default) to return a ggplot object, or \code{"plotly"}
#'   to return \code{plotly::ggplotly()}.
#'
#' @return A ggplot object if \code{return.type = "plot"}, or a plotly object
#'   if \code{return.type = "plotly"}.
#' @export
PlotFeatureMeanVariance <- function(object,
                                    assay = "RNA",
                                    mean_col = "mean_lognorm",
                                    var_col = "variance_lognorm",
                                    return.type = c("plot", "plotly")) {
  return.type <- match.arg(return.type)

  if (!assay %in% names(object)) {
    stop(sprintf("Assay '%s' not found in the provided object.", assay))
  }
  if (!all(c(mean_col, var_col) %in% colnames(object[[assay]][[1]]))) {
    stop(sprintf("Columns '%s' and/or '%s' not found in assay '%s'.",
                 mean_col, var_col, assay))
  }

  mean_var <- lapply(names(object[[assay]]), function(sample_name) {
    df <- object[[assay]][[sample_name]][, c("Feature", mean_col, var_col)]
    df <- data.frame(df, Sample = sample_name)
    colnames(df)[colnames(df) == mean_col] <- "mean"
    colnames(df)[colnames(df) == var_col] <- "variance"
    df
  })
  mean_var <- do.call(rbind, mean_var)

  p <- ggplot2::ggplot(mean_var,
                       ggplot2::aes(x = mean, y = variance, color = Sample)) +
    ggplot2::geom_smooth(method = "gam", se = FALSE, size = 1) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_color_manual(values = get_colors(length(unique(mean_var$Sample))))

  if (return.type == "plot") {
    return(p)
  } else {
    return(plotly::ggplotly(p))
  }
}

