

#' Plot Cell Metrics
#'
#' This function visualizes various metrics for samples.
#' It supports plotting violin plots of specified metrics, with options to log-transform the Y-axis,
#' combine plots, adjust aesthetics, and highlight outliers based on Median Absolute Deviation (MAD).
#'
#' @param object A Seurat object or a data frame containing metrics.
#' @param group.by A string indicating the column name in `object` to use for sample identification.
#' @param metrics A character vector specifying the metrics to plot. Each metric should correspond to
#'                a column in `object`.
#' @param color.by (Optional) A string indicating the column name to group data points in the plot.
#'                 Defaults to `color.by` if not specified.
#' @param split.by (Optional) A string indicating the column name to split the data into separate plots.
#' @param lwd.box Line width for the box in violin plots. Default is 0.2.
#' @param color Color to use for plotting. Can be a single color or a vector of colors corresponding to groups.
#' @param alpha Transparency level for the plots. Default is 0.6.
#' @param log.y Logical indicating whether to log-transform the Y-axis. Default is TRUE.
#' @param combine Logical indicating whether to combine plots into a single plot. Not implemented if TRUE. Default is FALSE.
#' @param ncol Number of columns to arrange the combined plots into, if `combine` is TRUE. Default is 1.
#' @param add.mad Logical indicating whether to add MAD-based outlier highlighting. Default is TRUE.
#' @param log.mad Logical indicating whether to log-transform the MAD calculation, applicable only when `add.mad` is TRUE. Default is TRUE.
#' @param nmads Number of MADs to use for determining outliers, applicable only when `add.mad` is TRUE. Default is 3.
#' @param type Specifies the type of threshold for outlier detection, applicable only when `add.mad` is TRUE. Can be "lower", "upper", or both.
#' @param guides A string specifying how guides should be treated in the layout. Details see \code{\link[patchwork]{wrap_plots}} function.
#' @param linewidth Line width for the MAD-based outlier highlighting, applicable only when `add.mad` is TRUE. Default is 0.5.
#' @param linetable An optional data frame containing predefined lines (thresholds) to add to the plots.
#'                  Must contain a column matching `split.by` and a column with threshold values.
#' @param ... Additional arguments.
#'
#' @return A list of ggplot objects, each representing a plot for a metric specified in `metrics`.
#'         If `combine` is TRUE, returns a single ggplot object combining all plots.
#'
#' @examples
#' \dontrun{
#' PlotCellMetrics(object = seurat_obj,
#'                    group.by = "orig.ident",
#'                    metrics = c("nFeature_RNA", "nCount_RNA"),
#'                    color = dittoSeq::dittoColors() )
#' }
#'
#' @export
#'

PlotCellMetrics <-
  function(object,
           group.by="orig.ident",
           metrics,
           color.by =group.by,
           split.by = NULL,
           lwd.box = 0.2,
           color=NULL,
           alpha = 0.6,
           log.y = T,
           combine = F,
           ncol = 1,
           log.mad = T,
           nmads = 3,
           type = "lower",
           guides = NULL,
           linewidth=0.5,
           add.mad=F,
           linetable=NULL,
           ...
  ){
    if("Seurat" %in% class(object)){
      object <- object@meta.data
    }else{
      object <- object
    }

    if(is.null(color)) {
      color <- get_colors(length(unique(object[, color.by])))
    }

    plot_list <- plotVln(object=object, x=group.by, y=metrics, log.y = log.y,combine = F,
                         group.by=color.by, split.by = split.by, color = color, lwd.box=lwd.box, guides = guides,...)
    list_name <- names(plot_list)
    if(add.mad){
      stat_list <- suppressMessages(CalculateMetricsPerSample.summary(object, sample.by = group.by,metrics = list_name))
      if(is.null(split.by)){
        object$temp <- "temp"
        split.by<- "temp"
      }

      plot_list <- mapply( function(h, data,i){
        h$group <- object[,match(split.by, colnames(object))][match(h[,"sample"], object[,group.by])]
        mad_result <- RunOutlier(h, nmads = nmads, log = log.mad,metrics = "Median", group.by = "group",method="mad",...)
        mad_coutoff <- attr(mad_result, "thresholds")
        mad_coutoff <- data.frame(rownames(mad_coutoff), mad_coutoff, check.names = F)
        colnames(mad_coutoff)[1] <- split.by
        if(length(type)==1){
          p1 <- data + ggplot2::geom_hline(ggplot2::aes(yintercept = .data[[type]]) , mad_coutoff,
                                  linetype = "dashed", color = "red", linewidth=linewidth)
        }else{
          p1 <- data + ggplot2::geom_hline(ggplot2::aes(yintercept = .data[[type[i]]]) , mad_coutoff,
                                  linetype = "dashed", color = "red", linewidth=linewidth)
        }
        return(p1)

      },h=stat_list, data=plot_list, i = 1:length(type) ,SIMPLIFY = F)
    }else{
      if(!is.null(linetable)){
        plot_list <- lapply(plot_list, function(pl){
          pl+ ggplot2::geom_hline(ggplot2::aes(yintercept = .data[[ colnames(linetable)[2] ]]) , linetable,
                         linetype = "dashed", color = "red", linewidth=linewidth)
        })
      }
    }
    if(length(plot_list)==1){
      plot_list <-  plot_list[[1]]
    }
    return(plot_list)

  }


plotVln <-
  function(object,
           x,
           y,
           group.by,
           split.by = NULL,
           lwd.box = 0.2,
           color,
           alpha = 0.8,
           log.y = T,
           combine = F,
           ncol = 1,
           guides = NULL,
           facet.cols=3,
           ...) {
    plot_list <- mapply(function(z,i){

      p1 <- suppressMessages(ggplot2::ggplot(object, mapping = ggplot2::aes( x=.data[[x]], y= .data[[z]])) +
                               ggplot2::geom_violin(scale="width", ggplot2::aes(fill=.data[[group.by]] ),
                                           color=NA,
                                           position = ggplot2::position_dodge(1),
                                           alpha=alpha)+
                               ggplot2::scale_fill_manual(values=color)+
                               ggplot2::scale_colour_manual(values =color)+
                               ggplot2::geom_boxplot(width=0.3,position = ggplot2::position_dodge(1),
                                                     ggplot2::aes(color=.data[[group.by]]),fill='white',outlier.size=0.5,
                                            alpha=1,lwd=lwd.box)+
                               ggplot2::scale_colour_manual(values=color)+
                               ggplot2::theme_classic(base_size = 13) +
                               ggplot2::theme(
                                 axis.line = ggplot2::element_line(color = "black", linewidth = 0.8),
                                 strip.text.x = ggplot2::element_text(size = 13),
                                 axis.text = ggplot2::element_text(color = "black", size = 13),
                                 axis.text.x = ggplot2::element_text(
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1
                                 ),
                                 ##angle
                                 axis.title.x = ggplot2::element_blank(),
                                 # axis.ticks.x = element_blank(),
                                 # legend.position = "bottom"
                               ))

      if(!is.null(split.by)){
        # p1 <- p1 + facet_grid(~get(split.by), scales = "free_x", space = "free_x",cols = facet.cols)
        p1 <- p1 + ggplot2::facet_wrap(~get(split.by), scales = "free_x", ncol  = facet.cols)

      }

      if( length(log.y)==1 ){
        if(log.y){
          p1<-p1+ggplot2::scale_y_log10()
        }
      }else{
        if(log.y[i]){
          p1<-p1+ggplot2::scale_y_log10()
        }
      }

      return(p1)
    }, z = y,i= 1: length(y) , SIMPLIFY = F)

    if(combine){
      plot_list <- patchwork::wrap_plots(plot_list, ncol = ncol, guides = guides)
    }
    return(plot_list)
  }


#' @title Plot cell filtration percentage results from quality control methods
#'
#' @description This function visualizes the percentage of cells filtered out based on quality control methods in a Seurat object.
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the metadata containing quality control information will be extracted for analysis.
#' @param type.detection Character string indicating the type of filtration method. Options are `"lq"` for low-quality cells or `"db"` for doublet cells. Default is `"lq"`.
#' @param sample.by Character string specifying the column name in the metadata used for grouping cells (e.g., `"orig.ident"` for different samples or conditions).
#' @param color Character string specifying the color used in the bar plot for representing filtered cells percentage. Default is `'#A6CEE3'`.
#' @param return.type Character string or vector indicating the type of output to return. Options are `"plot"`, `"interactive_table"`, or both. Default is `"plot"`.
#' @param color.shadow Character string specifying the shadow color for the bar plot text. Default is `"white"`.
#' @param color.text Character string specifying the color of the text in the bar plot. Default is `"black"`.
#'
#' @return The function returns either a `ggplot` object representing the bar plot or an interactive table, depending on the `return.type` specified.
#'   \itemize{
#'     \item If `"plot"` is specified, a bar plot showing the percentage of cells filtered out by the selected method for each sample is returned.
#'     \item If `"interactive_table"` is specified, an interactive table showing the count and percentage of filtered cells per sample is returned.
#'     \item If both are specified, a list containing both the plot and the interactive table is returned.
#'   }
#'
#' @export
PlotCellMethodFiltration <- function(object, type.detection ="lq", ##"lq", "db"
                                     sample.by="orig.ident",color='#A6CEE3', return.type="plot",
                                     color.shadow="white", color.text="black"){
  #
  if( length(setdiff(return.type, c("plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` or/and `interactive_table` ")
  }

  if( "Seurat" %in% is(object)){
    object <- object@meta.data
  }

  metaname <- colnames(object)

  switch(type.detection,
         "lq" = {
           grep_col <- grep("^lq_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           subtitle =  "Low quality cell filtered"
         },
         "db" = {
           grep_col <- grep("^db_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           subtitle =  "Doublet filtered"
         },
         stop("Invalid `type.detection`, only: `lq` or/and `db` ")
  )

  if(length(grep_col)==0){
    stop("No quality control metrics found for the specified type.detection!")
  }

  plot_table <- object[, match(c(sample.by, grep_col), metaname)]
  colnames(plot_table)[-1] <- stringr::str_sub(colnames(plot_table)[-1], 4, -1)
  out <- list()
  plot_count <- data.table(plot_table)[, lapply(.SD, function(x) sum(x=="Fail",na.rm = T) ), by=list(sample = get(sample.by)),
                                       .SDcols=setdiff(colnames(plot_table), sample.by) ]
  plot_pct <- data.table(plot_table)[, lapply(.SD, function(x) (sum(x=="Fail",na.rm = T)/length(x)) ), by=list(sample = get(sample.by)),
                                     .SDcols=setdiff(colnames(plot_table), sample.by)]

  if("interactive_table" %in% return.type){
    out$pct <- .re_table(plot_pct,csv.name = paste0(type.detection,"_pct.csv"),elementId = paste0(type.detection, "_pct-table"), right.sparkline = T , down.sparkline = T,first_name = "Sample", subtitle = paste0(subtitle, " percentage"), number_type = "custom", number_fmr = "paste0( round(value*100,2),'%')")
    out$count <- .re_table(plot_count,csv.name = paste0(type.detection,"_count.csv"),elementId = paste0(type.detection, "_count-table"),right.sparkline = T , down.sparkline = T,first_name = "Sample", subtitle = paste0(subtitle, " number"))
    return(out)
  }
  if("plot" %in% return.type){
    plot_pct <- melt(plot_pct, id.vars="sample" )
    plot_text <- data.table(plot_table)[, lapply(.SD, function(x){
      paste0(round((sum(x=="Fail",na.rm = T)/length(x))*100,1),"%","(",sum(x=="Fail",na.rm = T),")")
    }), by=list(sample = get(sample.by)), .SDcols=setdiff(colnames(plot_table), sample.by)]

    out$plot <- plotBar(plot_pct, sample_col = "sample",split.by = "variable",value_cols = "value",color.shadow = color.shadow, color.text = color.text,
                        data_text = melt(plot_text,id.vars = "sample"),
                        WidetoLong = F, reverse=T, color =  rep(color,  length(unique(plot_table[sample.by,])) *dim(plot_table)[2])  ,xlim = T)+
      ggplot2::theme(legend.position = "none")+ggplot2::labs(subtitle = paste0(subtitle, " percentage"))
  }
  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}

#' @title Plot overlap of filtered cells using upset plot
#'
#' @description This function creates an UpSet plot to visualize the overlap of cells filtered out by different quality control methods (e.g., low quality or doublet cells).
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the metadata containing quality control information will be extracted for analysis.
#' @param type.detection Character string indicating the type of filtration method. Options are `"lq"` for low-quality cells or `"db"` for doublet cells. Default is `"lq"`.
#' @param split.by Optional. Character string specifying the column in the metadata to split the data by (e.g., `"orig.ident"` to visualize separately for each sample). If `NULL`, all cells are considered together. Default is `NULL`.
#' @param nintersects Numeric value indicating the maximum number of intersections to display in the UpSet plot. Default is `20`.
#'
#' @return A list of UpSet plots if `split.by` is specified; otherwise, a single UpSet plot showing the overlap of filtered cells across different quality control metrics.
#'
#' @export
PlotCellMethodUpset <- function(object, type.detection="lq", split.by=NULL, nintersects=20){
  if( "Seurat" %in% is(object)){
    object <- object@meta.data
  }

  if( !is.null(split.by) ){
    object <- split(object, object[[split.by]])
  }else{
    object <- list(all=object)
    split.by="orig.ident"
  }


  plot_list <- lapply(object, function(x){
    metaname <- colnames(x)
    switch(type.detection,
           "lq" = {
             grep_col <- grep("^lq_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           },
           "db" = {
             grep_col <- grep("^db_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           },
           stop("Invalid `type.detection`, only: `lq` or/and `db` ")
    )
    plot_table <- x[, match(c(split.by, grep_col), metaname)]
    colnames(plot_table)[-1] <- stringr::str_sub(colnames(plot_table)[-1], 4, -1)
    data_binary <- as.data.frame(lapply(plot_table[,2:(length(colnames(plot_table))), drop=F], function(x) as.integer(x %in% 'Fail')))
    plot_upset <- UpSetR::upset(data_binary,sets = colnames(data_binary), order.by = "freq",show.numbers ="yes",nintersects=nintersects, mb.ratio = c(0.6, 0.4))
    return(plot_upset)
  })
  names(plot_list) <- names(object)
  if(length(plot_list)==1){
    plot_list <-  plot_list[[1]]
  }
  return(plot_list)
}

#' @title Plot quality control metrics with violin plots for filtered cells
#'
#' @description This function generates violin plots to visualize the distribution of selected quality control metrics for cells filtered by different quality control methods (e.g., low quality or doublet cells).
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the metadata containing quality control information will be extracted for visualization.
#' @param type.detection Character string indicating the type of filtration method. Options are `"lq"` for low-quality cells or `"db"` for doublet cells. Default is `"lq"`.
#' @param split.by Optional. Character string specifying the column in the metadata to split the data by (e.g., `"orig.ident"` to visualize separately for each sample). If `NULL`, all cells are considered together. Default is `NULL`.
#' @param metrics.by Character string specifying the quality control metric to plot. Default is `"nFeature_RNA"`.
#' @param color Vector of colors to use for plotting the different cell types (`"Fail"` vs. `"Pass"`). Default uses two colors.
#'
#' @return A list of `ggplot` violin plot objects if `split.by` is specified; otherwise, a single `ggplot` object visualizing the chosen quality control metric.
#'
#' @export
PlotCellMethodVln <- function(object, type.detection="lq", split.by=NULL, metrics.by ="nFeature_RNA", color = get_colors(2)[2:1] ){
  if( "Seurat" %in% is(object)){
    object <- object@meta.data
  }

  if( !is.null(split.by) ){
    object <- split(object, object[[split.by]])
  }else{
    object <- list(all=object)
  }

  plot_list <- lapply(object, function(x){
    metaname <- colnames(x)
    switch(type.detection,
           "lq" = {
             grep_col <- grep("^lq_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           },
           "db" = {
             grep_col <- grep("^db_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
           },
           stop("Invalid `type.detection`, only: `lq` or/and `db` ")
    )
    plot_table <- x[, match(c("orig.ident", grep_col), metaname)]
    colnames(plot_table)[-1] <- stringr::str_sub(colnames(plot_table)[-1], 4, -1)

    plot_table <- plot_table[match(rownames(x), rownames(plot_table)),]
    vln_metrics <- melt(data.table(plot_table), id.var= "orig.ident", value.name="Type")
    vln_metrics[[metrics.by]] <- rep(x[[metrics.by]], dim(plot_table)[2]-1 )

    vln_metrics <- na.omit(vln_metrics)
    plot_out <- plotVln(vln_metrics, x = "variable", y=metrics.by, group.by = "Type", color = color,
                        facet.cols=1 )[[1]]
    return(plot_out)
  })
  names(plot_list) <- names(object)
  if(length(plot_list)==1){
    plot_list <-  plot_list[[1]]
  }
  return(plot_list)
}


plotScatter <- function(object, x, y, group.by = NULL, color = NULL, ggside = TRUE,
                        log.x = TRUE, log.y = TRUE, split.by = NULL, ncol = 2, guide.nrow = 10,
                        raster.cutoff = 100000, size = 0.3, size.raster = 1,
                        label = NULL, size.label = 3) {
  # --- Color handling ---
  if (is.null(color)) {
    if (!is.null(group.by)) {
      color <- get_colors(length(unique(object[[group.by]])))
    } else {
      color <- get_colors(2)[2] # Default color if no grouping
    }
  }
  # --- geom_func (point layer) handling ---
  # Use scattermore for large datasets for performance
  if (length(object[[x]]) > raster.cutoff) {
    if (!requireNamespace("scattermore", quietly = TRUE)) {
      warning("Package 'scattermore' not found. Please install it for rasterized plotting of large datasets. Falling back to ggplot2::geom_point.")
      geom_func <- ggplot2::geom_point(size = size, alpha = 1)
    } else {
      geom_func <- scattermore::geom_scattermore(pointsize = size.raster, alpha = 1)
    }
  } else {
    if (!is.null(group.by)) {
      geom_func <- ggplot2::geom_point(size = size, alpha = 1)
    } else {
      geom_func <- ggplot2::geom_point(size = size, alpha = 1, color = color[1])
    }
  }
  # --- Create base ggplot object ---
  if (!is.null(group.by)) {
    p1 <- ggplot2::ggplot(data = object,
                          mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[group.by]]))
  } else {
    p1 <- ggplot2::ggplot(data = object,
                          mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  }
  # --- Add point layer ---
  p1 <- p1 + geom_func
  # --- Conditionally add label layer ---
  if (!is.null(label)) {
    # Ensure the label column exists in the data
    if (!label %in% colnames(object)) {
      warning(paste0("Label column '", label, "' not found in the data object. Skipping text labels."))
    } else {
      # Map the label column to the 'label' aesthetic
      p1 <- p1 +
        ggplot2::aes(label = .data[[label]]) +
        ggrepel::geom_text_repel(size = size.label, show.legend = FALSE,
                                 max.overlaps = Inf) # show.legend=F to prevent labels appearing in legend
    }
  }
  # --- Other ggplot layers and theme ---
  p1 <- p1 +
    # scientific_theme(base_size = 15)+ # Uncomment if you have this theme defined
    ggplot2::guides(color = ggplot2::guide_legend(nrow = guide.nrow, override.aes = list(size = 4))) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::scale_colour_manual(values = color)
  if (log.x) {
    p1 <- p1 + ggplot2::scale_x_log10()
  }
  if (log.y) {
    p1 <- p1 + ggplot2::scale_y_log10()
  }
  # --- ggside extension ---
  if (ggside) {
    if (!requireNamespace("ggside", quietly = TRUE)) {
      warning("Package 'ggside' not found. Skipping side density plots.")
    } else {
      p1 <- p1 +
        ggside::geom_xsidedensity(show.legend = FALSE, linewidth = 0.6) +
        ggside::scale_xsidey_continuous(expand = c(0, 0), labels = c(NULL, NULL, NULL, NULL)) + # Adjust side axis labels
        ggside::geom_ysidedensity(show.legend = FALSE, linewidth = 0.6) +
        ggside::scale_ysidex_continuous(expand = c(0, 0), labels = c(NULL, NULL, NULL, NULL)) + # Adjust side axis labels
        # ggside::theme_ggside_void() # Use this for a completely blank side theme
        ggplot2::theme(ggside.panel.scale = 0.25, ggside.axis.line = ggplot2::element_blank(),
                       ggside.axis.ticks = ggplot2::element_blank())
    }
  }
  # --- Faceting (facet_wrap) ---
  if (!is.null(split.by)) {
    p1 <- p1 + ggplot2::facet_wrap(~ get(split.by), ncol = ncol)
  }
  return(p1)
}

# plotScatter <- function(object, x, y, group.by=NULL, color=NULL, ggside=T,
#                         log.x=T, log.y=T, split.by =NULL, ncol=2,guide.nrow=10,
#                         raster.cutoff=100000, size =0.3, size.raster= 1){
#
#   if(is.null(color)){
#     if(!is.null(group.by)){
#       color <- get_colors( length(unique(object[[group.by]])) )
#       }else{
#       color <- get_colors(2)[2]
#       }
#   }
#
#   if (length(object[[x]]) > raster.cutoff) {
#     geom_func <- scattermore::geom_scattermore(pointsize = size.raster, alpha = 1)
#   } else {
#     if( !is.null(group.by) ){
#       geom_func <- ggplot2::geom_point(size = size, alpha = 1)
#     }else{
#       geom_func <- ggplot2::geom_point(size = size, alpha = 1, color=color[1])
#     }
#   }
#
#   if( !is.null(group.by) ){
#     p1 <- ggplot2::ggplot(data = object,
#                  mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], color =.data[[group.by]]))
#
#   }else{
#     p1 <- ggplot2::ggplot(data = object,
#                  mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]] ))
#   }
#
#   p1 <- p1 +
#     geom_func +
#     # scientific_theme(base_size = 15)+
#     ggplot2::guides(color=ggplot2::guide_legend(nrow = guide.nrow, override.aes = list(size=4)))+
#     ggplot2::theme_classic(base_size = 15)+
#     ggplot2::scale_colour_manual(values = color)
#
#   if(log.x){
#     p1 <- p1 + ggplot2::scale_x_log10()
#   }
#   if(log.y){
#     p1 <- p1 + ggplot2::scale_y_log10()
#   }
#   if(ggside){
#     p1 <- p1 +
#       ggside::geom_xsidedensity(show.legend	=F, linewidth=0.6)+
#       ggside::scale_xsidey_continuous(expand = c(0, 0),labels = c(NULL,NULL,NULL,NULL))+
#       ggside::geom_ysidedensity(show.legend	=F, linewidth=0.6)+
#       ggside::scale_ysidex_continuous(expand = c(0, 0),labels = c(NULL,NULL,NULL,NULL))+
#       # ggside::theme_ggside_void()
#       ggplot2::theme(ggside.panel.scale=0.25,ggside.axis.line=ggplot2::element_blank(),
#             ggside.axis.ticks=ggplot2::element_blank() )
#   }
#   if(!is.null(split.by)){
#     p1 <- p1 + ggplot2::facet_wrap(~ get(split.by), ncol = ncol)
#   }
#   return(p1)
# }

plotScatter2 <- function(object, x, y, group.by=NULL, color=NULL, ggside=T,
                        log.x=T, log.y=T, split.by =NULL, ncol=2,guide.nrow=10,
                        raster.cutoff=100000, size =0.3, size.raster= 1){

  if(is.null(color)){
    if(!is.null(group.by)){
      color <- get_colors( length(unique(object[[group.by]])) )
    }else{
      color <- get_colors(2)[2]
    }
  }

  if (length(object[[x]]) > raster.cutoff) {
    geom_func <- scattermore::geom_scattermore(pointsize = size.raster, alpha = 1)
  } else {
    if( !is.null(group.by) ){
      geom_func <- ggplot2::geom_point(ggplot2::aes(color =.data[[group.by]]),size = size, alpha = 1)
    }else{
      geom_func <- ggplot2::geom_point(size = size, alpha = 1, color=color[1])
    }
  }

  if( !is.null(group.by) ){
    p1 <- ggplot2::ggplot(data = object,
                          mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]]))

  }else{
    p1 <- ggplot2::ggplot(data = object,
                          mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]] ))
  }

  p1 <- p1 +
    geom_func +
    # scientific_theme(base_size = 15)+
    ggplot2::guides(color=ggplot2::guide_legend(nrow = guide.nrow, override.aes = list(size=4)))+
    ggplot2::theme_classic(base_size = 15)+
    ggplot2::scale_colour_manual(values = color)

  if(log.x){
    p1 <- p1 + ggplot2::scale_x_log10()
  }
  if(log.y){
    p1 <- p1 + ggplot2::scale_y_log10()
  }
  if(ggside){
    p1 <- p1 +
      ggside::geom_xsidedensity(show.legend	=F, linewidth=0.6)+
      ggside::scale_xsidey_continuous(expand = c(0, 0),labels = c(NULL,NULL,NULL,NULL))+
      ggside::geom_ysidedensity(show.legend	=F, linewidth=0.6)+
      ggside::scale_ysidex_continuous(expand = c(0, 0),labels = c(NULL,NULL,NULL,NULL))+
      ggplot2::theme(ggside.panel.scale=0.25,ggside.axis.line=ggplot2::element_blank(),
                     ggside.axis.ticks=ggplot2::element_blank() )
  }
  if(!is.null(split.by)){
    p1 <- p1 + ggplot2::facet_wrap(~ get(split.by), ncol = ncol)
  }
  return(p1)
}


plotScatter3D <- function(object, x, y, z , log.x=T, log.y=T, log.z=T, color=NULL,group.by=NULL,
                          ticks=NULL, text=NULL){
  automated_round <- function(x) {
    power <- floor(log10(x))
    base <- 10^power
    rounded_num <- base * floor(x / base)
    rounded_num
  }
  log_ticks <- function(data,length.out=5) {
    log_range <- log10(range(data, na.rm = TRUE))
    ticks <- 10^seq(from = ifelse(log_range[1]<0,0,log_range[1]) , to = log_range[2], length.out = length.out)
    ticks <- sapply(ticks[-1], automated_round)
    return(ticks)
  }

  if(!is.null(text)){
    p1 <- plotly::plot_ly(object, x = ~ get(x), y = ~ get(y), z = ~ get(z),  type = "scatter3d",mode="markers",
                          size = I(3), colors = color,text=~get(text) )
  }else{
    if(!is.null(group.by)){
      p1 <- plotly::plot_ly(object, x = ~ get(x), y = ~ get(y), z = ~ get(z), color = ~ get(group.by), type = "scatter3d",mode="markers",
                            size = I(3), colors = color )
    }else{
      p1 <- plotly::plot_ly(object, x = ~ get(x), y = ~ get(y), z = ~ get(z), type = "scatter3d",mode="markers",
                            size = I(3), colors = color )
    }

  }

  scene <- list()
  scene$xaxis$title <- x
  scene$yaxis$title <- y
  scene$zaxis$title <- z
  if(log.x){
    scene$xaxis$type <- "log"
    if( !is.null(ticks) ){
      scene$xaxis$tickvals <- log_ticks(object[[x]],length.out=ticks)
    }
    scene$xaxis$tickfont = list(size = 13)
  }
  if(log.y){
    scene$yaxis$type <- "log"
    if( !is.null(ticks) ){
      scene$yaxis$tickvals <- log_ticks(object[[y]],length.out=ticks)
    }
    scene$yaxis$tickfont = list(size = 13)

  }
  if(log.z){
    scene$zaxis$type <- "log"
    if( !is.null(ticks) ){
      scene$zaxis$tickvals <- log_ticks(object[[z]],length.out=ticks)
    }
    scene$zaxis$tickfont = list(size = 13)

  }
  scene$aspectratio = list(x = 1, y = 1, z = 1)

  p1 <- plotly::layout(p1, scene= scene)
  return(p1)
}


#' @title Plot cell metrics with scatter plots
#'
#' @description This function generates scatter plots to visualize the relationships between selected quality control metrics for cells in a Seurat object.
#'  The function supports both 2D and 3D scatter plots.
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the metadata containing the cell metrics will be extracted.
#' @param group.by Optional. Character string specifying the column name in the metadata to group cells by (e.g., `"cluster"` or `"cell_type"`). If specified, points in the scatter plot will be colored by the group.
#' @param metrics.by Character vector specifying the metrics to plot. Must contain either 2 or 3 metrics for generating 2D or 3D scatter plots, respectively.
#' @param split.by Optional. Character string specifying the column in the metadata to split the data by (e.g., `"orig.ident"` to visualize separately for each sample). If `NULL`, all cells are considered together. Default is `NULL`.
#' @param color Optional. A vector of colors to use for distinguishing groups in the scatter plot. If `NULL`, default colors are generated.
#' @param ggside Logical. If `TRUE`, includes marginal density plots (side plots) in the scatter plot for 2D visualizations. Default is `TRUE`.
#' @param log.x Logical. If `TRUE`, the x-axis will be log-transformed. Default is `FALSE`.
#' @param log.y Logical. If `TRUE`, the y-axis will be log-transformed. Default is `FALSE`.
#' @param log.z Logical. If `TRUE`, the z-axis will be log-transformed for 3D plots. Default is `FALSE`.
#' @param raster.cutoff Numeric value specifying the threshold for the number of points at which rasterization should occur, which can improve rendering performance for large datasets. Default is `100000`.
#' @param size Numeric value indicating the size of the points in the scatter plot for non-rasterized data. Default is `0.3`.
#' @param size.raster Numeric value indicating the size of the points in the scatter plot when rasterized. Default is `1`.
#' @param ticks Optional. Specifies the tick positions for the axes in the 3D plot. Default is `NULL`.
#'
#' @return A list of `ggplot` or `plotly` plot objects if `split.by` is specified; otherwise, a single plot object showing the relationships between the specified metrics.
#'
#' @export
PlotCellMetricsScatter <- function(object, group.by=NULL, metrics.by = NULL, split.by=NULL, color=NULL,
                                   ggside=T, log.x = F, log.y = F,log.z=F, raster.cutoff=100000, size =0.3, size.raster= 1, ticks=NULL ){
  if( "Seurat" %in% is(object)){
    object <- object@meta.data
  }

  if( !is.null(split.by) ){
    object <- split(object, object[[split.by]])
  }else{
    object <- list(all=object)
  }

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
        color <- get_colors(2)
      }
    }

    if(length(metrics.by)==2){
      plot_out <- plotScatter(object=x, x = metrics.by[1] ,y = metrics.by[2], group.by = group.by, color=color,
                              ggside=ggside, log.x = log.x, log.y = log.y, raster.cutoff=raster.cutoff, size =size, size.raster= size.raster)
    }else{
      plot_out <- plotScatter3D(object=x, x = metrics.by[1] ,y = metrics.by[2], z=metrics.by[3], group.by = group.by, color=color,
                                log.x = log.x, log.y = log.y, log.z=log.z, ticks=ticks)
      plot_out <- plotly::layout(plot_out, title =group.by
                         # ,
                         # height = 600
      )
    }

    return(plot_out)
  })
  names(plot_list) <- names(object)
  if(length(plot_list)==1){
    plot_list <-  plot_list[[1]]
  }
  return(plot_list)
}





#' Plot ADT Doublet Cutoff
#'
#' This function plots the ADT (Antibody-Derived Tag) cutoff values for specified features in a Seurat object.
#'
#' @param object A Seurat object containing the ADT data.
#' @param cutoff_table A data frame containing the cutoff values for the specified features obtained from `RunDB_ADT` function.
#' @param feature1 A character string specifying the first feature to plot.
#' @param feature2 A character string specifying the second feature to plot.
#' @param sample.by A character string specifying the metadata column to use for subsetting by sample. Default is "orig.ident".
#' @param group.by A character string specifying the metadata column to use for grouping. Default is NULL.
#' @param color A character string specifying the color to use for the plot. Default is NULL.
#' @param Sample A character vector specifying the samples to include in the plot. Default is NULL.
#' @param log.x A logical value indicating whether to log-transform the x-axis. Default is TRUE.
#' @param log.y A logical value indicating whether to log-transform the y-axis. Default is TRUE.
#' @param ggside A logical value indicating whether to use ggside for plotting. Default is TRUE.
#' @param size A numeric value specifying the size of the points in the plot. Default is 1.
#'
#' @return A ggplot object representing the scatter plot of the specified features with vertical and horizontal lines representing the cutoff values.
#'
#' @export
PlotADTCutoff <- function(object, cutoff_table, feature1=NULL, feature2=NULL,
                          sample.by="orig.ident",group.by=NULL,color=NULL, Sample=NULL, log.x=T, log.y=T,ggside=T,size=1) {
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }
  object$orig.ident <- object@meta.data[[sample.by]]
  object <- subset(object, subset = orig.ident %in% Sample)
  Seurat::DefaultAssay(object) <- "ADT"
  expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = "counts")

  # object <- Seurat::NormalizeData(object, normalization.method = 'CLR', margin = 2)
  # expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = "data")

  log_sums <- BPCells::colSums(log1p(expADT ), na.rm = TRUE)
  scale_factors <- exp(log_sums / nrow(expADT))
  expADT <- log1p(  BPCells::multiply_cols(expADT ,1/ scale_factors) )

  expADT1 <- as.numeric(as.matrix(expADT[feature1, ,drop=T]))
  expADT2 <- as.numeric(as.matrix(expADT[feature2, ,drop=T]))
  cut1 <- cutoff_table[cutoff_table$Sample %in% Sample & cutoff_table$feature1 %in% feature1, "cutoff1"]
  cut2 <- cutoff_table[cutoff_table$Sample %in% Sample & cutoff_table$feature2 %in% feature2, "cutoff2"]
  if(is.null(group.by)){
    expbind <- data.frame(expADT1, expADT2)
    colnames(expbind) <- c(feature1, feature2)
  }else{
    expbind <- data.frame(expADT1, expADT2, group=object@meta.data[[group.by]])
    colnames(expbind) <- c(feature1, feature2,group.by)
  }

  p1 <- plotScatter2(expbind, x=feature1, y=feature2, group.by=group.by, color =color,log.x = log.x, log.y = log.y, ggside = T,size = size )+
    ggplot2::geom_vline(xintercept = cut1, color="red", linewidth = 1)+
    ggplot2::geom_hline(yintercept = cut2, color="red", linewidth = 1)
  return(p1)
}





#' @title Get range of quality control metrics for different tissues
#'
#' @description This function retrieves and visualizes the range of specified quality control metrics for different tissue types. It generates bar plots for each metric, showing the distribution of values.
#'
#' @param tissue Character vector specifying the tissue type(s) to analyze.
#' @param metrics Character vector specifying the quality control metrics of interest (`"nCount_min"`, `"nCount_max"`, `"nFeature_min"`, `"nFeature_max"`, `"nCell_min"`, `"MT"`, `"RB"`, `"HB"`).
#'
#' @return A list of `ggplot` objects, each representing a bar plot for the specified metric.
#'
#' @export
GetMetricsRange <- function(tissue, metrics= c("nCount_min", "nCount_max",
                                               "nFeature_min", "nFeature_max",
                                               "nCell_min", "MT", "RB", "HB")) {
  # Filter data for the specified tissue
  data <- singlecellcutoff
  data <- data[data$Tissue %in% tissue, , drop = FALSE]

  # Function to create a bar plot
  create_plot <- function(plot_data, x_label, y_label = "Number") {
    col_name <- names(plot_data)[3]
    plot_data <- plot_data %>%
      dplyr::group_by(.data[[col_name]]) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::ungroup()

    ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(.data[[col_name]]), y = count, fill = factor(.data[[col_name]]))) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5, color = "black", size = 4) +
      ggplot2::labs(x = x_label, y = y_label) +
      ggplot2::scale_fill_manual(values = get_colors(length(unique(plot_data[[col_name]])))) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(fill = "Cutoff")
  }

  # Process each metric and generate plots
  results <- lapply(metrics, function(metric) {
    # Create a data frame for plotting
    plot_data <- data.frame(
      x = seq_along(data[[metric]]),
      y = data[[metric]],
      metric = data[[metric]]
    )

    # Generate the plot
    create_plot(plot_data, x_label = metric)
  })
  names(results) <- metrics

  if(length(results)==1){
    return(results[[1]])
  }

  # Return the list of plots
  return(results)
}
