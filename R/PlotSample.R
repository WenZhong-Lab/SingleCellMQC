
#' @title Plot Sample Cell Type Percentage
#' @description This function generates a bar plot or stacked bar plot showing the percentage of cell types in each sample.
#' @param object A Seurat object or a data frame containing the metadata information.
#' @param sample.by The column name in `object` that contains the sample information.
#' @param celltype.by The column name in `object` that contains the cell type information.
#' @param color.stackbar A vector of colors to use for the stacked bar plot. The length of the vector should be equal to the number of unique cell types. Default is NULL.
#' @param color.bar The color to use for the bar plot. Default is "#A6CEE3".
#' @param plot.type The type of plot to generate. Supported values include "stackbar" and "bar". Default is "stackbar".
#' @param nrow.bar The number of rows for the bar plot, if `plot.type` is "bar". Default is NULL.
#' @param ncol.bar The number of columns for the bar plot, if `plot.type` is "bar". Default is NULL.
#' @param levels The order of the samples to be displayed on the x-axis. If NULL, the order is determined by the order of the unique values in the `sample.by` column.
#' @param color.text.bar The color of the text annotations on the bar plot. Default is "black"
#' @param color.shadow.bar The color of the shadow behind the text annotations on the bar plot. Default is "white"
#' @param return.type The type of output to return. Supported values include "plot" and "interactive_table".
#' @param maxWidth The maximum width of the interactive table. Default is 85.
#'
#' @return Depending on the input parameters, this function can return a ggplot2::ggplot object for direct visualization in R, or an interactive data table that can be used for custom analyses or exported as CSV.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object with metadata information:
#' PlotSampleCellTypePCT(seurat_obj, sample.by = "orig.ident", celltype.by = "ScType")
#' }
#' @seealso \link{ShowSampleMetricsName}
#' @export
#'
PlotSampleCellTypePCT <- function(object, sample.by="orig.ident", celltype.by="ScType",
                                  color.stackbar=NULL,
                                  color.bar='#A6CEE3',
                                  plot.type="stackbar",nrow.bar=NULL,ncol.bar=NULL,
                                  levels=NULL,color.text.bar="black", color.shadow.bar="white",
                                  return.type = "plot", maxWidth=85){

  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }


  ##

  if( length(setdiff(return.type, c("table", "plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: 'table', `plot` or/and `interactive_table` ")
  }

  out <- list()

  metadata <- metadata[,c(sample.by, celltype.by)]
  colnames(metadata)[2] <- "CellType"
  plot_data <- data.table(metadata)[, .(count=.N), by=.(Sample=get(sample.by), CellType )]
  plot_data[, freq := round(count / sum(count) *100,2 ), by = Sample]
  if(!is.null(levels)){
    plot_data$Sample <- factor(plot_data$Sample, levels = levels)
  }


  if("table" %in% return.type ){
    freq_table <- plot_data
    freq_table <- data.table::dcast(freq_table, Sample ~ CellType, value.var="freq")
    freq_table <- as.data.frame(freq_table)
    out$table <- freq_table
  }


  if("interactive_table" %in% return.type){
    freq_table <- plot_data
    freq_table <- data.table::dcast(freq_table, Sample ~ CellType, value.var="freq")
    interactive_table <- .re_table(freq_table, csv.name = "freq_table", subtitle = "CellType% per Sample", number_type = "custom",number_fmr = "paste0(round(value,2), '%')",
                                   elementId = "freqtable", right.sparkline = T , down.sparkline = T, maxWidth=maxWidth)
    out$interactive_table <- interactive_table
  }

  if("plot" %in% return.type){
    if(is.null(color.stackbar)) {
      color.stackbar <- get_colors(length(unique(plot_data$CellType)))
    }
    plot_out <- switch (plot.type,
                        "stackbar" = ggplot2::ggplot(data = plot_data)+
                          ggplot2::geom_col (ggplot2::aes(x=Sample, y = freq, fill=CellType))+
                          ggplot2::theme_classic(base_size = 15) +
                          ggplot2::theme(
                            axis.text = ggplot2::element_text(color = "black"),
                            axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,vjust = 1),
                            axis.title.x = ggplot2::element_blank()
                          )+ggplot2::labs(y="Proportion")+
                          ggplot2::scale_fill_manual(values = color.stackbar ),
                        "bar" = {
                          color.bar <- rep(color.bar, length(unique(plot_data$Sample)))
                          ano_data <- plot_data
                          ano_data$freq <- paste0(ano_data$freq, "%(",ano_data$count,")")
                          return(plotBar(plot_data, sample_col = "Sample",split.by = "CellType",value_cols = "freq",
                                         reverse=F,WidetoLong = F, color =  color.bar,xlim = F,data_text = ano_data,
                                         nrow=nrow.bar,color.text = color.text.bar,color.shadow = color.shadow.bar,ncol = ncol.bar)+ggplot2::theme(legend.position = "none") )
                        },
                        stop("Invalid `plot.type`, only `stackbar`, `bar`. ")
    )
    out$plot <- plot_out

  }


  if(length(out)==1){
    out <- out[[1]]
  }

  return(out)
}




#' @title Plot Sample Metrics
#' @description This function generates a bar plot or interactive data table showing the specified metrics for each sample in the Seurat object.
#' @param object A Seurat object.
#' @param metrics A character vector specifying the metrics to plot.
#' @param type The type of metrics to show. Supported values include "count", "summary", and "Metrics_10x". Default is "count".
#' @param metrics_rename A character vector specifying the names to use for the metrics in the plot. Default is NULL.
#' @param color A character vector specifying the colors to use for the bar plot. Default is "#A6CEE3".
#' @param color.shadow The color of the shadow behind the text annotations on the bar plot. Default is "white".
#' @param reorder A character vector specifying the order of the samples to be displayed on the x-axis. Default is NULL.
#' @param color.text The color of the text annotations on the bar plot. Default is "black".
#' @param return.type The type of output to return. Supported values include "plot" and "interactive_table". Default is "plot".
#' @param csv.name The name of the CSV file to download. Default is "count".
#' @param elementId The ID of the element to use for the interactive table. Default is "count-table".
#' @param table.subtitle The subtitle to add to the interactive table. Default is "Metrics".
#' @param maxWidth The maximum width of the interactive table. Default is 85.
#'
#' @return Depending on the input parameters, this function can return a ggplot2::ggplot object for direct visualization in R, or an interactive data table that can be used for custom analyses or exported as CSV.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object:
#' # Visualisation of the count number of metrics
#' PlotSampleMetrics(seurat_obj, return.type = "plot")
#' PlotSampleMetrics(pbmc, type="count", metrics=c("nCell", "nGene_RNA"), return.type = "plot" )
#'
#' # Visualisation of the summary of percent.mt
#' PlotSampleMetrics(seurat_obj, type="summary", metrics=c("percent.mt"), return.type = "plot"  )
#'
#' # Visualisation  of the Metrics_10x
#' PlotSampleMetrics(seurat_obj, type="Metrics_10x",  return.type = "interactive_table"  )
#'
#' # return an interactive table
#' PlotSampleMetrics(seurat_obj, return.type = "interactive_table")
#'
#' }
#' @seealso \link{ShowSampleMetricsName}
#' @export
#'
PlotSampleMetrics <- function(object,
                              metrics= c("nCell", "nGene_RNA","nPro_ADT","nChain_TCR", "nCell_TCR", "TCR%", "nChain_BCR", "nCell_BCR", "BCR%"),
                              type= "count",
                              metrics_rename = NULL,
                              color= '#A6CEE3' ,
                              color.shadow ="white",
                              reorder=NULL,
                              color.text = "black",
                              return.type = "plot", ##"plot", "interactive_table"
                              csv.name='count',
                              elementId="count-table",
                              table.subtitle="Metrics", maxWidth=85){

  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  QC_misc <- GetSingleCellMQCData(object)

  if( length(names(QC_misc$perQCMetrics$perSample))==0 ){
    stop(" Please run CalculateMetrics() function first!! ")
  }

  if( length(setdiff(return.type, c("plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` or/and `interactive_table` ")
  }

  count_table <-  switch (type,
                          "count" = QC_misc$perQCMetrics$perSample$count,
                          "summary" = QC_misc$perQCMetrics$perSample$summary,
                          "Metrics_10x" = QC_misc$perQCMetrics$perSample$Metrics_10x,
                          stop("Invalid `type`, only `count`, `summary`, `Metrics_10x`. ")
  )

  if("nCell" %in% metrics){
    cat(paste0(
      paste0("> - ",nrow(count_table) , " samples contain ", sum(count_table[["nCell"]]), " cells, ",
             sum(count_table[["nCell_TCR"]]), " cells contain TCR information, ", sum(count_table[["nCell_BCR"]]), " cells contain BCR information.\n") ,
      paste0("> - ",sum(count_table[["nCell"]]<1000) , " samples less than 1000 cells : `" ,
             paste0(count_table$sample[count_table[["nCell"]]<1000], collapse = "`, `") ,"`\n") ,
      paste0("> - ",sum(count_table[["nCell"]]>20000) , " samples more than 20000 cells : `" ,
             paste0(count_table$sample[count_table[["nCell"]]>20000], collapse = "`, `") ,"`\n")
    ))
  }



  if(type=="summary"){
    inter_metrics <- intersect(metrics, names(count_table))


    if(length(inter_metrics) > 1){
      out_list <- lapply(inter_metrics, function(x){
        PlotSampleMetrics(object, type="summary", metrics=x,  metrics_rename=metrics_rename,  color= color ,
                          color.shadow =color.shadow,
                          reorder=reorder,
                          color.text = color.text,
                          return.type = return.type,
                          csv.name=paste0("summary","_",x ),
                          elementId=paste0("summary", "-",x , "-table"),
                          table.subtitle=x, maxWidth=maxWidth )

      })
      names(out_list) <- inter_metrics
      return(out_list)
    }

    if("nFeature_RNA" %in% inter_metrics){
      index <- count_table$nFeature_RNA$Median < 1000
      cat(paste0("> - ",sum(index) , " samples less than 1000 'Median genes per cell' : `" ,
                 paste0(count_table$nFeature_RNA[,1,drop=T][index], collapse = "`, `") ,"`\n"))
    }

    index <- match(inter_metrics, names(count_table))
    count_table <- count_table[[index]]
  }else{
    inter_metrics <- intersect(metrics, colnames(count_table))
    index <- match(inter_metrics, colnames(count_table))
    count_table <- count_table[, c(1, index)]
  }



  if(!is.null(metrics_rename)){
    colnames(count_table)[-1] <- metrics_rename
  }
  if(!is.null(reorder)){
    count_table[,1] <- factor(count_table[,1], levels = reorder)
  }

  out <- list()

  if( "interactive_table" %in% return.type){
    p1<-list()
    rownames(count_table) <-NULL

    p1$base <- .re_table(count_table,csv.name=csv.name, elementId=elementId,  subtitle = table.subtitle,first_name = "Sample",maxWidth=maxWidth,
                         down.sparkline = T,number_type = "auto"
    )
    results <- data.table(count_table[,-1])[, .(
      Max = sapply(.SD, max),
      Min = sapply(.SD, min),
      Mean = sapply(.SD, function(x) round(mean(x),3) ),
      Median = sapply(.SD, function(x) round(mean(x),3) )
    )]
    results <- t(results)
    colnames(results) <- colnames(count_table[,-1])
    results<- data.frame( Type= rownames(results), results, check.names = F )
    rownames(results)<-NULL
    p1$statistics <- .re_table(results,csv.name=paste0(csv.name,"_summary"), subtitle = paste0(table.subtitle," summary table"),
                               elementId=paste0("summary", elementId),  down.sparkline = T, number_type = "auto")
    out$interactive_table <- p1
  }


  if( "plot" %in% return.type){
    color <- rep(color, length(unique(count_table[,1])) )
    count_table<-suppressWarnings(data.table::melt(as.data.table(count_table), id.vars="sample", variable.name = "variable", value.name = "value"))
    p1 <- plotBar(count_table, sample_col = "sample",value_cols = "value",split.by = "variable",
                  WidetoLong = F, reverse=F, color = color ,color.shadow=color.shadow,color.text=color.text)+
      ggplot2::theme(legend.position = "none")
    out$plot <- p1

  }

  if(length(out)==1){
    out <- out[[1]]
  }

  return(out)
}




plotBar <-  function(data,
                     WidetoLong=T,
                     sample_col=NULL,
                     split.by=NULL,
                     value_cols=setdiff(colnames(data),sample_col),
                     color.by=sample_col,
                     data_text=data,
                     color= c("#377eb8", "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")
                     ,subtitle = "",
                     reverse=F,
                     nrow=1,
                     ncol=NULL,
                     color.shadow="white",
                     color.text = "black",
                     xlim=F) {
  # Define the scientific ggplot2::theme
  scientific_theme <- function(base_size = 12, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        text = ggplot2::element_text(color = "black"),
        # axis.title = ggplot2::element_text(size = rel(1.2), face = "bold"),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
        axis.text.x = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(hjust = 1),
        axis.line = ggplot2::element_line(color = "black"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "white", color = "black"),
        strip.text = ggplot2::element_text(face = "bold")
        # ,legend.position = "none"
      )
  }
  # Convert the data to long format
  if(WidetoLong){
    data <- data.table::setDT(data)
    data_long <- data.table::melt(data, measure.vars = value_cols, variable.name = "Variable", value.name = "Value")
    data_long$Value <- data_long$Value*100
    data_text <- data.table::setDT(data_text)
    data_text <- data.table::melt(data_text, measure.vars = value_cols, variable.name = "Variable", value.name = "Value")
    if(reverse){
      facet=sample_col
      sample_col="Variable"
    }else{
      sample_col=sample_col
      facet="Variable"
    }
    variable.name = "Variable"
    value.name="Value"
  }else{
    data_long<-data
    if(reverse){
      facet=sample_col
      sample_col=split.by
    }else{
      sample_col=sample_col
      facet=split.by
    }
    variable.name = color.by
    value.name=value_cols

  }


  if(color.by!=sample_col){
    color.by=variable.name
  }

  # Create the plot
  p <- ggplot2::ggplot(data_long, ggplot2::aes_string(y = sample_col, x = value.name, fill = color.by)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::facet_wrap(~ get(facet), scales = "free_x", nrow = nrow, ncol =ncol ) +
    shadowtext::geom_shadowtext(ggplot2::aes(y = data_text[[sample_col]], x = -Inf, label = paste(" ", data_text[[value.name]]) ), vjust = 0.5, hjust = 0,
                                bg.colour = color.shadow,color=color.text
    ) +
    ggplot2::scale_fill_manual(values = color) +
    scientific_theme() +
    ggplot2::labs(
      fill  = "Sample",
      y = sample_col,
      subtitle = subtitle
    )

  if(xlim){
    p <- p + ggplot2::coord_cartesian(xlim = c(0, max(data_long[[value.name]])))
  }

  return(p)
}


#' Plot V(D)J Sample Statistics
#'
#' Generates a plot or interactive data table for V(D)J metrics, subtype, CDR3 length, clonal quantification, clonal overlap, and gene usage.
#' This function is part of a suite of tools designed for exploring immunological data within single-cell datasets, focusing on V(D)J recombination data.
#'
#' @param object A Seurat object or a list containing the original VDJ table for each sample.
#' \itemize{
#' \item If type is one of "pct" or "subtype", a Seurat object processed with \code{\link{CalculateMetricsPerSample}} or \code{\link{CalculateMetrics}} is required as input.
#' \item If type is one of "CDR3", "clonalQuant", "clonalOverlap", or "geneUsage", lists obtained by the \code{\link{Read10XData}} function are required as input.
#' }
#' The input ideally contains V(D)J recombination data processed through the appropriate `scRepertoire` or similar analysis pipelines.
#'
#' @param color A character vector specifying colors to be used in the plot. The length of the vector should match the sample size.
#' @param type A character string indicating the type of V(D)J visualization or data to return. Supported types include:
#' \itemize{
#' \item "pct": for the percentage of chain type, Seurat object as input.
#' \item "subtype": for the percentage of chain subtype, Seurat object as input.
#' \item "CDR3": for the length distribution of the CDR3 sequences, raw list of "filtered_contig_annotations" data as input, details see \code{\link[scRepertoire]{clonalLength}} function from `scRepertoire` package.
#' \item "clonalQuant": for quantifying unique clonotypes, raw list of "filtered_contig_annotations" data as input, details see \code{\link[scRepertoire]{clonalQuant}} function from `scRepertoire` package.
#' \item "clonalOverlap": for measures of similarity between samples, raw list of "filtered_contig_annotations" data as input, details see \code{\link[scRepertoire]{clonalOverlap}} function from `scRepertoire` package.
#' \item "geneUsage": for the proportion of V gene usage, raw list of "filtered_contig_annotations" data as input, details see \code{\link[scRepertoire]{percentGenes}} function from `scRepertoire` package.
#' }
#' @param sample_name Optional; a character vector specifying sample names, applicable only for the `type` parameter values: "CDR3", "clonalQuant", "clonalOverlap", or "geneUsage".
#' @param size Numeric; the size of points in scatter plots. Defaults to 4.
#'
#' @param return.scRepertoire A logical indicating whether to return the list obtained from `scRepertoire::combineTCR` or `scRepertoire::combineBCR`, applicable only for the `type` parameter values: "CDR3", "clonalQuant", "clonalOverlap", or "geneUsage".
#' @param scRepertoire_data A list processed by `scRepertoire::combineTCR` or `scRepertoire::combineBCR`, applicable only for the `type` parameter values: "CDR3", "clonalQuant", "clonalOverlap", or "geneUsage". Default is NULL.
#' @param return.type Character string specifying the output type, either "plot" or "interactive_table". Defaults to "plot".
#'
#' @return Depending on the input parameters, this function returns either a ggplot2::ggplot object for direct visualization in R or an interactive data table that can be used for custom analyses or exported as a CSV. The specific output is determined by the `type` and `return.interactive_table` parameters.
#'
#' @export
#'
PlotSampleVDJ <- function(object, color=NULL,
                          type="pct", sample_name=names(object),
                          size=4,
                          return.type="plot",
                          return.scRepertoire=F,
                          scRepertoire_data=NULL
){
  if(type %in% c("pct", "subtype")){
    QC_misc <- GetSingleCellMQCData(object)

  }
  if( length(setdiff(type, c("pct", "subtype", "CDR3", "clonalQuant", "clonalOverlap", "geneUsage")))!=0 ){
    stop("Invalid `return.type`, only: `pct`, `subtype`, `CDR3`, `clonalQuant`, `clonalOverlap`, `geneUsage`")
  }

  if( length(setdiff(return.type, c("plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` or/and `interactive_table` ")
  }

  ##

  out <- list()

  if(type == "pct"){
    if(is.null(color)) {
      color <- get_colors(length(unique(QC_misc$perQCMetrics$perSample$count[,1])))
    }
    intersect_type <- intersect(colnames(QC_misc$perQCMetrics$perSample$count), c("sample", "TRA%", "TRB%", "IGH%", "IGK%", "IGL%", "TCR%", "BCR%") )
    point_data <- QC_misc$perQCMetrics$perSample$count[, intersect_type]
    point_data$sample <- as.character(point_data$sample)

    text <- c()
    if("TCR%" %in% intersect_type){
      text <-c(text,paste0("> - ", sum(point_data[["TRA%"]] > point_data[["TRB%"]]), " TRA% > TRB% samples, ", sum(point_data[["TRA%"]] <= point_data[["TRB%"]]), " TRA% <= TRB% samples." ,"\n" ))
    }
    if("BCR%" %in% intersect_type){
      text <-c(text,paste0("> - ", sum(point_data[["IGH%"]] > (point_data[["IGK%"]]+ point_data[["IGL%"]]) ), " IGH% > (IGK+IGL)% samples, ", sum(point_data[["IGH%"]] <= (point_data[["IGK%"]]+ point_data[["IGL%"]])), " IGH% <= (IGK+IGL)% samples." ,"\n" ))
    }

    if("TCR%" %in% intersect_type & "BCR%" %in% intersect_type){
      text <-c(text, paste0("> - ", sum(point_data[["BCR%"]] > point_data[["TCR%"]]), " BCR% > TCR% samples : `", paste0(point_data$sample[point_data[["BCR%"]] > point_data[["TCR%"]]], collapse = "`, `") ,"`\n"   ) )
    }

    cat(paste0(text,collapse = ""))


    if("interactive_table" %in% return.type){
      out$interactive_table <- PlotSampleMetrics(object, return.type ="interactive_table" ,metrics = c("TRA%", "TRB%", "IGH%", "IGK%", "IGL%", "TCR%", "BCR%",
                                                                                                       "nCell_TRA", "nCell_TRB","nCell_IGH","nCell_IGK","nCell_IGL","nCell_TCR","nCell_BCR",
                                                                                                       "nChain_TRA", "nChain_TRB","nChain_IGH","nChain_IGK","nChain_IGL","nChain_TCR","nChain_BCR"
      ), csv.name = "VDJ_pct", elementId = "vdj_pct-table", table.subtitle="V(D)J chain")
    }


    if("plot" %in% return.type){

      if("TCR%" %in%  intersect_type){
        out$plot$TCR <- ggplot2::ggplot(point_data, ggplot2::aes(x = `TRA%`, y = `TRB%`, label=sample)) +
          ggplot2::geom_point(ggplot2::aes(color =sample), shape = 19, size = size) +
          ggplot2::guides(color=ggplot2::guide_legend(title="Sample")) +
          ggplot2::scale_color_manual(values = color) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
            legend.position="none",
            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
            axis.text = ggplot2::element_text(color = "black"),
            axis.text.x = ggplot2::element_text(angle = 0,hjust = 1,vjust = 1)##
          )+
          ggrepel::geom_text_repel(show.legend = F) +
          ggplot2::labs(subtitle = "TCR proportion " , x= "TRA%",
               y="TRB%") +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "black")
      }

      if("BCR%" %in%  intersect_type){
        out$plot$BCR <- ggplot2::ggplot(point_data, ggplot2::aes(x = `IGH%` , y =`IGK%` + `IGL%`, label=sample)) +
          ggplot2::geom_point(ggplot2::aes(color =sample), shape = 19, size = size) +
          ggplot2::guides(color=ggplot2::guide_legend(title="Sample")) +
          ggplot2::scale_color_manual(values = color) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
            legend.position="none",
            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
            axis.text = ggplot2::element_text(color = "black"),
            axis.text.x = ggplot2::element_text(angle = 0,hjust = 1,vjust = 1)##
          )+
          ggrepel::geom_text_repel(show.legend = F) +
          ggplot2::labs(subtitle = "BCR proportion " , x= "IGH%",
               y="IGLC%") +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "black")
      }

      if(!is.null(out$plot$TCR) & !is.null(out$plot$BCR) ){
        out$plot$VDJ <- ggplot2::ggplot(point_data, ggplot2::aes(x =`TCR%` , y = `BCR%` , label=sample)) +
          ggplot2::geom_point(ggplot2::aes(color =sample), shape = 19, size = size) +
          ggplot2::guides(color=ggplot2::guide_legend(title="Sample")) +
          ggplot2::scale_color_manual(values = color) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
            legend.position="none",
            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
            axis.text = ggplot2::element_text(color = "black"),
            axis.text.x = ggplot2::element_text(angle = 0,hjust = 1,vjust = 1)##
          )+
          ggrepel::geom_text_repel(show.legend = F) +
          ggplot2::labs(subtitle = "Chain proportion " , x= "TCR%",
               y="BCR%") +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "black")
      }
    }
  }

  if(type == "subtype"){


    if("interactive_table" %in% return.type){
      index <- na.omit(unique(object@meta.data[, c("receptor_subtype")]))
      index <- c( paste0(index, "%"), index)
      out$interactive_table <- PlotSampleMetrics(object, return.type = "interactive_table" ,metrics = index, csv.name = "VDJ_subtype", elementId = "vdj_subtype-table",
                                                 table.subtitle="VDJ subtype")
    }

    if("plot" %in% return.type){
      if(is.null(color)) {
        color <- get_colors(length(unique(object@meta.data[,"receptor_subtype"])))
      }
      subtype_data <- object@meta.data[, c("orig.ident","receptor_subtype")]
      chain <- data.table(subtype_data)[, .(count = .N), by=c("orig.ident", "receptor_subtype")][, proportion := count / sum(count), by = orig.ident][!is.na(receptor_subtype)]
      group.by = "orig.ident"
      x = "proportion"
      chain_name = "receptor_subtype"
      out$plot <- ggplot2::ggplot(data = chain)+
        ggplot2::geom_col (ggplot2::aes(x=.data[[group.by]], y = .data[[x]], fill=.data[[chain_name]]))+
        ggplot2::theme_classic(base_size = 15) +
        ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
          # panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
          axis.text = ggplot2::element_text(color = "black"),
          axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,vjust = 1),
          axis.title.x = ggplot2::element_blank()
        )+ggplot2::labs(y="Proportion")+
        ggplot2::scale_fill_manual(values = color )
    }
  }

  if(type %in% c("CDR3", "clonalQuant", "clonalOverlap", "geneUsage") ){
    if(is.null(color)) {
      color <- get_colors(length(sample_name))
    }
    if( is.null(scRepertoire_data) ){
      if (packageVersion("scRepertoire") > "1.99.0") {
        colname <- colnames(object[[1]])
        index <- match("cdr3", colname)
        object <- lapply(object, function(x){
          x <- x[, c(1:12, index:length(colname))]
          return(x)
        })

        #
        if(sum(grepl("TR",unique(object[[1]]$chain))) > 0){
          combined <- scRepertoire::combineTCR(object, samples = sample_name)
        }else{
          combined <- scRepertoire::combineBCR(object, samples = sample_name)
        }
        combined <- lapply(combined, function(x){
          x[[3]] <- dplyr::na_if(x[[3]] ,"NA")
          return(x)
        })
        if(return.scRepertoire){
          return(combined)
        }
      }else{
        stop("Error: `scRepertoire` version must > 1.99.0")
      }
    }else{
      combined = scRepertoire_data
    }
  }

  if(type == "CDR3"){
    if("plot" %in% return.type){
      aa<-scRepertoire::clonalLength(combined, cloneCall="aa", chain = "both")
      data_len<-aa$data
      p1 <- suppressMessages(ggplot2::ggplot(data_len, ggplot2::aes(x = length)) +
                               ggplot2::stat_density(
                                 ggplot2::aes(color = values),
                                 linewidth = 0.7,
                                 geom = "line",
                                 position = "identity"
                               ) +
                               ggplot2::labs(subtitle = "Distribution of CDR3 lengths" , x="CDR3 length") +
                               ggplot2::scale_color_manual(values =color )+
                               ggplot2::theme_bw(base_size = 15)+
                               ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
                                 panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
                                 axis.text = ggplot2::element_text(color = "black"),
                                 axis.text.x = ggplot2::element_text(angle = 0,hjust = 1,vjust = 1)##
                               )+
                               ggplot2::scale_x_log10())
    }
    out$plot <- p1
  }

  if(type == "clonalQuant"){
    if(sum(grepl("TR",unique(object[[1]]$chain))) > 0){
      csv.name="TCR_clonalQuant"
      elementId= "TCR_clonalQuant-table"
    }else{
      csv.name="BCR_clonalQuant"
      elementId= "BCR_clonalQuant-table"
    }

    p1 <- suppressMessages(scRepertoire::clonalQuant(combined,
                                                     cloneCall="strict",
                                                     chain = "both",
                                                     scale = TRUE)+
                             ggplot2::scale_fill_manual(values = color))

    if("plot" %in% return.type){
      out$plot <- p1
    }

    if("interactive_table" %in% return.type){
      p1 <- p1$data[, c(2,1,3,4)]
      colnames(p1) <- c("Sample", "contigs", "total", "percentage%")
      out$interactive_table <- .re_table(p1,csv.name=csv.name, elementId= elementId, down.sparkline = T,subtitle = csv.name,number_type = "auto")
    }
  }

  if(type == "clonalOverlap"){
    if(sum(grepl("TR",unique(object[[1]]$chain))) > 0){
      csv.name="TCR_morisita"
      elementId= "TCR_morisita-table"
    }else{
      csv.name="BCR_morisita"
      elementId= "BCR_morisita-table"
    }

    p1 <- scRepertoire::clonalOverlap(combined,
                                      cloneCall = "strict",
                                      method = "morisita")
    p1$data <- na.omit(p1$data )
    if("plot" %in% return.type){
      out$plot <- p1
    }

    if("interactive_table" %in% return.type){
      out$interactive_table <- .re_table(p1$data,csv.name=csv.name, elementId= elementId,
                                         down.sparkline = T,subtitle = "Morisita")
    }

  }

  if(type == "geneUsage"){
    if(sum(grepl("TR",unique(object[[1]]$chain))) > 0){
      df.genes1 <- scRepertoire::percentGenes(combined,
                                              chain = "TRA",
                                              gene = "Vgene",
                                              exportTable = TRUE)
      df.genes2 <- scRepertoire::percentGenes(combined,
                                              chain = "TRB",
                                              gene = "Vgene",
                                              exportTable = TRUE)
      csv.name="TCR_geneUsage"
      elementId= "TCR_geneUsage-table"
    }else{
      df.genes1 <- scRepertoire::percentGenes(combined,
                                              chain = "IGH",
                                              gene = "Vgene",
                                              exportTable = TRUE)
      df.genes2 <- scRepertoire::percentGenes(combined,
                                              chain = "IGL",
                                              gene = "Vgene",
                                              exportTable = TRUE)
      csv.name="BCR_geneUsage"
      elementId= "BCR_geneUsage-table"
    }
    df.genes <- cbind(df.genes1, df.genes2)

    if("interactive_table" %in% return.type){
      re_data <- data.frame(Sample=rownames(df.genes), df.genes, check.names = F, row.names = NULL)
      out$interactive_table <- .re_table(re_data,csv.name=csv.name, elementId= elementId,right.sparkline = T,
                                         down.sparkline = T,subtitle = "geneUsage", number_type = "custom", number_fmr = "paste0( round(value*100,2),'%' )")
    }
    if("plot" %in% return.type){
      distance.matrix <- dist(scale(df.genes, center=TRUE, scale=TRUE),
                              method="euclidean")
      mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)# eig=TRUE willl return eigen values
      mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
      mds.values <- mds.stuff$points
      df <- data.frame(mds.values, V3=rownames(mds.values))

      out$plot <- ggplot2::ggplot(df, ggplot2::aes(x = X1, y = X2,label=df$V3)) +
        ggplot2::geom_point(ggplot2::aes(color =df$V3), shape = 19, size = 4) +
        ggplot2::guides(color=ggplot2::guide_legend(title="Samples")) +
        ggplot2::scale_color_manual(values = color) +
        ggplot2::theme_bw(base_size = 15)+
        ggplot2::theme(#text = ggplot2::element_text(face = "bold"),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
          axis.text = ggplot2::element_text(color = "black"),
          axis.text.x = ggplot2::element_text(angle = 0,hjust = 1,vjust = 1),##
          legend.position = "none"
        )+
        ggrepel::geom_text_repel(show.legend = F) +
        ggplot2::labs(subtitle = "Multi-dimensional Scaling " , x=paste0("Dim1 (", mds.var.per[1], "%)"),
             y=paste0("Dim2 (", mds.var.per[2], "%)"))
    }
  }

  if(length(out)==1){
    out <- out[[1]]
  }
}


CalculateFeaturePCT <- function(object,  split.by="orig.ident", assay="RNA", feature_list= list(Female="XIST", Male=c("DDX3Y", "UTY", "RPS4Y1") ) ){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  exp <- Seurat::GetAssayData(object, slot = "counts", assay = assay)
  row_indices <- split(1:ncol(exp), object@meta.data[[split.by]])

  out_list <- lapply( row_indices, function(y){

    detect_ratio <- lapply(feature_list, function(x){
      index <- intersect(x, rownames(exp))
      index <- match(index, rownames(exp))
      index_colsum <- Matrix::colSums(exp[index, y, drop=F])
      index_ratio <- sum(index_colsum>0) / length(index_colsum)
      return(index_ratio)
    })
    detect_ratio <- do.call(cbind, detect_ratio)
    return(detect_ratio)
  })

  out_list <- do.call(rbind, out_list)
  rownames(out_list) <- names(row_indices)
  return(out_list)
}



#' @title Plot Sample Label
#' @description The function is used to plot sample labels based on the percentage of features in each sample.
#' @param object A Seurat object or a data frame containing the percentage of features in each sample.
#' @param sample.by The column name in the meta.data slot of the Seurat object to split the data by.
#' @param assay The assay to use for the feature calculation.
#' @param feature_list A list of features to calculate the percentage of. The list should be named, with each element containing a vector of feature names.
#' @param color A vector of colors to use for the bar fill, corresponding to the unique values in `color.by`.
#' @param return.type A character vector indicating the type of output to return. Supported values include "plot" and "interactive_table".
#'
#' @return A ggplot2::ggplot object representing the bar plot.
#' @export
#'
PlotSampleLabel <- function(object,  sample.by="orig.ident", assay="RNA", feature_list= list(Female="XIST", Male=c("DDX3Y", "UTY", "RPS4Y1") ),
                            color=c("#A6CEE3", "#B2DF8A"),
                            return.type="plot"){
  if( ("Seurat" %in% is(object)) ){
    object <- CalculateFeaturePCT(object, split.by=sample.by, assay=assay, feature_list=feature_list)
  }

  ##

  if( length(setdiff(return.type, c("plot", "interactive_table", "table")))!=0 ){
    stop("Invalid `return.type`, only: `plot` , `table`, `interactive_table` ")
  }


  object <- data.frame(split=rownames(object), object, check.names = F)
  plot_data <- data.table::melt( data.table(object), id.vars=c("split"))
  out<-list()
  if("plot" %in% return.type ){
    out$plot <- ggplot2::ggplot(data = plot_data) +
      ggplot2::geom_col (ggplot2::aes(x=split, y = value, fill=variable),
               color = "white",alpha=0.8)+
      ggplot2::scale_fill_manual(values = color ) +
      ggplot2::labs(y = "Percent")+
      ggplot2::theme_bw(base_size = 14)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            axis.title.x = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
            axis.text = ggplot2::element_text(color = "black")
      )+
      ggplot2::labs(fill = "Group", subtitle = "The pct of features")
  }


  rownames(object) <- NULL
  index <- apply(object[,-1], 1, function(x){
    max_index <- which.max(x)
  })
  object$Predict <- colnames(object)[-1][index]
  colnames(object)[1] <- "Sample"

  if("interactive_table" %in% return.type){
    out$interactive_table <- .re_table(object, csv.name = "SampleLabel", elementId = "SampleLabel", right.sparkline = F, down.sparkline = F, subtitle = "The pct of features", maxWidth = NULL)
  }

  if("table" %in% return.type){
    out$table <- object
  }


  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}
