
#' Run Cell Annotation with ScType
#'
#' This function perfroms cell annotation using ScType method on a Seurat object.
#'
#' @param object A Seurat object.
#' @param use_internal_data Logical; if TRUE, uses an internal sctype dataset of immune cell markers for prediction. Default: TRUE.
#' @param group.by A character string specifying the metadata column to group cells by before prediction. Default: "seurat_clusters".
#' @param preprocess A character string specifying the dimensionality reduction method. Supported values include:
#' \itemize{
#' \item "rna.umap": Perform the analysis to UMAP step on RNA assay.
#' \item "adt.umap": Perform the analysis to UMAP step on ADT assay.
#' \item "wnn.umap": Perform the analysis to UMAP step using weighted nearest neighbor (WNN) method on RNA+ADT assays.
#' }
#' @param cutoff The cutoff value used for deciding when to label a cell type as 'Unknown'. Default:4
#' @param type.tissue A character string specifying the tissue type to use for marker selection. Default: "Immune system".
#' @param type.condition A character string specifying the condition type to use for marker selection. Default: NULL.
#' @param type.cell A character string specifying the cell type to use for marker selection. Default: "Normal".
#' @param data_source A character string specifying the data source to use for marker selection. Default: "ScType".
#' \itemize{
#' \item "ScType": Use the ScType database for marker selection.
#' \item "Cell_Taxonomy": Use the Cell_Taxonomy database for marker selection.
#' }
#' @param ntop An integer specifying the number of top markers to use for prediction. Default: 30.
#' @param split.by A character string specifying the metadata column to split the data by before prediction. Default: NULL.
#' @param return.name  A character string specifying the name of the new column to add to the Seurat object. Default: "ScType".
#' @param genelist A list of gene lists to use for prediction. Default: NULL.
#' @param resolution The resolution to use for clustering. Default is 1.
#' @param ... Additional arguments to pass to the `RunPipeline` function.
#'
#' @return A Seurat object with a new column "ScType" in the metadata slot. The column contains the predicted cell types.
#'
#' @details
#' The function first preprocesses the data if requested, and then uses the provided or internal gene list to score cells based on the ScType method.
#'
#' @examples
#' \dontrun{
#' # Assuming `seuratObj` is a Seurat object
#' seuratObj <- RunScType(seuratObj, group.by="seurat_clusters", type.tissue = "Immune system")
#' }
#'
#' @export
#' @references
#' Ianevski, Aleksandr et al. “Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data.” Nature communications vol. 13,1 1246. 10 Mar. 2022, doi:10.1038/s41467-022-28803-w
#' @seealso \code{\link{ShowDatabaseTissue}} for displaying the available tissue types in the ScType database.
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#' @importFrom methods as
#' @importFrom stats median na.omit

RunScType <- function(object, split.by= NULL, return.name="ScType",
                      genelist= NULL, use_internal_data =TRUE,
                      group.by="seurat_clusters",
                      preprocess=NULL,
                      resolution=1,
                      cutoff=3,
                      type.tissue = "Immune system",
                      type.condition = NULL, type.cell="Normal",  data_source="c", ntop=30,
                      ...
){

  if( !("Seurat" %in% class(object))){
    stop("Invalid object, must be a Seurat object")
  }
  Seurat::DefaultAssay(object) <- "RNA"
  if(use_internal_data & is.null(genelist) ){
    if(data_source=="Main"){
      genelist <- marker_all$Main
    }else{
      genelist = .marker_database(type.tissue=type.tissue, type.condition =type.condition, type.cell = type.cell, data_source=data_source,  ntop=ntop)
    }
  }else{
    genelist = genelist
  }
  data <- object
  if(!is.null(split.by)){
    cat(">>>>>> Split by: ", split.by, " ","\n")
    object_list <- splitObject(object, split.by = split.by, tmpdir="./temp/SingleCellMQC_tempScType/" )
    out <- lapply(names(object_list), function(x){
      cat(">>>>>>>>>>>>>>> ", "\n")
      cat( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "--------- Cell anotation: ", x, " ","\n") )
      cat(">>>>>>>>>>>>>>> ", "\n")
      sc_data <- object_list[[x]]

      if ( "BPCells" %in% attr(class(Seurat::GetAssayData(sc_data, assay = "RNA", slot = "counts")), "package") ) {
        sc_data <- SeuratObject::SetAssayData(object = sc_data,
                                              assay = "RNA",
                                              slot = "counts",
                                              new.data = as(Seurat::GetAssayData(sc_data, assay = "RNA", slot = "counts"), "dgCMatrix") )
      }

      if(!is.null(preprocess)){

        sc_data <- suppressMessages(RunPipeline(sc_data, preprocess =preprocess, resolution=resolution, ...))
      }else{
        sc_data <- suppressMessages(RunPipeline(sc_data, preprocess ="rna.umap", resolution=resolution, ...))
      }
      ano_result <- .runScType(object = sc_data, genelist = genelist, group.by=group.by, cutoff=cutoff)
      rm(sc_data)
      gc()
      return(ano_result)
    })
    names(out) <- NULL
    out <- do.call(rbind, out)
  }else{
    if(!is.null(preprocess)){
      object <- suppressMessages(RunPipeline(object, preprocess =preprocess, resolution=resolution, ...))
    }
    out <- .runScType(object = object, genelist = genelist, group.by=group.by, cutoff=cutoff)
  }

  colnames(out) <- return.name
  if(data_source=="Main"){
    out[,1][out[,1] == "NKT"] <- "T cell"
    out[,1][out[,1] %in% c("Mast", "Neutrophil") ] <- "Granulocyte"
    out[,1][out[,1] %in% c("Epithelium", "Unknown", "RBC", "Platelet") ] <- "Other"
    out$Immune <- out[[return.name]]
    out$Immune[out$Immune %in% c("T cell", "B cell", "NK", "DC", "Mon/Mac", "Granulocyte") ] <- "Immune"
    out$Immune[out$Immune %in% setdiff(unique(out[[return.name]]), c("T cell", "B cell", "NK", "DC", "Mon/Mac", "Granulocyte") ) ] <- "Other"
  }

  data <- Seurat::AddMetaData(data, metadata =out)
  return(data)
}

.runScType <- function(object, genelist, group.by, cutoff){
  ##
  gene_name <- unique( as.character(c(do.call(c, genelist[[1]] ) , do.call(c, genelist[[2]]) )) )
  object <- Seurat::ScaleData(object, features=c(gene_name),verbose = F)


  ##sketch
  group.by_final = group.by

  assay <- "RNA"

  if( "RenameDims" %in% class(SeuratObject::GetAssayData(object, assay = assay, slot = "scale.data")) ){
    scale_data <- methods::as(object = SeuratObject::GetAssayData(object, assay = assay, slot = "scale.data"), Class = "matrix")
  }else{
    scale_data <- SeuratObject::GetAssayData(object, assay=assay, slot = "scale.data")
  }

  es.max = .sctype_score(scRNAseqData = scale_data, scaled = TRUE, gs = genelist[[1]], gs2 = genelist[[2]])
  ###
  cL_resutls = do.call("rbind", lapply( unique(na.omit(object@meta.data[[group.by]])), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames( object@meta.data[object@meta.data[[group.by]] %in% cl, ] )]) , decreasing = !0)
    utils::head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(object@meta.data[[group.by]]==cl, na.rm=T)), 10)
  }))
  sctype_scores = cL_resutls %>% dplyr::group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/cutoff] = "Unknown"
  out <- sctype_scores$type[match(object@meta.data[[group.by_final]], sctype_scores$cluster)]
  out <- data.frame(ScType=out)
  rownames(out) <- rownames(object@meta.data)
  return(out)

}



#' @importFrom data.table .N .I :=
.marker_database <- function(type.tissue = "Immune system",  type.condition = NULL, type.cell=NULL, data_source="ScType", ntop=10,
                             Cell_Taxonomy.min_count=1, min.feature=2){
  database <- switch (data_source,
                      "ScType" = marker_all$ScType,
                      "Cell_Taxonomy" = marker_all$Cell_Taxonomy,
                      stop("Invalid data_source")
  )
  marker_data <- database$pos
  cell_all <- marker_data[, .(count=.N), by=.(Cell_standard, Cell_Marker)]

  select_data <- marker_data
  if(!is.null(type.tissue)){
    if(!(type.tissue %in% select_data$Tissue_standard ) ){
      stop("Invalid type.tissue")
    }
    select_data <- select_data[select_data$Tissue_standard %in% type.tissue, ]
  }
  if(!is.null(type.condition) & data_source=="Cell_Taxonomy" ){
    select_data <- select_data[select_data$Condition %in% type.condition, ]
  }

  if (!is.null(type.cell) ) {
    select_data <- select_data[select_data$Cell_type %in% type.cell, ]
  }

  select_data <- select_data[, .(count=.N), by=.(Cell_standard, Cell_Marker)]

  cell_count <- select_data[cell_all, on = .(Cell_standard, Cell_Marker), nomatch = 0][, .(Cell_standard, Cell_Marker, num.select=count, num.all = i.count)]

  if(data_source=="Cell_Taxonomy"){
    cell_count <- cell_count[cell_count$num.select>=Cell_Taxonomy.min_count,]
  }

  cell_order <- cell_count[order(num.select, num.all, decreasing = T), .SD[1:ntop], by = Cell_standard]
  cell_order <- stats::na.omit(cell_order)
  pos <- split(cell_order$Cell_Marker, cell_order$Cell_standard )



  if(data_source=="ScType"){
    st_neg <- database$neg
    if(!is.null(type.tissue)){
      st_neg <-st_neg[st_neg$Tissue_standard %in% type.tissue, ]
    }
    if(!is.null(type.cell)){
      st_neg <-st_neg[st_neg$Cell_type %in% type.cell, ]
    }
    st_neg <- split(st_neg$Cell_Marker,st_neg$Cell_standard)

  }
  neg <- switch (data_source,
                 "ScType" =   lapply(st_neg, function(x){
                   if( sum(is.na(x))>0 ){
                     x=character()
                   }
                   return(x)
                 }),
                 "Cell_Taxonomy" = stats::setNames(replicate(length(pos), character(), simplify = FALSE), names(pos) ),
                 stop("Invalid data_source")
  )
  lengths(pos)
  pos <- pos[names(pos)[lengths(pos) > min.feature]]
  neg <- neg[names(pos)[lengths(pos) > min.feature]]

  genelist <- list(pos=pos, neg=neg)
  return(genelist)
}
.sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T);
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }

  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  es.max
}





.MetricsOutlier <- function(object,
                            sample.by="orig.ident",
                            metrics.by = NULL,
                            type= NULL,
                            log= T,
                            nmads =rep(3, length(metrics.by)),
                            aggregate_type="median"){
  metrics_summary <- switch (aggregate_type,
                             "median" = data.table(object)[, lapply(.SD, function(x) stats::median( as.numeric(x), na.rm = TRUE) ),
                                                           by = sample.by, .SDcols = metrics.by],
                             "mean" = data.table(object)[, lapply(.SD, function(x) mean( as.numeric(x), na.rm = TRUE)),
                                                         by = sample.by, .SDcols = metrics.by],
                             "sum" = data.table(object)[, lapply(.SD, function(x) sum( as.numeric(x), na.rm = TRUE)),
                                                        by = sample.by, .SDcols = metrics.by],
                             stop("Invalid `aggregate_type`, only `median`, `mean`, `sum`. ")
  )

  metrics_out <- lapply(seq_along(metrics.by), function(x){
    x_mad <- scater::isOutlier(metrics_summary[, get(metrics.by[x])], type=type[x], log=log[x] ,nmads = nmads[x])
    x_out <- data.frame(sample=metrics_summary[,get(sample.by)], metrics_name=metrics.by[x], metrics_value=metrics_summary[, get(metrics.by[x])],
                        cutoff =attributes(x_mad)$thresholds[[type[x]]], type=type[x], log=log[x],isFlagged = x_mad)
    return(x_out)
  })
  metrics_out <- do.call(rbind, metrics_out)

  metrics_list <- lapply(seq_along(metrics.by), function(x){
    out <- as.character(metrics_out$sample[(metrics_out$metrics_name %in% metrics.by[x]) & metrics_out$isFlagged==TRUE ])
    if( !S4Vectors::isEmpty(out) & sum(is.na(out))==0 ){
      cat(" ", metrics.by[x], "warning samples:", paste0(out,collapse = "," ),"\n")
    }
    return(out)
  })
  names(metrics_list) <- metrics.by
  return(list(result_table=metrics_out, outlier_list=metrics_list) )
}

.outlier_reactable <- function(object, subtitle= paste0("Metrics warning results (MAD Statistics)"), csv_name="isFlagged",
                               header_font_size=14, font_size=14, cell_padding=4, maxWidth=70){
  colnames(object)[1] <- "Sample"

  plot_data <- object

  colnames(plot_data)[1:4] <- c("Sample", "Metrics", "Value", "Cutoff")
  plot_data <- plot_data[plot_data$isFlagged, ]
  if(  "split.by" %in% colnames(plot_data)){
    index <- "split.by"
  }else{
    index <- NULL
  }

  plot_data <- plot_data[,c(index, "Sample", "Metrics", "isFlagged",  "Value", "Cutoff", "type", "log")]

  plot_data$Value <- round(plot_data$Value, 3)
  plot_data$Cutoff <- round(plot_data$Cutoff, 3)
  out <- htmltools::browsable(
    htmltools::tagList(
      htmltools::div(htmltools::tags$label("Group by", `for` = "max-grouping-select")),
      htmltools::tags$select(
        id = "max-grouping-select",
        onchange = "Reactable.setGroupBy('max-grouping-table-metrics', this.value ? [this.value] : [])",
        htmltools::tags$option("None", value = ""),
        lapply(c("Metrics",  "Sample"), htmltools::tags$option)
      ),
      htmltools::tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('max-grouping-table-metrics', '", csv_name,".csv')")),
      htmltools::tags$hr("aria-hidden" = "true"),
      reactable::reactable(
        plot_data,
        columns = list(
          Sample = reactable::colDef(aggregate = "unique"),
          Metrics = reactable::colDef(aggregate = "unique"),
          isFlagged = reactable::colDef(aggregate = "count")
        ),
        defaultPageSize = 10,
        minRows = 5,
        elementId = "max-grouping-table-metrics",
        showPageSizeOptions = TRUE,
        theme = reactablefmtr::cosmo(header_font_size =header_font_size, font_size =font_size,  cell_padding =cell_padding),
        showPagination = TRUE,
        filterable = T
      )
      %>% .add_custom_text(text=subtitle, font_size = header_font_size)
    )
  )

  return(out)
}





.sampleScatterPlot <- function(object, x, y, color.by,color, label.by, cutoff, y.title="Value", point.size=2){
  set.seed(123)
  p1 <- ggplot2::ggplot(object, mapping = ggplot2::aes( x=.data[[x]], y= .data[[y]])) +
    ggplot2::geom_violin(linewidth=0, color="white",ggplot2::aes(fill = .data[[x]]),show.legend = F )+
    ggplot2::scale_fill_manual(values =  color)+
    ggplot2::geom_jitter(width = 0.2, height = 0,  size = point.size, ggplot2::aes( color=.data[[color.by]])) +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")  )+
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.text = ggplot2::element_text(color = "black", size = 13),
      ##angle
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )+  ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[[label.by]]),
      box.padding = 0.3
    )+
    ggplot2::geom_hline(yintercept = cutoff , color = "red", linewidth=1)+
    ggplot2::labs(y=y.title)

  return(p1)
}

#' @title Find Metrics Outlier Samples
#' @description This function identifies outlier samples based on the median absolute deviation (MAD) of the metrics in a Seurat object.
#'
#' @param object A Seurat object or a data frame.
#' @param sample.by A character string specifying the metadata column to group cells by before prediction. Default: "orig.ident".
#' @param metrics.by A character vector specifying the metrics to use for outlier detection. Default: c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent.mt", "percent.isotype").
#' @param split.by A character string specifying the metadata column to split the data by before prediction. Default: NULL.
#' @param type A character vector specifying the type of outlier detection for each metric. Supported values include:
#' \itemize{
#' \item "lower": Detect lower outliers.
#' \item "higher": Detect higher outliers.
#' }
#' @param log A logical vector specifying whether to log-transform the metrics before outlier detection. Default: c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE).
#' @param nmads A numeric vector specifying the number of median absolute deviations (MADs) to use for outlier detection. Default: c(3, 3, 3, 3, 3, 3).
#' @param aggregate_type A character string specifying the type of aggregation to use for outlier detection. Supported values include: "median", "mean", "sum". Default:"median"
#' @param color A character vector specifying the colors to use for the plot. Default: c('#ed556a', 'white').
#' @param return.type A character vector specifying the type of output to return. Supported values include:
#' \itemize{
#' \item "table": Return a table of outlier samples.
#' \item "interactive_table": Return an interactive table of outlier samples.
#' \item "plot": Return a plot of outlier samples.
#' }
#' @return A list containing the outlier samples and an interactive table of outlier samples.
#' @export
#' @importFrom ggplot2 .data
#'

FindSampleMetricsWarning <- function(object,
                                     sample.by="orig.ident",
                                     metrics.by = NULL,
                                     split.by = NULL,
                                     type= NULL,
                                     log= NULL,
                                     nmads = rep(3, length(metrics.by)),
                                     aggregate_type="median",
                                     color=c('#A6CEE3'),
                                     return.type=c("table")
){
  if("Seurat" %in% class(object)){
    QC_misc <- GetSingleCellMQCData(object)
    object <- object@meta.data
  }else{
    object <- object
    QC_misc <- NULL
  }

  if( length(setdiff(return.type, c("plot", "interactive_table", "table")))!=0 ){
    stop("Invalid `return.type`, only: `plot`, `interactive_table`, or/and `table`. ")
  }

  ##

  if(is.null(metrics.by)){
    metrics.by <- if (is.null(metrics.by)) {
      c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent.mt", "percent.isotype")
    } else {
      metrics.by
    }

    type <- c("lower", "lower", "lower", "lower", "higher", "higher")
    log <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
    nmads <- c(3, 3, 3, 3, 3, 3)
  }else{
    if(is.null(type)  ){
      stop("Invalid `type`, only `lower`, `higher` ")
    }else{
      if( sum(type %in% "lower")==0 & sum(type %in% "higher" )==0  ){
        stop("Invalid `type`, only `lower`, `higher` ")
      }
    }

  }

  if(is.null(log) | length(log)!=length(metrics.by) ){
    stop("Invalid `log`, length must equal to `metrics.by` ")
  }


  if(!is.null(QC_misc)){
    metrics_count <- intersect(metrics.by, colnames(QC_misc$perQCMetrics$perSample$count))
    ###
    ###
    object_count <- QC_misc$perQCMetrics$perSample$count
  }
  metrics_temp <- metrics.by
  metrics.by<- intersect(metrics.by, colnames(object))

  if( !is.null(split.by) ){
    object_list <- split(object, object[[split.by]] )
    out <- lapply(names(object_list), function(x){
      cat(">>>>>>", x, "\n")
      out <- .MetricsOutlier(object = object_list[[x]], sample.by=sample.by, metrics.by=metrics.by, type=type, log=log, nmads = nmads, aggregate_type=aggregate_type)
      out$result_table <- data.frame(split.by=x, out$result_table)
      return(out)
    })

    result_table <- lapply(out, function(x){
      x$result_table
    })
    names(result_table) <- names(object_list)
    rbind_table <- do.call(rbind, result_table)
    rownames(rbind_table) <- NULL

    outlier_list <- lapply(out, function(x){
      x$outlier_list
    })
    names(outlier_list) <- names(object_list)

    out <- list()
    if("table" %in% return.type){
      out$table <- list(table=rbind_table, list=outlier_list)
    }


    if("interactive_table" %in% return.type){
      out$interactive_table <- .outlier_reactable(rbind_table[,c(2:8,1)], subtitle = paste0("Metrics outlier results (MAD Statistics)"), csv_name = paste0("outlier_metrics_results") )
    }

    if("plot" %in% return.type){
      plot_list <- lapply(names(result_table), function(x){
        plot_data <-  result_table[[x]]
        plot_data$label <- ifelse(plot_data$isFlagged, as.character(plot_data$sample), "")
        plot_data <- split(plot_data, plot_data$metrics_name)

        plot_out <- lapply(names(plot_data), function(y){
          plot_data[[x]]$isFlagged <- factor(plot_data[[y]]$isFlagged, levels = c("TRUE", "FALSE") )
          .sampleScatterPlot(plot_data[[y]], x="metrics_name", y="metrics_value", color.by = "isFlagged", label.by = "label", cutoff= unique(plot_data[[y]]$cutoff),
                             y.title = paste0(aggregate_type, " ", y),color=color )
        })
        plot_out <- patchwork::wrap_plots(plot_out, guides = "collect")
        return(plot_out)
      })
      names(plot_list) <- names(result_table)
      out$plot <- plot_list

    }

    if(length(out)==1){
      out <- out[[1]]
    }

    return(out)
  }else{
    out <- .MetricsOutlier(object = object, sample.by=sample.by, metrics.by=metrics.by, type=type[match(metrics.by, metrics_temp)],
                           log=log[match(metrics.by, metrics_temp)], nmads = nmads[match(metrics.by, metrics_temp)], aggregate_type=aggregate_type)


    if(!is.null(metrics_count)){
      out_count <- .MetricsOutlier(object = object_count, sample.by="sample", metrics.by=metrics_count,
                                   type=type[match(metrics_count, metrics_temp)], log=log[match(metrics_count, metrics_temp)],
                                   nmads = nmads[match(metrics_count, metrics_temp)], aggregate_type=aggregate_type)
      out$result_table <- rbind(out$result_table, out_count$result_table )
      out$outlier_list <- c(out$outlier_list, out_count$outlier_list)
    }


    result <- list()
    if("table" %in% return.type ){
      result$table <- list(table=out$result_table, list=out$outlier_list)
    }

    if("interactive_table" %in% return.type ){
      mydata <- out$result_table
      interactive_table <- .outlier_reactable(mydata)
      result$interactive_table <- interactive_table
    }


    if("plot" %in% return.type){
      plot_data <- out$result_table
      plot_data$label <- ifelse(plot_data$isFlagged, as.character(plot_data$sample), "")
      plot_data <- split(plot_data, plot_data$metrics_name)
      plot_out <- lapply(names(plot_data), function(x){
        plot_data[[x]]$isFlagged <- factor(plot_data[[x]]$isFlagged, levels = c("TRUE", "FALSE") )
        .sampleScatterPlot(plot_data[[x]], x="metrics_name", y="metrics_value", color.by = "isFlagged", label.by = "label", cutoff= unique(plot_data[[x]]$cutoff), y.title = paste0(aggregate_type, " ", x),color=color )
      })
      plot_out <- patchwork::wrap_plots(plot_out, guides = "collect")
      result$plot <- plot_out
    }

    if(length(result)==1){
      result <- result[[1]]
    }

    return(result)
  }
}





#' @title Find Inter-Sample Cell Type Proportion Outliers
#'
#' @description
#' This function detects outliers in cell type composition across samples using PCA (Principal Component Analysis) and DBSCAN (Density-Based Spatial Clustering of Applications with Noise). It can be applied to a Seurat object or a data frame containing sample and cell type metadata.
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the metadata will be extracted. If a data frame is provided, it should contain the necessary columns for sample and cell type information.
#' @param sample.by A character string specifying the metadata column that identifies samples. Default: `"orig.ident"`.
#' @param celltype.by A character string specifying the metadata column that contains cell type information. Default: `"ScType"`.
#' @param split.by A character string specifying an optional metadata column to split the analysis by. Default: `NULL`.
#' @param color.pca A character vector specifying the colors to use for the PCA plot. If `NULL`, default colors will be used. Default: `NULL`.
#' @param color.bar A character string specifying the color to use for the bar plot. Default: `"#56B4E9"`.
#' @param confidence_level A numeric value specifying the confidence level for outlier detection (e.g., `0.95` for 95% confidence interval). Default: `0.95`.
#' @param labelOnlyOutlier A logical value specifying whether to label only the outlier samples in the PCA plot. Default: `TRUE`.
#' @param minPts An integer specifying the minimum number of points required to form a cluster in DBSCAN. Default: `3`.
#' @param method A character string specifying the method for outlier detection. Supported methods include:
#' - `"ellipse"`: Use Euclidean distance for PCA-based outlier detection.
#' - `"dbscan"`: Use DBSCAN for density-based outlier detection.
#' Default: `"ellipse"`.
#' @param top An integer specifying the number of top principal components to display in the contribution plots. Default: `Inf` (all components).
#' @param return.type A character vector specifying the type of output to return. Supported values include:
#' - `"table"`: Return a table of outlier results.
#' - `"plot"`: Return a plot of outlier results.
#' - `"interactive_table"`: Return an interactive table of outlier results.
#' Default: `"table"`.
#'
#' @return A list containing the following components based on the `return.type` parameter:
#' - If `"table"` is specified: A list with two elements:
#'   - `outlier`: A data frame containing the outlier detection results.
#'   - `contribution`: A data frame containing the contribution of each cell type to the principal components.
#' - If `"plot"` is specified: A list with four elements:
#'   - `pca`: A ggplot object showing the PCA plot with outliers labeled.
#'   - `contrib_pc1`: A ggplot object showing the contribution of cell types to PC1.
#'   - `contrib_pc2`: A ggplot object showing the contribution of cell types to PC2.
#'   - `contrib_pc1_2`: A ggplot object showing the combined contribution of cell types to PC1 and PC2.
#' - If `"interactive_table"` is specified: A list with two elements:
#'   - `outlier`: An interactive table of outlier results.
#'   - `contribution`: An interactive table of cell type contributions to the principal components.
#'
#' @seealso [RunScType] for cell type annotation.
#'
#' @export
FindInterSamplePCTOutlier <- function(object,
                                      sample.by = "orig.ident",
                                      celltype.by = "ScType",
                                      split.by=NULL,
                                      color.pca =NULL,
                                      color.bar= "#56B4E9",
                                      confidence_level = 0.95,
                                      labelOnlyOutlier = T,
                                      minPts=3,
                                      method="ellipse",
                                      top=Inf,
                                      return.type=c("table")

) {

  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }


  if( length(setdiff(return.type, c("table", "plot", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `table` or/and `plot` or/and `interactive_table` ")
  }

  out <- list()

  ##

  metadata <- metadata[,c(sample.by, celltype.by, split.by)]
  colnames(metadata)[2] <- "CellType"
  plot_data <- data.table(metadata)[, .(count=.N), by=.(Sample=get(sample.by), CellType )]
  plot_data[, freq := round(count / sum(count) *100,2 ), by = Sample]
  plot_data <- data.table::dcast(plot_data, Sample~ CellType, value.var = "freq")
  plot_data <- replace(plot_data, is.na(plot_data), 0)

  if(dim(plot_data)[1]==1){
    return(NA)
  }

  pca_res <- stats::prcomp(plot_data[,-1],scale.=T)
  pca_out <- data.frame(Sample=plot_data$Sample ,pca_res$x )
  percentage <- round(pca_res$sdev/sum(pca_res$sdev) * 100,  2)
  percentage <- paste0(colnames(pca_out)[2:3], " (", paste(as.character(percentage), "%", ")", sep = ""))
  pca_out$label <- pca_out$Sample
  rownames(pca_out) <- pca_out$Sample

  if(!is.null(split.by)){
    pca_out[[split.by]] <- metadata[[split.by]][match(pca_out$Sample, metadata[[sample.by]])]
  }

  contrib_pc1 <- factoextra::fviz_contrib(pca_res,choice="var",top = top, axes=1, fill = color.bar, color = color.bar)
  contrib_pc2 <- factoextra::fviz_contrib(pca_res,choice="var",top = top, axes=2, fill = color.bar, color = color.bar)
  contrib_pc1_2 <- factoextra::fviz_contrib(pca_res,choice="var",top = top, axes=1:2, fill = color.bar, color = color.bar)


  contrib_data <- data.frame(Group=contrib_pc1$data$name,
                             PC1_contrib = contrib_pc1$data$contrib,
                             PC2_contrib = contrib_pc2$data$contrib,
                             PC1_2_contrib = contrib_pc1_2$data$contrib
  )

  if("plot" %in% return.type ){
    if(is.null(color.pca)){
      if(is.null(split.by)){
        color.pca <- get_colors(2)[2]
      }else{
        color.pca <- get_colors(length(unique(metadata[[split.by]])))
      }
    }
    utils::capture.output(plot_PCA <- pcaCLOutlier(pca_out, confidence_level = confidence_level, return_type = "plot", color = color.pca, custom_group = split.by, custom_label = "label", method = method, minPts = minPts, labelOnlyOutlier=labelOnlyOutlier))
    plot_PCA <- plot_PCA+
      ggplot2::labs(subtitle = "PCA of cell type proportions" , x= percentage[1],
                    y=percentage[2])
    out$plot <- list(pca=plot_PCA, contrib_pc1=contrib_pc1, contrib_pc2=contrib_pc2, contrib_pc1_2 = contrib_pc1_2)

  }
  outlier=pcaCLOutlier(pca_out, confidence_level = confidence_level, return_type = "table", color = color.pca, custom_group = split.by, custom_label = "label",method = method, minPts = minPts)

  if("interactive_table" %in% return.type){
    temp <- outlier[, c(1:3, 7:8)]
    rownames(temp) <- NULL
    outlier_int_table <- .re_table_fmtr(temp,csv.name="CellPCToutlier", maxWidth = NULL,subtitle = "Cell PCT Outlier results",
                                        elementId="PCT-table",  down.sparkline = T,
                                        bar_type = "reactablefmtr::data_bars(
                             count,
                             text_position = 'inside-base',
                             number_fmt = scales::number,
                             fill_color 	='#A6CEE3'
                           )")

    inter_contri_table <- .re_table(contrib_data,csv.name="PCA_Contribution", maxWidth = NULL,subtitle="Contribution of groups to PCs",
                                    elementId="PCA_Contribution-table",  down.sparkline = T,other_sticky_column = 2, number_fmr = "paste0(round(value,2), '%')",
                                    number_type = "custom",
    )
    out$interactive_table <- list(outlier=outlier_int_table, contribution=inter_contri_table)
  }

  if("table" %in% return.type){
    out$table <- list(outlier=outlier, contribution=contrib_data)
  }


  if(length(out)==1){
    out <- out[[1]]
  }

  return(out)
}

pcaCLOutlier <- function(data, confidence_level = 0.95, return_type = "plot", color = get_colors(2), custom_group = NULL, custom_label = "ID", method = "ellipse", minPts = 3, labelOnlyOutlier = TRUE) {

  find_outliers <- function(df, level) {
    if (nrow(df) <= 3) {
      df$isOutlier <- FALSE
      return(df)
    }
    ellipse_points <- ggplot2::ggplot_build(
      ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::stat_ellipse(level = level)
    )$data[[1]]
    polygon <- cbind(ellipse_points$x, ellipse_points$y)

    points <- cbind(df$PC1, df$PC2)
    in_ellipse <- sp::point.in.polygon(points[,1], points[,2], polygon[,1], polygon[,2])
    df$isOutlier <- in_ellipse == 0
    return(df)
  }

  CompDBSCAN <- function(df, minPts) {
    pct_data <- df[, c("PC1", "PC2")]
    if(dim(pct_data)[1] < minPts) {
      df$isOutlier <- FALSE
      return(df)
    }
    kNN_distances <- dbscan::kNNdist(pct_data, k = minPts - 1)
    sorted_distances <- sort(kNN_distances)
    second_diff <- diff(sorted_distances, differences = 2)
    elbow_index <- which.max(second_diff) + 1
    eps_value <- sorted_distances[elbow_index]
    out <- dbscan::dbscan(pct_data, eps = eps_value, minPts = minPts)
    outliers <- rownames(pct_data)[out$cluster == 0]
    df$isOutlier <- rownames(df) %in% outliers
    return(df)
  }

  if (is.null(custom_group)) {
    if (method == "ellipse") {
      data_outlier <- find_outliers(data, confidence_level)
    } else if (method == "dbscan") {
      data_outlier <- CompDBSCAN(data, minPts)
    } else {
      stop("Invalid method. Use 'ellipse' or 'dbscan'.")
    }
    data_outlier["Group"] <- "All"
    custom_group <- "Group"
  } else {
    if (method == "ellipse") {
      data_outlier <- data %>%
        dplyr::group_by(.data[[custom_group]]) %>%
        dplyr::group_split() %>%
        lapply(function(group) find_outliers(group, confidence_level)) %>%
        dplyr::bind_rows()
    } else if (method == "dbscan") {
      data_outlier <- data %>%
        dplyr::group_by(.data[[custom_group]]) %>%
        dplyr::group_split() %>%
        lapply(function(group) CompDBSCAN(group, minPts)) %>%
        dplyr::bind_rows()
    } else {
      stop("Invalid method. Use 'ellipse' or 'dbscan'.")
    }
  }

  cat("Outlier samples:", paste0(data_outlier$label[data_outlier$isOutlier], collapse = ", "), "\n")

  if (return_type == "table") {
    return(data_outlier)
  } else if (return_type == "plot") {
    label_data <- if (labelOnlyOutlier) dplyr::filter(data_outlier, isOutlier) else data_outlier
    p <- ggplot2::ggplot(data_outlier, ggplot2::aes(x = PC1, y = PC2, color = .data[[custom_group]])) +
      ggplot2::geom_point(ggplot2::aes(shape = isOutlier), size = 3, alpha = 1) +
      ggplot2::scale_shape_manual(values = c(16, 17)) +
      ggplot2::stat_ellipse(level = confidence_level,show.legend = F) +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = .data[[custom_label]]),
        size = 4,
        show.legend = FALSE
      ) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        # panel.grid.major = ggplot2::element_blank(),
        # panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.text = ggplot2::element_text(color = "black"),
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, vjust = 1),
        legend.position = if (is.null(custom_group)) "none" else "right"
      ) +
      ggplot2::scale_color_manual(values = color)
    return(p)
  } else {
    stop("Invalid return_type. Use 'table' or 'plot'.")
  }
}




#' Generate Cell Ranger Alerts for QC Metrics
#'
#' This function generates alerts for 10X Genomics Single Cell QC metrics.
#'
#' @param object The input object, either a Seurat object or a data frame containing QC metrics.
#' @param return.type A character string indicating the return type: either "table" or "interactive_table".
#'                    Default is "table". Accepted values are "table" or "interactive_table".
#'
#' @return Depending on the `return.type` parameter, this function returns:
#'         - A data frame if `return.type` is "table".
#'         - An interactive HTML table if `return.type` is "interactive_table".
#'
#' @export
#'
#' @examples
#'\dontrun{
#' # Example 1: Using a Seurat object
#' alerts <- CellRangerAlerts(seurat_object, return.type = "table")
#'
#' # Example 2: Using a data frame directly
#' alerts <- CellRangerAlerts(qc_metrics_df, return.type = "interactive_table")
#'}
#'
CellRangerAlerts <- function(object, return.type="table"){
  if( !("Seurat" %in% methods::is(object)) ){
    object <- object
  }else{
    object <- GetSingleCellMQCData(object)[["perQCMetrics"]][["perSample"]][["Metrics_10x"]]
  }

  if(is.null(object)){
    return(NULL)
  }

  if( length(setdiff(return.type, c("table", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `table` or/and `interactive_table` ")
  }

  CandidateMetrics <- intersect(colnames(object), unique(c(mertrics10x_table$Metrics_name1, mertrics10x_table$Metrics_name2)) )

  ##pre
  out <- lapply(CandidateMetrics, function(x){
    mertrics_value <- object[[x]]
    alert_temp <- mertrics10x_table[mertrics10x_table$Metrics_name1 %in% x | mertrics10x_table$Metrics_name2 %in% x,]
    alert_data <- switch ( as.character(alert_temp$index),
                           "0" = {data.table(Sample=object[[1]], Error = !(mertrics_value < alert_temp$acceptable) , Warning= !(mertrics_value < alert_temp$targeted) )},
                           "1" = {data.table(Sample=object[[1]], Error = !(mertrics_value > alert_temp$acceptable) , Warning= !(mertrics_value > alert_temp$targeted) )},
                           "2" = {data.table(Sample=object[[1]], Error = !(mertrics_value >= alert_temp$acceptable) , Warning= !(mertrics_value >= alert_temp$targeted) )}
    )
    alert_data <- data.table::melt(alert_data, id.vars="Sample")
    alert_data$Metrics <- x

    alert_data$`Error cutoff` <- rep(alert_temp$acceptable, length(object[[1]])*2 )
    alert_data$`Warning cutoff` <- rep(alert_temp$targeted, length(object[[1]])*2)
    alert_data$Value <- rep(object[[x]], 2)
    return(alert_data)
  })
  names(out) <- CandidateMetrics

  out_table <- list()

  ##Alerts_table
  Alerts_table <- lapply(out, function(x){
    if( sum(x$value, na.rm = T)==0 ){
      return(NULL)
    }else{
      alerts_data <- x[x$value,]
      alerts_data <- alerts_data[, .SD[which.max(variable == "Error")], by = setdiff(colnames(alerts_data), "variable") ]
      return(alerts_data)
    }

  })
  Alerts_table <- do.call(rbind, Alerts_table)
  Alerts_table <- as.data.frame(Alerts_table)
  Alerts_table <- Alerts_table[, -grep("value", colnames(Alerts_table))]
  Alerts_table <- Alerts_table[, c("Sample", "variable", "Metrics", "Value", "Error cutoff", "Warning cutoff")]
  colnames(Alerts_table)[1:2] <- c("Sample", "Alerts type")
  if( "table" %in% return.type){
    out_table$table <- Alerts_table
  }
  if( "interactive_table" %in% return.type){
    out_table$interactive_table <-  htmltools::browsable(
      htmltools::tagList(
        htmltools::div(htmltools::tags$label("Group by", `for` = "max-grouping-select")),
        htmltools::tags$select(
          id = "max-grouping-select",
          onchange = "Reactable.setGroupBy('max-grouping-table-cellranger', this.value ? [this.value] : [])",
          htmltools::tags$option("None", value = ""),
          lapply(c("Alerts type", "Metrics", "Sample"), htmltools::tags$option)
        ),
        htmltools::tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('max-grouping-table-cellranger', 'Metrics_10X_Alerts.csv')")),
        htmltools::tags$hr("aria-hidden" = "true"),

        reactable::reactable(
          Alerts_table,
          columns = list(
            Sample = reactable::colDef(aggregate = "unique"),
            Metrics = reactable::colDef(aggregate = "count"),
            `Alerts type` = reactable::colDef(aggregate = "unique")
          ),
          defaultPageSize = 10,
          minRows = 5,
          elementId = "max-grouping-table-cellranger",
          showPageSizeOptions = TRUE,
          theme = reactablefmtr::cosmo(header_font_size =14, font_size =14,  cell_padding =4),
          showPagination = TRUE,
          filterable = T
        ) %>% .add_custom_text(text="10X Metrics Alerts", font_size = 16)
      )
    )
  }


  if(length(out_table)==1){
    out_table <- out_table[[1]]
  }
  return(out_table)
}



#' Find Common PCT Outliers
#'
#' This function identifies common pct (percentage) outliers based on a built-in database of cell proportions.
#'
#' @param object A Seurat object or a data frame containing the metadata.
#' @param sample.by A character string specifying the metadata column to use for identifying samples. Default is "orig.ident".
#' @param celltype.by A character string specifying the metadata column to use for identifying cell types. Default is "ScType".
#' @param return.type A character vector specifying the type of output to return. Options are "table" and/or "interactive_table". Default is "table".
#' @param tissue A character string specifying the tissue type to filter by. If NULL, the function will stop and prompt the user to specify a tissue.
#'
#' @return Depending on the `return.type` parameter, the function returns:
#' \itemize{
#'   \item If "table" is specified, a data frame containing the outlier results.
#'   \item If "interactive_table" is specified, an interactive table generated using the `reactable` package.
#'   \item If both are specified, a list containing both the table and the interactive table.
#' }
#'
#' @export
#'
#' @details
#' To detect outlier samples based on cell type composition, one approach involves compiling and integrating data from established tissue atlases,
#' such as those from the DISCO database and additional single-cell reference datasets,
#' to define the common range of inter-individual variability for major typical cell types within each tissue.
#' The function then identifies samples exhibiting deviations from these established reference proportions, marking them as potential outliers.
#'
#'
FindCommonPCTOutlier <- function(object,
                                 sample.by = "orig.ident",
                                 celltype.by = "ScType",
                                 return.type=c("table"),
                                 tissue = NULL

) {

  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }
  if( length(setdiff(return.type, c("table", "interactive_table")))!=0 ){
    stop("Invalid `return.type`, only: `table` or/and `interactive_table` ")
  }

  if(is.null(tissue)){
    stop("Error: The tissue must be specified!!, please run `ShowCommonPCTTissue` to get the available tissue. ")
  }

  out <- list()
  out_table <- .perSamplePCTOutlier(metadata, tissue = tissue, celltype.by = celltype.by, sample.by = sample.by)
  if("table" %in% return.type){
    out$table <- out_table
  }

  if("interactive_table" %in% return.type){
    out$interactive_table <- .outlier_reactable_pct(out_table$table)
  }

  if(length(out)==1){
    out <- out[[1]]
  }

  return(out)

}

.outlier_reactable_pct <- function(object, subtitle="Common range outlier results", csv_name="outlier_pct_range_persample",
                                   header_font_size=14, font_size=14, cell_padding=4, maxWidth=70){
  colnames(object)[1] <- "Sample"

  plot_data <- object
  plot_data <- plot_data[plot_data$isOutlier, ]


  out <- htmltools::browsable(
    htmltools::tagList(
      htmltools::div(htmltools::tags$label("Group by", `for` = "pct-grouping-select")),
      htmltools::tags$select(
        id = "pct-grouping-select",
        onchange = "Reactable.setGroupBy('pct-grouping-table', this.value ? [this.value] : [])",
        htmltools::tags$option("All", value = ""),
        lapply(c("isOutlier", "CellType", "Sample"), htmltools::tags$option)
      ),
      htmltools::tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('pct-grouping-table', '", csv_name,".csv')")),
      htmltools::tags$hr("aria-hidden" = "true"),
      reactable::reactable(
        plot_data,
        columns = list(
          Sample = reactable::colDef(aggregate = "unique"),
          CellType = reactable::colDef(aggregate = "unique"),
          isOutlier = reactable::colDef(aggregate = "count")
        ),
        defaultPageSize = 10,
        minRows = 5,
        elementId = "pct-grouping-table",
        showPageSizeOptions = TRUE,
        theme = reactablefmtr::cosmo(header_font_size =13, font_size =13,  cell_padding =4),
        showPagination = TRUE,
        filterable = T
      ) %>% .add_custom_text(text=subtitle, font_size = 16)
    )
  )
  return(out)
}


.perSamplePCTOutlier <- function(object, tissue, celltype.by, sample.by, type="min-max"){
  metadata <- object[, c(sample.by, celltype.by)]
  pct_table <- data.table::data.table(metadata)[, list(.N), by=.(Sample=get(sample.by), CellType=get(celltype.by))]
  pct_table[, total := sum(N), by = Sample]
  pct_table[, proportion := N / total]

  tissue_table <- pct_stat_out[pct_stat_out$tissue==tissue, ]
  index <- match(pct_table$CellType, tissue_table$CellType)
  if(type == "min-max"){
    pct_table$min_cutoff <- round(tissue_table$min[index]*100,3)
    pct_table$max_cutoff <- round(tissue_table$max[index]*100,3)
    pct_table$isOutlier <- pct_table$proportion < tissue_table$min[index] | pct_table$proportion > tissue_table$max[index]
  }else if(type == "Q1-Q3"){
    pct_table$Q1_cutoff <- round(tissue_table$Q1[index]*100,3)
    pct_table$Q3_cutoff <- round(tissue_table$Q3[index]*100,3)
    pct_table$isOutlier <- pct_table$proportion < tissue_table$Q1[index] | pct_table$proportion > tissue_table$Q3[index]
  }else{
    stop("Error: The type must be 'min-max' or 'Q1-Q3'!!")
  }
  pct_table$proportion <- round(pct_table$proportion*100,3)
  colnames(pct_table)[5:7] <- paste0(colnames(pct_table)[5:7], "%")

  out <-  list(table=pct_table)

  metrics.by <- unique(pct_table$CellType)
  pct_list <- lapply(seq_along(metrics.by), function(x){
    out <- pct_table$Sample[(pct_table$CellType %in% metrics.by[x]) & pct_table$isOutlier==TRUE ]
    if( !S4Vectors::isEmpty(out) & sum(is.na(out))==0 ){
      cat(" ", metrics.by[x], "outlier samples:", paste0(out,collapse = "," ),"\n")
    }
    return(out)
  })
  names(pct_list) <- metrics.by
  pct_list <- pct_list[!(sapply(pct_list, is.null) | sapply(pct_list, function(x) length(x) == 0))]
  out$list <- pct_list

  ##
  if(length(pct_list)!=0){

    # subset freq >= 3
    freq_table <- as.data.frame(table( do.call(c, pct_list) ))
    freq_table <- freq_table[freq_table$Freq >=3,]
    freq_table <- freq_table[ order(freq_table$Freq, decreasing =T) , ]
    cat_freq_word3ormore <- paste0( "Note: ", length(freq_table[,1]), " warning samples (celltype% outliers >= 3 ): ", paste0(freq_table[,1], collapse = ", ")," \n")
    cat(cat_freq_word3ormore)
    out$warining = as.character(freq_table[,1])
  }else{
    out$warining = NA
  }
  return(out )
}



