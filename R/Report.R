

#' @title Generate Comprehensive Single-Cell QC HTML Report
#'
#' @description This function orchestrates the generation of an interactive HTML report for single-cell data quality control (QC).
#' It consolidates various QC modules, including sample-level, cell-level, feature-level, and batch effect assessments,
#' into a unified and easily digestible document. Users can customize which QC modules are included and specify
#' parameters relevant to their data and analysis goals.
#'
#' @param object A Seurat or similar single-cell object containing preprocessed single-cell data.
#'               It is expected that basic preprocessing steps (e.g., using `RunPreprocess` or equivalent)
#'               have already been applied.
#' @param VDJ_data An optional list containing processed VDJ (T-cell receptor or B-cell receptor) data.
#'                 If provided, relevant VDJ-specific QC metrics will be integrated into the report.
#'                 Defaults to `NULL`.
#' @param sample.by The name of the column in metadata that uniquely identifies different samples.
#'                  This is crucial for grouping and analyzing data across different samples.
#'                  Defaults to `"orig.ident"`.
#' @param outputFile The full path, including filename, where the generated HTML report will be saved.
#'                   The directory will be created if it does not exist.
#'                   Defaults to `"./SingleCellMQC/SingleCellMQC.html"`.
#' @param RNA_cluster_name The name of the column in metadata that stores RNA-based cluster assignments.
#'                         Used for visualizing RNA clustering and batch effects.
#'                         Defaults to `"rna_cluster"`.
#' @param ADT_cluster_name The name of the column in metadata that stores ADT-based cluster assignments.
#'                         Used for visualizing ADT clustering and batch effects in multi-modal data.
#'                         Defaults to NULL.
#' @param RNA.batch.by A character vector specifying one or more column names in metadata
#'                     that represent categorical variables used for RNA batch identification or for
#'                     grouping cells in PCA/UMAP plots related to RNA data.
#'                     Defaults to `"orig.ident"`.
#' @param ADT.batch.by A character vector specifying one or more column names in metadata
#'                     that represent categorical variables used for ADT batch identification or for
#'                     grouping cells in PCA/UMAP plots related to ADT data.
#'                     Defaults to NULL.
#' @param celltype.by The name of the column in metadata that contains cell type annotations.
#'                    This is used for cell type-specific QC metrics and visualizations.
#'                    Defaults to `"ScType"`.
#' @param tissue An optional character string specifying the tissue source of the samples.
#'               This information can be incorporated into sample-level QC plots for context.
#'               Defaults to `NULL`.
#' @param do.sample A logical value indicating whether the sample-level QC module should be included
#'                  in the report. Defaults to `TRUE`.
#' @param do.cell A logical value indicating whether the cell-level QC module should be included
#'                in the report. Defaults to `TRUE`.
#' @param do.feature A logical value indicating whether the feature-level (gene/ADT) QC module
#'                   should be included in the report. Defaults to `TRUE`.
#' @param do.batch A logical value indicating whether the batch effect assessment module should be
#'                 included in the report. Defaults to `TRUE`.
#' @param RNA.covariate.formula An optional formula (e.g., `~ (1|batch1) + (1|batch2)`) specifying covariates for
#'                              RNA data used in batch effects analysis.
#'                              Defaults to `NULL`.
#' @param ADT.covariate.formula An optional formula (e.g., `~ (1|batch1) + (1|batch2)`) specifying covariates for
#'                              ADT data used in batch effects analysis.
#'                              Defaults to `NULL`.
#' @param RNA.covariate.other A character vector of additional covariates (column names from `meta.data`)
#'                            relevant for RNA data analysis, but not necessarily part of a formula.
#'                            Defaults to `NULL`.
#' @param ADT.covariate.other A character vector of additional covariates (column names from `meta.data`)
#'                            relevant for ADT data analysis.
#'                            Defaults to `NULL`.
#' @param do.covariate.split.celltype A logical value. If `TRUE`, covariate analysis plots (if applicable)
#'                                    will be split or faceted by cell type.
#'                                    Defaults to `FALSE`.
#' @param marker_genes A character vector of well-known marker genes to be specifically highlighted
#'                     in the feature QC section, such as those used for cell type identification.
#'                     Defaults to a predefined list of immune cell markers.
#'
#' @return This function does not explicitly return a value. Instead, it generates and saves
#'         a comprehensive HTML report to the path specified by `outputFile`.
#' @export
RunReport <- function(object=NULL, VDJ_data=NULL, sample.by="orig.ident", outputFile="./SingleCellMQC/SingleCellMQC.html",
                      do.sample = TRUE,
                      do.cell = TRUE,
                      do.feature=TRUE,
                      do.batch=TRUE,
                      RNA_cluster_name="rna_cluster",
                      ADT_cluster_name=NULL,
                      RNA.batch.by= "orig.ident",
                      ADT.batch.by=NULL,
                      RNA.covariate.formula = NULL,
                      ADT.covariate.formula = NULL,
                      RNA.covariate.other=NULL,
                      ADT.covariate.other=NULL,
                      celltype.by="ScType",
                      do.covariate.split.celltype=F,
                      tissue=NULL,
                      marker_genes=c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                     "CD14","FCGR3A","CD68","FCN1","ITGAX")
){


  outputFile <- R.utils::getAbsolutePath(outputFile)
  plot_out = dirname(outputFile)
  output_file = basename(outputFile)
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")


  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/"), recursive = TRUE)
  }

  metadata <- object@meta.data
  labels <- unique(object@meta.data[, sample.by, drop=T])

  nsample <- length(unique(metadata[,sample.by]))
  SCMQC <- GetSingleCellMQCData(object)
  out_list <- list()
  out_list$nsample <- nsample
  out_list$SCMQC <- SCMQC

  # return(out_list)
  cat("\n")
  message("Making QC Report!")
  cat("Output report to:", outputFile, "\n")

  #sample qc
  if( do.sample ){
    out_list$Sample <- report.sample(object=object,
                                     VDJ_data=VDJ_data,
                                     sample.by=sample.by,
                                     color=NULL,
                                     plot_out=plot_out,
                                     section= 1:3,
                                     celltype.by = celltype.by,
                                     tissue=tissue)
  }

  #cell qc
  if( do.cell){
    out_list$Cell <- report.cell(object=object, plot_out=plot_out,sample.by=sample.by, section = c(1:2))
  }

  # feature qc
  if(do.feature){
    out_list$Feature <- report.feature(object = object ,sample.by = sample.by,plot_out = plot_out,
                                       celltype.by = celltype.by,
                                       marker_name = marker_genes)
  }


  #batch qc
  if( do.batch){
    out_list$Batch <- report.batch(object, RNA_cluster_name=RNA_cluster_name, ADT_cluster_name=ADT_cluster_name,
                                   sample.by=sample.by,
                                   RNA.batch.by= RNA.batch.by, ADT.batch.by=ADT.batch.by,
                                   celltype.by=celltype.by,
                                   RNA.covariate.other=RNA.covariate.other,
                                   ADT.covariate.other=ADT.covariate.other,
                                   RNA.covariate.formula=RNA.covariate.formula,
                                   ADT.covariate.formula=ADT.covariate.formula,
                                   do.covariate.split.celltype=do.covariate.split.celltype,
                                   plot_out, section=1:3)
  }

  cat("\n")
  message("Making QC Report!")
  cat("Output report to:", outputFile, "\n")
  reportRmd <- system.file("report/0_main.Rmd",
                           package = "SingleCellMQC")
  css_file <- system.file("report", "bootstrap.css", package = "SingleCellMQC")
  target_dir <- plot_out
  file.copy(css_file, target_dir, overwrite = TRUE)

  rmarkdown::render(reportRmd,
                    output_dir = plot_out,
                    output_file = output_file,
                    params = list(data = out_list, cell_dir=plot_dir) )

  return(NULL)
}


report.sample <- function(object, VDJ_data=NULL,
                          sample.by="orig.ident",
                          color=NULL,
                          plot_out=NULL,
                          tissue=NULL,
                          celltype.by=NULL,
                          split.by=NULL,
                          section= c(1:3)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Sample <- list()
  SCMQC <- GetSingleCellMQCData(object)
  labels <- unique(object@meta.data[, sample.by, drop=T])

  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Summary key tasks"))
  Sample$summary <- SummarySample(object = object, sample.by=sample.by, tissue = tissue, celltype.by=celltype.by)
  Sample$summary <- .re_table(Sample$summary, csv.name = "SummarySample", elementId = "SummarySample"  )

  if(1 %in% section){
    Sample$Section1  <- .sampleSection1(object, VDJ_data=VDJ_data, plot_out=plot_out,sample.by=sample.by, split.by=split.by)
  }

  if(2 %in% section & !is.null(celltype.by)){
    Sample$Section2  <- .sampleSection2(object,  plot_out=plot_out, celltype.by = celltype.by, tissue=tissue,sample.by=sample.by, split.by=split.by)
  }

  if(3 %in% section){
    Sample$Section3  <- .sampleSection3(object,sample.by=sample.by)
  }
  return(Sample)
}

.sampleSection1 <- function(object, VDJ_data=NULL, plot_out=NULL, split.by=NULL, sample.by="orig.ident"){
  outlist <- list()
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment"))

  ##Alerts
  Alert_table <- CellRangerAlerts(object, return.type = "table" )
  if(!is.null(Alert_table)){
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment --- Checking 10X metrics"))
    Alert_sample_list <- split(Alert_table$Sample, Alert_table$`Alerts type`)
    Alert_sample_list$Error <- unique(Alert_sample_list$Error)
    Alert_sample_list$Warning <- unique(Alert_sample_list$Warning)

    outlist$Alert$table <- Alert_sample_list
    outlist$Alert$interactive_table <- CellRangerAlerts(object, return.type = "interactive_table")

    outlist$Metrics$Metrics_10X <- PlotSampleMetrics(object, type = "Metrics_10x",
                                                     metrics = ShowSampleMetricsName(object, type = "Metrics_10x"),
                                                     return.type = "interactive_table",
                                                     csv.name = "metrics_10x",
                                                     elementId ="metrics_10x-table" ,
                                                     table.subtitle = "10x Genomics metrics")
  }


  ##Metrics information
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment --- Checking RNA & ADT Metrics"))
  outlist$Metrics$note <- utils::capture.output(outlist$Metrics$count <- PlotSampleMetrics(object, type = "count",  return.type = "interactive_table",metrics= c("nCell", "nGene_RNA","nPro_ADT"),
                                                                                           csv.name = "GEX_metrics", elementId ="GEX_metrics-table",maxWidth = NULL))
  outlist$Metrics$note <- c(outlist$Metrics$note , utils::capture.output(outlist$Metrics$summary <- PlotSampleMetrics(object, type = "summary",  return.type = "interactive_table", metrics = ShowSampleMetricsName(object, type = "summary"), maxWidth = NULL)))

  ##metrics outlier
  temp <- FindSampleMetricsWarning(object, return.type =c("table", "interactive_table"), split.by = split.by, sample.by = sample.by )

  if(is.null(split.by)){
    outlist$outlier_metrics$table <- temp$table$list[ unlist(lapply(temp$table$list, function(x) {length(x)!=0} )) ]
  }else{
    result <- lapply(names(temp$table$list[[1]]),
                     function(field) {
                       Reduce(union, lapply(temp$table$list, `[[`, field))
                     })
    names(result) <- names(temp$table$list[[1]])
    result <- Filter(function(x) length(x) > 0, result)
    outlist$outlier_metrics$table <- result
  }

  outlist$outlier_metrics$interactive_table <- temp$interactive_table


  ##VDJ

  if( "receptor_subtype" %in% colnames(object@meta.data) ){
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment --- V(D)J"))
    outlist$VDJ$note <- utils::capture.output(outlist$VDJ$`V(D)J chain` <- PlotSampleVDJ(object, type="pct", return.type = "interactive_table"))
    outlist$VDJ$`V(D)J subtype` <- PlotSampleVDJ(object, type="subtype", return.type = "interactive_table")
  }

  if(!is.null(VDJ_data)){

    ##TCR
    if( length(VDJ_data) !=0 ){
      if( "TCR" %in% names(VDJ_data) ){
        scRepertoire_data <- PlotSampleVDJ(VDJ_data$TCR, color =NULL, type = "CDR3", return.scRepertoire = T)

        ##clonalQuant
        cat("-----TCR clonalQuant \n")
        outlist$VDJ$clonalQuant$TCR <- PlotSampleVDJ(VDJ_data$TCR, color = NULL, type = "clonalQuant", scRepertoire_data =scRepertoire_data , return.type = "interactive_table")

        ##
        cat("-----TCR geneUsage \n")
        outlist$VDJ$geneUsage$TCR <- PlotSampleVDJ(VDJ_data$TCR, color = NULL, type = "geneUsage", scRepertoire_data =scRepertoire_data , return.type = "interactive_table")

        ##CDR3
        cat("-----TCR CDR3 \n")
        if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Sample_QC/VDJ/CDR3"))) {
          dir.create(paste0(plot_out,"/SingleCellMQC/plot/Sample_QC/VDJ/CDR3"), recursive = TRUE)
        }


        outlist$VDJ$CDR3$TCR <- PlotSampleVDJ(VDJ_data$TCR, color = NULL, type = "CDR3", scRepertoire_data =scRepertoire_data , return.type = "plot")

        ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "TCR_CDR3", ".png"),
                        plot= outlist$VDJ$CDR3$TCR ,
                        width= 10,
                        height = 8, dpi=150)
        # suppressWarnings(htmlwidgets::saveWidget( plotly::partial_bundle(plotly::ggplotly(outlist$VDJ$CDR3$TCR)), paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "TCR_CDR3", ".html"), selfcontained = F) )
        suppressWarnings(htmlwidgets::saveWidget( (plotly::ggplotly(outlist$VDJ$CDR3$TCR)), paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "TCR_CDR3", ".html"), selfcontained = F) )

        outlist$VDJ$CDR3$TCR <- "TCR_CDR3"
      }
    }

    ##BCR
    if( length(VDJ_data) !=0 ){
      if( "BCR" %in% names(VDJ_data) ){
        scRepertoire_data <- PlotSampleVDJ(VDJ_data$BCR, color =NULL, type = "CDR3", return.scRepertoire = T)

        ##clonalQuant
        cat("-----BCR clonalQuant \n")
        outlist$VDJ$clonalQuant$BCR <- PlotSampleVDJ(VDJ_data$BCR, color = NULL, type = "clonalQuant", scRepertoire_data =scRepertoire_data , return.type = "interactive_table")

        ##geneUsage
        cat("-----BCR geneUsage \n")
        outlist$VDJ$geneUsage$BCR <- PlotSampleVDJ(VDJ_data$BCR, color = NULL, type = "geneUsage", scRepertoire_data =scRepertoire_data , return.type = "interactive_table")

        ##CDR3
        cat("-----BCR CDR3 \n")
        if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Sample_QC/VDJ/CDR3"))) {
          dir.create(paste0(plot_out,"/SingleCellMQC/plot/Sample_QC/VDJ/CDR3"), recursive = TRUE)
        }
        outlist$VDJ$CDR3$BCR <- PlotSampleVDJ(VDJ_data$BCR, color = NULL, type = "CDR3", scRepertoire_data =scRepertoire_data , return.type = "plot")


        ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "BCR_CDR3", ".png"),
                        plot= outlist$VDJ$CDR3$BCR ,
                        width= 10,
                        height = 8, dpi=150)
        # suppressWarnings(htmlwidgets::saveWidget( plotly::partial_bundle(plotly::ggplotly(outlist$VDJ$CDR3$BCR)), paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "BCR_CDR3", ".html"), selfcontained = F) )
        suppressWarnings(htmlwidgets::saveWidget( (plotly::ggplotly(outlist$VDJ$CDR3$BCR)), paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "BCR_CDR3", ".html"), selfcontained = F) )

        outlist$VDJ$CDR3$BCR <- "BCR_CDR3"

      }
    }
  }
  return(outlist)
}

.sampleSection2 <- function(object, plot_out=NULL, celltype.by=NULL, tissue=NULL, split.by=NULL, sample.by="orig.ident"){
  outlist <- list()
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Outlier sample detection by celltype"))


  if( is.null(celltype.by)){
    if("ScType" %in% colnames(object@meta.data)  ){
      celltype.by <- "ScType"
    }
  }

  if( !is.null(celltype.by)){
    if( !is.null(tissue) ){
      Common_result <- FindCommonPCTOutlier(object,sample.by=sample.by, tissue = tissue, return.type = c("table", "interactive_table"), celltype.by = celltype.by )
      outlist$outlier_comp$com_table <- Common_result$table$list[ unlist(lapply(Common_result$table$list, function(x) {length(x)!=0} )) ]
      outlist$outlier_comp$com_interactive_table <- Common_result$interactive_table
      outlist$outlier_comp$warning <- Common_result$table$warining
    }

    outlist$outlier_comp$CellTypePCT <- PlotSampleCellTypePCT(object,sample.by=sample.by,  celltype.by = celltype.by, return.type = "interactive_table", maxWidth = NULL)

    utils::capture.output(ellipse <- FindInterSamplePCTOutlier(object, sample.by=sample.by, return.type =c("table", "plot", "interactive_table"), celltype.by = celltype.by,  method = "ellipse", split.by = split.by))
    utils::capture.output(dbscan <- FindInterSamplePCTOutlier(object,sample.by=sample.by,  return.type =c("table", "plot"), celltype.by = celltype.by,  method = "dbscan", split.by = split.by))

      outlist$outlier_comp$pca_ellipse <- suppressWarnings(plotly::ggplotly(ellipse[["plot"]][["pca"]]+  ggplot2::geom_point(ggplot2::aes(shape = isOutlier, text=Sample), size = 3, alpha = 0.7)   ))
      outlist$outlier_comp$pca_dbscan <- suppressWarnings(plotly::ggplotly(dbscan[["plot"]][["pca"]]+  ggplot2::geom_point(ggplot2::aes(shape = isOutlier, text=Sample), size = 3, alpha = 0.7)))

      outlist$outlier_comp$interactive_table <- ellipse[["interactive_table"]][["contribution"]]
      outlist$outlier_comp$name_ellipse <- ellipse$table$outlier$Sample[ellipse$table$outlier$isOutlier]
      outlist$outlier_comp$name_dbscan <- dbscan$table$outlier$Sample[dbscan$table$outlier$isOutlier]

  }

  return(outlist)

}

.sampleSection3 <- function(object, sample.by=sample.by){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample identity verification "))
  outlist <- list()
  outlist$identity <- PlotSampleLabel(object, sample.by=sample.by, return.type = "interactive_table" )
  outlist$num <- PlotSampleLabel(object,sample.by=sample.by, return.type = "table")$Predict
  return(outlist)
}


report.cell <- function(object, sample.by="orig.ident", color=NULL,  plot_out=NULL, section= c(1:2)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Cell <- list()

  if(1 %in% section){
    Cell$Section1  <- .cellSection1(object,plot_out=plot_out, sample.by=sample.by)
  }

  if(2 %in% section){
    Cell$Section2  <-  .cellSection2(object,plot_out=plot_out, sample.by=sample.by)
  }


  return(Cell)
}

.cellSection1 <- function(object, plot_out=NULL,sample.by="orig.ident"){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Cell QC: Low quality cell information"))
  outlist <- list()

  #lq
  ##Filteration
  outlist$lq_Filtration <- PlotCellMethodFiltration(object, return.type = "interactive_table", type.detection = "lq")

  ##upset
  outlist$lq_upset <- PlotCellMethodUpset(object,  type.detection = "lq")
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/LowCell/upset"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/LowCell/upset"), recursive = TRUE)
  }
  pdf(paste0(plot_out, "/SingleCellMQC/plot/Cell_QC/LowCell/upset/", "low_upset_all", ".pdf"),
      width = 8,
      height = 5)
  print( outlist$lq_upset)
  dev.off()


  ##pre post
  outlist$lq_pre_post <- StatsCellFilter(object, return.type = "interactive_table", type.detection = "lq")

  #scater
  gene_mt <- PlotCellMetricsScatter(object, split.by = sample.by, metrics.by = c("nFeature_RNA", "percent.mt"), ggside = T, size = 1 , color = "black")
  gene_count <- PlotCellMetricsScatter(object, split.by = sample.by, metrics.by = c("nFeature_RNA", "nCount_RNA"), ggside = T, size = 1 , color = "black")
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/LowCell/scatter"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/LowCell/scatter"), recursive = TRUE)
  }

  temp <- lapply( names(gene_mt), function(y){
    p1 = patchwork::wrap_plots( list(gene_mt[[y]],gene_count[[y]]) , ncol = 2)
    ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Cell_QC/LowCell/scatter/", y, ".png"),
                    plot=p1,
                    width= 11,
                    height = 5.5, dpi=150)
    return(NULL)
  } )
  outlist$scatter <- names(gene_mt)
  return(outlist)

}
.cellSection2 <- function(object, plot_out=NULL,sample.by="orig.ident"){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Cell QC: Doublet information"))
  outlist <- list()

  #db
  ##Filteration

  outlist$db_Filtration <- PlotCellMethodFiltration(object, return.type = "interactive_table", type.detection = "db")

  ##upset
  metadata <- object@meta.data
  grep_col <- grep("^db_(?!.*score$)", colnames(metadata), perl = TRUE, value = TRUE)

  if(length(grep_col)>1 ){
    outlist$db_upset <- PlotCellMethodUpset(object,  type.detection = "db")
    if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/Doublet/upset"))) {
      dir.create(paste0(plot_out,"/SingleCellMQC/plot/Cell_QC/Doublet/upset"), recursive = TRUE)
    }
    pdf(paste0(plot_out, "/SingleCellMQC/plot/Cell_QC/Doublet/upset/", "db_upset_all", ".pdf"),
        width = 8,
        height = 5)
    print( outlist$db_upset)
    dev.off()
  }


  ##pre post
  outlist$db_pre_post <- StatsCellFilter(object, return.type = "interactive_table", type.detection = "db")


  return(outlist)
}



report.feature <- function(object, marker_name= c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                                  "CD14","FCGR3A","CD68","FCN1","ITGAX"),
                           celltype.by=NULL,
                           sample.by="orig.ident",
                           plot_out=NULL,
                           section= c(1:1)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Feature <- list()

  if(1 %in% section){
    Feature$Section1  <- .featureSection1(object, marker_name, celltype.by=celltype.by, sample.by=sample.by, plot_out=plot_out)
  }

  return(Feature)
}

.featureSection1 <- function(object, marker_name=c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                                   "CD14","FCGR3A","CD68","FCN1","ITGAX"),
                             celltype.by=NULL,
                             sample.by="orig.ident",
                             plot_out=NULL,
                             split.by="orig.ident"){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Feature QC: Feature Metrics"))

  outlist <- list()
  marker_name <-marker_name

  if(!is.null(celltype.by) & !is.null(marker_name) ){
    #marker
    outlist$markerPCT <-StatsFeaturePCT(object, feature=marker_name, group.by = celltype.by,
                                        split.by = sample.by, assay = "RNA",
                                        slot="counts",
                                        return.type="interactive_table")
  }

  ##plot
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Feature_QC/FeatureMetrics/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Feature_QC/FeatureMetrics/"), recursive = TRUE)
  }

  out <- CalculateMetricsPerFeature(object, sample.by = sample.by, add.Seurat = F)
  for_name = names(out$RNA)

  temp <- lapply(for_name, function(x){
    p1 <- PlotFeatureMetrics(out, metric ="pct", assay ="RNA", sample =x)
    p2 <- PlotFeatureMetricsScatter(out, metrics.by = c("mean_lognorm", "variance_lognorm"), assay = "RNA",sample = x, ggside = F )+
      ggplot2::labs(title = "Mean-variance plot")
    suppressWarnings(p2 <- AddFeaturePlotLabel(p2, metric_order_col = "mean_lognorm"))
    suppressWarnings(p2 <- AddFeaturePlotLabel(p2, metric_order_col = "variance_lognorm"))

    bottom_row <- patchwork::wrap_plots(p2, patchwork::plot_spacer(), nrow = 1, widths = c(3, 1))
    final_p3 <- patchwork::wrap_plots(p1, bottom_row, ncol = 1, heights = c(1, 1.8))
    suppressWarnings(ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Feature_QC/FeatureMetrics/", x, "_RNA.png"),
                                     plot = final_p3,
                                     width= 10,
                                     height = 10,
                                     dpi=150 ))

    # if object have ADT assay
    if( "ADT" %in% names(object@assays) ){
      p1 <- PlotFeatureMetrics(out, metric ="pct", assay ="ADT", sample =x)
      p2 <- PlotFeatureMetricsScatter(out, metrics.by = c("mean_clr", "variance_clr"), assay = "ADT", sample = x, ggside = F )+
        ggplot2::labs(title = "Mean-variance plot")
      suppressWarnings(p2 <- AddFeaturePlotLabel(p2, metric_order_col = "mean_clr"))
      suppressWarnings(p2 <- AddFeaturePlotLabel(p2, metric_order_col = "variance_clr"))
      bottom_row <- patchwork::wrap_plots(p2, patchwork::plot_spacer(), nrow = 1, widths = c(3, 1))
      final_p3 <- patchwork::wrap_plots(p1, bottom_row, ncol = 1, heights = c(1, 1.8))
      suppressWarnings(ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Feature_QC/FeatureMetrics/", x, "_ADT.png"),
                                       plot = final_p3,
                                       width= 10,
                                       height = 10,
                                       dpi=150 ))
    }

  })
  outlist$base <- for_name
  if( "ADT" %in% names(object@assays) ){
    outlist$assay <- c("RNA", "ADT")
  }else{
    outlist$assay <- c("RNA")
  }
  return(outlist)

}





report.batch <- function(object, RNA_cluster_name=NULL, ADT_cluster_name=NULL,
                         sample.by="orig.ident",
                         RNA.batch.by= "orig.ident",
                         ADT.batch.by="orig.ident",
                         celltype.by=NULL,
                         RNA.covariate.other=NULL,
                         ADT.covariate.other=NULL,
                         RNA.covariate.formula=NULL,
                         ADT.covariate.formula=NULL,
                         do.covariate.split.celltype=F,
                         plot_out, section=1:3){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Batch <- list()

  if(1 %in% section){
    Batch$Section1  <- .batchSection1(object,
                                      RNA.batch.by=RNA.batch.by,
                                      ADT.batch.by=ADT.batch.by,
                                      plot_out=plot_out,
                                      sample.by=sample.by,
                                      RNA.covariate.other=RNA.covariate.other,
                                      ADT.covariate.other=ADT.covariate.other,
                                      RNA.covariate.formula=RNA.covariate.formula,
                                      ADT.covariate.formula=ADT.covariate.formula,
                                      celltype.by=celltype.by,
                                      do.covariate.split.celltype=do.covariate.split.celltype
    )
  }

  if(2 %in% section){
    Batch$Section2  <- .batchSection2(object, RNA_cluster_name=RNA_cluster_name,
                                      ADT_cluster_name=ADT_cluster_name,
                                      RNA.batch.by=RNA.batch.by,
                                      ADT.batch.by=ADT.batch.by,
                                      plot_out=plot_out,
                                      sample.by=sample.by)
  }

  if(3 %in% section){
    Batch$Section3  <- .batchSection3(object,
                                      celltype.by=celltype.by,
                                      RNA.batch.by=RNA.batch.by,
                                      ADT.batch.by=ADT.batch.by,
                                      sample.by=sample.by)
  }


  return(Batch)
}

.batchSection1 <- function(object,
                           sample.by="orig.ident",
                           RNA.batch.by= "orig.ident",
                           ADT.batch.by="orig.ident", plot_out,
                           RNA.covariate.other=NULL,
                           ADT.covariate.other=NULL,
                           RNA.covariate.formula =NULL,
                           ADT.covariate.formula=NULL,
                           do.covariate.split.celltype=F,
                           celltype.by=NULL){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Sample-level"))
  metadata <- getMetaData(object)
  outlist <- list()

    if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/pseudobulk/"))) {
      dir.create(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/pseudobulk/"), recursive = TRUE)
    }

    if( checkAssay(object, "RNA") ){
      message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    "------ Run Batch QC: Pseudobulk PCA (RNA)"))
      index = union(RNA.batch.by, RNA.covariate.other)
      PseudobulkRNA <- RunPseudobulkData(object, assay = "RNA", slot = "counts",
                                         sample.by = sample.by,
                                         other_cols_contain =  index)
      pseudobulk_metadata <- PseudobulkRNA[["pseudobulk_metadata"]]
      suppressWarnings(suppressMessages(p1_base <- PlotReducedDim(PseudobulkRNA, group.by =sample.by, size = 4,
                                                               pca.label = sample.by)))

      if(!is.null(RNA.covariate.formula)){
        suppressMessages( varpart_pca<- RunVarPartPseudobulkPCA(PseudobulkRNA, formula = RNA.covariate.formula))
        p_var <- PlotVarPartStackBar(varpart_pca)
      }


      temp <- lapply(index, function(x){
        # 每次循环都从原始的 p1_base 开始，避免前一个循环的颜色设置影响下一个
        p1 <- p1_base
        tempdata <- p1$data
        tempdata[[ x ]] <- pseudobulk_metadata[[x]] # 直接将 x 列添加到 tempdata
        p1$data <- tempdata
        # 判断 pseudobulk_metadata[[x]] 的数据类型
        if (is.numeric(pseudobulk_metadata[[x]])) {
          # 如果是数值型，使用连续渐变色
          suppressMessages(p1 <- p1 +ggplot2::guides(color = "none", fill = "none")+
                             ggplot2::aes_string(fill = x, color = x) +
                             ggplot2::scale_fill_viridis_c(option = "D", direction = -1) +
                             ggplot2::scale_color_viridis_c(option = "D", direction = -1) +
                             ggplot2::labs(fill = x, color = x)+
                             ggplot2::guides(fill = ggplot2::guide_colorbar(title = x),
                                             color = ggplot2::guide_colorbar(title = x))
      )
        } else {
          # 如果是非数值型，使用离散色板
          if(x == sample.by){
            suppressMessages(p1 <- p1 +
                               ggplot2::aes_string(fill = x, color = x) +
                               ggplot2::scale_fill_manual(values=get_colors( length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::scale_color_manual(values=get_colors(length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::labs(fill = x, color = x)+
                               ggplot2::theme(legend.position = "none"))
          }else{
            suppressMessages(p1 <- p1 +
                               ggplot2::aes_string(fill = x, color = x) +
                               ggplot2::scale_fill_manual(values=get_colors( length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::scale_color_manual(values=get_colors(length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::labs(fill = x, color = x))
          }
        }
        suppressWarnings(htmlwidgets::saveWidget( (plotly::ggplotly(p1)),
                                                  paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_RNA.html"),
                                                  selfcontained = F) )
        if(!is.null(RNA.covariate.formula)){
          p2 <- patchwork::wrap_plots(list(p1, p_var), ncol=1, heights = c(3,0.5))
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_RNA.png"),
                          plot=p2,
                          width=10,
                          height = 8, dpi=150)
        }else{
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_RNA.png"),
                          plot=p1,
                          width=10,
                          height = 8, dpi=150)
        }

      })
      outlist$pseudobulk_RNA <-paste0(index, "_RNA")


      if(!is.null(RNA.covariate.formula)){
        message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      "------ Run Batch QC: Covariate Impact (RNA)"))
        varpart_pca<- RunVarPartPseudobulk(PseudobulkRNA, formula = RNA.covariate.formula)
        if( !is.null(varpart_pca) ){
          p_var <- PlotVarPartVln(varpart_pca)+ggplot2::labs(title = "All cells")
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/VarPart_all_RNA.png"),
                          plot=p_var,
                          width=10,
                          height = 8, dpi=150)
          outlist$pseudobulk_RNA_VP = "VarPart_all_RNA.png"
        }


        if(!is.null(celltype.by) & do.covariate.split.celltype ){
          message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        "------ Run Batch QC: Covariate Impact per celltype (RNA)"))
          PseudobulkRNA <- RunPseudobulkData(object, assay = "RNA", slot = "counts", cluster.by = celltype.by,
                                             sample.by = sample.by,
                                             other_cols_contain =  index)
          varpart_pca<- RunVarPartPseudobulk(PseudobulkRNA, formula = RNA.covariate.formula)
          if( !is.null(varpart_pca) ){
            p_var <- PlotVarPartVln(varpart_pca, do.split = T)
            ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/VarPart_celltype_RNA.png"),
                            plot=p_var,
                            width=15,
                            height = 15, dpi=150)
            outlist$pseudobulk_RNA_VP_celltype= "VarPart_celltype_RNA.png"
          }
        }

      }
    }

    if( checkAssay(object, "ADT") ){
      message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    "------ Run Batch QC: Pseudobulk PCA (ADT)"))
      index = union(ADT.batch.by, ADT.covariate.other)
      PseudobulkADT <- RunPseudobulkData(object, assay = "ADT", slot = "counts",
                                         sample.by = sample.by,
                                         other_cols_contain =  index)
      pseudobulk_metadata <- PseudobulkADT[["pseudobulk_metadata"]]
      suppressWarnings(suppressMessages(p1_base <- PlotReducedDim(PseudobulkADT, group.by =sample.by, size = 4,
                                                                  pca.label = sample.by)))

      if(!is.null(ADT.covariate.formula)){
        suppressMessages(varpart_pca<- RunVarPartPseudobulkPCA(PseudobulkADT, formula = ADT.covariate.formula))
        p_var <- PlotVarPartStackBar(varpart_pca)
      }


      temp <- lapply(index, function(x){
        # 每次循环都从原始的 p1_base 开始，避免前一个循环的颜色设置影响下一个
        p1 <- p1_base
        tempdata <- p1$data
        tempdata[[ x ]] <- pseudobulk_metadata[[x]] # 直接将 x 列添加到 tempdata
        p1$data <- tempdata
        # 判断 pseudobulk_metadata[[x]] 的数据类型
        if (is.numeric(pseudobulk_metadata[[x]])) {
          # 如果是数值型，使用连续渐变色
          suppressMessages(p1 <- p1 +ggplot2::guides(color = "none", fill = "none")+
                             ggplot2::aes_string(fill = x, color = x) +
                             ggplot2::scale_fill_viridis_c(option = "D", direction = -1) +
                             ggplot2::scale_color_viridis_c(option = "D", direction = -1) +
                             ggplot2::labs(fill = x, color = x)+
                             ggplot2::guides(fill = ggplot2::guide_colorbar(title = x),
                                             color = ggplot2::guide_colorbar(title = x))
          )
        } else {
          # 如果是非数值型，使用离散色板
          if(x == sample.by){
            suppressMessages(p1 <- p1 +
                               ggplot2::aes_string(fill = x, color = x) +
                               ggplot2::scale_fill_manual(values=get_colors( length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::scale_color_manual(values=get_colors(length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::labs(fill = x, color = x)+
                               ggplot2::theme(legend.position = "none"))
          }else{
            suppressMessages(p1 <- p1 +
                               ggplot2::aes_string(fill = x, color = x) +
                               ggplot2::scale_fill_manual(values=get_colors( length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::scale_color_manual(values=get_colors(length(unique(pseudobulk_metadata[[x]]))  ))+
                               ggplot2::labs(fill = x, color = x))
          }
        }
        suppressWarnings(htmlwidgets::saveWidget( (plotly::ggplotly(p1)),
                                                  paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_ADT.html"),
                                                  selfcontained = F) )
        if(!is.null(ADT.covariate.formula)){
          p2 <- patchwork::wrap_plots(list(p1, p_var), ncol=1, heights = c(3,0.5))
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_ADT.png"),
                          plot=p2,
                          width=10,
                          height = 8, dpi=150)
        }else{
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/", x, "_ADT.png"),
                          plot=p1,
                          width=10,
                          height = 8, dpi=150)
        }

      })
      outlist$pseudobulk_ADT <-paste0(index, "_ADT")


      if(!is.null(ADT.covariate.formula)){
        message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      "------ Run Batch QC: Covariate Impact (ADT)"))
        varpart_pca<- RunVarPartPseudobulk(PseudobulkADT, formula = ADT.covariate.formula)
        if( !is.null(varpart_pca) ){
          p_var <- PlotVarPartVln(varpart_pca)+ggplot2::labs(title = "All cells")
          ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/VarPart_all_ADT.png"),
                          plot=p_var,
                          width=10,
                          height = 8, dpi=150)
          outlist$pseudobulk_ADT_VP = "VarPart_all_ADT.png"
        }


        if(!is.null(celltype.by) & do.covariate.split.celltype){
          message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        "------ Run Batch QC: Covariate Impact per celltype (ADT)"))
          PseudobulkADT <- RunPseudobulkData(object, assay = "ADT", slot = "counts", cluster.by = celltype.by,
                                             sample.by = sample.by,
                                             other_cols_contain =  index)
          varpart_pca<- RunVarPartPseudobulk(PseudobulkADT, formula = ADT.covariate.formula)
          if( !is.null(varpart_pca) ){
            p_var <- PlotVarPartVln(varpart_pca, do.split = T)
            ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/pseudobulk/VarPart_celltype_ADT.png"),
                            plot=p_var,
                            width=15,
                            height = 15, dpi=150)
            outlist$pseudobulk_ADT_VP_celltype="VarPart_celltype_ADT.png"
          }
        }

      }
    }

  return(outlist)
}



.batchSection2 <- function(object, RNA_cluster_name=NULL, ADT_cluster_name=NULL,
                           sample.by="orig.ident",
                           RNA.batch.by= "orig.ident",
                           ADT.batch.by="orig.ident", plot_out){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Cell-level"))
  metadata <- object@meta.data
  outlist <- list()
  ## batch test
  if(!is.null(RNA_cluster_name)){
    if(RNA_cluster_name %in% colnames(object@meta.data) ){
      outlist$test_RNA <- utils::capture.output(temp <- RunBatchTest(object, cluster.by = RNA_cluster_name),  type = "message")
    }
  }

  if(!is.null(ADT_cluster_name)){
    if(ADT_cluster_name %in% colnames(object@meta.data)){
      outlist$test_ADT <- utils::capture.output(temp <- RunBatchTest(object, cluster.by = ADT_cluster_name),  type = "message")
    }
  }

  ##cell
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Cell clustering"))
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/cell/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/cell/"), recursive = TRUE)
  }
  batch_cellplot <- function(object, batch.by, name, reduction, cluster){
    labels <- unique(object@meta.data[[batch.by]])
    groups <- ceiling(seq_along(labels) / 50)
    grouped_data <- split(labels, groups)
    names(grouped_data) <- seq_along(grouped_data)
    cluster_data <- SeuratObject::Embeddings(object = object, reduction= reduction)
    index <- union(cluster, batch.by)
    cluster_data <- data.frame(cluster_data, object@meta.data[, index, drop=F])
    cluster_data[[cluster]] <- as.character(cluster_data[[cluster]])

    plotbatch_list <- lapply(names(grouped_data), function(x){
      cluster_data <- cluster_data[cluster_data[[batch.by]] %in% grouped_data[[x]],]
      color <- get_colors( length(unique(cluster_data[[cluster]])))
      p <- plotScatter(cluster_data, x= colnames(cluster_data)[1] , y= colnames(cluster_data)[2], group.by = cluster, color = color,size=0.5,
                       log.x = F,log.y = F,ggside = F,split.by = batch.by, ncol = 5,guide.nrow=10, raster.cutoff=100000)+ggplot2::theme(legend.position = "none")

      p[["facet"]][["params"]][["ncol"]] <- 5
      p[["facet"]][["params"]][["nrow"]] <- ceiling(length(grouped_data[[x]]) / 5)

      ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/cell/",name, "_",batch.by,"_", x,".png"),
                      plot=p,
                      width=15,
                      height = 0.5+3*( ceiling(length(grouped_data[[x]])/5) ), dpi=100)
    })
    return( paste0( name,"_", batch.by, "_", names(grouped_data))  )
  }

  if( "RNA" %in% names(object@assays) && !is.null(RNA_cluster_name) ){
    if(RNA_cluster_name %in% colnames(object@meta.data) && "rna.umap" %in% Seurat::Reductions(object)  ){
      SeuratObject::DefaultAssay(object) <- "RNA"
      temp <- lapply(RNA.batch.by, function(y){
        batch_cellplot(object, batch.by = y, name = "RNA", reduction = "rna.umap", cluster = RNA_cluster_name)
      })
      names(temp) <- paste0("RNA_",RNA.batch.by)
      outlist$cell_RNA <- temp

    }
  }

  if( "ADT" %in% names(object@assays) && !is.null(ADT_cluster_name) ){
    if(ADT_cluster_name %in% colnames(object@meta.data) && "adt.umap" %in% Seurat::Reductions(object)  ){
      SeuratObject::DefaultAssay(object) <- "ADT"
      temp <- lapply(ADT.batch.by, function(y){
        batch_cellplot(object, batch.by = y, name = "ADT", reduction = "adt.umap", cluster = ADT_cluster_name)
      })
      names(temp) <- paste0("ADT_",ADT.batch.by)
      outlist$cell_ADT <- temp
    }
  }
  return(outlist)
}


.batchSection3 <- function(object,
                           sample.by="orig.ident",
                           RNA.batch.by= "orig.ident",
                           ADT.batch.by="orig.ident",
                           celltype.by=NULL){
  outlist <- list()
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Feature-level"))
  if( !is.null(sample.by) & !is.null(RNA.batch.by) & !is.null(celltype.by) ){
    GTE_table <- lapply(RNA.batch.by, function(x){
      GTE_out <- RunBatchGTE(object, assay = "RNA", batch.by = x, group.by = celltype.by)
      rownames(GTE_out) <- NULL
      GTE_out_50 <- GTE_out[order(-GTE_out$GTE_Overall)[1:50], 1:2]
      .re_table(GTE_out_50, maxWidth = NULL, csv.name = paste0("GTE_RNA_", x), elementId = paste0("GTE_RNA_", x),subtitle = "Top 50" )
    })
    names(GTE_table) <- RNA.batch.by
    outlist$RNA = GTE_table
  }

  if( !is.null(sample.by) & !is.null(ADT.batch.by) & !is.null(celltype.by) ){
    GTE_table <- lapply(ADT.batch.by, function(x){
      GTE_out <- RunBatchGTE(object, assay = "ADT", batch.by = x, group.by = celltype.by)
      rownames(GTE_out) <- NULL
      GTE_out_50 <- GTE_out[order(-GTE_out$GTE_Overall)[1:50], 1:2]
      .re_table(GTE_out_50, maxWidth = NULL, csv.name = paste0("GTE_ADT_", x), elementId = paste0("GTE_ADT_", x),subtitle = "Top 50"  )
    })
    names(GTE_table) <- ADT.batch.by
    outlist$ADT = GTE_table

  }
  return(outlist)

}



#' @title Run Preprocessing Pipeline for Single-Cell Data
#'
#' @description
#' This function runs a comprehensive preprocessing pipeline for single-cell RNA-seq data. It includes steps for calculating quality metrics, cell type annotation, doublet detection, low-quality cell filtering, and dimensionality reduction (e.g., UMAP). The function is designed to work with Seurat objects and supports customization of doublet detection and low-quality cell filtering methods.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param sample.by A character string specifying the column in the metadata that identifies different samples or batches. Default is `"orig.ident"`.
#' @param db.method A character vector specifying the method(s) to use for doublet detection. Options include `"scDblFinder"`. Default is `"scDblFinder"`.
#' @param lq.method A character vector specifying the method(s) to use for low-quality cell filtering. Options include `"MAD"`. Default is `"MAD"`.
#' @param preprocess The method of preprocessing and clustering to use. Options include:
#' \itemize{
#' \item "rna.pca": Perform the analysis to PCA step on RNA assay.
#' \item "rna.umap": Perform the analysis to UMAP step on RNA assay.
#' \item "rna.umap.har": Perform the analysis to UMAP step on RNA assay with Harmony integration.
#' \item "adt.pca": Perform the analysis to PCA step on ADT assay.
#' \item "adt.umap": Perform the analysis to UMAP step on ADT assay.
#' \item "adt.umap.har": Perform the analysis to UMAP step on ADT assay with Harmony integration.
#' \item "wnn.umap": Perform the analysis to UMAP step using weighted nearest neighbor (WNN) method on RNA+ADT assays.
#' \item "wnn.umap.har": Perform the analysis to UMAP step using weighted nearest neighbor (WNN) method on RNA+ADT assays with Harmony integration.
#' }
#' The default setting: "rna.umap".
#' @param ... Additional arguments to be passed to \code{`RunPipeline`} function.
#'
#' @return
#' A Seurat object with the following updates:
#' - Quality metrics added to the metadata.
#' - Cell type annotations added to the metadata.
#' - Doublet scores and classifications added to the metadata.
#' - Low-quality cell classifications added to the metadata.
#' - Dimensionality reduction results (e.g., UMAP) added to the object.
#'
#' @details
#' The function performs the following steps:
#' 1. **Calculate Quality Metrics**: Adds quality metrics (e.g., number of genes, mitochondrial percentage) to the metadata.
#' 2. **Cell Type Annotation**: Runs cell type annotation using the `RunScType` function.
#' 3. **Doublet Detection**: Detects doublets using the specified method(s) (e.g., `scDblFinder`).
#' 4. **Low-Quality Cell Detection**: Filters low-quality cells using the specified method(s) (e.g., `MAD`).
#' 5. **Dimensionality Reduction**: Runs a preprocessing pipeline (e.g., UMAP) for visualization and downstream analysis.
#'
#' @export
RunPreprocess <- function(object, sample.by="orig.ident", db.method=c("scDblFinder", "hybrid"), lq.method=c("MAD") , preprocess = "rna.umap", ... ){
  object <- CalculateMetrics(object, sample.by = sample.by)
  object <- RunScType(object, split.by = sample.by)
  object <- RunDbt(object, methods = db.method, add.Seurat = T,split.by = sample.by)
  object <- RunLQ(object, methods = lq.method, add.Seurat = T,split.by = sample.by)
  object <- RunPipeline(object, preprocess = preprocess, har.group.by = sample.by, ...)
  return(object)
}
