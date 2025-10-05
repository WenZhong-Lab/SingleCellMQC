RunReport1 <- function(object=NULL, VDJ_data=NULL, sample.by="orig.ident", outputFile="./SingleCellMQC/SingleCellMQC.html",
                      section.sample = c(1:3), section.cell = c(1:2), section.feature=1,section.batch=c(1:3),
                      RNA_cluster_name="rna_cluster",
                      ADT_cluster_name="adt_cluster",
                      RNA.batch.by= "orig.ident",
                      ADT.batch.by="orig.ident",
                      RNA.batch.covariate=c( "orig.ident" ,"percent.mt", "nFeature_RNA", "nCount_RNA"),
                      ADT.batch.covariate=c( "orig.ident" , "nFeature_ADT", "nCount_ADT"),
                      celltype.by="ScType",
                      tissue=NULL
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

  if( !is.null(section.sample) ){
    out_list$Sample <- report.sample(object=object,
                                     VDJ_data=VDJ_data,
                                     sample.by=sample.by,
                                     color=NULL,
                                     plot_out=plot_out,
                                     section= section.sample,
                                     celltype.by = celltype.by,
                                     tissue=tissue)
  }

  print(outputFile)
  print(output_file)
  print(plot_dir)


  # quarto::quarto_render(
  #   input = paste0(outputFile, "/0_main.qmd"), # .qmd
  #   output_file = "SingleCellMQC.html",
  #   execute_params = list(data = NULL, cell_dir = plot_dir)
  # )
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
                          tissue,celltype.by=NULL,
                          split.by=NULL,
                          section= c(1:3)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Sample <- list()
  SCMQC <- GetSingleCellMQCData(object)
  labels <- unique(object@meta.data[, sample.by, drop=T])

  if(1 %in% section){
    Sample$Section1  <- .sampleSection1(object, VDJ_data=VDJ_data, plot_out=plot_out,sample.by=sample.by, split.by=split.by)
  }

  if(2 %in% section){
    Sample$Section2  <- .sampleSection2(object,  plot_out=plot_out, celltype.by = celltype.by, tissue=tissue,sample.by=sample.by, split.by=split.by)
  }

  if(3 %in% section){
    Sample$Section3  <- .sampleSection3(object,sample.by=sample.by)
  }
  return(Sample)
}

.sampleSection1 <- function(object, VDJ_data=NULL, plot_out=NULL, split.by=NULL, sample.by="orig.ident"){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment"))
  outlist <- list()

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
        pdf(paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "TCR_CDR3", ".pdf"),
            width = 10,
            height = 8)
        print( outlist$VDJ$CDR3$TCR )
        dev.off()

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

        pdf(paste0(plot_out, "/SingleCellMQC/plot/Sample_QC/VDJ/CDR3/", "BCR_CDR3", ".pdf"),
            width = 10,
            height = 8)
        print( outlist$VDJ$CDR3$BCR )
        dev.off()

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
