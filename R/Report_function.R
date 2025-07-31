report.sample <- function(object, VDJ_data=NULL, sample.by="orig.ident", color=NULL,  plot_out=NULL,tissue,celltype.by=NULL, section= c(1:3)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Sample <- list()
  SCMQC <- GetSingleCellMQCData(object)
  labels <- unique(object@meta.data[, sample.by, drop=T])

  if(1 %in% section){
    Sample$Section1  <- .sampleSection1(object, VDJ_data=VDJ_data, plot_out=plot_out)
  }


  if(2 %in% section){
    Sample$Section2  <- .sampleSection2(object,  plot_out=plot_out, celltype.by = celltype.by, tissue=tissue)
  }

  if(3 %in% section){
    Sample$Section3  <- .sampleSection3(object)
  }
  return(Sample)
}

.sampleSection1 <- function(object, VDJ_data=NULL, plot_out=NULL, split.by=NULL){
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

    outlist$Metrics$Metrics_10X <- PlotSampleMetrics(object, type = "Metrics_10x", metrics = ShowSampleMetricsName(object, type = "Metrics_10x"),
                                                     return.type = "interactive_table", csv.name = "metrics_10x", elementId ="metrics_10x-table" , table.subtitle = "10x Genomics metrics")
  }


  ##Metrics information
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample quality assessment --- Checking RNA & ADT Metrics"))
  outlist$Metrics$note <- utils::capture.output(outlist$Metrics$count <- PlotSampleMetrics(object, type = "count",  return.type = "interactive_table",metrics= c("nCell", "nGene_RNA","nPro_ADT"),
                                                                                           csv.name = "GEX_metrics", elementId ="GEX_metrics-table",maxWidth = NULL))
  outlist$Metrics$note <- c(outlist$Metrics$note , utils::capture.output(outlist$Metrics$summary <- PlotSampleMetrics(object, type = "summary",  return.type = "interactive_table", metrics = ShowSampleMetricsName(object, type = "summary"), maxWidth = NULL)))

  ##metrics outlier
  temp <- FindSampleMetricsWarning(object, return.type =c("table", "interactive_table"), split.by = split.by)

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

.sampleSection2 <- function(object, plot_out=NULL, celltype.by=NULL, tissue=NULL, split.by=NULL){
  outlist <- list()
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Outlier sample detection by celltype"))


  if( is.null(celltype.by)){
    if("ScType" %in% colnames(object@meta.data)  ){
      celltype.by <- "ScType"
    }
  }

  if( !is.null(celltype.by)){
    if( !is.null(tissue) ){
      Common_result <- FindCommonPCTOutlier(object, tissue = tissue, return.type = c("table", "interactive_table"), celltype.by = celltype.by )
      outlist$outlier_comp$com_table <- Common_result$table$list[ unlist(lapply(Common_result$table$list, function(x) {length(x)!=0} )) ]
      outlist$outlier_comp$com_interactive_table <- Common_result$interactive_table
    }

    outlist$outlier_comp$CellTypePCT <- PlotSampleCellTypePCT(object, celltype.by = celltype.by, return.type = "interactive_table", maxWidth = NULL)

    utils::capture.output(ellipse <- FindInterSamplePCTOutlier(object, return.type =c("table", "plot", "interactive_table"), celltype.by = celltype.by,  method = "ellipse", split.by = split.by))
    utils::capture.output(dbscan <- FindInterSamplePCTOutlier(object, return.type =c("table", "plot"), celltype.by = celltype.by,  method = "dbscan", split.by = split.by))

    if(!is.na(ellipse)){
      outlist$outlier_comp$pca_ellipse <- suppressWarnings(plotly::ggplotly(ellipse[["plot"]][["pca"]]+  ggplot2::geom_point(ggplot2::aes(shape = isOutlier, text=Sample), size = 3, alpha = 0.7)   ))
      outlist$outlier_comp$pca_dbscan <- suppressWarnings(plotly::ggplotly(dbscan[["plot"]][["pca"]]+  ggplot2::geom_point(ggplot2::aes(shape = isOutlier, text=Sample), size = 3, alpha = 0.7)))

      outlist$outlier_comp$interactive_table <- ellipse[["interactive_table"]][["contribution"]]
      outlist$outlier_comp$name_ellipse <- ellipse$table$outlier$Sample[ellipse$table$outlier$isOutlier]
      outlist$outlier_comp$name_dbscan <- dbscan$table$outlier$Sample[dbscan$table$outlier$isOutlier]

    }

  }

  return(outlist)

}

.sampleSection3 <- function(object){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),  "------ Run Sample QC: Sample identity verification "))
  outlist <- list()
  outlist$identity <- PlotSampleLabel(object, return.type = "interactive_table" )
  outlist$num <- PlotSampleLabel(object, return.type = "table")$Predict
  return(outlist)
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


.featureSection1 <- function(object, marker_name=c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                                   "CD14","FCGR3A","CD68","FCN1","ITGAX"), celltype.by=NULL, sample.by="orig.ident", plot_out=NULL,split.by="orig.ident"){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Feature QC: Feature quality assessment"))

  outlist <- list()
  marker_name <-marker_name
  if( is.null(celltype.by)){
    if("ScType" %in% colnames(object@meta.data)  ){
      celltype.by <- "ScType"
    }
  }
  #marker
  outlist$markerPCT <-StatsFeaturePCT(object, feature=marker_name, group.by = celltype.by, split.by = sample.by, assay = "RNA", slot="counts", return.type="interactive_table")

  ##plot
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Feature_QC/FeatureMetrics/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Feature_QC/FeatureMetrics/"), recursive = TRUE)
  }
  for_name <- unique(object@meta.data[[sample.by]])
  out <- CalculateMetricsPerFeature(object, add.Seurat = F)

  temp <- lapply(for_name, function(x){
    p1 <- PlotFeatureMetrics(out, type="pct", assay ="RNA", sample =x)
    p2 <- PlotFeatureMetrics(out, type="mean", assay ="RNA", sample =x)
    p3 <- PlotFeatureMetrics(out, type="variance.standardized", assay ="RNA", sample =x)

    p4 <- patchwork::wrap_plots( list(p1,p2,p3) , ncol = 1)
    ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Feature_QC/FeatureMetrics/", x, "_RNA.png"),
                    plot=p4,
                    width= 10,
                    height = 12, dpi=150)
    # if object have ADT assay
    if( "ADT" %in% names(object@assays) ){
      p5 <- PlotFeatureMetrics(out, type="pct", assay ="ADT", sample =x)
      p6 <- PlotFeatureMetrics(out, type="mean", assay ="ADT", sample =x)
      p7 <- PlotFeatureMetrics(out, type="variance.standardized", assay ="ADT", sample =x)
      p8 <- patchwork::wrap_plots( list(p5,p6,p7) , ncol = 1)
      ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Feature_QC/FeatureMetrics/", x, "_ADT.png"),
                      plot=p8,
                      width= 10,
                      height = 12, dpi=150)
    }

  })
  outlist$base <- for_name
  if( "ADT" %in% names(object@assays) ){
    outlist$assay <- c("RNA", "ADT")
  }else{
    outlist$assay <- c("RNA")
  }

  ##sample
  if(!is.null(split.by)){
    ##RNA
    if( "RNA" %in% names(object@assays) ){
      SeuratObject::DefaultAssay(object) <- "RNA"
      out <- RunVarExplained(object, assay = "RNA", variables = split.by)
      out <- data.frame(Feature= rownames(out), out)
      top_20 <- out[order(out[,2], decreasing = TRUE), ][1:20, ]
      rownames(top_20) <-NULL
      outlist$VE_RNA <- .re_table(top_20, maxWidth = NULL, csv.name = "VarExplained_RNA", elementId = "VarExplained_RNA", subtitle = "Top 20 features variance explained in RNA assay")
    }

    ##ADT
    if( "ADT" %in% names(object@assays) ){
      SeuratObject::DefaultAssay(object) <- "ADT"
      out <- RunVarExplained(object, assay = "ADT", variables = split.by)
      out <- data.frame(Feature= rownames(out), out)
      top_20 <- out[order(out[,2], decreasing = TRUE), ][1:20, ]
      rownames(top_20) <-NULL
      outlist$VE_ADT <- .re_table(top_20, maxWidth = NULL, csv.name = "VarExplained_ADT", elementId = "VarExplained_ADT", subtitle = "Top 20 features variance explained in ADT assay")
    }
  }


  return(outlist)

}




report.feature <- function(object, marker_name= c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                                  "CD14","FCGR3A","CD68","FCN1","ITGAX"), celltype.by=NULL, sample.by="orig.ident", plot_out=NULL, section= c(1:1)){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Feature <- list()

  if(1 %in% section){
    Feature$Section1  <- .featureSection1(object, marker_name, celltype.by=celltype.by, sample.by=sample.by, plot_out=plot_out)
  }

  return(Feature)
}


.batchSection1 <- function(object,
                           sample.by="orig.ident",
                           RNA.batch.by= "orig.ident",
                           ADT.batch.by="orig.ident", plot_out,
                           RNA.batch.covariate=c( "orig.ident" ,"percent.mt", "nFeature_RNA", "nCount_RNA"),
                           ADT.batch.covariate=c( "orig.ident" , "nFeature_ADT", "nCount_ADT"),
                           celltype.by=NULL){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Sample-level"))
  metadata <- object@meta.data
  outlist <- list()
  {
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Pseudobulk PCA"))

    if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/"))) {
      dir.create(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/"), recursive = TRUE)
    }

    if( "RNA" %in% names(object@assays) ){
      SeuratObject::DefaultAssay(object) <- "RNA"
      suppressMessages(p1 <- PlotReducedDim(object , plot.type = "pseudobulk", assay="RNA", group.by.pseudobulk = RNA.batch.by[1], sample.by.pseudobulk = sample.by))
      temp <- lapply(RNA.batch.by, function(x){
        tempdata <- p1$data
        tempdata[[ RNA.batch.by[1] ]] <- object@meta.data[[x]][match(rownames(tempdata), object@meta.data[[sample.by]])]
        p1$data <- tempdata
        if(x==sample.by){
          suppressMessages(p1 <- p1 + ggplot2::scale_fill_manual(values=get_colors( length(unique(metadata[[x]]))  ))+ ggplot2::scale_color_manual(values=get_colors(length(unique(metadata[[x]]))  ))+
                             ggplot2::labs(fill = x, color = x)+ ggplot2::theme(legend.position = "none"))
        }else{
          suppressMessages(p1 <- p1 + ggplot2::scale_fill_manual(values=get_colors( length(unique(metadata[[x]]))  ))+ ggplot2::scale_color_manual(values=get_colors(length(unique(metadata[[x]]))  ))+
                             ggplot2::labs(fill = x, color = x))
        }

        ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_RNA.png"),
                        plot=p1,
                        width=10,
                        height = 8, dpi=150)

        # suppressWarnings(htmlwidgets::saveWidget( plotly::partial_bundle(plotly::ggplotly(p1)), paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_RNA.html"), selfcontained = F) )
        suppressWarnings(htmlwidgets::saveWidget( (plotly::ggplotly(p1)), paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_RNA.html"), selfcontained = F) )

      })
      outlist$pseudobulk_RNA <-paste0(RNA.batch.by, "_RNA")
    }

    if ("ADT" %in% names(object@assays)) {
      SeuratObject::DefaultAssay(object) <- "ADT"
      suppressMessages(p1 <- PlotReducedDim(object, plot.type = "pseudobulk", assay = "ADT", group.by.pseudobulk = ADT.batch.by[1], sample.by.pseudobulk = sample.by))
      temp <- lapply(ADT.batch.by, function(x) {
        tempdata <- p1$data
        tempdata[[ADT.batch.by[1]]] <- object@meta.data[[x]][match(rownames(tempdata), object@meta.data[[sample.by]])]
        p1$data <- tempdata
        if (x == sample.by) {
          suppressMessages(p1 <- p1 +
                             ggplot2::scale_fill_manual(values = get_colors(length(unique(metadata[[x]])))) +
                             ggplot2::scale_color_manual(values = get_colors(length(unique(metadata[[x]])))) +
                             ggplot2::labs(fill = x, color = x) +
                             ggplot2::theme(legend.position = "none"))
        } else {
          suppressMessages(p1 <- p1 +
                             ggplot2::scale_fill_manual(values = get_colors(length(unique(metadata[[x]])))) +
                             ggplot2::scale_color_manual(values = get_colors(length(unique(metadata[[x]])))) +
                             ggplot2::labs(fill = x, color = x))
        }

        ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_ADT.png"),
                        plot = p1,
                        width = 10,
                        height = 8, dpi = 150)

        # suppressWarnings(htmlwidgets::saveWidget(plotly::partial_bundle(plotly::ggplotly(p1)),
        #                                          paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_ADT.html"),
        #                                          selfcontained = FALSE))
        suppressWarnings(htmlwidgets::saveWidget((plotly::ggplotly(p1)),
                                                 paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Clustering/pseudobulk/", x, "_ADT.html"),
                                                 selfcontained = FALSE))
      })
      outlist$pseudobulk_ADT <- paste0(ADT.batch.by, "_ADT")
    }

  }



  {
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Covariate Impact"))
    ## CovariateImpact
    if( "RNA" %in% names(object@assays) ){
      message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    "------ Run Batch QC: Covariate Impact RNA data"))
      SeuratObject::DefaultAssay(object) <- "RNA"
      outlist$CV_RNA <- PlotCovariateImpact(object, assay="RNA" ,variables =RNA.batch.covariate, pseudobulk.celltype.by=celltype.by,  return.type = "interactive_table")
    }

    if( "ADT" %in% names(object@assays) ){
      SeuratObject::DefaultAssay(object) <- "ADT"
      if(dim(object)[1] < 10){
        message("--------nADT < 10")
      }else{
        message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      "------ Run Batch QC: Covariate Impact ADT data"))
        outlist$CV_ADT <- PlotCovariateImpact(object, assay="ADT" ,variables =ADT.batch.covariate, pseudobulk.celltype.by=celltype.by,  return.type = "interactive_table")
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
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/cell/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Clustering/cell/"), recursive = TRUE)
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

      ggplot2::ggsave(paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Clustering/cell/",name, "_",batch.by,"_", x,".png"),
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
                           RNA.batch.covariate=c( "orig.ident" ,"percent.mt", "nFeature_RNA", "nCount_RNA"),
                           ADT.batch.covariate=c( "orig.ident" , "nFeature_ADT", "nCount_ADT"),
                           plot_out){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------ Run Batch QC: Feature-level"))
  outlist <- list()
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/"), recursive = TRUE)
  }
  if (!dir.exists(paste0(plot_out,"/SingleCellMQC/excel/Batch_effect/"))) {
    dir.create(paste0(plot_out,"/SingleCellMQC/excel/Batch_effect/"), recursive = TRUE)
  }

  ##
  if( "RNA" %in% names(object@assays) ){

    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Variance explained  per feature (RNA data) "))

    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Pseudobulk RNA"))

    SeuratObject::DefaultAssay(object) <- "RNA"
    out <- RunVarExplainedPerFeature(object, assay = "RNA", variables = RNA.batch.covariate, pseudobulk.sample.by = sample.by, type = "pseudobulk")
    p1 <- PlotVEPerFeature(out, plot.type="density")
    ggplot2::ggsave(
      paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Feature/","RNA_pseudobulk.png"),
      plot = p1,
      width = 8,
      height = 6,
      dpi = 150
    )
    write.csv(out, file = paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/RNA_pseudobulk_VE_perFeature.csv") )
    outlist$CV_RNA <- PlotVEPerFeature(out, plot.type="density", return.type = "interactive_table", csv.name = "RNA_pseudobulk_VE_perFeature")

    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Cell RNA"))
    out <- RunVarExplainedPerFeature(object, assay = "RNA", variables = RNA.batch.covariate, type = "cell")
    p1 <- PlotVEPerFeature(out, plot.type="density")
    ggplot2::ggsave(
      paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Feature/","RNA_cell.png"),
      plot = p1,
      width = 8,
      height = 6,
      dpi = 150
    )
    write.csv(out, file = paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/RNA_cell_VE_perFeature.csv") )
    outlist$Cell_CV_RNA <- PlotVEPerFeature(out, plot.type="density", return.type = "interactive_table", csv.name = "RNA_cell_VE_perFeature")
  }

  if( "ADT" %in% names(object@assays) ){
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Variance explained  per feature (ADT data) "))
    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Pseudobulk ADT"))
    SeuratObject::DefaultAssay(object) <- "ADT"
    out <- RunVarExplainedPerFeature(object, assay = "ADT", variables = ADT.batch.covariate, pseudobulk.sample.by = sample.by, type = "pseudobulk")
    p1 <- PlotVEPerFeature(out, plot.type="density")
    ggplot2::ggsave(
      paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Feature/","ADT_pseudobulk.png"),
      plot = p1,
      width = 8,
      height = 6,
      dpi = 150
    )
    write.csv(out, file = paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/ADT_pseudobulk_VE_perFeature.csv") )
    outlist$CV_ADT <- PlotVEPerFeature(out, plot.type="density", return.type = "interactive_table",csv.name = "ADT_pseudobulk_VE_perFeature")

    message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  "------ Run Batch QC: Cell ADT"))
    out <- RunVarExplainedPerFeature(object, assay = "ADT", variables = ADT.batch.covariate, type = "cell")
    p1 <- PlotVEPerFeature(out, plot.type="density")
    ggplot2::ggsave(
      paste0(plot_out, "/SingleCellMQC/plot/Batch_effect/Feature/","ADT_cell.png"),
      plot = p1,
      width = 8,
      height = 6,
      dpi = 150
    )
    write.csv(out, file = paste0(plot_out,"/SingleCellMQC/plot/Batch_effect/Feature/ADT_cell_VE_perFeature.csv") )
    outlist$Cell_CV_ADT <- PlotVEPerFeature(out, plot.type="density", return.type = "interactive_table", csv.name = "ADT_cell_VE_perFeature")
  }

  return(outlist)
}




report.batch <- function(object, RNA_cluster_name=NULL, ADT_cluster_name=NULL,
                         sample.by="orig.ident",
                         RNA.batch.by= "orig.ident",
                         ADT.batch.by="orig.ident",
                         celltype.by=NULL,
                         RNA.batch.covariate=c( "orig.ident" ,"percent.mt", "nFeature_RNA", "nCount_RNA"),
                         ADT.batch.covariate=c( "orig.ident" , "nFeature_ADT", "nCount_ADT"),
                         plot_out, section=1){
  plot_dir <- paste0(plot_out,"/SingleCellMQC/plot/")
  Batch <- list()

  if(1 %in% section){
    Batch$Section1  <- .batchSection1(object, RNA.batch.by=RNA.batch.by, ADT.batch.by=ADT.batch.by, plot_out=plot_out, sample.by=sample.by,
                                      RNA.batch.covariate=RNA.batch.covariate, ADT.batch.covariate=ADT.batch.covariate, celltype.by=celltype.by
    )
  }

  if(2 %in% section){
    Batch$Section2  <- .batchSection2(object, RNA_cluster_name=RNA_cluster_name, ADT_cluster_name=ADT_cluster_name, RNA.batch.by=RNA.batch.by, ADT.batch.by=ADT.batch.by, plot_out=plot_out, sample.by=sample.by)
  }

  if(3 %in% section){
    Batch$Section3  <- .batchSection3(object,sample.by=sample.by, RNA.batch.covariate=RNA.batch.covariate, ADT.batch.covariate=ADT.batch.covariate,  plot_out=plot_out)
  }
  return(Batch)
}


#' @title Generate QC HTML Report
#'
#' @description This function generates a comprehensive HTML report for single-cell data quality control (QC).
#' The report includes modules for sample QC, cell QC, feature QC, and batch QC, with customizable sections for each module.
#'
#' @param object Seurat object containing single-cell data. This object should be preprocessed using the `RunPreprocess` function or other appropriate preprocessing functions.
#' @param VDJ_data List of VDJ data. Default is NULL.
#' @param sample.by Column name in the metadata indicating sample information. Default is "orig.ident".
#' @param outputFile Path to the output HTML file. Default is "./SingleCellMQC/SingleCellMQC.html".
#' @param color Color palette for plots. Default is NULL.
#' @param section.sample Sections to include in the Sample QC module. Default is `c(1:3)`.
#'   - **`1`**: Run Sample QC: Sample quality assessment.
#'   - **`2`**: Run Sample QC: Outlier sample detection by cell type.
#'   - **`3`**: Run Sample QC: Sample identity verification.
#' @param section.cell Sections to include in the Cell QC module. Default is `c(1:2)`.
#'   - **`1`**: Run Cell QC: Low-quality cell information.
#'   - **`2`**: Run Cell QC: Doublet information.
#' @param section.feature Sections to include in the Feature QC module. Default is `1`.
#'   - **`1`**: Run Feature QC: Feature quality assessment.
#' @param section.batch Sections to include in the Batch QC module. Default is `1:3`.
#'   - **`1`**: Run Batch QC: Sample level.
#'   - **`2`**: Run Batch QC: Cell level.
#'   - **`3`**: Run Batch QC: Feature level.
#' @param RNA_cluster_name Column name in the metadata indicating RNA cluster information. Default is "rna_cluster". Only for: section.batch (2).
#' @param ADT_cluster_name Column name in the metadata indicating ADT cluster information. Default is "adt_cluster". Only for: section.batch (2).
#' @param RNA.batch.by Character vector. Column name(s) in the metadata indicating RNA batch information (categorical variable) or other categorical variables for PCA or UMAP plot. Default is "orig.ident". Only for: section.batch (1, 2).
#' @param ADT.batch.by Character vector. Column name(s) in the metadata indicating ADT batch information (categorical variable) or other categorical variables for PCA or UMAP plot. Default is "orig.ident". Only for: section.batch (1, 2).
#' @param RNA.batch.covariate Character vector. Covariates (can be categorical or continuous) for RNA `SVD` and `Variance Explained` analysis. Default is `c("orig.ident", "percent.mt", "nFeature_RNA", "nCount_RNA")`. section.batch (1, 3).
#' @param ADT.batch.covariate Character vector. Covariates (can be categorical or continuous) for ADT `SVD` and `Variance Explained` analysis. Default is `c("orig.ident", "nFeature_ADT", "nCount_ADT")`. Only for: section.batch (1, 3).
#' @param celltype.by Column name in the metadata indicating cell type information. Default is "ScType". Only for: section.sample (2), section.batch (1).
#' @param tissue Tissue name. Default is NULL. Only for: section.sample (1, 2).
#'
#' @return A HTML report is generated and saved to the specified output file path. No value is returned.
#' @export
#'
RunReport <- function(object=NULL, VDJ_data=NULL, sample.by="orig.ident", outputFile="./SingleCellMQC/SingleCellMQC.html", color=NULL,
                      section.sample = c(1:3), section.cell = c(1:2), section.feature=1,section.batch=c(1:3),
                      RNA_cluster_name="rna_cluster", ADT_cluster_name="adt_cluster", RNA.batch.by= "orig.ident", ADT.batch.by="orig.ident",
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

  if( !is.null(section.sample) ){
    out_list$Sample <- report.sample(object=object, VDJ_data=VDJ_data, sample.by=sample.by, color=color,  plot_out=plot_out, section= section.sample,celltype.by = celltype.by, tissue=tissue)
  }

  if( !is.null(section.cell) ){
    out_list$Cell <- report.cell(object=object, plot_out=plot_out,sample.by=sample.by, section = section.cell)
  }

  if( !is.null(section.feature) ){
    out_list$Feature <- report.feature(object=object, marker_name= c("CD3D","CD3E","CD19", "MS4A1","CD79A",
                                                                     "CD14","FCGR3A","CD68","FCN1","ITGAX"), celltype.by=celltype.by, sample.by=sample.by, plot_out=plot_out, section=1)
  }

  if( !is.null(section.batch) ){
    out_list$Batch <- report.batch(object, RNA_cluster_name=RNA_cluster_name, ADT_cluster_name=ADT_cluster_name,
                                   sample.by=sample.by,
                                   RNA.batch.by= RNA.batch.by, ADT.batch.by=ADT.batch.by,
                                   celltype.by=celltype.by,
                                   RNA.batch.covariate=RNA.batch.covariate,
                                   ADT.batch.covariate=ADT.batch.covariate,
                                   plot_out, section=section.batch)
  }

  # return(out_list)
  cat("\n")
  message("Making QC Report!")
  cat("Output report to:", outputFile, "\n")
  reportRmd <- system.file("report/0_main.Rmd",
                           package = "SingleCellMQC")
  css_file <- system.file("report", "bootstrap.css", package = "SingleCellMQC")
  target_dir <- plot_out
  file.copy(css_file, target_dir, overwrite = TRUE)

  rmarkdown::render(reportRmd, output_dir = plot_out, output_file = output_file, params = list(data = out_list, cell_dir=plot_dir))
  return(NULL)
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
#' @examples
#' \dontrun{
#' # Load a Seurat object
#' data("pbmc_small", package = "Seurat")
#'
#' # Run the preprocessing pipeline
#' pbmc_small <- RunPreprocess(pbmc_small, sample.by = "orig.ident")
#'
#' # View the updated Seurat object
#' print(pbmc_small)
#' }
#'
#' @export
RunPreprocess <- function(object, sample.by="orig.ident", db.method=c("scDblFinder", "hybrid"), lq.method=c("MAD") , preprocess = "rna.umap", ... ){
  object <- CalculateMetrics(object)
  object <- RunScType(object, split.by = sample.by)
  object <- RunDbt(object, methods = db.method, add.Seurat = T,split.by = sample.by)
  object <- RunLQ(object, methods = lq.method, add.Seurat = T,split.by = sample.by)
  object <- RunPipeline(object, preprocess = preprocess,...)
  return(object)
}
