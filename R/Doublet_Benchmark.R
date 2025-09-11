
clusterBench <- function(pca_value, clusters,  filter_columns){
  sil_value <- bluster::approxSilhouette(pca_value, clusters = clusters)
  sil_value <- mean(sil_value$width)
  ch_value <- fpc::calinhara(x=pca_value, clustering= as.numeric(clusters), cn=length(unique(as.numeric(clusters))) )
  db_value <- clusterSim::index.DB(x=pca_value, cl= as.numeric(clusters) )$DB
  return(data.frame(method=filter_columns, asw=sil_value, ch= ch_value, db=db_value))
}

benchmark <- function(object, filter_columns, resolution=1, type="db"){
  if(filter_columns!="no_filter"){
    utils::capture.output(object <- FilterCells(object, filter_columns =filter_columns ))
  }

  # pipeline
  object <- Seurat::NormalizeData(object = object, verbose = F)
  object <- Seurat::FindVariableFeatures(object = object, verbose = F)
  object <- Seurat::ScaleData(object = object, verbose = F)
  object <- Seurat::RunPCA(object = object, verbose = F)
  object <- Seurat::FindNeighbors(object = object, dims = 1:20, verbose = F)
  object <- Seurat::FindClusters(object = object, resolution =resolution, cluster.name = paste0("cluster_", resolution), verbose = F )

  out <- lapply(paste0("cluster_", resolution), function(cluster_name){
    cluste_label =object@meta.data[[cluster_name]]
    pca_value = SeuratObject::Embeddings(object, reduction = "pca")
    value <- clusterBench(pca_value, cluste_label, filter_columns=filter_columns )
    if(type=="db"){
      if(filter_columns!="no_filter"){
        object <- suppressWarnings(RunScType(object, group.by=cluster_name, data_source = "Main", return.name="filter_ano"))
        ari <- aricode::ARI(c1 = object$ScType, c2 = object$filter_ano)
        nmi <- aricode::NMI(c1 = object$ScType, c2 = object$filter_ano)
        value$ari <- ari
        value$nmi <- nmi
      }else{
        value$ari <- NA
        value$nmi <- NA
      }
    }
    return(value)
  })
  out <- do.call(rbind, out)
  out <- data.frame(resolution=resolution,out )
  return(out)
}

#' Benchmarking of doublet detection methods
#'
#' Benchmarking of RNA-based doublet detection methods
#'
#' @param object A Seurat object containing single-cell data.
#' @param split.by String specifying metadata column to split object by. Default is "orig.ident".
#' @param method_columns Character vector of method names (column names) for doublet detection.
#' @param resolution A numeric vector of cluster resolutions to benchmark. Default is \code{seq(1, 1.5, 0.1)}.
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellBenchmark/".
#'
#' @returns A data.frame with one row per sample/method/resolution, containing clustering metric results.
#' @export
#'
RunBenchmarkDoublet <- function(object, split.by="orig.ident", method_columns,resolution=seq(1,1.5,0.1),
                                   BPtmpdir= "./temp/SingleCellMQC_BPCellBenchmark/"){
  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  method_columns <- c("no_filter", method_columns)
  p <- progressr::progressor(along = 1:length(split_object))

  Benchmark_out <- smart_lapply(split_object, function(x){
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(x, assay = "RNA", slot = "counts")), "package") ) {
      x <- SeuratObject::SetAssayData(object = x,
                                      assay = "RNA",
                                      slot = "counts",
                                      new.data = as(Seurat::GetAssayData(x, assay = "RNA", slot = "counts"), "dgCMatrix") )
    }

    out  <- lapply(method_columns, function(y){
      benchmark(x, filter_columns = y, resolution = resolution, type="db")
    } )
    out <- do.call(rbind, out)
    p()
    gc()
    return(out)
  })
  Benchmark_out <- do.call(rbind, Benchmark_out)
  Benchmark_out <- data.frame( sample=rep(names(split_object), each=length(method_columns)*length(resolution) ),  Benchmark_out)

  return(Benchmark_out)
}

rankBenchValue <- function(object) {
  # Validate input
  if (!is.data.frame(object)) {
    stop("Input must be a data.frame")
  }

  required_cols <- c("sample", "method", "asw", "ch", "db", "ari", "nmi")
  if (!all(required_cols %in% colnames(object))) {
    stop("Input data.frame must contain all required columns: ",
         paste(required_cols, collapse = ", "))
  }

  # Step 1: Calculate mean values by sample and method
  mean_values <- object %>%
    dplyr::group_by(sample, method) %>%
    dplyr::summarise(
      dplyr::across(c(asw, ch, db, ari, nmi),
                    ~mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Step 2: Rank metrics within each sample
  ranked <- mean_values %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      # Higher values get better rank (1 = best)
      asw_rank = rank(-asw, ties.method = "min"),
      ch_rank = rank(-ch, ties.method = "min"),
      ari_rank = rank(-ari, ties.method = "min"),
      nmi_rank = rank(-nmi, ties.method = "min"),
      # For DB, lower values get better rank (1 = best)
      db_rank = rank(db, ties.method = "min"),
      # Calculate average rank
      overall = rowMeans(dplyr::across(dplyr::ends_with("_rank")))
    )
  ranked$overall[ranked$method %in% "no_filter"]<-NA
  ranked <- ranked  %>%
    dplyr::mutate(
      # Rank overall within each sample group
      overall_rank = rank(overall, ties.method = "min")
    ) %>%
    dplyr::ungroup()

  # Final output
  final_ranks <- ranked %>%
    dplyr::arrange(sample, overall)
  final_ranks$ari[final_ranks$method %in% "no_filter"]<-NA
  final_ranks$nmi[final_ranks$method %in% "no_filter"]<-NA
  final_ranks$ari_rank[final_ranks$method %in% "no_filter"]<-NA
  final_ranks$nmi_rank[final_ranks$method %in% "no_filter"]<-NA
  final_ranks$overall[final_ranks$method %in% "no_filter"]<-NA
  final_ranks$overall_rank[final_ranks$method %in% "no_filter"]<-NA


  return(final_ranks)
}


#' Visualize Benchmarking Rank Distribution
#'
#' Plots the distribution of overall ranks for different doublet detection methods,
#' showing how often each rank is achieved per method. Returns both the barplot
#' and the summarized mean metric table.
#'
#' @param object A data.frame, as produced by \code{RunBenchmarkDoublet}.
#' @param plot.ncol Number of columns in the facet plot grid. Default is 5.
#'
#' @return A list with two elements:
#'   \item{plot}{A ggplot2 object of the barplot.}
#'   \item{mean_value}{A data.frame of summarized metrics/ranks.}
#' @export
#'
PlotBenchmarkDoublet <- function(object, plot.ncol=5){
  ave_value  <- rankBenchValue(object)
  df_counts <- ave_value %>%
    dplyr::group_by(method, overall_rank) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  df_counts <- df_counts[df_counts$method!="no_filter",]
  df_counts$overall_rank <- as.character(df_counts$overall_rank)
  p1 <- ggplot2::ggplot(df_counts, ggplot2::aes(x = overall_rank, y = count, fill=overall_rank)) +
    ggplot2::theme_classic(base_size = 14)+
    ggplot2::scale_fill_manual(values = get_colors(length( unique(df_counts$overall_rank)) ))+
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~ method, ncol = plot.ncol) +
    ggplot2::labs(x = "Overal Rank", y = "Count")+
    ggplot2::theme(text = ggplot2::element_text(color = "black"),
                   axis.text =  ggplot2::element_text(color = "black"),
                   legend.position = "none")
  return(list(plot=p1, mean_value = ave_value))
}



