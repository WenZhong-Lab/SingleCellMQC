
#' Low-Quality Cells Detection
#'
#' This function handles multiple quality assessment methods to identify low-quality cells in single-cell data.
#'
#' @param object A Seurat object.
#' @param methods A character vector specifying the quality assessment methods to apply (one or more). Available method options are:
#' \itemize{
#'   \item "MAD": Outliers based on the median absolute deviation. Suitable for RNA and ADT.
#'   \item "ddqc": Only RNA. An unsupervised adaptive QC framework to perform flexible and data-driven QC at the level of cell types while retaining critical biological insights and improved power for downstream analysis. (Subramanian, Ayshwarya et al., Genome biology, 2022, PMID: 36575523)
#'   \item "miQC": Only RNA. An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data. (Hippen, Ariel A et al., PLoS computational biology, 2021, PMID: 34428202)
#'   \item "fixed": Fixed cutoffs for RNA and ADT metrics.
#'   }
#'   The default is all methods.
#' @param add.Seurat Logical; if TRUE, adds the low-quality cell detection results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the low-quality cell detection results.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param ... Additional arguments.
#'
#' @references
#' 1.Hippen, Ariel A et al. “miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.” PLoS computational biology vol. 17,8 e1009290. 24 Aug. 2021, doi:10.1371/journal.pcbi.1009290
#'
#' 2.Subramanian, Ayshwarya et al. “Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics.” Genome biology vol. 23,1 267. 27 Dec. 2022, doi:10.1186/s13059-022-02820-w
#'
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with the low-quality cell detection results.
#' @details The 'Fail' value represents low-quality cells.
#'
#' @examples
#' \dontrun{
#' # Assuming 'seuratObj' is a Seurat object prepared for low-quality cell detection:
#' seuratObj <- RunLQ(object = seuratObj, methods = c("MAD", "ddqc", "miQC"), add.Seurat = TRUE)
#' }
#'
#' @seealso \code{\link{RunLQ_MAD}}, \code{\link{RunLQ_miQC}}, \code{\link{RunLQ_ddqc}}
#' @export
RunLQ <- function(object, methods=c("ddqc", "miQC", "MAD"), add.Seurat = T,split.by=NULL, ...) {
  if (!is.character(methods)) {
    stop("methods should be a character vector.")
  }

  valid_methods <- c("ddqc", "miQC", "MAD", "fixed")
  if (!all(methods %in% valid_methods)) {
    stop("Unknown method(s) found. Valid methods are: ", paste(valid_methods, collapse = ", "))
  }

  results <- lapply(methods, function(method) {
    func_name <- paste0("RunLQ_", method)
    func <- match.fun(func_name)
    if(method == "fixed"){
      out <- func(object, add.Seurat=F, ...)
    }else{
      out <- func(object ,add.Seurat=F, split.by=split.by, ...)
    }
    out <- data.table(cell=rownames(out), out)
    return(out)
  }
  )

  if (length(methods) > 1) {
    out <- Reduce(function(x, y){
      merge_data <- merge(x, y, by = "cell", all = TRUE)
    }, results)
  } else {
    out <- results[[1]]
  }
  out <- data.frame(out)
  rownames(out) <- out[, 1]
  out <- out[, -1, drop=F]

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}


#' Low quality cells detection based on fixed cutoff
#'
#' Identifies low-quality cells based on the fixed cutoff on a Seurat object.
#'
#' @param object A Seurat object to be filtered.
#' @param min.nFeature_RNA (Optional) Minimum number of RNA features per cell.
#' @param min.nCount_RNA (Optional) Minimum RNA count per cell.
#' @param percent.mt (Optional) Maximum percentage of mitochondrial genes per cell.
#' @param percent.rb (Optional) Maximum percentage of ribosomal genes per cell.
#' @param percent.hb (Optional) Maximum percentage of hemoglobin genes per cell.
#' @param min.nFeature_ADT (Optional) Minimum number of ADT features per cell.
#' @param min.nCount_ADT (Optional) Minimum ADT count per cell.
#' @param add.Seurat (Optional) Logical. If TRUE, the final filtering result is added
#'                   to the Seurat object as a new metadata column named `db_fixed`.
#'                   Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return If `add.Seurat` is FALSE, returns a vector of "Pass" or "Fail" for each cell
#'         based on the filtering criteria. If `add.Seurat` is TRUE, returns the modified
#'         Seurat object with the final filtering result added as a new metadata column
#'         named `db_fixed`.
#'
#' @export
RunLQ_fixed <- function(object, min.nFeature_RNA=NULL,  min.nCount_RNA=NULL,
                        percent.mt=NULL, percent.rb=NULL, percent.hb=NULL,
                        min.nFeature_ADT=NULL,  min.nCount_ADT=NULL,
                        add.Seurat = TRUE,...
){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  # Initialize the output data frame
  out <- data.frame(row.names = colnames(object))

  # Filter based on nFeature_RNA
  if(!is.null(min.nFeature_RNA)){
    out$min.nFeature_RNA <- ifelse(object$nFeature_RNA >= min.nFeature_RNA, "Pass", "Fail")
  }

  # Filter based on nCount_RNA
  if(!is.null(min.nCount_RNA)){
    out$min.nCount_RNA <- ifelse(object$nCount_RNA >= min.nCount_RNA, "Pass", "Fail")
  }


  # Filter based on percent.mt
  if(!is.null(percent.mt)){
    out$percent.mt <- ifelse(object$percent.mt <= percent.mt, "Pass", "Fail")
  }

  # Filter based on percent.rb
  if(!is.null(percent.rb)){
    out$percent.rb <- ifelse(object$percent.rb <= percent.rb, "Pass", "Fail")
  }

  # Filter based on percent.hb
  if(!is.null(percent.hb)){
    out$percent.hb <- ifelse(object$percent.hb <= percent.hb, "Pass", "Fail")
  }

  # Filter based on nFeature_ADT
  if(!is.null(min.nFeature_ADT)){
    out$min.nFeature_ADT <- ifelse(object$nFeature_ADT >= min.nFeature_ADT, "Pass", "Fail")
  }

  # Filter based on nCount_ADT
  if(!is.null(min.nCount_ADT)){
    out$min.nCount_ADT <- ifelse(object$nCount_ADT >= min.nCount_ADT, "Pass", "Fail")
  }

  # Combine all results into a single "Pass" or "Fail" column
  out$lq_fixed <- ifelse(rowSums(out == "Fail") == 0, "Pass", "Fail")

  # If add.Seurat is TRUE, add the lq_fixed result to the Seurat object

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out[,"lq_fixed", drop = FALSE])
    return(object)
  } else {
    return(out[,"lq_fixed", drop = FALSE])
  }
}





RunOutlier <-
  function(object,
           method = "mad",#sd
           nmads  = 3,
           nsds = 3 ,
           constant = 1.4826,
           log = F,
           group.by = NULL,
           share.median = F,
           share.mean=F,
           share.mad = F,
           share.sd = F,
           metrics = NULL,
           type = c("both", "lower", "higher"),
           ...) {

    if (!method %in% c("mad", "sd")) {
      stop("Invalid method. Please choose either 'mad' or 'sd'.")
    }
    median_dif=mean_dif=NULL
    if (is.null(metrics)) {
      object <- data.frame(metrics = object)
      metrics <- "metrics"
    }
    if (log) {
      object[[metrics]] <- log(object[[metrics]])
    }

    if(is.null(group.by)){
      data <- split(object,as.character( rep("metrics", length(object[[metrics]]))))
    }else{
      data <- split(object,as.character( object[[group.by]]))
    }

    out <- lapply(data, function(x){
      if("mad" %in% method){
        if( !is.null(group.by) &  share.median){
          stat_value <- median(object[[metrics]],na.rm =T)
        }else{
          stat_value <- median(x[[metrics]],na.rm =T)
        }
        if(!is.null(group.by) & share.mad){
          stat_dif <- data.table(object[,c(metrics, group.by),drop=F])[, median_dif:=  get(metrics)- median(get(metrics)),by=list(get(group.by)) ]$median_dif
        }else{
          stat_dif <- x[[metrics]]-stat_value
        }
        cal_value <- median(abs( na.omit(stat_dif) )) * constant

        nds = nmads
      }else{
        if(!is.null(group.by) & share.mean){
          stat_value <- mean(object[[metrics]],na.rm=T)
        }else{
          stat_value <- mean(x[[metrics]],na.rm =T)
        }
        if(!is.null(group.by) & share.sd){
          stat_dif <- data.table(object[,c(metrics, group.by),drop=F])[, mean_dif:=  (get(metrics)- mean(get(metrics)))^2 ,by=list(get(group.by)) ]$mean_diff
        }else{
          stat_dif <- na.omit(x[[metrics]]-stat_value)
        }
        cal_value <- sqrt(sum(stat_dif)/ (length(stat_dif)-1))

        nds = nsds
      }
      cutoff <- nds * cal_value
      cut_low <- stat_value - cutoff
      cut_high <-  stat_value + cutoff


      if(type[1] == "lower"){
        filtered_x <- x[[metrics]] < cut_low
        if(log){
          cut_low <- exp(cut_low)
        }
        cut_high <- Inf
      }else if(type[1] == "higher"){
        filtered_x <- x[[metrics]] > cut_high
        cut_low <- Inf
        if (log) {
          cut_high <- exp(cut_high)
        }
      }else{
        filtered_x <- abs(x[[metrics]] - stat_value) < cutoff
        if (log) {
          cut_low <- exp(cut_low)
          cut_high <- exp(cut_high)
        }
      }

      value <- c("lower"=cut_low, "higher"=cut_high)
      attr(filtered_x, "thresholds") <- value
      names(filtered_x) <- rownames(x)
      return(filtered_x)
    })

    thresholds <- do.call(rbind, lapply(out, function(x){
      attr(x, "thresholds")
    }))
    names(out) <- NULL
    out  <- do.call(c, out)
    attr(out, "thresholds") <- thresholds
    return(out)
  }



.get_named_column <- function(df, col_name) {
  result <- df[[col_name]]
  names(result) <- rownames(df)
  return(result)
}

#' Low quality cells detection with MAD
#'
#' Identifies low-quality cells based on the median absolute deviation (MAD) method on a Seurat object.
#'
#' @param object Seurat object or data frame containing metrics for quality control. RNA metrics should contain "nCount_RNA", "nFeature_RNA" and "percent.mt". If you have ADT data, it should contain "nCount_ADT", "nFeature_ADT" and "percent.isotype". These metrics can all be generated via the \code{\link{CalculateMetrics}} function.
#' @param nmads Numeric. The number of MADs to use for identifying outliers. Default is 3.
#' @param split.by Optional; Column name by which to group the data before calculating statistics.
#' @param share.medians Logical; if TRUE, forces all groups to have the same median. Default is FALSE.
#' @param share.mads Logical; if TRUE, forces all groups to have the same MAD. Default is FALSE.
#' @param add.Seurat Logical; if TRUE, adds the MAD results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the MAD results.
#' @param ... Additional arguments to be passed to \code{\link[scater]{isOutlier}} function from the `scater` package.
#'
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with outlier results.
#'
#' @export
RunLQ_MAD <- function(object, nmads  = 3, split.by=NULL, add.Seurat=T, share.medians=F, share.mads=F,  ...) {
  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }

  known_params <- names(formals(scater::isOutlier))
  filtered_args <- list(...)
  filtered_args <- setdiff(filtered_args[names(filtered_args) %in% known_params],c("log", "batch", "type", "nmads") )

  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "--------Low-quality method: MAD --------"))
  cat("--------RNA --------")
  if( length(intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(metadata)) )!=3 ){
    stop("metadata colnames not obtain `nCount_RNA`, `nFeature_RNA`, `percent.mt`")
  }

  nCount_RNA <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nCount_RNA"),
                                                  "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads, type="lower"), filtered_args))
  nFeature_RNA <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nFeature_RNA"),
                                                    "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads, type="lower"), filtered_args))
  mt <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="percent.mt"),
                                          "batch"=metadata[[split.by]], "log"=F, "nmads"=nmads, type="higher"), filtered_args))
  mad_result <- data.frame(nCount_RNA, nFeature_RNA, mt)
  colnames(mad_result) <- c( paste0(c("lq_RNA_umi_","lq_RNA_gene_","lq_RNA_mt_"), nmads, c("mad")) )

  if( length(intersect(c("nFeature_ADT", "nCount_ADT", "percent.isotype"), colnames(metadata))==2) ){
    cat("--------ADT --------")
    nCount_ADT <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nCount_ADT"),
                                                    "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads, type="lower"), filtered_args))
    nFeature_ADT <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nFeature_ADT"),
                                                      "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads, type="lower"), filtered_args))
    isotype <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="percent.isotype"),
                                                 "batch"=metadata[[split.by]], "log"=F, "nmads"=nmads, type="higher"), filtered_args))
    mad_result <- cbind(mad_result, nCount_ADT, nFeature_ADT, isotype)
    colnames(mad_result) <- c( paste0(c("lq_RNA_umi_","lq_RNA_gene_","lq_RNA_mt_", "lq_ADT_umi_","lq_ADT_pro_", "lq_ADT_isotype_"), nmads, c("mad")) )
  }

  # mad_result[[paste0("lq_RNA_", nmads, "mad")]] <- nCount_RNA | nFeature_RNA | mt
  mad_result <- as.data.frame(lapply(mad_result, function(x) {
    ifelse(x == TRUE, "Fail", "Pass")
  }))
  rownames(mad_result) <- names(nCount_RNA)
  mad_result <- mad_result[match(rownames(metadata), names(nCount_RNA) ), ]

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = mad_result)
    return(object)
  }

  return(mad_result)
}


#' Low quality cells detection with miQC
#'
#' Performs quality control using the miQC method on a Seurat object, an adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.
#'
#' @param object A Seurat object.
#' @param model_type (character) What type of model to generate. A linear mixture model ("linear") based on mitochondrial percentage and library complexity is recommended. B-spline ("spline") and two-degree polynomial ("polynomial") models are also supported. For a simpler model, a one-dimensional gaussian mixture model ("one_dimensional") based on mitochondrial percentage only is available. Default = "linear".
#' @param name.nFeature_RNA The name of the column in metadata representing the number of detected features per cell.
#' @param name.percent.mt The name of the column in metadata representing MT pct per cell.
#' @param add.Seurat Logical; if TRUE, adds the miQC results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the miQC results.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param random.seed Random seed for reproducibility.
#' @param ... Additional arguments to be passed to \code{\link[miQC]{filterCells}} function from the `miQC` package.
#' @param BPtmpdir Temporary directory for BPCells matrix processing.
#'
#' @return A Seurat object with appended miQC metadata if `add.Seurat` is TRUE; otherwise, a data frame with the miQC results.
#' @details The 'Fail' value represents low-quality cells.
#' @examples
#' \dontrun{
#' # Assuming `seuratObject` is a Seurat object
#' seuratObject <- RunLQ_miQC(seuratObject, add.Seurat = TRUE)
#' }
#' @seealso \code{\link{RunLQ}}
#' @export
#' @references
#' Hippen, Ariel A et al. “miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.” PLoS computational biology vol. 17,8 e1009290. 24 Aug. 2021, doi:10.1371/journal.pcbi.1009290
#'
RunLQ_miQC <- function(object, model_type = "linear",  split.by= "orig.ident",
                       name.nFeature_RNA = "nFeature_RNA", name.percent.mt = "percent.mt", add.Seurat = TRUE,
                       random.seed=1 ,BPtmpdir="./temp/SingleCellMQC_miQC/", ...) {
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("miQC", quietly = TRUE)){
    stop("Error: Please install the miQC package first.")
  }

  known_params <- names(formals(miQC::filterCells))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]


  metadata <- object@meta.data[, c(name.nFeature_RNA, name.percent.mt), drop = FALSE]
  colnames(metadata)[1:2] <- c("detected", "subsets_mito_percent")
  object <- Seurat::AddMetaData(object, metadata = metadata)

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Low-quality detection (miQC)"))


  p <- progressr::progressor(along = 1:length(split_object))

  out <- lapply( split_object, function(x){
    set.seed(random.seed)
    sce <- suppressWarnings(Seurat::as.SingleCellExperiment(x, assay = "RNA"))
    model <- miQC::mixtureModel(sce,model_type=model_type)

    out <- do.call(miQC::filterCells, c(list("sce" = sce, "model"=model, "verbose"=F), filtered_args))
    index <- colnames(out)

    lq_miQC <- data.frame( lq_miQC = rep("Fail", dim(x)[2] ) )
    rownames(lq_miQC) <- rownames(x@meta.data)
    lq_miQC$lq_miQC[ match(index,rownames(lq_miQC))  ] <- "Pass"

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(x, sample.by = split.by), ": ", sum( lq_miQC$lq_miQC=="Fail" ), " fail cells, ", sum( lq_miQC$lq_miQC=="Pass" ), " pass cells" ))
    rm(x)
    gc()
    p()

    return(lq_miQC)
  }
  )
  names(out) <- NULL
  out <- do.call(rbind, out)

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}



#' Low quality cells detection with ddqc
#'
#' Performs quality control using the ddqc method on a Seurat object, an unsupervised adaptive QC framework.
#'
#' @param object A Seurat object.
#' @param add.Seurat Logical; if TRUE, adds the ddqc results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the ddqc results.
#' @param ... Additional arguments to be passed to ddqcR::ddqc.metrics function.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param BPtmpdir Temporary directory for BPCells matrix processing.
#'
#' @return A Seurat object with appended ddqc metadata if `add.Seurat` is TRUE; otherwise, a data frame with the ddqc results.
#' @details The 'Fail' value represents low-quality cells.
#'
#' @examples
#' \dontrun{
#' # Assuming `seuratObject` is a Seurat object
#' seuratObject <- RunLQ_ddqc(seuratObject, add.Seurat = TRUE)
#' }
#' @seealso \code{\link{RunLQ}}
#' @export
#' @references
#' Subramanian, Ayshwarya et al. “Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics.” Genome biology vol. 23,1 267. 27 Dec. 2022, doi:10.1186/s13059-022-02820-w
#'

RunLQ_ddqc <- function(object,add.Seurat=TRUE,split.by="orig.ident", BPtmpdir="./temp/SingleCellMQC_ddqc/", ...) {
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("ddqcR", quietly = TRUE)){
    stop("Error: Please install the ddqcR package first.")
  }


  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)

  # split_object <- Seurat::SplitObject(object = object, split.by = split.by)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Low-quality detection (ddqc)"))


  known_params <- formals(ddqcR::ddqc.metrics)
  known_params$data <-NULL
  filtered_args <- list(...)
  updated_params <- .replace_params(known_params, filtered_args)

  p <- progressr::progressor(along = 1:length(split_object))

  out <- lapply( split_object , function(x){
    object_sample <- x
    SeuratObject::DefaultAssay(object_sample) <- "RNA"

    if( is(SeuratObject::GetAssayData(object_sample, assay = "RNA", slot = "counts"), "IterableMatrix" ) ){
      object_sample <- SeuratObject::SetAssayData(object_sample, assay = "RNA", slot = "counts",
                                                  new.data = as(object = SeuratObject::GetAssayData(object_sample, assay = "RNA", slot = "counts"), Class = "dgCMatrix") )
    }

    result <- data.frame(lq_ddqc=rep(FALSE, length(rownames(object_sample@meta.data)) ))
    rownames(result) <- rownames(object_sample@meta.data)
    object_sample <- ddqcR::initialQC(object_sample)
    # out <- !(do.call(ddqcR::ddqc.metrics, c(list("data" = object_sample), filtered_args))[, "passed.qc", drop = F])
    out <- !(ddqcR::ddqc.metrics(data = object_sample, res = known_params$res, threshold = known_params$threshold,
                                 do.counts = known_params$do.counts, do.genes = known_params$do.genes, do.mito = known_params$do.mito,
                                 do.ribo = known_params$do.ribo,n.genes.lower.bound = known_params$n.genes.lower.bound,
                                 percent.mito.upper.bound = known_params$percent.mito.upper.bound, random.state = known_params$random.state
    )[, "passed.qc", drop = F])
    out <- as.data.frame(out)
    colnames(out) <- "lq_ddqc"

    result[match(rownames(out),rownames(result)), 1] <- out$lq_ddqc

    result$lq_ddqc <- ifelse(result$lq_ddqc == TRUE , "Fail", "Pass")

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_sample, sample.by = split.by), ": ", sum( result$lq_ddqc=="Fail" ), " fail cells, ", sum( result$lq_ddqc=="Pass" ), " pass cells" ))
    rm(object_sample)
    gc()
    p()

    return(result)

  } )

  names(out) <- NULL
  out <- do.call(rbind, out)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}


# DropletUtils ------------------------------------------------------------

#The NA value represents an empty droplet with count < 100, which can be customized with the lower value
#FDR>0.05 is considered to be the captured droplet
#' Run empty droplet identification
#'
#' This function identifies empty droplets from single-cell RNA-seq data using the DropletUtils package.
#'
#'
#' @param dir_name Character. The path to the directory containing 10X count data.
#' @param lower Numeric. The lower UMI count threshold to consider a droplet as non-empty. Default is 100.
#' @param FDR Numeric. The false discovery rate threshold for identifying retained droplets. Default is 0.01.
#'
#' @return A SingleCellExperiment object containing only the retained droplets that passed the thresholding criteria.
#' @details
#' The NA value represents an empty droplet with count < 100, which can be customized with the `lower` parameter.
#' FDR > 0.05 is considered to be the captured droplet.
#'
#' @export
RunEmptyDroplet<-function(dir_name=NULL,lower=100,FDR=0.01){

  if(!requireNamespace("DropletUtils", quietly = TRUE)){
    stop("Error: Please install the DropletUtils package first.")
  }

  sce <- DropletUtils::read10xCounts(dir_name)
  bcrank1 <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(sce))

  # Only showing unique points for plotting speed.
  uniq1 <- !duplicated(bcrank1$rank)
  plot(bcrank1$rank[uniq1], bcrank1$total[uniq1], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  graphics::abline(h=S4Vectors::metadata(bcrank1)$inflection, col="darkgreen", lty=2)
  graphics::abline(h=S4Vectors::metadata(bcrank1)$knee, col="dodgerblue", lty=2)
  graphics::legend("bottomleft", legend=c("Inflection", "Knee"),
                   col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

  set.seed(100)
  e.out <- DropletUtils::emptyDrops(SingleCellExperiment::counts(sce),lower = lower)
  summary_table<-summary(e.out$FDR <= FDR)


  retained <- sce[,which(e.out$FDR <= FDR)]
  bcrank2 <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(retained))
  uniq2 <- !duplicated(bcrank2$rank)
  plot(bcrank2$rank[uniq2], bcrank2$total[uniq2], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  print(summary_table)
  return(retained)
}

